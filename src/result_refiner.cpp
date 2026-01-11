#include "refinement/result_refiner.h"
#include "core/differentiation.h"
#include <cmath>
#include <algorithm>
#include <set>
#include <limits>

#ifdef ENABLE_HIGH_PRECISION
#include "hp/result_refiner_hp.h"
#include "hp/polynomial_hp.h"
#include "hp/precision_context.h"
#include "hp/high_precision_types.h"
#include "hp/precision_conversion.h"
#endif

namespace polynomial_solver {

ResultRefiner::ResultRefiner() {
}

ResultRefiner::~ResultRefiner() {
}

RefinementResult ResultRefiner::refine(
    const SubdivisionSolverResult& solver_result,
    const PolynomialSystem& original_system,
    const RefinementConfig& config) const
{
    RefinementResult result;

    const std::size_t dim = original_system.dimension();
    const std::size_t num_resolved = solver_result.num_resolved;

    // Only support 1D for now
    if (dim != 1) {
        // For 2D+, just return empty result
        // (2D may have infinite roots in degenerate regions)
        return result;
    }

    const Polynomial& poly = original_system.equation(0);

    // If no resolved boxes, try to handle unresolved boxes
    if (num_resolved == 0) {
        // Add all unresolved boxes to unverified list
        for (std::size_t i = 0; i < solver_result.boxes.size(); ++i) {
            result.unverified_boxes.push_back(i);
        }

        // Merge unverified boxes into problematic regions (1D only)
        if (!result.unverified_boxes.empty()) {
            result.problematic_regions = mergeUnverifiedBoxes1D(
                result.unverified_boxes, solver_result, poly, config);
        }

        return result;
    }

    // Step 1: Refine each resolved box using Newton's method
    std::vector<std::size_t> verified_indices;
    std::vector<RefinedRoot> candidate_roots;

    for (std::size_t i = 0; i < num_resolved; ++i) {
        const SubdivisionBoxResult& box = solver_result.boxes[i];

        // Refine using Newton's method
        double refined_location;
        double residual;

        if (!refineRoot1D(box, poly, config, refined_location, residual)) {
            result.unverified_boxes.push_back(i);
            continue;
        }

        // Check residual tolerance
        if (std::abs(residual) >= config.residual_tolerance) {
            result.unverified_boxes.push_back(i);
            continue;
        }

        // Verify root and determine multiplicity at refined location
        std::vector<double> point{refined_location};
        double first_nonzero_deriv = 0.0;
        unsigned int mult = estimateMultiplicity(
            point, original_system, config.max_multiplicity, 1e-10, first_nonzero_deriv);

        // Compute exclusion radius based on multiplicity and derivative
        double exclusion_radius = computeExclusionRadiusFromDerivative(
            mult, first_nonzero_deriv, config.target_tolerance, config.exclusion_multiplier);

        // Estimate condition number
        double condition_estimate = estimateConditionNumber1D(
            refined_location, poly, first_nonzero_deriv);

        // Estimate actual error from condition number and residual
        // |error| ≈ κ * |residual| / |f'(x)|
        double estimated_error = condition_estimate * std::abs(residual) / std::max(std::abs(first_nonzero_deriv), 1e-14);

        // Flag if higher precision is needed
        // Threshold: if estimated error > 1e-10, double precision may be insufficient
        bool needs_higher_precision = (estimated_error > 1e-10);

        // Create refined root
        RefinedRoot refined;
        refined.location = point;
        refined.residual = std::vector<double>{residual};
        refined.max_error = std::vector<double>{estimated_error};
        refined.multiplicity = mult;
        refined.first_nonzero_derivative = first_nonzero_deriv;
        refined.exclusion_radius = exclusion_radius;
        refined.source_boxes.push_back(i);
        refined.verified = true;
        refined.depth = box.depth;
        refined.condition_estimate = condition_estimate;
        refined.needs_higher_precision = needs_higher_precision;

        verified_indices.push_back(i);
        candidate_roots.push_back(refined);
    }

    // Step 2: Merge nearby roots and cancel duplicates
    std::set<std::size_t> cancelled;
    std::vector<bool> merged(candidate_roots.size(), false);

    for (std::size_t i = 0; i < candidate_roots.size(); ++i) {
        if (merged[i]) continue;

        RefinedRoot& root = candidate_roots[i];
        double radius = root.exclusion_radius;  // Use pre-computed exclusion radius

        // Check for nearby candidate roots to merge
        for (std::size_t j = i + 1; j < candidate_roots.size(); ++j) {
            if (merged[j]) continue;

            double dist = computeDistance(root.location, candidate_roots[j].location);
            if (dist < radius) {
                // Merge j into i
                root.source_boxes.insert(
                    root.source_boxes.end(),
                    candidate_roots[j].source_boxes.begin(),
                    candidate_roots[j].source_boxes.end());
                merged[j] = true;
                cancelled.insert(verified_indices[j]);
            }
        }
        
        // Check all other boxes (including unresolved) for cancellation
        for (std::size_t j = 0; j < solver_result.boxes.size(); ++j) {
            if (j < num_resolved && std::find(verified_indices.begin(), 
                                              verified_indices.end(), j) != verified_indices.end()) {
                continue;  // Skip verified roots (handled above)
            }
            
            const SubdivisionBoxResult& other_box = solver_result.boxes[j];
            double dist = computeDistance(root.location, other_box.center);
            
            if (dist < radius) {
                cancelled.insert(j);
            }
        }
        
        // Add to final results
        result.roots.push_back(root);
    }
    
    // Step 3: Collect cancelled box indices
    result.cancelled_boxes.assign(cancelled.begin(), cancelled.end());
    std::sort(result.cancelled_boxes.begin(), result.cancelled_boxes.end());

    // Step 4: Add unresolved boxes from solver to unverified list
    // These are boxes that the solver couldn't resolve (e.g., due to degeneracy)
    for (std::size_t i = num_resolved; i < solver_result.boxes.size(); ++i) {
        result.unverified_boxes.push_back(i);
    }

    // Step 5: Merge unverified boxes into problematic regions (1D only)
    if (!result.unverified_boxes.empty()) {
        result.problematic_regions = mergeUnverifiedBoxes1D(
            result.unverified_boxes, solver_result, poly, config);
    }

    return result;
}

bool ResultRefiner::refineRoot1D(
    const SubdivisionBoxResult& box,
    const Polynomial& poly,
    const RefinementConfig& config,
    double& refined_location,
    double& residual) const
{
    // Get derivative
    Polynomial dpoly = Differentiation::derivative(poly, 0, 1);

    // Start from box center
    double x = box.center[0];
    double lower = box.lower[0];
    double upper = box.upper[0];

    // Newton iteration
    for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
        double f = poly.evaluate(x);
        double df = dpoly.evaluate(x);

        // Check for convergence with condition-aware criterion
        if (std::abs(f) < config.residual_tolerance) {
            // Residual is small - but is the error also small?
            // For ill-conditioned problems, small residual doesn't guarantee small error
            // Check condition number to estimate actual error

            // Estimate actual error using Newton's method error formula
            // For a simple root: |x - r| ≈ |f(x)| / |f'(r)| ≈ |f(x)| / |f'(x)|
            //
            // NOTE: This error bound is still not rigorous and can be overly conservative.
            // Issues:
            // 1. The formula assumes f'(x) ≈ f'(r), which may not hold far from the root
            // 2. Does not account for rounding errors in polynomial evaluation
            // 3. Can be 10-1000x more conservative than actual error
            //
            // TODO: Consider replacing with interval Newton method for rigorous error bounds
            // that properly account for:
            // - Rounding errors in coefficient representation
            // - Rounding errors in polynomial evaluation
            // - Distance from current iterate to true root
            //
            // For now, this provides a reasonable heuristic for convergence checking.
            double estimated_error = std::abs(f) / std::max(std::abs(df), 1e-14);

            // Accept root only if estimated error is within target tolerance
            if (estimated_error <= config.target_tolerance) {
                refined_location = x;
                residual = f;
                return true;
            }

            // Residual is small but estimated error is too large
            // This indicates an ill-conditioned problem requiring higher precision
            // Continue iterating (though unlikely to improve with double precision)
            // Will eventually hit max iterations and return false
        }

        // Check for zero derivative (multiplicity or failure)
        if (std::abs(df) < 1e-14) {
            // Try to subdivide the box and check sign changes
            double mid = 0.5 * (lower + upper);
            double f_lower = poly.evaluate(lower);
            double f_mid = poly.evaluate(mid);
            double f_upper = poly.evaluate(upper);

            // Check which half has sign change
            if (f_lower * f_mid < 0.0) {
                // Root in [lower, mid]
                upper = mid;
                x = 0.5 * (lower + upper);
            } else if (f_mid * f_upper < 0.0) {
                // Root in [mid, upper]
                lower = mid;
                x = 0.5 * (lower + upper);
            } else {
                // No clear sign change, might be multiple root
                // Use current best estimate
                refined_location = x;
                residual = f;

                // For multiple roots, derivative is near zero, so condition check
                // would give unreliable results. Just check residual.
                // The multiplicity detection later will handle this case.
                return std::abs(f) < config.residual_tolerance;
            }
            continue;
        }

        // Newton step
        double x_new = x - f / df;

        // Check if step is valid
        if (!isValidNewtonStep(x, x_new, lower, upper, f, poly.evaluate(x_new))) {
            // Step went outside box or increased residual
            // Subdivide and check sign changes
            double mid = 0.5 * (lower + upper);
            double f_lower = poly.evaluate(lower);
            double f_mid = poly.evaluate(mid);
            double f_upper = poly.evaluate(upper);

            if (f_lower * f_mid < 0.0) {
                upper = mid;
                x = 0.5 * (lower + upper);
            } else if (f_mid * f_upper < 0.0) {
                lower = mid;
                x = 0.5 * (lower + upper);
            } else {
                // No sign change, use bisection toward current best
                if (x < mid) {
                    upper = mid;
                } else {
                    lower = mid;
                }
                x = 0.5 * (lower + upper);
            }
        } else {
            // Accept Newton step
            x = x_new;

            // Update box bounds to maintain sign change
            double f_lower = poly.evaluate(lower);
            double f_upper = poly.evaluate(upper);
            double f_x = poly.evaluate(x);

            if (f_lower * f_x < 0.0) {
                upper = x;
            } else if (f_x * f_upper < 0.0) {
                lower = x;
            }
            // If no sign change, keep bounds as is
        }
    }

    // Max iterations reached, check if we're close enough with condition-aware criterion
    double f = poly.evaluate(x);
    refined_location = x;
    residual = f;

    if (std::abs(f) < config.residual_tolerance) {
        // Check condition number to verify error is acceptable
        double df = dpoly.evaluate(x);
        Polynomial ddpoly = Differentiation::derivative(dpoly, 0, 1);
        double ddf = ddpoly.evaluate(x);

        double kappa = std::abs(ddf) / (std::abs(df) * std::abs(df) + 1e-100);
        double estimated_error = kappa * std::abs(f) / std::max(std::abs(df), 1e-14);

        return estimated_error <= config.target_tolerance;
    }

    return false;
}

bool ResultRefiner::isValidNewtonStep(
    double x_old, double x_new,
    double lower, double upper,
    double f_old, double f_new) const
{
    // Check if new position is within bounds
    if (x_new < lower || x_new > upper) {
        return false;
    }

    // Check if residual decreased
    if (std::abs(f_new) > std::abs(f_old)) {
        return false;
    }

    return true;
}

bool ResultRefiner::refineRoot1DWithMultiplicity(
    double initial_guess,
    double lower,
    double upper,
    const Polynomial& poly,
    unsigned int multiplicity,
    const RefinementConfig& config,
    double& refined_location,
    double& residual) const
{
    // Get derivative
    Polynomial dpoly = Differentiation::derivative(poly, 0, 1);

    double x = initial_guess;

    // Modified Newton iteration: x_{n+1} = x_n - m * f(x_n) / f'(x_n)
    for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
        double f = poly.evaluate(x);
        double df = dpoly.evaluate(x);

        // Check for convergence
        if (std::abs(f) < config.residual_tolerance) {
            refined_location = x;
            residual = f;
            return true;
        }

        // Check for zero derivative
        if (std::abs(df) < 1e-14) {
            // For multiple roots, try using higher-order derivatives
            // Compute f^(m)(x) where m is the multiplicity
            Polynomial deriv_m = poly;
            for (unsigned int k = 0; k < multiplicity; ++k) {
                deriv_m = Differentiation::derivative(deriv_m, 0, 1);
            }
            double f_m = deriv_m.evaluate(x);

            // If f^(m)(x) is non-zero, we can estimate the step size
            if (std::abs(f_m) > 1e-14) {
                // For a root with multiplicity m: f(x) ≈ c * (x - r)^m
                // So: (x - r) ≈ (f(x) * m! / f^(m)(x))^(1/m)
                double factorial_m = 1.0;
                for (unsigned int k = 2; k <= multiplicity; ++k) {
                    factorial_m *= k;
                }

                double ratio = f * factorial_m / f_m;
                double step = std::copysign(std::pow(std::abs(ratio), 1.0 / multiplicity), ratio);

                double x_new = x - step;

                // Clamp to bounds
                if (x_new < lower) x_new = 0.5 * (x + lower);
                if (x_new > upper) x_new = 0.5 * (x + upper);

                x = x_new;
                continue;
            }

            // If we can't make progress, fail
            return false;
        }

        // Modified Newton step
        double x_new = x - multiplicity * f / df;

        // Clamp to bounds
        if (x_new < lower) {
            x_new = 0.5 * (x + lower);
        }
        if (x_new > upper) {
            x_new = 0.5 * (x + upper);
        }

        x = x_new;
    }

    // Check final residual
    double f = poly.evaluate(x);
    if (std::abs(f) < config.residual_tolerance) {
        refined_location = x;
        residual = f;
        return true;
    }

    return false;
}

double ResultRefiner::estimateConditionNumber1D(
    double location,
    const Polynomial& poly,
    double derivative_value) const
{
    // For a root-finding problem, the condition number relates the error in
    // the root to the residual. A practical estimate uses the ratio of
    // second derivative to first derivative squared.
    //
    // Theory: For f(x) near a root r, we have:
    //   f(x) ≈ f'(r)*(x-r) + f''(r)*(x-r)^2/2
    //
    // If we have residual ε = f(x), then:
    //   |x-r| ≈ |ε/f'(r)| * (1 + |f''(r)|*|x-r|/(2|f'(r)|))
    //
    // The condition number is approximately:
    //   κ ≈ 1 + |f''(r)|/(2|f'(r)|) * typical_error
    //
    // For a more robust estimate, we use:
    //   κ ≈ max(1, |f''(x)| / |f'(x)|^2)

    if (std::abs(derivative_value) < 1e-14) {
        // Derivative too small, condition number is very large
        return 1e16;
    }

    // Compute second derivative
    Polynomial dpoly = Differentiation::derivative(poly, 0, 1);
    Polynomial ddpoly = Differentiation::derivative(dpoly, 0, 1);
    double f_prime = derivative_value;
    double f_double_prime = ddpoly.evaluate(location);

    // Condition number estimate: |f''| / |f'|^2
    // This gives the sensitivity of the root to perturbations
    double kappa = std::abs(f_double_prime) / (std::abs(f_prime) * std::abs(f_prime));

    // For very ill-conditioned problems, also check higher derivatives
    // If f'''(x) is large relative to f'(x), the problem is even worse
    Polynomial dddpoly = Differentiation::derivative(ddpoly, 0, 1);
    double f_triple_prime = dddpoly.evaluate(location);
    double kappa_higher = std::abs(f_triple_prime) / (std::abs(f_prime) * std::abs(f_prime) * std::abs(f_prime));

    // Take the maximum as a conservative estimate
    kappa = std::max(kappa, kappa_higher);

    return std::max(1.0, kappa);
}

unsigned int ResultRefiner::estimateMultiplicityOstrowski(
    double x1, double x2, double x3) const
{
    // Ostrowski's method (1973) for multiplicity estimation
    // For a root of multiplicity m, Newton iteration satisfies:
    //   e_{n+1} ≈ ((m-1)/m) * e_n
    // where e_n = x_n - r is the error at iteration n
    //
    // From 3 consecutive iterates, we can estimate:
    //   p = 1/2 + (x₁ - x₂) / (x₃ - 2x₂ + x₁)
    // Then m = floor(p)

    double numerator = x1 - x2;
    double denominator = x3 - 2.0 * x2 + x1;

    // Avoid division by zero
    if (std::abs(denominator) < 1e-15) {
        return 1;  // Can't estimate, assume simple root
    }

    double p_est = 0.5 + numerator / denominator;

    // Apply floor with minimum value 1
    int multiplicity = static_cast<int>(std::floor(p_est));
    if (multiplicity < 1) {
        multiplicity = 1;
    }

    // Sanity check: cap at reasonable maximum
    if (multiplicity > 20) {
        multiplicity = 20;
    }

    return static_cast<unsigned int>(multiplicity);
}

unsigned int ResultRefiner::estimateMultiplicityOstrowskiFromPoint(
    double x0, const Polynomial& poly) const
{
    // Perform 3 Newton iterations to get x1, x2, x3
    Polynomial dpoly = Differentiation::derivative(poly, 0, 1);

    double x = x0;
    std::vector<double> iterates;
    iterates.push_back(x);

    for (int i = 0; i < 3; ++i) {
        double f = poly.evaluate(x);
        double df = dpoly.evaluate(x);

        if (std::abs(df) < 1e-15) {
            // Derivative too small, can't continue Newton
            return 1;
        }

        x = x - f / df;
        iterates.push_back(x);
    }

    // Now we have x0, x1, x2, x3
    return estimateMultiplicityOstrowski(iterates[1], iterates[2], iterates[3]);
}

unsigned int ResultRefiner::estimateMultiplicity(
    const std::vector<double>& point,
    const PolynomialSystem& system,
    unsigned int max_order,
    double derivative_threshold,
    double& first_nonzero_deriv) const
{
    const std::size_t dim = system.dimension();
    first_nonzero_deriv = 0.0;

    // For 1D case, check derivatives directly
    if (dim == 1) {
        // Check each equation (though typically only one for 1D)
        for (const Polynomial& eq : system.equations()) {
            // Check derivatives from order 1 to max_order
            for (unsigned int order = 1; order <= max_order; ++order) {
                Polynomial deriv = Differentiation::derivative(eq, 0, order);
                double deriv_val = deriv.evaluate(point);

                if (std::abs(deriv_val) > derivative_threshold) {
                    // Found first non-zero derivative at order 'order'
                    // This means multiplicity = order
                    first_nonzero_deriv = deriv_val;
                    return order;
                }
            }
        }

        // All derivatives zero up to max_order
        // Multiplicity is at least max_order + 1
        first_nonzero_deriv = 0.0;
        return max_order + 1;
    }

    // For multi-dimensional case, check gradient and higher-order derivatives
    for (unsigned int order = 1; order <= max_order; ++order) {
        bool has_nonzero = false;
        double max_deriv = 0.0;

        // Check all equations
        for (const Polynomial& eq : system.equations()) {
            // For first order, check gradient (all partial derivatives)
            if (order == 1) {
                std::vector<Polynomial> grad = Differentiation::gradient(eq);
                for (std::size_t axis = 0; axis < dim; ++axis) {
                    double deriv_val = grad[axis].evaluate(point);
                    if (std::abs(deriv_val) > derivative_threshold) {
                        has_nonzero = true;
                        if (std::abs(deriv_val) > std::abs(max_deriv)) {
                            max_deriv = deriv_val;
                        }
                    }
                }
            } else {
                // For higher orders, check partial derivatives along each axis
                // Note: For true multiplicity, we should check all mixed partials,
                // but for simplicity we check diagonal derivatives
                for (std::size_t axis = 0; axis < dim; ++axis) {
                    Polynomial deriv = Differentiation::derivative(eq, axis, order);
                    double deriv_val = deriv.evaluate(point);
                    if (std::abs(deriv_val) > derivative_threshold) {
                        has_nonzero = true;
                        if (std::abs(deriv_val) > std::abs(max_deriv)) {
                            max_deriv = deriv_val;
                        }
                    }
                }
            }
        }

        if (has_nonzero) {
            // First non-zero derivative at order 'order'
            first_nonzero_deriv = max_deriv;
            return order;
        }
    }

    // All derivatives zero up to max_order
    // Multiplicity is at least max_order + 1
    first_nonzero_deriv = 0.0;
    return max_order + 1;
}

double ResultRefiner::computeExclusionRadiusFromDerivative(
    unsigned int multiplicity,
    double first_nonzero_deriv,
    double tolerance,
    double multiplier) const
{
    // Avoid division by zero
    double abs_deriv = std::abs(first_nonzero_deriv);
    if (abs_deriv < 1e-100) {
        // Derivative is essentially zero, use large exclusion radius
        return multiplier * 1e-3;  // 0.1% of domain
    }

    if (multiplicity == 1) {
        // Simple root: r ≈ tolerance / |f'(x)|
        // The root is approximately at distance r where |f(x+r)| ≈ |f'(x)| * r = tolerance
        double radius = multiplier * tolerance / abs_deriv;

        // Clamp to reasonable range
        radius = std::max(radius, multiplier * tolerance);  // At least tolerance-based
        radius = std::min(radius, 0.1);  // At most 10% of domain

        return radius;
    } else {
        // Multiple root: r ≈ (tolerance / |f^(m)(x)|)^(1/m)
        // For a root of multiplicity m: f(x+r) ≈ f^(m)(x) * r^m / m!
        // Setting |f(x+r)| = tolerance: r ≈ (tolerance * m! / |f^(m)(x)|)^(1/m)

        // Compute factorial (for small m)
        double factorial = 1.0;
        for (unsigned int i = 2; i <= multiplicity && i <= 10; ++i) {
            factorial *= static_cast<double>(i);
        }

        double exponent = 1.0 / static_cast<double>(multiplicity);
        double radius = multiplier * std::pow(tolerance * factorial / abs_deriv, exponent);

        // Clamp to reasonable range
        radius = std::max(radius, multiplier * std::pow(tolerance, exponent));
        radius = std::min(radius, 0.1);

        return radius;
    }
}

double ResultRefiner::computeExclusionRadius(
    unsigned int multiplicity,
    double tolerance,
    double multiplier) const
{
    if (multiplicity == 1) {
        // Simple root: linear exclusion
        return multiplier * tolerance;
    } else {
        // Multiple root: scale by tolerance^(1/m)
        double exponent = 1.0 / static_cast<double>(multiplicity);
        return multiplier * std::pow(tolerance, exponent);
    }
}

double ResultRefiner::computeDistance(
    const std::vector<double>& p1,
    const std::vector<double>& p2) const
{
    if (p1.size() != p2.size()) {
        return std::numeric_limits<double>::infinity();
    }

    double sum_sq = 0.0;
    for (std::size_t i = 0; i < p1.size(); ++i) {
        double diff = p1[i] - p2[i];
        sum_sq += diff * diff;
    }

    return std::sqrt(sum_sq);
}

std::vector<ProblematicRegion> ResultRefiner::mergeUnverifiedBoxes1D(
    const std::vector<std::size_t>& unverified_indices,
    const SubdivisionSolverResult& solver_result,
    const Polynomial& poly,
    const RefinementConfig& config) const
{
    std::vector<ProblematicRegion> regions;

    if (unverified_indices.empty()) {
        return regions;
    }

    // Sort by center position
    std::vector<std::size_t> sorted_indices = unverified_indices;
    std::sort(sorted_indices.begin(), sorted_indices.end(),
        [&solver_result](std::size_t a, std::size_t b) {
            return solver_result.boxes[a].center[0] < solver_result.boxes[b].center[0];
        });

    // Merge nearby boxes into regions
    // Use a threshold based on box sizes - adjacent boxes should be merged
    double merge_threshold = 1e-6;  // Boxes closer than this are merged

    ProblematicRegion current;
    current.lower = solver_result.boxes[sorted_indices[0]].lower[0];
    current.upper = solver_result.boxes[sorted_indices[0]].upper[0];
    current.box_indices.push_back(sorted_indices[0]);

    for (std::size_t i = 1; i < sorted_indices.size(); ++i) {
        std::size_t idx = sorted_indices[i];
        const auto& box = solver_result.boxes[idx];

        // Check if this box is adjacent or close to current region
        if (box.lower[0] - current.upper <= merge_threshold) {
            // Merge into current region
            current.upper = std::max(current.upper, box.upper[0]);
            current.box_indices.push_back(idx);
        } else {
            // Finish current region and start new one
            regions.push_back(current);

            current = ProblematicRegion();
            current.lower = box.lower[0];
            current.upper = box.upper[0];
            current.box_indices.push_back(idx);
        }
    }
    regions.push_back(current);

    // Try to refine each region
    for (auto& region : regions) {
        refineProblematicRegion1D(region, poly, config);
    }

    return regions;
}

void ResultRefiner::refineProblematicRegion1D(
    ProblematicRegion& region,
    const Polynomial& poly,
    const RefinementConfig& config) const
{
    region.refinement_attempted = true;

    // Start from center of region
    double x0 = (region.lower + region.upper) / 2.0;

    // Estimate multiplicity at starting point
    PolynomialSystem system(std::vector<Polynomial>{poly});
    std::vector<double> point{x0};
    double first_nonzero_deriv = 0.0;
    unsigned int mult = estimateMultiplicity(
        point, system, config.max_multiplicity, 1e-10, first_nonzero_deriv);

    region.multiplicity = mult;

    // Try Newton refinement
    double refined_x;
    double residual;

    bool converged = refineRoot1D_fromPoint(x0, poly, config, refined_x, residual);

    if (!converged) {
        region.refinement_succeeded = false;
        region.refined_root = x0;  // Use starting point
        region.residual = poly.evaluate(x0);
        region.condition_estimate = estimateConditionNumber1D(x0, poly, first_nonzero_deriv);
        region.needs_higher_precision = true;
        return;
    }

    // Update multiplicity estimate at refined location
    point[0] = refined_x;
    mult = estimateMultiplicity(
        point, system, config.max_multiplicity, 1e-10, first_nonzero_deriv);

    double condition = estimateConditionNumber1D(refined_x, poly, first_nonzero_deriv);

    // Estimate actual error using condition number
    double est_error = condition * std::abs(residual) / std::max(std::abs(first_nonzero_deriv), 1e-14);

    // Check condition-aware convergence
    region.refinement_succeeded = (est_error <= config.target_tolerance);
    region.refined_root = refined_x;
    region.multiplicity = mult;
    region.residual = residual;
    region.condition_estimate = condition;
    region.needs_higher_precision = (est_error > 1e-10);
}

bool ResultRefiner::refineRoot1D_fromPoint(
    double x0,
    const Polynomial& poly,
    const RefinementConfig& config,
    double& refined_location,
    double& residual) const
{
    // New workflow based on detection limits analysis:
    // 1. Do 3 standard Newton iterations to get close to root
    // 2. Use Ostrowski to estimate multiplicity from the 3 iterates
    // 3. Use modified Newton with detected multiplicity to converge

    // Use power-basis differentiation if polynomial has power coefficients
    // This is more efficient (simpler formula) and avoids conversion
    Polynomial dpoly;
    if (poly.hasPowerCoefficients()) {
        dpoly = Differentiation::differentiateAxisPower(poly, 0);
    } else {
        dpoly = Differentiation::derivative(poly, 0, 1);
    }

    double x = x0;
    std::vector<double> iterates;
    iterates.push_back(x);

    unsigned int estimated_mult = 1;
    bool multiplicity_detected = false;

    for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
        double f = poly.evaluate(x);
        double df = dpoly.evaluate(x);

        // After 3 iterations, estimate multiplicity using Ostrowski
        if (iter == 2 && iterates.size() == 4 && !multiplicity_detected) {
            estimated_mult = estimateMultiplicityOstrowski(
                iterates[1], iterates[2], iterates[3]);
            multiplicity_detected = true;
        }

        // Check if derivative is too small (multiple root)
        if (std::abs(df) < 1e-14) {
            // For multiple roots, check convergence based on residual only
            if (std::abs(f) < config.residual_tolerance) {
                refined_location = x;
                residual = f;
                return true;
            }

            // If we haven't detected multiplicity yet, try Ostrowski from current point
            if (!multiplicity_detected) {
                estimated_mult = estimateMultiplicityOstrowskiFromPoint(x, poly);
                multiplicity_detected = true;
            }

            // Use modified Newton: x_new = x - m * f / f'
            if (estimated_mult > 1 && std::abs(df) > 1e-15) {
                double step = static_cast<double>(estimated_mult) * f / df;
                x = x - step;
                continue;
            }

            // Can't make progress
            return false;
        }

        // Check for convergence with condition-aware criterion
        if (std::abs(f) < config.residual_tolerance) {
            // Compute second derivative for condition estimation
            Polynomial ddpoly = Differentiation::derivative(dpoly, 0, 1);
            double ddf = ddpoly.evaluate(x);

            // Estimate condition number: κ ≈ |f''| / |f'|²
            double kappa = std::abs(ddf) / (std::abs(df) * std::abs(df) + 1e-100);

            // Estimate actual error: |error| ≈ κ × |residual| / |f'|
            double estimated_error = kappa * std::abs(f) / std::max(std::abs(df), 1e-14);

            // Accept root only if estimated error is within target tolerance
            if (estimated_error <= config.target_tolerance) {
                refined_location = x;
                residual = f;
                return true;
            }
        }

        // Choose Newton step based on detected multiplicity
        double step;
        if (multiplicity_detected && estimated_mult > 1) {
            // Modified Newton for multiple roots
            step = static_cast<double>(estimated_mult) * f / df;
        } else {
            // Standard Newton for simple roots
            step = f / df;
        }

        double x_new = x - step;
        x = x_new;

        // Store iterate for Ostrowski (only first 3)
        if (iter < 3) {
            iterates.push_back(x);
        }
    }

    // Max iterations reached - check final residual
    double f = poly.evaluate(x);
    residual = f;

    if (std::abs(f) < config.residual_tolerance) {
        // Check condition-aware convergence one more time
        double df = dpoly.evaluate(x);
        Polynomial ddpoly = Differentiation::derivative(dpoly, 0, 1);
        double ddf = ddpoly.evaluate(x);

        double kappa = std::abs(ddf) / (std::abs(df) * std::abs(df) + 1e-100);
        double estimated_error = kappa * std::abs(f) / std::max(std::abs(df), 1e-14);

        if (estimated_error <= config.target_tolerance) {
            refined_location = x;
            return true;
        }
    }

    return false;
}

#ifdef ENABLE_HIGH_PRECISION
bool ResultRefiner::refineRoot1DWithPrecisionEscalation(
    double initial_guess,
    const Polynomial& poly,
    const RefinementConfig& config,
    RefinedRoot& refined_root) const
{
    // Step 1: Try double precision refinement first
    double refined_x;
    double residual;
    bool converged = refineRoot1D_fromPoint(initial_guess, poly, config, refined_x, residual);

    if (!converged) {
        // Refinement failed in double precision
        refined_root.verified = false;
        refined_root.needs_higher_precision = true;
        refined_root.location = std::vector<double>(1, initial_guess);
        return false;
    }

    // Step 2: Estimate multiplicity and condition number
    PolynomialSystem system(std::vector<Polynomial>{poly});
    std::vector<double> point{refined_x};
    double first_nonzero_deriv = 0.0;
    unsigned int mult = estimateMultiplicity(
        point, system, config.max_multiplicity, 1e-10, first_nonzero_deriv);

    double condition = estimateConditionNumber1D(refined_x, poly, first_nonzero_deriv);
    double est_error = condition * std::abs(residual) / std::max(std::abs(first_nonzero_deriv), 1e-14);

    // Step 3: Check if double precision is sufficient
    if (est_error <= config.target_tolerance && condition < 1e5) {
        // Double precision is sufficient
        refined_root.location = std::vector<double>(1, refined_x);
        refined_root.residual = std::vector<double>(1, residual);
        refined_root.multiplicity = mult;
        refined_root.first_nonzero_derivative = first_nonzero_deriv;
        refined_root.condition_estimate = condition;
        refined_root.verified = true;
        refined_root.needs_higher_precision = false;
        return true;
    }

    // Step 4: Switch to high precision
    // Select precision based on condition number
    unsigned int precision_bits;
    if (condition < 1e5) {
        precision_bits = 256;  // ~77 decimal digits
    } else if (condition < 1e10) {
        precision_bits = 512;  // ~154 decimal digits
    } else {
        precision_bits = 1024; // ~308 decimal digits
    }

    std::cout << "Switching to " << precision_bits << "-bit precision (condition="
              << std::scientific << std::setprecision(2) << condition
              << ", est_error=" << est_error << ")" << std::endl;

    // Set precision context
    PrecisionContext ctx(precision_bits);

    // Convert polynomial to high precision
    PolynomialHP poly_hp(poly);

    // Configure high-precision refinement
    RefinementConfigHP config_hp;
    config_hp.max_newton_iters = config.max_newton_iters * 2; // Allow more iterations
    config_hp.max_multiplicity = config.max_multiplicity;

    // Only pass multiplicity hint if condition number is reasonable AND multiplicity is low
    // For very ill-conditioned problems or high multiplicities (>3), double-precision
    // multiplicity detection may be unreliable, so let HP refiner detect it from scratch
    if (condition < 1e8 && mult <= 3) {
        config_hp.multiplicity_hint = mult; // Pass the multiplicity detected in double precision
    } else {
        config_hp.multiplicity_hint = 0; // Let HP refiner detect multiplicity from scratch
    }

    // Set tolerance based on precision
    // Use more conservative tolerances that are achievable
    // For multiple roots, convergence is slower, so we need realistic targets
    if (precision_bits >= 1024) {
        config_hp.target_tolerance_str = "1e-100";  // ~1/3 of available precision
        config_hp.residual_tolerance_str = "1e-100";
    } else if (precision_bits >= 512) {
        config_hp.target_tolerance_str = "1e-50";   // ~1/3 of available precision
        config_hp.residual_tolerance_str = "1e-50";
    } else {
        config_hp.target_tolerance_str = "1e-25";   // ~1/3 of available precision
        config_hp.residual_tolerance_str = "1e-25";
    }

    // Refine with high precision
    RefinedRootHP result_hp = ResultRefinerHP::refineRoot1D(refined_x, poly_hp, config_hp);



    if (!result_hp.converged) {
        // HP refiner didn't fully converge, but may have improved the result
        // For multiple roots, full convergence is difficult, so accept partial results
        // if we have error bounds
        if (!result_hp.has_guaranteed_bounds) {
            // No bounds available - refinement truly failed
            refined_root.verified = false;
            refined_root.needs_higher_precision = true;
            refined_root.location = std::vector<double>(1, refined_x);
            refined_root.condition_estimate = condition;
            return false;
        }

        // We have bounds but didn't meet tolerance - accept as best effort
        std::cout << "  Note: HP refiner achieved partial convergence (error bounds available)" << std::endl;
    }

    // Step 5: Convert high-precision result back to double precision structure
    refined_root.location = std::vector<double>(1, toDouble(result_hp.location));
    refined_root.residual = std::vector<double>(1, toDouble(result_hp.residual));
    refined_root.multiplicity = result_hp.multiplicity;
    refined_root.first_nonzero_derivative = toDouble(result_hp.first_nonzero_derivative);
    refined_root.condition_estimate = toDouble(result_hp.condition_estimate);
    refined_root.verified = true;
    refined_root.needs_higher_precision = false; // Successfully refined with HP

    // Store error bounds if available
    if (result_hp.has_guaranteed_bounds) {
        refined_root.max_error = std::vector<double>(1, toDouble(result_hp.max_error));
    }

    std::cout << "High-precision refinement succeeded: root=" << refined_root.location[0]
              << ", multiplicity=" << refined_root.multiplicity
              << ", error=" << (result_hp.has_guaranteed_bounds ? toDouble(result_hp.max_error) : 0.0)
              << std::endl;

    return true;
}
#endif

//=============================================================================
// CurveRefiner Implementation
//=============================================================================

CurveRefiner::CurveRefiner(const Polynomial& curve_poly)
    : curve_poly_(curve_poly)
{
    if (curve_poly.dimension() != 2) {
        throw std::invalid_argument(
            "CurveRefiner requires a 2-dimensional polynomial, got dimension " +
            std::to_string(curve_poly.dimension()));
    }

    // Precompute partial derivatives
    dg_dx_ = Differentiation::derivative(curve_poly_, 0, 1);
    dg_dy_ = Differentiation::derivative(curve_poly_, 1, 1);
}

CurveRefinedPoint CurveRefiner::refine(
    double x0, double y0,
    const CurveRefinementConfig& config) const
{
    CurveRefinedPoint result;
    result.x = x0;
    result.y = y0;
    result.converged = false;
    result.iterations = 0;

    double x = x0, y = y0;

    for (unsigned int iter = 0; iter < config.max_iterations; ++iter) {
        result.iterations = iter + 1;

        // Evaluate g and its gradient at current point
        double g = curve_poly_.evaluate({x, y});
        double gx = dg_dx_.evaluate({x, y});
        double gy = dg_dy_.evaluate({x, y});

        // Check convergence
        result.residual = std::abs(g);
        if (result.residual < config.residual_tolerance) {
            result.x = x;
            result.y = y;
            result.converged = true;
            return result;
        }

        // Gradient magnitude squared
        double grad_sq = gx * gx + gy * gy;

        if (grad_sq < config.min_gradient_norm) {
            // At or near a singular point (gradient too small)
            result.x = x;
            result.y = y;
            return result;
        }

        // Newton step along gradient: move by -g/|∇g|² × ∇g
        double t = -g / grad_sq;
        x += t * gx;
        y += t * gy;
    }

    // Did not converge within max iterations
    result.x = x;
    result.y = y;
    result.residual = std::abs(curve_poly_.evaluate({x, y}));
    return result;
}

std::vector<CurveRefinedPoint> CurveRefiner::refineMultiple(
    const std::vector<std::pair<double, double>>& points,
    const CurveRefinementConfig& config) const
{
    std::vector<CurveRefinedPoint> results;
    results.reserve(points.size());

    for (const auto& pt : points) {
        results.push_back(refine(pt.first, pt.second, config));
    }

    return results;
}

//=============================================================================
// Function-based curve refinement (numerical gradient)
//=============================================================================

CurveRefinedPoint refineCurveNumerical(
    const std::function<double(double, double)>& g,
    double x0, double y0,
    const CurveRefinementConfig& config)
{
    CurveRefinedPoint result;
    result.x = x0;
    result.y = y0;
    result.converged = false;
    result.iterations = 0;

    double x = x0, y = y0;
    const double h = 1e-8;  // Step size for numerical gradient

    for (unsigned int iter = 0; iter < config.max_iterations; ++iter) {
        result.iterations = iter + 1;

        // Evaluate function at current point
        double val = g(x, y);
        result.residual = std::abs(val);

        // Check convergence
        if (result.residual < config.residual_tolerance) {
            result.x = x;
            result.y = y;
            result.converged = true;
            return result;
        }

        // Compute numerical gradient using central differences
        double gx = (g(x + h, y) - g(x - h, y)) / (2 * h);
        double gy = (g(x, y + h) - g(x, y - h)) / (2 * h);

        // Gradient squared norm
        double grad_sq = gx * gx + gy * gy;

        // Check for singular point (gradient too small)
        if (grad_sq < config.min_gradient_norm * config.min_gradient_norm) {
            result.x = x;
            result.y = y;
            return result;  // converged = false
        }

        // Newton step along gradient direction
        // Project point onto zero set: (x,y) -= (g / |∇g|²) * ∇g
        double t = -val / grad_sq;
        x += t * gx;
        y += t * gy;
    }

    // Did not converge within max iterations
    result.x = x;
    result.y = y;
    result.residual = std::abs(g(x, y));
    return result;
}

std::vector<CurveRefinedPoint> refineCurveNumericalMultiple(
    const std::function<double(double, double)>& g,
    const std::vector<std::pair<double, double>>& points,
    const CurveRefinementConfig& config)
{
    std::vector<CurveRefinedPoint> results;
    results.reserve(points.size());

    for (const auto& pt : points) {
        results.push_back(refineCurveNumerical(g, pt.first, pt.second, config));
    }

    return results;
}

} // namespace polynomial_solver

