#include "result_refiner.h"
#include "differentiation.h"
#include <cmath>
#include <algorithm>
#include <set>
#include <limits>

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

    if (num_resolved == 0) {
        return result;
    }

    // Only support 1D for now
    if (dim != 1) {
        // For 2D+, just return empty result
        // (2D may have infinite roots in degenerate regions)
        return result;
    }

    const Polynomial& poly = original_system.equation(0);

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
        refined.max_error = std::vector<double>{config.target_tolerance};
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

        // Check for convergence
        if (std::abs(f) < config.residual_tolerance) {
            refined_location = x;
            residual = f;
            return true;
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

    // Max iterations reached, check if we're close enough
    double f = poly.evaluate(x);
    refined_location = x;
    residual = f;
    return std::abs(f) < config.residual_tolerance;
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

} // namespace polynomial_solver

