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

        // Estimate multiplicity at refined location
        std::vector<double> point{refined_location};
        unsigned int mult = estimateMultiplicity(
            point, original_system, config.max_multiplicity);

        // Create refined root
        RefinedRoot refined;
        refined.location = point;
        refined.residual = std::vector<double>{residual};
        refined.max_error = std::vector<double>{config.target_tolerance};
        refined.multiplicity = mult;
        refined.source_boxes.push_back(i);
        refined.verified = true;
        refined.depth = box.depth;

        verified_indices.push_back(i);
        candidate_roots.push_back(refined);
    }
    
    // Step 2: Merge nearby roots and cancel duplicates
    std::set<std::size_t> cancelled;
    std::vector<bool> merged(candidate_roots.size(), false);
    
    for (std::size_t i = 0; i < candidate_roots.size(); ++i) {
        if (merged[i]) continue;
        
        RefinedRoot& root = candidate_roots[i];
        double radius = computeExclusionRadius(
            root.multiplicity,
            config.target_tolerance,
            config.exclusion_multiplier);
        
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

unsigned int ResultRefiner::estimateMultiplicity(
    const std::vector<double>& point,
    const PolynomialSystem& system,
    unsigned int max_order) const
{
    const std::size_t dim = system.dimension();
    const double derivative_threshold = 1e-10;

    // Check derivatives up to max_order
    for (unsigned int order = 1; order <= max_order; ++order) {
        bool has_nonzero = false;

        // Check all equations
        for (const Polynomial& eq : system.equations()) {
            // For first order, check gradient
            if (order == 1) {
                std::vector<Polynomial> grad = Differentiation::gradient(eq);
                for (std::size_t axis = 0; axis < dim; ++axis) {
                    double deriv_val = grad[axis].evaluate(point);
                    if (std::abs(deriv_val) > derivative_threshold) {
                        has_nonzero = true;
                        break;
                    }
                }
            } else {
                // For higher orders, check all partial derivatives
                // For simplicity, check diagonal second derivatives
                for (std::size_t axis = 0; axis < dim; ++axis) {
                    Polynomial deriv = Differentiation::derivative(eq, axis, order);
                    double deriv_val = deriv.evaluate(point);
                    if (std::abs(deriv_val) > derivative_threshold) {
                        has_nonzero = true;
                        break;
                    }
                }
            }

            if (has_nonzero) break;
        }

        if (has_nonzero) {
            return order;  // First non-zero derivative order
        }
    }

    // All derivatives zero up to max_order
    return max_order;
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

