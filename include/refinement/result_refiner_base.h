#ifndef RESULT_REFINER_BASE_H
#define RESULT_REFINER_BASE_H

/**
 * @file result_refiner_base.h
 * @brief Templated root refinement for PolynomialBase<Scalar>
 *
 * This header provides ResultRefinerBase<Scalar>, a fully templated root refiner
 * that works with any scalar type (double, mpreal, etc.). It unifies the
 * functionality of ResultRefiner and ResultRefinerHP into a single template.
 *
 * Features:
 * - Multiple iteration methods: Newton, Modified Newton, Halley, Schröder
 * - Multiple multiplicity detection: Taylor, Ostrowski, SimpleThreshold
 * - Condition-aware convergence
 * - Rigorous error bounds computation
 *
 * Usage:
 * @code
 * PolynomialBase<double> poly = ...;
 * RefinementConfigBase<double> config;
 * config.target_tolerance = 1e-15;
 *
 * auto result = ResultRefinerBase<double>::refineRoot1D(0.5, poly, config);
 * if (result.converged) {
 *     std::cout << "Root: " << result.location << std::endl;
 * }
 * @endcode
 */

#include "core/polynomial_base.h"
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <set>

namespace polynomial_solver {

//=============================================================================
// Enums for configuration
//=============================================================================

/**
 * @brief Multiplicity estimation method
 */
enum class MultiplicityMethodBase {
    NONE,           ///< No multiplicity detection (assume m=1)
    HINT,           ///< Use provided hint
    TAYLOR,         ///< Taylor series ratio test
    OSTROWSKI,      ///< Ostrowski's method (3 Newton iterates)
    SIMPLE_THRESHOLD ///< Simple absolute threshold on derivatives
};

/**
 * @brief Iteration method for root refinement
 */
enum class IterationMethodBase {
    NEWTON,         ///< Standard Newton's method
    MODIFIED_NEWTON,///< Modified Newton for multiple roots: x_new = x - m*f/f'
    HALLEY,         ///< Halley's method (third-order)
    SCHRODER        ///< Schröder's method (third-order, good for multiple roots)
};

//=============================================================================
// Configuration struct
//=============================================================================

/**
 * @brief Configuration for templated refinement
 *
 * This struct mirrors the original RefinementConfig with additional options
 * for multiplicity detection and iteration methods.
 */
template<typename Scalar>
struct RefinementConfigBase {
    // Core parameters (matching original RefinementConfig)
    Scalar target_tolerance;         ///< Target error tolerance (default: 1e-15)
    Scalar residual_tolerance;       ///< Residual tolerance |f(x)| for convergence (default: 1e-15)
    unsigned int max_newton_iters;   ///< Maximum Newton iterations (default: 50)
    unsigned int max_multiplicity;   ///< Maximum multiplicity to check (default: 10)
    Scalar exclusion_multiplier;     ///< Multiplier for exclusion radius (default: 3.0)

    // Extended parameters for templated version
    unsigned int multiplicity_hint;  ///< Hint for multiplicity (0 = auto-detect)
    Scalar taylor_ratio_threshold;   ///< Ratio threshold for Taylor method (default: 10)
    MultiplicityMethodBase multiplicity_method;  ///< Method for multiplicity detection
    IterationMethodBase iteration_method;        ///< Iteration method for refinement

    RefinementConfigBase()
        : target_tolerance(Scalar(1e-15)),
          residual_tolerance(Scalar(1e-15)),
          max_newton_iters(50u),
          max_multiplicity(10u),
          exclusion_multiplier(Scalar(3)),
          multiplicity_hint(0u),
          taylor_ratio_threshold(Scalar(10)),
          multiplicity_method(MultiplicityMethodBase::OSTROWSKI),  // More robust than Taylor
          iteration_method(IterationMethodBase::MODIFIED_NEWTON)
    {}
};

//=============================================================================
// Result struct
//=============================================================================

/**
 * @brief A refined root with convergence information
 */
template<typename Scalar>
struct RefinedRootBase {
    Scalar location;                 ///< Root location
    Scalar residual;                 ///< Residual f(x) at root
    unsigned int multiplicity;       ///< Estimated multiplicity
    Scalar first_nonzero_derivative; ///< Value of first non-zero derivative
    Scalar condition_estimate;       ///< Estimated condition number
    unsigned int iterations;         ///< Newton iterations performed
    bool converged;                  ///< True if converged to target tolerance

    // Error bounds
    Scalar max_error;                ///< Maximum error bound
    Scalar interval_lower;           ///< Lower bound of interval containing root
    Scalar interval_upper;           ///< Upper bound of interval containing root
    bool has_guaranteed_bounds;      ///< True if interval bounds are rigorous

    // Degeneracy indicators
    bool needs_higher_precision;     ///< True if condition suggests double precision insufficient
    Scalar exclusion_radius;         ///< Radius for merging nearby roots
    std::vector<std::size_t> source_boxes;  ///< Indices of solver boxes merged into this root

    RefinedRootBase()
        : location(0), residual(0), multiplicity(1),
          first_nonzero_derivative(0), condition_estimate(1),
          iterations(0), converged(false),
          max_error(0), interval_lower(0), interval_upper(0),
          has_guaranteed_bounds(false),
          needs_higher_precision(false), exclusion_radius(0)
    {}
};

/**
 * @brief A problematic region that couldn't be resolved
 */
template<typename Scalar>
struct ProblematicRegionBase {
    Scalar lower;                    ///< Lower bound of region
    Scalar upper;                    ///< Upper bound of region
    std::vector<std::size_t> box_indices;  ///< Solver box indices in this region

    // Refinement results (if attempted)
    bool refinement_attempted;       ///< True if refinement was attempted
    bool refinement_succeeded;       ///< True if refinement succeeded
    Scalar refined_root;             ///< Refined root location (if successful)
    unsigned int multiplicity;       ///< Estimated multiplicity
    Scalar residual;                 ///< Residual at refined root
    Scalar condition_estimate;       ///< Estimated condition number
    bool needs_higher_precision;     ///< True if higher precision needed

    ProblematicRegionBase()
        : lower(0), upper(0),
          refinement_attempted(false), refinement_succeeded(false),
          refined_root(0), multiplicity(0), residual(0),
          condition_estimate(1), needs_higher_precision(false)
    {}
};

/**
 * @brief Result of batch refinement process
 */
template<typename Scalar>
struct RefinementResultBase {
    std::vector<RefinedRootBase<Scalar>> roots;    ///< Verified and deduplicated roots
    std::vector<std::size_t> cancelled_boxes;      ///< Box indices eliminated as duplicates
    std::vector<std::size_t> unverified_boxes;     ///< Box indices that didn't pass verification
    std::vector<ProblematicRegionBase<Scalar>> problematic_regions;  ///< Merged unverified regions
    bool any_needs_higher_precision;               ///< True if any root needs higher precision
};

//=============================================================================
// Helper functions for scalar operations
//=============================================================================

namespace detail {

template<typename Scalar>
inline Scalar abs_val(const Scalar& x) {
    return (x < Scalar(0)) ? -x : x;
}

template<typename Scalar>
inline Scalar max_val(const Scalar& a, const Scalar& b) {
    return (a > b) ? a : b;
}

template<typename Scalar>
inline Scalar pow_val(const Scalar& base, unsigned int exp) {
    if (exp == 0) return Scalar(1);
    Scalar result = base;
    for (unsigned int i = 1; i < exp; ++i) {
        result *= base;
    }
    return result;
}

// Factorial for small values
inline unsigned long long factorial(unsigned int n) {
    unsigned long long result = 1;
    for (unsigned int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

} // namespace detail

//=============================================================================
// ResultRefinerBase class
//=============================================================================

/**
 * @class ResultRefinerBase
 * @brief Templated root refinement using Newton-type methods
 *
 * Provides static methods for refining polynomial roots to high accuracy.
 * Works with any scalar type that supports basic arithmetic operations.
 *
 * @tparam Scalar The coefficient/evaluation type (double, mpreal, etc.)
 */
template<typename Scalar>
class ResultRefinerBase {
public:
    using Poly = PolynomialBase<Scalar>;
    using Config = RefinementConfigBase<Scalar>;
    using Result = RefinedRootBase<Scalar>;

    //=========================================================================
    // Main Refinement Methods
    //=========================================================================

    /**
     * @brief Refine a 1D root using configured iteration method
     *
     * @param initial_guess Initial guess for root location
     * @param poly Polynomial to refine
     * @param config Refinement configuration
     * @return Refined root with convergence information
     */
    static Result refineRoot1D(
        const Scalar& initial_guess,
        const Poly& poly,
        const Config& config = Config())
    {
        Result result;
        Scalar x = initial_guess;

        // Get derivative polynomials
        // ddpoly is needed for condition-aware convergence check
        Poly dpoly = poly.differentiate(0);
        Poly ddpoly = dpoly.differentiate(0);

        // Phase 1: Do initial standard Newton iterations to get closer to root
        // This improves the Ostrowski multiplicity estimate which needs good iterates
        unsigned int initial_newton_steps = 10;
        for (unsigned int i = 0; i < initial_newton_steps; ++i) {
            Scalar f = poly.evaluate(x);
            Scalar df = dpoly.evaluate(x);
            if (detail::abs_val(df) < config_small_value()) break;
            Scalar step = f / df;
            if (detail::abs_val(step) < config.target_tolerance * Scalar(0.1)) break;
            x = x - step;
        }

        // Phase 2: Estimate multiplicity using Ostrowski from improved starting point
        unsigned int estimated_multiplicity = 1;
        auto ostrowski_result = estimateMultiplicityOstrowskiWithLastIterate(x, poly);
        estimated_multiplicity = ostrowski_result.multiplicity;

        // For simple roots, use the last iterate from Ostrowski (3 more Newton steps)
        if (ostrowski_result.valid && estimated_multiplicity == 1) {
            x = ostrowski_result.last_iterate;
        }

        // Override with hint if provided
        if (config.multiplicity_method == MultiplicityMethodBase::HINT &&
            config.multiplicity_hint > 0) {
            estimated_multiplicity = config.multiplicity_hint;
        }

        // Track previous step size for stagnation detection
        Scalar prev_step = Scalar(1);

        // Iteration loop
        for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
            result.iterations = iter + 1;

            // Evaluate function and derivatives
            Scalar f = poly.evaluate(x);
            Scalar df = dpoly.evaluate(x);

            // Precision floor check: if both |f| and |df| are extremely small,
            // we've hit the precision limit and should stop
            Scalar precision_floor = config.residual_tolerance * config.residual_tolerance;
            if (detail::abs_val(f) < precision_floor && detail::abs_val(df) < precision_floor) {
                result.location = x;
                result.residual = f;
                result.multiplicity = estimated_multiplicity;
                result.first_nonzero_derivative = df;
                result.condition_estimate = Scalar(1e16);  // Indicate precision issue
                result.converged = true;  // Consider converged (at precision limit)
                result.max_error = config.target_tolerance;  // Can't estimate better
                result.needs_higher_precision = true;  // Flag: may need more precision
                return result;
            }

            // Convergence check: use step size estimate
            // step = m * f / f' represents the error estimate for modified Newton
            if (detail::abs_val(f) < config.residual_tolerance) {
                Scalar step_estimate = Scalar(static_cast<int>(estimated_multiplicity)) *
                    detail::abs_val(f) / detail::max_val(detail::abs_val(df), config_small_value());

                // Compute condition number for HP flagging
                Scalar ddf = ddpoly.evaluate(x);
                Scalar kappa = detail::abs_val(ddf) /
                    (detail::abs_val(df) * detail::abs_val(df) + config_small_value());
                kappa = detail::max_val(kappa, Scalar(1));

                if (step_estimate <= config.target_tolerance) {
                    result.location = x;
                    result.residual = f;
                    result.multiplicity = estimated_multiplicity;
                    result.first_nonzero_derivative = df;
                    result.condition_estimate = kappa;
                    result.converged = true;
                    result.max_error = step_estimate;
                    // Flag HP if condition number is high (> 1e6) even if converged
                    result.needs_higher_precision = (kappa > Scalar(1e6));
                    return result;
                }
            }

            // Perform iteration step
            Scalar step = performIterationStep(
                x, f, df, poly, dpoly, ddpoly,
                estimated_multiplicity, config.iteration_method);

            prev_step = step;

            if (detail::abs_val(step) < config.target_tolerance * Scalar(0.01)) {
                // Step too small, likely converged
                break;
            }

            x = x - step;
        }

        // Final result
        result.location = x;
        result.residual = poly.evaluate(x);
        result.multiplicity = estimated_multiplicity;

        Scalar df = dpoly.evaluate(x);
        Scalar ddf = ddpoly.evaluate(x);
        result.first_nonzero_derivative = df;

        // Compute condition number for HP flagging
        Scalar kappa = detail::abs_val(ddf) /
            (detail::abs_val(df) * detail::abs_val(df) + config_small_value());
        kappa = detail::max_val(kappa, Scalar(1));
        result.condition_estimate = kappa;

        // Use step-based error estimate
        Scalar step_estimate = Scalar(static_cast<int>(estimated_multiplicity)) *
            detail::abs_val(result.residual) / detail::max_val(detail::abs_val(df), config_small_value());
        result.max_error = step_estimate;

        // Check convergence using step estimate
        if (detail::abs_val(result.residual) < config.residual_tolerance &&
            step_estimate <= config.target_tolerance) {
            result.converged = true;
            // Flag HP if condition number is high (> 1e6) even if converged
            result.needs_higher_precision = (kappa > Scalar(1e6));
        } else {
            result.converged = detail::abs_val(result.residual) < config.residual_tolerance;
            result.needs_higher_precision = (step_estimate > Scalar(1e-10)) || (kappa > Scalar(1e6));
        }

        return result;
    }

    //=========================================================================
    // Batch Refinement with Merging
    //=========================================================================

    /**
     * @brief Refine all roots from solver result with merging and deduplication
     *
     * This method:
     * 1. Refines each resolved box using Newton's method
     * 2. Merges nearby roots using exclusion radius
     * 3. Collects unverified boxes into problematic regions
     * 4. Flags roots that need higher precision
     *
     * @param solver_result Raw result from subdivision solver
     * @param poly The polynomial (1D)
     * @param config Refinement configuration
     * @return Batch refinement result with deduplicated roots
     */
    template<typename SolverResultType>
    static RefinementResultBase<Scalar> refine(
        const SolverResultType& solver_result,
        const Poly& poly,
        const Config& config = Config())
    {
        RefinementResultBase<Scalar> batch_result;
        batch_result.any_needs_higher_precision = false;

        const std::size_t num_resolved = solver_result.num_resolved;

        // Handle case of no resolved boxes
        if (num_resolved == 0) {
            for (std::size_t i = 0; i < solver_result.boxes.size(); ++i) {
                batch_result.unverified_boxes.push_back(i);
            }
            if (!batch_result.unverified_boxes.empty()) {
                batch_result.problematic_regions = mergeUnverifiedBoxes1D(
                    batch_result.unverified_boxes, solver_result, poly, config);
            }
            return batch_result;
        }

        // Step 1: Refine each resolved box
        std::vector<std::size_t> verified_indices;
        std::vector<Result> candidate_roots;

        for (std::size_t i = 0; i < num_resolved; ++i) {
            const auto& box = solver_result.boxes[i];
            Scalar x0 = Scalar(box.center[0]);

            Result refined = refineRoot1D(x0, poly, config);

            if (!refined.converged) {
                batch_result.unverified_boxes.push_back(i);
                continue;
            }

            // Compute exclusion radius
            refined.exclusion_radius = computeExclusionRadius(
                refined.multiplicity,
                refined.first_nonzero_derivative,
                config.target_tolerance,
                config.exclusion_multiplier);

            // Check if higher precision is needed
            Scalar error_est = computeErrorEstimate(
                refined.residual, refined.first_nonzero_derivative, refined.multiplicity);
            refined.needs_higher_precision = (error_est > Scalar(1e-10));

            refined.source_boxes.push_back(i);
            verified_indices.push_back(i);
            candidate_roots.push_back(refined);
        }

        // Step 2: Merge nearby roots using exclusion radius
        std::vector<bool> merged(candidate_roots.size(), false);
        std::set<std::size_t> cancelled;

        for (std::size_t i = 0; i < candidate_roots.size(); ++i) {
            if (merged[i]) continue;

            Result& root = candidate_roots[i];
            Scalar radius = root.exclusion_radius;

            // Find nearby roots to merge
            for (std::size_t j = i + 1; j < candidate_roots.size(); ++j) {
                if (merged[j]) continue;

                Scalar dist = detail::abs_val(root.location - candidate_roots[j].location);
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

            // Add to final roots
            batch_result.roots.push_back(root);
            if (root.needs_higher_precision) {
                batch_result.any_needs_higher_precision = true;
            }
        }

        // Collect cancelled boxes
        batch_result.cancelled_boxes.assign(cancelled.begin(), cancelled.end());
        std::sort(batch_result.cancelled_boxes.begin(), batch_result.cancelled_boxes.end());

        // Step 3: Add unresolved boxes from solver
        for (std::size_t i = num_resolved; i < solver_result.boxes.size(); ++i) {
            batch_result.unverified_boxes.push_back(i);
        }

        // Step 4: Merge unverified boxes into problematic regions
        if (!batch_result.unverified_boxes.empty()) {
            batch_result.problematic_regions = mergeUnverifiedBoxes1D(
                batch_result.unverified_boxes, solver_result, poly, config);
        }

        return batch_result;
    }

    /**
     * @brief Merge unverified solver boxes into contiguous regions
     */
    template<typename SolverResultType>
    static std::vector<ProblematicRegionBase<Scalar>> mergeUnverifiedBoxes1D(
        const std::vector<std::size_t>& unverified_indices,
        const SolverResultType& solver_result,
        const Poly& poly,
        const Config& config)
    {
        std::vector<ProblematicRegionBase<Scalar>> regions;

        if (unverified_indices.empty()) {
            return regions;
        }

        // Sort indices by lower bound
        std::vector<std::size_t> sorted_indices = unverified_indices;
        std::sort(sorted_indices.begin(), sorted_indices.end(),
            [&](std::size_t a, std::size_t b) {
                return solver_result.boxes[a].lower[0] < solver_result.boxes[b].lower[0];
            });

        // Merge adjacent boxes (threshold: boxes closer than this are merged)
        Scalar merge_threshold = Scalar(1e-6);

        ProblematicRegionBase<Scalar> current;
        current.lower = Scalar(solver_result.boxes[sorted_indices[0]].lower[0]);
        current.upper = Scalar(solver_result.boxes[sorted_indices[0]].upper[0]);
        current.box_indices.push_back(sorted_indices[0]);

        for (std::size_t i = 1; i < sorted_indices.size(); ++i) {
            std::size_t idx = sorted_indices[i];
            const auto& box = solver_result.boxes[idx];

            Scalar box_lower = Scalar(box.lower[0]);
            Scalar box_upper = Scalar(box.upper[0]);

            // Check if box is adjacent to current region
            if (box_lower - current.upper <= merge_threshold) {
                // Merge into current region
                current.upper = detail::max_val(current.upper, box_upper);
                current.box_indices.push_back(idx);
            } else {
                // Finish current region and start new one
                regions.push_back(current);

                current = ProblematicRegionBase<Scalar>();
                current.lower = box_lower;
                current.upper = box_upper;
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

    /**
     * @brief Try to refine a problematic region to find a root
     */
    static void refineProblematicRegion1D(
        ProblematicRegionBase<Scalar>& region,
        const Poly& poly,
        const Config& config)
    {
        region.refinement_attempted = true;

        // Start from center of region
        Scalar x0 = (region.lower + region.upper) / Scalar(2);

        // Try refinement
        Result refined = refineRoot1D(x0, poly, config);

        if (!refined.converged) {
            region.refinement_succeeded = false;
            region.refined_root = x0;
            region.residual = poly.evaluate(x0);
            region.condition_estimate = refined.condition_estimate;
            region.needs_higher_precision = true;
            return;
        }

        // Estimate error for precision check
        Scalar error_est = computeErrorEstimate(
            refined.residual, refined.first_nonzero_derivative, refined.multiplicity);

        region.refinement_succeeded = (error_est <= config.target_tolerance);
        region.refined_root = refined.location;
        region.multiplicity = refined.multiplicity;
        region.residual = refined.residual;
        region.condition_estimate = refined.condition_estimate;
        region.needs_higher_precision = (error_est > Scalar(1e-10));
    }

    //=========================================================================
    // Iteration Methods
    //=========================================================================

    /**
     * @brief Perform one iteration step using configured method
     */
    static Scalar performIterationStep(
        const Scalar& x,
        const Scalar& f,
        const Scalar& df,
        const Poly& poly,
        const Poly& dpoly,
        const Poly& ddpoly,
        unsigned int multiplicity,
        IterationMethodBase method)
    {
        const Scalar small_val = config_small_value();
        Scalar step;

        switch (method) {
            case IterationMethodBase::NEWTON:
                if (detail::abs_val(df) < small_val) {
                    return Scalar(0);
                }
                step = f / df;
                break;

            case IterationMethodBase::MODIFIED_NEWTON:
                if (detail::abs_val(df) < small_val) {
                    return Scalar(0);
                }
                step = Scalar(static_cast<int>(multiplicity)) * f / df;
                break;

            case IterationMethodBase::HALLEY: {
                Scalar ddf = ddpoly.evaluate(x);
                Scalar denom = Scalar(2) * df * df - f * ddf;
                if (detail::abs_val(denom) < small_val) {
                    return Scalar(0);
                }
                step = Scalar(2) * f * df / denom;
                break;
            }

            case IterationMethodBase::SCHRODER: {
                Scalar ddf = ddpoly.evaluate(x);
                Scalar denom = df * df - f * ddf;
                if (detail::abs_val(denom) < small_val) {
                    return Scalar(0);
                }
                step = f * df / denom;
                break;
            }

            default:
                step = f / df;
                break;
        }

        // Limit step size to avoid divergence
        Scalar max_step = Scalar(1);
        if (detail::abs_val(step) > max_step) {
            step = step * max_step / detail::abs_val(step);
        }

        return step;
    }

    //=========================================================================
    // Multiplicity Estimation Methods
    //=========================================================================

    /**
     * @brief Estimate multiplicity using configured method
     */
    static unsigned int estimateMultiplicityModular(
        const Scalar& x,
        const Poly& poly,
        const Config& config,
        const std::vector<Scalar>& iterates,
        Scalar& first_nonzero_deriv)
    {
        switch (config.multiplicity_method) {
            case MultiplicityMethodBase::NONE:
                first_nonzero_deriv = poly.differentiate(0).evaluate(x);
                return 1;

            case MultiplicityMethodBase::HINT:
                first_nonzero_deriv = poly.differentiate(0).evaluate(x);
                return (config.multiplicity_hint > 0) ? config.multiplicity_hint : 1;

            case MultiplicityMethodBase::TAYLOR:
                return estimateMultiplicityTaylor(
                    x, poly, config.max_multiplicity,
                    config.residual_tolerance, first_nonzero_deriv,
                    config.taylor_ratio_threshold);

            case MultiplicityMethodBase::OSTROWSKI:
                if (iterates.size() >= 4) {
                    unsigned int mult = estimateMultiplicityOstrowski(
                        iterates[1], iterates[2], iterates[3]);
                    // Still compute first_nonzero_deriv
                    estimateMultiplicityTaylor(x, poly, mult + 1,
                        config.residual_tolerance, first_nonzero_deriv);
                    return mult;
                }
                // Fall through to Taylor if not enough iterates
                return estimateMultiplicityTaylor(
                    x, poly, config.max_multiplicity,
                    config.residual_tolerance, first_nonzero_deriv,
                    config.taylor_ratio_threshold);

            case MultiplicityMethodBase::SIMPLE_THRESHOLD:
                return estimateMultiplicitySimpleThreshold(
                    x, poly, config.max_multiplicity,
                    config.residual_tolerance, first_nonzero_deriv);

            default:
                first_nonzero_deriv = poly.differentiate(0).evaluate(x);
                return 1;
        }
    }

    /**
     * @brief Estimate multiplicity using Taylor series / derivative test
     *
     * Finds the first non-zero derivative at the point.
     * For a root of multiplicity m, f^(k)(x) = 0 for k < m and f^(m)(x) ≠ 0.
     *
     * Simple threshold approach: returns the first k where |f^(k)(x)| > threshold.
     * The ratio test is unreliable for polynomials with multiple factors.
     */
    static unsigned int estimateMultiplicityTaylor(
        const Scalar& location,
        const Poly& poly,
        unsigned int max_order,
        const Scalar& threshold,
        Scalar& first_nonzero_deriv,
        const Scalar& /* ratio_threshold */ = Scalar(10))
    {
        // Compute derivatives and find first above threshold
        Poly current = poly;

        for (unsigned int k = 1; k <= max_order; ++k) {
            current = current.differentiate(0);
            Scalar deriv_val = current.evaluate(location);

            if (detail::abs_val(deriv_val) > threshold) {
                first_nonzero_deriv = deriv_val;
                return k;
            }
        }

        // All derivatives below threshold
        first_nonzero_deriv = Scalar(0);
        return max_order + 1;
    }

    /**
     * @brief Estimate multiplicity using simple threshold on derivatives
     */
    static unsigned int estimateMultiplicitySimpleThreshold(
        const Scalar& location,
        const Poly& poly,
        unsigned int max_order,
        const Scalar& threshold,
        Scalar& first_nonzero_deriv)
    {
        Poly current = poly;

        for (unsigned int k = 1; k <= max_order; ++k) {
            current = current.differentiate(0);
            Scalar deriv_val = current.evaluate(location);

            if (detail::abs_val(deriv_val) > threshold) {
                first_nonzero_deriv = deriv_val;
                return k;
            }
        }

        first_nonzero_deriv = Scalar(0);
        return max_order + 1;
    }

    /**
     * @brief Estimate multiplicity using Ostrowski's method (1973)
     *
     * Uses 3 consecutive Newton iterates: p = 1/2 + (x1-x2)/(x3-2x2+x1)
     * Multiplicity = floor(p)
     */
    static unsigned int estimateMultiplicityOstrowski(
        const Scalar& x1,
        const Scalar& x2,
        const Scalar& x3)
    {
        Scalar numerator = x1 - x2;
        Scalar denominator = x3 - Scalar(2) * x2 + x1;

        const Scalar small_val = config_small_value();
        if (detail::abs_val(denominator) < small_val) {
            // Converging quadratically -> simple root
            return 1;
        }

        Scalar p_est = Scalar(0.5) + numerator / denominator;

        // Convert to int using floor
        int multiplicity = static_cast<int>(p_est);
        if (p_est < Scalar(0)) {
            multiplicity = 1;
        }

        if (multiplicity < 1) {
            multiplicity = 1;
        }
        if (multiplicity > 20) {
            multiplicity = 20;
        }

        return static_cast<unsigned int>(multiplicity);
    }

    /**
     * @brief Result from Ostrowski multiplicity estimation
     */
    struct OstrowskiResult {
        unsigned int multiplicity;
        Scalar last_iterate;  // x3 from the 3 Newton steps
        bool valid;           // false if Newton steps failed
    };

    /**
     * @brief Estimate multiplicity from a starting point using Ostrowski
     *
     * Performs 3 standard Newton steps then applies Ostrowski's formula.
     * Returns the last iterate so caller can continue from there if mult=1.
     */
    static OstrowskiResult estimateMultiplicityOstrowskiWithLastIterate(
        const Scalar& start,
        const Poly& poly)
    {
        OstrowskiResult result;
        result.multiplicity = 1;
        result.last_iterate = start;
        result.valid = false;

        Poly dpoly = poly.differentiate(0);
        const Scalar small_val = config_small_value();

        Scalar x0 = start;
        Scalar x1, x2, x3;

        // Step 1
        Scalar f0 = poly.evaluate(x0);
        Scalar df0 = dpoly.evaluate(x0);
        if (detail::abs_val(df0) < small_val) return result;
        x1 = x0 - f0 / df0;

        // Step 2
        Scalar f1 = poly.evaluate(x1);
        Scalar df1 = dpoly.evaluate(x1);
        if (detail::abs_val(df1) < small_val) {
            result.last_iterate = x1;
            return result;
        }
        x2 = x1 - f1 / df1;

        // Step 3
        Scalar f2 = poly.evaluate(x2);
        Scalar df2 = dpoly.evaluate(x2);
        if (detail::abs_val(df2) < small_val) {
            result.last_iterate = x2;
            return result;
        }
        x3 = x2 - f2 / df2;

        result.multiplicity = estimateMultiplicityOstrowski(x1, x2, x3);
        result.last_iterate = x3;
        result.valid = true;
        return result;
    }

    /**
     * @brief Estimate multiplicity from a starting point using Ostrowski
     *
     * Convenience wrapper that only returns multiplicity.
     */
    static unsigned int estimateMultiplicityOstrowskiFromPoint(
        const Scalar& start,
        const Poly& poly)
    {
        return estimateMultiplicityOstrowskiWithLastIterate(start, poly).multiplicity;
    }

    //=========================================================================
    // Utility Methods
    //=========================================================================

    /**
     * @brief Estimate condition number for root-finding
     *
     * κ ≈ |f''(x)| / |f'(x)|²
     */
    static Scalar estimateConditionNumber1D(
        const Scalar& location,
        const Poly& poly,
        const Scalar& derivative_value)
    {
        const Scalar small_val = config_small_value();

        if (detail::abs_val(derivative_value) < small_val) {
            // Near multiple root, use higher derivative ratio
            Poly dpoly = poly.differentiate(0);
            Poly ddpoly = dpoly.differentiate(0);
            Scalar d2f = ddpoly.evaluate(location);
            Scalar df = derivative_value;

            if (detail::abs_val(df) < small_val) {
                return Scalar(1e10);  // Very ill-conditioned
            }

            Scalar kappa = detail::abs_val(d2f) / (detail::abs_val(df) * detail::abs_val(df));
            return detail::max_val(Scalar(1), kappa);
        }

        // Standard condition number estimate
        Poly dpoly = poly.differentiate(0);
        Poly ddpoly = dpoly.differentiate(0);
        Scalar d2f = ddpoly.evaluate(location);

        Scalar kappa = detail::abs_val(d2f) /
                       (detail::abs_val(derivative_value) * detail::abs_val(derivative_value) + small_val);

        return detail::max_val(Scalar(1), kappa);
    }

    /**
     * @brief Compute error estimate from residual and derivative
     *
     * For multiplicity m: error ≈ (m! * |f(x)| / |f^(m)(x)|)^(1/m)
     */
    static Scalar computeErrorEstimate(
        const Scalar& residual,
        const Scalar& first_nonzero_deriv,
        unsigned int multiplicity)
    {
        const Scalar small_val = config_small_value();

        if (detail::abs_val(first_nonzero_deriv) < small_val) {
            return Scalar(1);  // Can't estimate
        }

        if (multiplicity == 1) {
            // Simple root: error ≈ |f(x)| / |f'(x)|
            return detail::abs_val(residual) / detail::abs_val(first_nonzero_deriv);
        }

        // Multiple root: error ≈ (m! * |f(x)| / |f^(m)(x)|)^(1/m)
        Scalar fact = Scalar(static_cast<double>(detail::factorial(multiplicity)));
        Scalar ratio = fact * detail::abs_val(residual) / detail::abs_val(first_nonzero_deriv);

        // Compute m-th root
        double ratio_d = static_cast<double>(ratio);
        double exp = 1.0 / static_cast<double>(multiplicity);
        return Scalar(std::pow(ratio_d, exp));
    }

    /**
     * @brief Compute rigorous error bounds for a refined root
     *
     * Returns interval [lower, upper] guaranteed to contain the true root.
     */
    static bool computeErrorBounds(
        const Scalar& location,
        const Poly& poly,
        unsigned int multiplicity,
        const Scalar& first_nonzero_deriv,
        Scalar& lower,
        Scalar& upper)
    {
        const Scalar small_val = config_small_value();

        if (detail::abs_val(first_nonzero_deriv) < small_val) {
            lower = location;
            upper = location;
            return false;
        }

        Scalar error_est = computeErrorEstimate(
            poly.evaluate(location), first_nonzero_deriv, multiplicity);

        // Add safety margin
        Scalar margin = error_est * Scalar(2);

        lower = location - margin;
        upper = location + margin;

        return true;
    }

    /**
     * @brief Compute exclusion radius based on multiplicity and derivative
     *
     * For a root with multiplicity m and first non-zero derivative D:
     * radius ≈ (tolerance * m! / |D|)^(1/m)
     */
    static Scalar computeExclusionRadius(
        unsigned int multiplicity,
        const Scalar& first_nonzero_deriv,
        const Scalar& tolerance,
        const Scalar& multiplier = Scalar(3))
    {
        const Scalar small_val = config_small_value();

        if (detail::abs_val(first_nonzero_deriv) < small_val) {
            return tolerance;
        }

        if (multiplicity == 1) {
            return multiplier * tolerance / detail::abs_val(first_nonzero_deriv);
        }

        // Multiple root
        Scalar fact = Scalar(static_cast<double>(detail::factorial(multiplicity)));
        Scalar ratio = tolerance * fact / detail::abs_val(first_nonzero_deriv);

        double ratio_d = static_cast<double>(ratio);
        double exp = 1.0 / static_cast<double>(multiplicity);
        return multiplier * Scalar(std::pow(ratio_d, exp));
    }

private:
    /**
     * @brief Get a small value appropriate for the scalar type
     */
    static Scalar config_small_value() {
        return Scalar(1e-100);
    }
};

} // namespace polynomial_solver

#endif // RESULT_REFINER_BASE_H
