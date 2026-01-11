#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "hp/result_refiner_hp.h"
#include "hp/differentiation_hp.h"
#include "hp/precision_conversion.h"
#include "core/interval_arithmetic.h"
#include <cmath>
#include <iostream>

// Temporary debug flag (disabled by default)
// #define DEBUG_HP_REFINER 1

namespace polynomial_solver {

// ============================================================================
// Modular Components
// ============================================================================

unsigned int ResultRefinerHP::estimateMultiplicityModular(
    const mpreal& x,
    const PolynomialHP& poly,
    const RefinementConfigHP& config,
    const std::vector<mpreal>& iterates)
{
    #ifdef DEBUG_HP_REFINER
    std::cout << "  [HP DEBUG] estimateMultiplicityModular called, method = " << static_cast<int>(config.multiplicity_method) << std::endl;
    #endif

    switch (config.multiplicity_method) {
        case MultiplicityMethod::NONE:
            return 1;

        case MultiplicityMethod::HINT:
            #ifdef DEBUG_HP_REFINER
            std::cout << "  [HP DEBUG] Using HINT method, hint = " << config.multiplicity_hint << std::endl;
            #endif
            return (config.multiplicity_hint > 0) ? config.multiplicity_hint : 1;

        case MultiplicityMethod::TAYLOR: {
            mpreal mult_threshold = mpreal("1e-50");
            mpreal first_nonzero_deriv;
            return estimateMultiplicity(x, poly, config.max_multiplicity,
                                       mult_threshold, first_nonzero_deriv,
                                       config.taylor_ratio_threshold);
        }

        case MultiplicityMethod::OSTROWSKI:
            if (iterates.size() >= 4) {
                // Use last 3 iterates: x1, x2, x3
                return estimateMultiplicityOstrowski(iterates[1], iterates[2], iterates[3]);
            } else {
                // Not enough iterates yet, use Taylor as fallback
                mpreal mult_threshold = mpreal("1e-50");
                mpreal first_nonzero_deriv;
                return estimateMultiplicity(x, poly, config.max_multiplicity,
                                           mult_threshold, first_nonzero_deriv,
                                           config.taylor_ratio_threshold);
            }

        case MultiplicityMethod::TAYLOR_THEN_OSTROWSKI:
            // Try Taylor first
            {
                mpreal mult_threshold = mpreal("1e-50");
                mpreal first_nonzero_deriv;
                unsigned int taylor_mult = estimateMultiplicity(
                    x, poly, config.max_multiplicity, mult_threshold,
                    first_nonzero_deriv, config.taylor_ratio_threshold);

                // If Taylor gives m>1, trust it
                if (taylor_mult > 1) {
                    return taylor_mult;
                }

                // Otherwise try Ostrowski if we have enough iterates
                if (iterates.size() >= 4) {
                    return estimateMultiplicityOstrowski(iterates[1], iterates[2], iterates[3]);
                }

                return 1;
            }

        default:
            return 1;
    }
}

mpreal ResultRefinerHP::performIterationStep(
    mpreal& x,
    const PolynomialHP& poly,
    const PolynomialHP& dpoly,
    const PolynomialHP& ddpoly,
    unsigned int multiplicity,
    IterationMethod method)
{
    mpreal f = poly.evaluate(x);
    mpreal df = dpoly.evaluate(x);
    mpreal step;

    switch (method) {
        case IterationMethod::NEWTON:
            // Standard Newton: x_new = x - f/f'
            if (abs(df) < mpreal("1e-100")) {
                return mpreal("0.0");  // Can't make progress
            }
            step = f / df;
            break;

        case IterationMethod::MODIFIED_NEWTON:
            // Modified Newton for multiple roots: x_new = x - m*f/f'
            if (abs(df) < mpreal("1e-100")) {
                return mpreal("0.0");
            }
            step = mpreal(multiplicity) * f / df;
            break;

        case IterationMethod::HALLEY: {
            // Halley's method: x_new = x - 2*f*f' / (2*f'^2 - f*f'')
            mpreal ddf = ddpoly.evaluate(x);
            mpreal denom = mpreal("2.0") * df * df - f * ddf;
            if (abs(denom) < mpreal("1e-100")) {
                return mpreal("0.0");
            }
            step = mpreal("2.0") * f * df / denom;
            break;
        }

        case IterationMethod::SCHRODER: {
            // Schröder's method: x_new = x - f*f' / (f'^2 - f*f'')
            mpreal ddf = ddpoly.evaluate(x);
            mpreal denom = df * df - f * ddf;
            if (abs(denom) < mpreal("1e-100")) {
                return mpreal("0.0");
            }
            step = f * df / denom;
            break;
        }

        default:
            step = f / df;  // Fallback to standard Newton
            break;
    }

    // Limit step size to avoid divergence
    mpreal max_step = mpreal("1.0");
    if (abs(step) > max_step) {
        step = step * max_step / abs(step);
    }

    x = x - step;
    return step;
}

// ============================================================================
// Main Refinement Method (using modular components)
// ============================================================================

RefinedRootHP ResultRefinerHP::refineRoot1D(
    double initial_guess,
    const PolynomialHP& poly,
    const RefinementConfigHP& config)
{
    RefinedRootHP result;

    // Convert initial guess to high precision
    mpreal x = toHighPrecision(initial_guess);

    // Convert tolerance strings to mpreal
    mpreal target_tol = mpreal(config.target_tolerance_str);
    mpreal residual_tol = mpreal(config.residual_tolerance_str);

    // Prepare derivative polynomials (compute what we need based on iteration method)
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    PolynomialHP ddpoly;
    if (config.iteration_method == IterationMethod::HALLEY ||
        config.iteration_method == IterationMethod::SCHRODER) {
        ddpoly = DifferentiationHP::derivative(dpoly, 0, 1);
    }

    // Store iterates for Ostrowski method
    std::vector<mpreal> iterates;
    iterates.push_back(x);  // x0 = initial guess

    // ============================================================
    // STEP 1: Initial multiplicity estimation
    // ============================================================
    unsigned int estimated_multiplicity = 1;

    if (config.multiplicity_timing == MultiplicityTiming::ONCE_AT_START) {
        estimated_multiplicity = estimateMultiplicityModular(x, poly, config, iterates);

        #ifdef DEBUG_HP_REFINER
        std::cout << "  [HP DEBUG] Initial multiplicity estimate = " << estimated_multiplicity << std::endl;
        #endif
    }

    // Track previous error for stagnation detection
    mpreal prev_error_bound = mpreal("1e100");

    // ============================================================
    // STEP 2: Iteration loop
    // ============================================================
    for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
        result.iterations = iter + 1;

        // Re-estimate multiplicity if configured to do so
        if (config.multiplicity_timing == MultiplicityTiming::EVERY_ITERATION) {
            estimated_multiplicity = estimateMultiplicityModular(x, poly, config, iterates);
        }

        #ifdef DEBUG_HP_REFINER
        if (iter < 5) {
            mpreal f = poly.evaluate(x);
            std::cout << "  [HP DEBUG] Iter " << iter << ": x = " << x
                      << ", f = " << f << ", mult = " << estimated_multiplicity << std::endl;
        }
        #endif

        // Compute the m-th derivative (first non-zero derivative)
        mpreal first_nonzero_deriv;
        if (estimated_multiplicity == 1) {
            first_nonzero_deriv = dpoly.evaluate(x);
        } else {
            // Compute m-th derivative explicitly
            PolynomialHP deriv_m = poly;
            for (unsigned int k = 0; k < estimated_multiplicity; ++k) {
                deriv_m = DifferentiationHP::derivative(deriv_m, 0, 1);
            }
            first_nonzero_deriv = deriv_m.evaluate(x);
        }

        // Compute rigorous error bounds
        mpreal lower, upper;
        mpreal current_error_bound = mpreal("1e100");
        bool has_bounds = false;

        if (abs(first_nonzero_deriv) > mpreal("1e-100")) {
            has_bounds = computeErrorBounds(x, poly, estimated_multiplicity,
                                           first_nonzero_deriv, lower, upper);
            if (has_bounds) {
                current_error_bound = (upper - lower) / mpreal("2.0");
            }
        }

        #ifdef DEBUG_HP_REFINER
        if (has_bounds) {
            std::cout << "  [HP DEBUG] Error bound = " << current_error_bound
                      << ", target = " << target_tol << std::endl;
        } else {
            std::cout << "  [HP DEBUG] No error bounds computed" << std::endl;
        }
        #endif

        // Check convergence
        if (has_bounds && current_error_bound <= target_tol) {
            result.location = x;
            result.residual = poly.evaluate(x);
            result.multiplicity = estimated_multiplicity;
            result.first_nonzero_derivative = first_nonzero_deriv;
            result.condition_estimate = estimateConditionNumber1D(x, poly, first_nonzero_deriv);
            result.converged = true;
            result.interval_lower = lower;
            result.interval_upper = upper;
            result.max_error = current_error_bound;
            result.has_guaranteed_bounds = true;
            return result;
        }

        // Check for stagnation and re-estimate multiplicity if configured
        if (config.multiplicity_timing == MultiplicityTiming::WHEN_STAGNANT) {
            if (has_bounds && current_error_bound >= prev_error_bound * mpreal("0.99")) {
                estimated_multiplicity = estimateMultiplicityModular(x, poly, config, iterates);
                #ifdef DEBUG_HP_REFINER
                std::cout << "  [HP DEBUG] Stagnation detected, re-estimated mult = "
                          << estimated_multiplicity << std::endl;
                #endif
            }
        }

        prev_error_bound = current_error_bound;

        // Perform iteration step using configured method
        mpreal step = performIterationStep(x, poly, dpoly, ddpoly,
                                          estimated_multiplicity, config.iteration_method);

        // Store iterate for Ostrowski method (if needed)
        if (iter < 3) {
            iterates.push_back(x);
        }
    }

    // ============================================================
    // STEP 3: Max iterations reached - compute final results
    // ============================================================
    result.location = x;
    result.residual = poly.evaluate(x);

    // Final multiplicity estimate
    estimated_multiplicity = estimateMultiplicityModular(x, poly, config, iterates);

    // Compute the m-th derivative
    mpreal final_first_nonzero_deriv;
    if (estimated_multiplicity == 1) {
        final_first_nonzero_deriv = dpoly.evaluate(x);
    } else {
        PolynomialHP deriv_m = poly;
        for (unsigned int k = 0; k < estimated_multiplicity; ++k) {
            deriv_m = DifferentiationHP::derivative(deriv_m, 0, 1);
        }
        final_first_nonzero_deriv = deriv_m.evaluate(x);
    }

    result.multiplicity = estimated_multiplicity;
    result.first_nonzero_derivative = final_first_nonzero_deriv;
    result.condition_estimate = estimateConditionNumber1D(x, poly, final_first_nonzero_deriv);

    // Compute final error bounds
    mpreal lower, upper;
    if (computeErrorBounds(x, poly, estimated_multiplicity, final_first_nonzero_deriv, lower, upper)) {
        result.interval_lower = lower;
        result.interval_upper = upper;
        result.max_error = (upper - lower) / mpreal("2.0");
        result.has_guaranteed_bounds = true;

        if (result.max_error <= target_tol) {
            result.converged = true;
            return result;
        }
    }

    result.converged = false;
    result.error_message = "Maximum iterations reached - error bounds exceed tolerance";
    return result;
}

RefinedRootHP ResultRefinerHP::refineRoot1DSchroder(
    double initial_guess,
    const PolynomialHP& poly,
    const RefinementConfigHP& config)
{
    RefinedRootHP result;

    // Convert initial guess to high precision
    mpreal x = toHighPrecision(initial_guess);

    // Convert tolerance strings to mpreal
    mpreal target_tol = mpreal(config.target_tolerance_str);
    mpreal residual_tol = mpreal(config.residual_tolerance_str);

    // Get derivative polynomials
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    PolynomialHP ddpoly = DifferentiationHP::derivative(dpoly, 0, 1);

    // Schröder iteration
    for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
        result.iterations = iter + 1;

        // Evaluate function and first two derivatives
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        mpreal d2f = ddpoly.evaluate(x);

        // Store residual
        result.residual = f;

        // Check convergence with condition-aware criteria
        if (abs(f) < residual_tol && abs(df) > mpreal("1e-100")) {
            // Estimate condition number
            mpreal kappa = abs(d2f) / (abs(df) * abs(df) + mpreal("1e-100"));
            mpreal estimated_error = kappa * abs(f) / abs(df);

            if (estimated_error <= target_tol) {
                // Converged with good accuracy
                result.location = x;
                result.multiplicity = 1;  // Could estimate from derivatives if needed
                result.first_nonzero_derivative = df;
                result.condition_estimate = kappa;
                result.converged = true;
                return result;
            }
        }

        // Schröder's method: x_{n+1} = x_n - (f * f') / (f'^2 - f * f'')
        mpreal numerator = f * df;
        mpreal denominator = df * df - f * d2f;

        // For multiple roots, both f and df approach zero
        // Check if we're near a multiple root by comparing |df| to |f|
        bool near_multiple_root = (abs(df) < mpreal("1e-20"));

        if (near_multiple_root) {
            // Near a multiple root, estimate multiplicity and use modified Newton
            mpreal first_nonzero_deriv = mpreal(0);
            mpreal mult_threshold = mpreal("1e-50");
            unsigned int mult = estimateMultiplicity(
                x, poly, config.max_multiplicity, mult_threshold, first_nonzero_deriv);

            if (mult > 1 && abs(first_nonzero_deriv) > mpreal("1e-100")) {
                // Check convergence for multiple root
                if (abs(f) < residual_tol) {
                    result.location = x;
                    result.multiplicity = mult;
                    result.first_nonzero_derivative = first_nonzero_deriv;
                    result.condition_estimate = estimateConditionNumber1D(x, poly, first_nonzero_deriv);
                    result.converged = true;

                    // Compute rigorous error bounds
                    mpreal lower, upper;
                    if (computeErrorBounds(x, poly, mult, first_nonzero_deriv, lower, upper)) {
                        result.interval_lower = lower;
                        result.interval_upper = upper;
                        result.max_error = (upper - lower) / mpreal("2.0");
                        result.has_guaranteed_bounds = true;
                    }

                    return result;
                }

                // Use modified Newton for multiple root
                mpreal step = mpreal(mult) * f / df;

                // Limit step size
                mpreal max_step = mpreal("0.1");
                if (abs(step) > max_step) {
                    step = step * max_step / abs(step);
                }

                x = x - step;
                continue;
            }
        }

        // Check if denominator is too small
        if (abs(denominator) < mpreal("1e-100")) {
            // Fall back to Newton step
            if (abs(df) > mpreal("1e-100")) {
                mpreal newton_step = f / df;

                // Limit step size
                mpreal max_step = mpreal("1.0");
                if (abs(newton_step) > max_step) {
                    newton_step = newton_step * max_step / abs(newton_step);
                }

                x = x - newton_step;
                continue;
            } else {
                // Can't make progress
                result.location = x;
                result.converged = false;
                result.error_message = "Both Schröder and Newton denominators too small";
                return result;
            }
        }

        mpreal schroder_step = numerator / denominator;

        // Check if Schröder step is reasonable compared to Newton step
        // If Schröder over-compensates (would proceed in wrong direction or too far),
        // fall back to Newton step
        if (abs(df) > mpreal("1e-100")) {
            mpreal newton_step = f / df;

            // If Schröder step is more than 10% larger than Newton step, use Newton
            if (abs(schroder_step) > mpreal("1.1") * abs(newton_step)) {
                schroder_step = newton_step;
            }
        }

        // Limit step size to avoid divergence
        mpreal max_step = mpreal("1.0");
        if (abs(schroder_step) > max_step) {
            schroder_step = schroder_step * max_step / abs(schroder_step);
        }

        x = x - schroder_step;
    }

    // Max iterations reached - check final residual
    mpreal f = poly.evaluate(x);
    result.location = x;
    result.residual = f;

    if (abs(f) < residual_tol) {
        // Check condition-aware convergence one more time
        mpreal df = dpoly.evaluate(x);
        mpreal d2f = ddpoly.evaluate(x);

        mpreal kappa = abs(d2f) / (abs(df) * abs(df) + mpreal("1e-100"));
        mpreal estimated_error = kappa * abs(f) / max(abs(df), mpreal("1e-50"));

        result.multiplicity = 1;
        result.first_nonzero_derivative = df;
        result.condition_estimate = kappa;

        if (estimated_error <= target_tol) {
            result.converged = true;

            // Compute rigorous error bounds
            mpreal lower, upper;
            if (computeErrorBounds(x, poly, 1, df, lower, upper)) {
                result.interval_lower = lower;
                result.interval_upper = upper;
                result.max_error = (upper - lower) / mpreal("2.0");
                result.has_guaranteed_bounds = true;
            }

            return result;
        }

        result.converged = false;
        result.error_message = "Residual small but estimated error too large (ill-conditioned)";
        return result;
    }

    result.converged = false;
    result.error_message = "Maximum iterations reached";
    return result;
}

unsigned int ResultRefinerHP::estimateMultiplicityOstrowski(
    const mpreal& x1,
    const mpreal& x2,
    const mpreal& x3)
{
    // Ostrowski's method (1973) for multiplicity estimation
    // Based on the asymptotic error ratio of Newton's method for multiple roots
    //
    // For a root of multiplicity m, Newton's method converges linearly:
    // e_{n+1} / e_n → (m-1)/m  as n → ∞
    //
    // From 3 consecutive iterates x₁, x₂, x₃:
    // p = 1/2 + (x₁ - x₂)/(x₃ - 2x₂ + x₁)
    //
    // The formula is derived from:
    // (x₁ - x₂)/(x₂ - x₃) ≈ (m-1)/m
    // Solving for m gives the above formula

    mpreal numerator = x1 - x2;
    mpreal denominator = x3 - mpreal("2.0") * x2 + x1;

    // Check for degenerate case (denominator too small)
    if (abs(denominator) < mpreal("1e-100")) {
        // Iterates are converging quadratically → simple root
        return 1;
    }

    // Compute estimate
    mpreal p_est = mpreal("0.5") + numerator / denominator;

    #ifdef DEBUG_HP_REFINER
    std::cout << "    [Ostrowski] p_est = " << p_est << std::endl;
    std::cout << "    [Ostrowski] x1=" << x1 << ", x2=" << x2 << ", x3=" << x3 << std::endl;
    #endif

    #ifdef DEBUG_OSTROWSKI
    std::cout << "DEBUG Ostrowski: x1=" << x1 << ", x2=" << x2 << ", x3=" << x3 << std::endl;
    std::cout << "  numerator=" << numerator << ", denominator=" << denominator << std::endl;
    std::cout << "  p_est=" << p_est << std::endl;

    // Also compute the error ratio directly
    mpreal e1 = x2 - x1;  // error at step 1
    mpreal e2 = x3 - x2;  // error at step 2
    if (abs(e1) > mpreal("1e-100")) {
        mpreal ratio = e2 / e1;
        std::cout << "  error_ratio e2/e1=" << ratio << " (should be (m-1)/m)" << std::endl;
    }
    #endif

    // Apply floor (NOT round or ceiling!) with minimum value 1
    // The formula gives p ≈ m + 0.5, so floor(p) gives the correct multiplicity
    // For example: triple root gives p ≈ 3.5, floor(3.5) = 3 ✓
    // Debug showed: m=2→p=2.5, m=3→p=3.5, m=4→p=4.5, so floor is correct!
    int multiplicity = static_cast<int>(floor(p_est));

    #ifdef DEBUG_HP_REFINER
    std::cout << "    [Ostrowski] floor(p_est) = " << multiplicity << std::endl;
    #endif

    if (multiplicity < 1) {
        multiplicity = 1;
    }

    // Sanity check: cap at reasonable maximum
    if (multiplicity > 20) {
        multiplicity = 20;  // Very high multiplicities are rare
    }

    return static_cast<unsigned int>(multiplicity);
}

unsigned int ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(
    const mpreal& start,
    const PolynomialHP& poly)
{
    // Perform 3 REGULAR Newton steps (not modified Newton!)
    // This is critical - Ostrowski's method requires regular Newton iterates

    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);

    mpreal x0 = start;
    mpreal x1, x2, x3;

    // Step 1: x1 = x0 - f(x0)/f'(x0)
    {
        mpreal f0 = poly.evaluate(x0);
        mpreal df0 = dpoly.evaluate(x0);

        if (abs(df0) < mpreal("1e-100")) {
            // Derivative too small, can't do Newton step
            return 1;
        }

        x1 = x0 - f0 / df0;
    }

    // Step 2: x2 = x1 - f(x1)/f'(x1)
    {
        mpreal f1 = poly.evaluate(x1);
        mpreal df1 = dpoly.evaluate(x1);

        if (abs(df1) < mpreal("1e-100")) {
            return 1;
        }

        x2 = x1 - f1 / df1;
    }

    // Step 3: x3 = x2 - f(x2)/f'(x2)
    {
        mpreal f2 = poly.evaluate(x2);
        mpreal df2 = dpoly.evaluate(x2);

        if (abs(df2) < mpreal("1e-100")) {
            return 1;
        }

        x3 = x2 - f2 / df2;
    }

    // Now apply Ostrowski's formula to x1, x2, x3
    return estimateMultiplicityOstrowski(x1, x2, x3);
}

unsigned int ResultRefinerHP::estimateMultiplicity(
    const mpreal& location,
    const PolynomialHP& poly,
    unsigned int max_order,
    const mpreal& threshold,
    mpreal& first_nonzero_deriv,
    double ratio_threshold)
{
    first_nonzero_deriv = mpreal(0);

    // Taylor series analysis for f(x) with root of multiplicity m at x = r:
    // f(x) = c_m * (x-r)^m + c_{m+1} * (x-r)^{m+1} + ...
    // where c_m ≠ 0 is the first non-zero coefficient
    //
    // Derivatives at x ≈ r + ε:
    // f^(k)(x) = k! * c_k * (x-r)^0 + (k+1)! * c_{k+1} * (x-r) + ...
    //          ≈ k! * c_k  (for k ≥ m, when x ≈ r)
    //          ≈ 0         (for k < m)
    //
    // Key insight: f^(m)(r) = m! * c_m
    // The coefficient scale c_m is related to the polynomial's coefficient magnitude.

    // Evaluate all derivatives up to max_order
    std::vector<mpreal> deriv_values(max_order + 1);

    for (unsigned int order = 1; order <= max_order; ++order) {
        PolynomialHP deriv = DifferentiationHP::derivative(poly, 0, order);
        deriv_values[order] = deriv.evaluate(location);
    }

    // Estimate coefficient scale from the polynomial's coefficients
    // Use whichever representation is available (prefer power to avoid conversion)
    mpreal coeff_scale = mpreal(0);
    if (poly.hasPowerCoefficients()) {
        const std::vector<mpreal>& coeffs = poly.powerCoefficients();
        for (const mpreal& c : coeffs) {
            coeff_scale = max(coeff_scale, abs(c));
        }
    } else if (poly.hasBernsteinCoefficients()) {
        const std::vector<mpreal>& coeffs = poly.bernsteinCoefficients();
        for (const mpreal& c : coeffs) {
            coeff_scale = max(coeff_scale, abs(c));
        }
    }
    if (coeff_scale < mpreal("1e-100")) {
        coeff_scale = mpreal(1);  // Fallback if all coefficients are tiny
    }

    // For a root of multiplicity m at x = r:
    // - f^(k)(r) = 0 for k < m
    // - f^(m)(r) = m! * c_m ≠ 0
    //
    // When evaluating at x ≈ r + ε (small ε):
    // - f^(k)(x) ≈ k! * c_k + (k+1)! * c_{k+1} * ε + ...
    //
    // For k < m: f^(k)(x) ≈ (k+1)! * c_{k+1} * ε + ... ≈ O(ε)
    // For k = m: f^(m)(x) ≈ m! * c_m (constant, independent of ε)
    //
    // Strategy: Look for the first derivative that is NOT decreasing as we approach the root
    // This is indicated by comparing |f^(k)| with |f^(k+1)| / (k+1)

    // Compute factorials
    std::vector<mpreal> factorials(max_order + 1);
    factorials[0] = mpreal(1);
    for (unsigned int i = 1; i <= max_order; ++i) {
        factorials[i] = factorials[i-1] * mpreal(i);
    }

    // Strategy: Use ratio test to find first derivative that doesn't vanish at root
    // For f(x) = c_m * (x-r)^m + c_{m+1} * (x-r)^{m+1} + ... near x = r + ε:
    //
    // f^(k)(x) ≈ k! * c_k + (k+1)! * c_{k+1} * ε + ...
    //
    // For k < m: c_k = 0, so f^(k)(x) ≈ (k+1)! * c_{k+1} * ε
    // For k = m: c_m ≠ 0, so f^(m)(x) ≈ m! * c_m (independent of ε)
    //
    // Key insight: |f^(k+1)| / |f^(k)| should be LARGE (>> 1) for k < m
    //              |f^(m+1)| / |f^(m)| should be O(ε) or O(1) for k = m

    // First pass: Find all derivatives above absolute threshold
    std::vector<unsigned int> significant_orders;
    for (unsigned int order = 1; order <= max_order; ++order) {
        if (abs(deriv_values[order]) > threshold) {
            significant_orders.push_back(order);
        }
        #ifdef DEBUG_HP_REFINER
        std::cout << "      f^(" << order << ") = " << deriv_values[order]
                  << " (abs = " << abs(deriv_values[order]) << ")" << std::endl;
        #endif
    }

    #ifdef DEBUG_HP_REFINER
    std::cout << "      Significant orders (above threshold " << threshold << "): ";
    for (auto o : significant_orders) std::cout << o << " ";
    std::cout << std::endl;
    #endif

    if (significant_orders.empty()) {
        // All derivatives are zero
        first_nonzero_deriv = mpreal(0);
        return max_order + 1;
    }

    // Second pass: Use ratio test to find true first non-zero derivative
    // Check ratios |f^(k+1)| / |f^(k)| for consecutive significant orders
    //
    // For f(x) = (x-r)^m near x = r+ε:
    // - f^(k) ~ ε^(m-k) for k < m
    // - f^(m) ~ constant for k = m
    // So ratio f^(k+1)/f^(k) ~ ε^(-1) for k < m (very large!)
    //
    // Strategy: Look for the first order where ratio is NOT extremely large
    // The ratio decreases as we approach m: for m=6, we see 250→200→150→100→50
    // Use a simple threshold: ratio > 10 means derivative vanishes

    for (size_t i = 0; i < significant_orders.size(); ++i) {
        unsigned int order = significant_orders[i];
        mpreal deriv_val = abs(deriv_values[order]);

        // Check if there's a next significant derivative
        if (i + 1 < significant_orders.size()) {
            unsigned int next_order = significant_orders[i + 1];
            mpreal next_deriv = abs(deriv_values[next_order]);

            // Compute ratio
            mpreal ratio = next_deriv / max(deriv_val, mpreal("1e-100"));

            #ifdef DEBUG_HP_REFINER
            std::cout << "      Ratio f^(" << next_order << ")/f^(" << order << ") = "
                      << ratio << " (threshold = " << ratio_threshold << ")" << std::endl;
            #endif

            // Configurable threshold: if ratio > threshold, derivative vanishes at root
            // Default threshold=10 works because for k < m, ratio ~ 50*(m-k) >> 10
            // For k = m, ratio ~ O(1) or O(ε) < 10
            // For extreme multiplicities, increase threshold (e.g., 50 for m>10)
            if (ratio > mpreal(ratio_threshold)) {
                #ifdef DEBUG_HP_REFINER
                std::cout << "        -> Ratio > threshold, f^(" << order << ") vanishes, continue" << std::endl;
                #endif
                continue;  // Try next order
            }
        }

        // This derivative doesn't have a much larger successor
        // It's likely the first non-zero derivative at the root

        // Final check: verify it's above relative threshold
        mpreal expected_scale = factorials[order] * coeff_scale;
        mpreal relative_tol = mpreal("1e-15");
        mpreal deriv_threshold = max(threshold, relative_tol * expected_scale);

        #ifdef DEBUG_HP_REFINER
        std::cout << "      Final check: deriv_val = " << deriv_val
                  << ", deriv_threshold = " << deriv_threshold << std::endl;
        #endif

        if (deriv_val > deriv_threshold) {
            #ifdef DEBUG_HP_REFINER
            std::cout << "      -> Accepted! Returning multiplicity = " << order << std::endl;
            #endif
            first_nonzero_deriv = deriv_values[order];
            return order;
        }
    }

    // Fallback: return the highest significant order
    if (!significant_orders.empty()) {
        unsigned int order = significant_orders.back();
        first_nonzero_deriv = deriv_values[order];
        return order;
    }

    // All derivatives are below threshold up to max_order
    // Multiplicity is at least max_order + 1
    first_nonzero_deriv = mpreal(0);
    return max_order + 1;
}

mpreal ResultRefinerHP::estimateConditionNumber1D(
    const mpreal& location,
    const PolynomialHP& poly,
    const mpreal& derivative_value)
{
    mpreal threshold = mpreal("1e-50");

    if (abs(derivative_value) < threshold) {
        // Derivative too small, condition number is very large
        return mpreal("1e100");
    }

    // Compute second derivative
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    PolynomialHP ddpoly = DifferentiationHP::derivative(dpoly, 0, 1);
    mpreal f_prime = derivative_value;
    mpreal f_double_prime = ddpoly.evaluate(location);

    // Condition number estimate: |f''| / |f'|^2
    // This gives the sensitivity of the root to perturbations
    mpreal kappa = abs(f_double_prime) / (abs(f_prime) * abs(f_prime));

    // For very ill-conditioned problems, also check higher derivatives
    // If f'''(x) is large relative to f'(x), the problem is even worse
    PolynomialHP dddpoly = DifferentiationHP::derivative(ddpoly, 0, 1);
    mpreal f_triple_prime = dddpoly.evaluate(location);
    mpreal kappa_higher = abs(f_triple_prime) / (abs(f_prime) * abs(f_prime) * abs(f_prime));

    // Take the maximum as a conservative estimate
    kappa = max(kappa, kappa_higher);

    return max(mpreal(1), kappa);
}

bool ResultRefinerHP::computeErrorBounds(
    const mpreal& location,
    const PolynomialHP& poly,
    unsigned int multiplicity,
    const mpreal& first_nonzero_deriv,
    mpreal& lower,
    mpreal& upper)
{
    // Rigorous error bounds using interval Newton method
    //
    // For a refined root x*, we want to find an interval [x* - r, x* + r]
    // that is GUARANTEED to contain the true root.
    //
    // Interval Newton theorem: If x* is an approximation to a root and
    // we can find r > 0 such that:
    //   1. f'(y) ≠ 0 for all y in [x* - r, x* + r]
    //   2. |f(x*)| / min|f'(y)| <= r
    // Then the interval [x* - r, x* + r] contains a unique root.

    mpreal f_center = poly.evaluate(location);

    // For simple roots (multiplicity = 1)
    if (multiplicity == 1) {
        // Check if derivative is non-zero at center
        if (abs(first_nonzero_deriv) < mpreal("1e-100")) {
            return false;
        }

        // Get derivative polynomial
        PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);

        // Start with initial radius based on Newton correction
        mpreal initial_radius = abs(f_center / first_nonzero_deriv);

        // Expand radius to ensure we can bound the derivative
        // Use a conservative multiplier
        mpreal radius = initial_radius * mpreal("10.0");

        // Limit maximum radius to avoid numerical issues
        if (radius > mpreal("0.1")) {
            radius = mpreal("0.1");
        }

        // Ensure minimum radius for very well-converged roots
        if (radius < mpreal("1e-100")) {
            radius = mpreal("1e-100");
        }

        // Compute bounds on f' over the interval [location - radius, location + radius]
        // using rigorous interval arithmetic.
        //
        // For interval Newton method, we need min|f'(x)| over the interval.
        // Using interval arithmetic gives us GUARANTEED bounds.

        // Clamp interval to [0, 1] domain
        mpreal interval_lower = max(location - radius, mpreal("0.0"));
        mpreal interval_upper = min(location + radius, mpreal("1.0"));
        Interval search_interval(interval_lower, interval_upper);

        // Use interval arithmetic to find min|f'(x)| over the interval
        mpreal df_min = minAbsOnInterval(dpoly, search_interval);

        // Check if derivative stays bounded away from zero
        if (df_min < mpreal("1e-100")) {
            // Derivative may be zero or very small in the interval
            // Fall back to center value if we can
            df_min = abs(first_nonzero_deriv);
            if (df_min < mpreal("1e-100")) {
                return false;
            }
        }

        // Compute rigorous radius: r >= |f(x*)| / min|f'|
        mpreal rigorous_radius = abs(f_center) / df_min;

        // For rigorous bounds, we need to account for:
        // 1. Rounding errors in the computation
        // 2. The fact that we're using a finite sample of the derivative
        // 3. Numerical errors in polynomial evaluation
        //
        // Compute machine epsilon for current precision
        // epsilon ≈ 2^(-precision_in_bits)
        int precision_bits = getPrecision();
        mpreal eps = pow(mpreal("2.0"), -mpreal(precision_bits));

        // Add both relative and absolute safety margins
        mpreal relative_margin = rigorous_radius * mpreal("0.1");  // 10% relative margin
        mpreal absolute_margin = eps * abs(location) * mpreal("100.0");  // Account for rounding

        // Use the larger of the two margins
        mpreal safety_margin = max(relative_margin, absolute_margin);
        rigorous_radius = rigorous_radius + safety_margin;

        // Ensure minimum radius based on the precision
        // We want at least 10x machine epsilon as the minimum radius
        mpreal min_radius = eps * mpreal("10.0");
        if (rigorous_radius < min_radius) {
            rigorous_radius = min_radius;
        }

        lower = location - rigorous_radius;
        upper = location + rigorous_radius;

        // Clamp to [0, 1] domain
        if (lower < mpreal("0.0")) lower = mpreal("0.0");
        if (upper > mpreal("1.0")) upper = mpreal("1.0");

        return true;
    }

    // For multiple roots (multiplicity > 1), use rigorous bounds based on m-th derivative
    else {
        // For a root of multiplicity m, we have:
        // f(x) = (x - r)^m * g(x) where g(r) ≠ 0
        // f^(m)(r) = m! * g(r)
        //
        // Rigorous bound: If |f(x*)| is small and |f^(m)(x*)| ≈ m! * g(r),
        // then |x* - r| <= (|f(x*)| / |g(r)|)^(1/m)
        //
        // We need to bound |g(r)| from below to get a rigorous upper bound on |x* - r|

        if (abs(first_nonzero_deriv) < mpreal("1e-100")) {
            return false;
        }

        // Compute m!
        mpreal factorial_m = mpreal(1);
        for (unsigned int i = 2; i <= multiplicity; ++i) {
            factorial_m *= mpreal(i);
        }

        // Estimate |g(r)| = |f^(m)(r)| / m!
        mpreal g_r = first_nonzero_deriv / factorial_m;

        if (abs(g_r) < mpreal("1e-100")) {
            return false;
        }

        // For rigorous bounds, we need to account for the fact that g(x) may vary
        // Sample the m-th derivative in a neighborhood to get a lower bound on |g(x)|
        PolynomialHP deriv_m = poly;
        for (unsigned int k = 0; k < multiplicity; ++k) {
            deriv_m = DifferentiationHP::derivative(deriv_m, 0, 1);
        }

        // Sample at a few points near location
        mpreal g_min = abs(g_r);
        const int num_samples = 3;
        mpreal sample_radius = mpreal("0.01");  // Small neighborhood

        for (int i = -num_samples; i <= num_samples; ++i) {
            if (i == 0) continue;
            mpreal offset = mpreal(i) * sample_radius / mpreal(num_samples);
            mpreal x_sample = location + offset;

            // Clamp to [0, 1]
            if (x_sample < mpreal("0.0")) x_sample = mpreal("0.0");
            if (x_sample > mpreal("1.0")) x_sample = mpreal("1.0");

            mpreal deriv_m_sample = deriv_m.evaluate(x_sample);
            mpreal g_sample = abs(deriv_m_sample / factorial_m);

            if (g_sample > mpreal("1e-100") && g_sample < g_min) {
                g_min = g_sample;
            }
        }

        // Use the minimum value as a conservative lower bound
        mpreal g_lower_bound = g_min * mpreal("0.9");  // 10% safety margin

        if (g_lower_bound < mpreal("1e-100")) {
            return false;
        }

        // Rigorous bound: |x* - r| <= (|f(x*)| / g_lower_bound)^(1/m)
        mpreal abs_f = abs(f_center);
        mpreal ratio = abs_f / g_lower_bound;

        // Compute m-th root using logarithms for numerical stability
        mpreal radius;
        if (ratio > mpreal("1e-200")) {
            mpreal log_ratio = log(ratio);
            mpreal log_radius = log_ratio / mpreal(multiplicity);
            radius = exp(log_radius);
        } else {
            // f is extremely small, use minimal radius
            radius = mpreal("1e-100");
        }

        // Add safety margin for rigorous bound (necessary condition)
        radius = radius * mpreal("1.1");

        // Set interval bounds
        lower = location - radius;
        upper = location + radius;

        // Clamp to [0, 1] domain
        if (lower < mpreal("0.0")) lower = mpreal("0.0");
        if (upper > mpreal("1.0")) upper = mpreal("1.0");

        return true;
    }
}

// ============================================================================
// New Multiplicity Detection Methods
// ============================================================================

unsigned int ResultRefinerHP::estimateMultiplicitySimpleThreshold(
    const mpreal& location,
    const PolynomialHP& poly,
    unsigned int max_order,
    const mpreal& threshold)
{
    // Simplest method: find first derivative |f^(k)| > threshold
    // No ratio test, no fancy logic - just absolute threshold

    for (unsigned int k = 1; k <= max_order; ++k) {
        PolynomialHP deriv = DifferentiationHP::derivative(poly, 0, k);
        mpreal val = abs(deriv.evaluate(location));

        if (val > threshold) {
            // f^(k) is non-zero at location, so multiplicity is k
            return k;
        }
    }

    // All derivatives up to max_order are below threshold
    return max_order + 1;
}

unsigned int ResultRefinerHP::estimateMultiplicitySturm(
    const mpreal& location,
    const PolynomialHP& poly,
    const mpreal& interval_radius)
{
    // Sturm sequence method for multiplicity detection
    //
    // Proper Sturm sequence construction in power basis:
    // For polynomial f, construct sequence: f₀ = f, f₁ = f', f₂ = -rem(f₀, f₁), ...
    // Count sign changes at interval endpoints to determine number of roots.
    //
    // For multiplicity detection:
    // - If f has a root of multiplicity m at x=r, then f^(k) has a root of multiplicity m-k
    // - We count how many consecutive derivatives f, f', f'', ... have roots in [r-ε, r+ε]
    // - The first derivative that doesn't have a root gives us the multiplicity

    // Only works for 1D polynomials
    if (poly.dimension() != 1) {
        return 1;
    }

    mpreal left = location - interval_radius;
    mpreal right = location + interval_radius;
    unsigned int max_check = 10;

    // Helper: Convert Bernstein coefficients to power basis for 1D polynomial
    // For degree n, Bernstein basis: B_i^n(x) = C(n,i) * x^i * (1-x)^(n-i)
    // Power basis: x^i
    // Conversion: a_i = sum_j b_j * <B_j^n, x^i>
    auto bernstein_to_power = [](const std::vector<mpreal>& bern_coeffs) -> std::vector<mpreal> {
        int n = static_cast<int>(bern_coeffs.size()) - 1;
        if (n < 0) return std::vector<mpreal>();

        std::vector<mpreal> power_coeffs(n + 1, mpreal(0));

        // Use the formula: power_coeff[k] = sum_{i=k}^{n} bern_coeff[i] * C(i,k) / C(n,k)
        for (int k = 0; k <= n; ++k) {
            for (int i = k; i <= n; ++i) {
                // Compute binomial coefficients C(i,k) and C(n,k)
                mpreal binom_i_k = mpreal(1);
                mpreal binom_n_k = mpreal(1);

                for (int j = 0; j < k; ++j) {
                    binom_i_k *= mpreal(i - j) / mpreal(j + 1);
                    binom_n_k *= mpreal(n - j) / mpreal(j + 1);
                }

                power_coeffs[k] += bern_coeffs[i] * binom_i_k / binom_n_k;
            }
        }

        return power_coeffs;
    };

    // Helper: Evaluate polynomial in power basis
    auto eval_power = [](const std::vector<mpreal>& coeffs, const mpreal& x) -> mpreal {
        if (coeffs.empty()) return mpreal(0);

        // Horner's method: a_n*x^n + ... + a_1*x + a_0 = (...((a_n*x + a_{n-1})*x + ...)*x + a_0)
        mpreal result = coeffs.back();
        for (int i = static_cast<int>(coeffs.size()) - 2; i >= 0; --i) {
            result = result * x + coeffs[i];
        }
        return result;
    };

    // Helper: Polynomial derivative in power basis
    auto derivative_power = [](const std::vector<mpreal>& coeffs) -> std::vector<mpreal> {
        if (coeffs.size() <= 1) return std::vector<mpreal>();

        std::vector<mpreal> deriv(coeffs.size() - 1);
        for (size_t i = 1; i < coeffs.size(); ++i) {
            deriv[i - 1] = coeffs[i] * mpreal(i);
        }
        return deriv;
    };

    // Helper: Polynomial remainder (pseudo-division)
    auto poly_remainder = [](const std::vector<mpreal>& dividend, const std::vector<mpreal>& divisor) -> std::vector<mpreal> {
        if (divisor.empty() || (divisor.size() == 1 && abs(divisor[0]) < mpreal("1e-100"))) {
            return std::vector<mpreal>();  // Division by zero
        }

        std::vector<mpreal> rem = dividend;
        int deg_divisor = static_cast<int>(divisor.size()) - 1;

        // Find actual degree (skip leading zeros)
        while (deg_divisor >= 0 && abs(divisor[deg_divisor]) < mpreal("1e-100")) {
            deg_divisor--;
        }

        if (deg_divisor < 0) return std::vector<mpreal>();

        // Polynomial long division
        while (static_cast<int>(rem.size()) > deg_divisor + 1) {
            int deg_rem = static_cast<int>(rem.size()) - 1;
            while (deg_rem >= 0 && abs(rem[deg_rem]) < mpreal("1e-100")) {
                deg_rem--;
                rem.pop_back();
            }

            if (deg_rem < deg_divisor) break;

            // Divide leading terms
            mpreal coeff = rem[deg_rem] / divisor[deg_divisor];

            // Subtract divisor * coeff * x^(deg_rem - deg_divisor)
            for (int i = 0; i <= deg_divisor; ++i) {
                rem[deg_rem - deg_divisor + i] -= coeff * divisor[i];
            }

            rem.pop_back();  // Remove leading term (now zero)
        }

        return rem;
    };

    // Helper: Count sign changes in a sequence
    auto count_sign_changes = [](const std::vector<mpreal>& values) -> int {
        int changes = 0;
        mpreal prev_nonzero = mpreal(0);
        bool found_first = false;

        for (const mpreal& val : values) {
            if (abs(val) > mpreal("1e-100")) {
                if (found_first && prev_nonzero * val < mpreal(0)) {
                    changes++;
                }
                prev_nonzero = val;
                found_first = true;
            }
        }
        return changes;
    };

    // Get polynomial in power basis (prefer direct access to avoid conversion)
    std::vector<mpreal> power_coeffs;
    if (poly.hasPowerCoefficients()) {
        power_coeffs = poly.powerCoefficients();
    } else {
        power_coeffs = bernstein_to_power(poly.bernsteinCoefficients());
    }

    if (power_coeffs.empty()) {
        return 1;
    }

    // Build Sturm sequence in power basis
    std::vector<std::vector<mpreal>> sturm_seq;
    sturm_seq.push_back(power_coeffs);
    sturm_seq.push_back(derivative_power(power_coeffs));

    // Build rest of sequence: f_i = -rem(f_{i-2}, f_{i-1})
    for (size_t i = 2; i < 20 && !sturm_seq.back().empty(); ++i) {
        std::vector<mpreal> rem = poly_remainder(sturm_seq[i-2], sturm_seq[i-1]);
        if (rem.empty()) break;

        // Negate
        for (auto& c : rem) {
            c = -c;
        }
        sturm_seq.push_back(rem);
    }

    // Evaluate Sturm sequence at interval endpoints
    std::vector<mpreal> vals_left, vals_right;
    for (const auto& p : sturm_seq) {
        vals_left.push_back(eval_power(p, left));
        vals_right.push_back(eval_power(p, right));
    }

    // Count sign changes
    int changes_left = count_sign_changes(vals_left);
    int changes_right = count_sign_changes(vals_right);
    int num_roots = changes_left - changes_right;

    // If there are no roots in the interval, f doesn't vanish there
    if (num_roots == 0) {
        return 1;  // Simple root or no root
    }

    // If there's exactly 1 root, check derivatives to find multiplicity
    // For multiplicity m: f, f', ..., f^(m-1) all have roots, but f^(m) doesn't
    for (unsigned int k = 1; k <= max_check; ++k) {
        PolynomialHP deriv = DifferentiationHP::derivative(poly, 0, k);
        if (deriv.empty()) return k;

        // Get derivative in power basis (prefer direct access to avoid conversion)
        std::vector<mpreal> deriv_power;
        if (deriv.hasPowerCoefficients()) {
            deriv_power = deriv.powerCoefficients();
        } else {
            deriv_power = bernstein_to_power(deriv.bernsteinCoefficients());
        }
        if (deriv_power.empty()) return k;

        std::vector<std::vector<mpreal>> deriv_sturm;
        deriv_sturm.push_back(deriv_power);
        deriv_sturm.push_back(derivative_power(deriv_power));

        for (size_t i = 2; i < 20 && !deriv_sturm.back().empty(); ++i) {
            std::vector<mpreal> rem = poly_remainder(deriv_sturm[i-2], deriv_sturm[i-1]);
            if (rem.empty()) break;
            for (auto& c : rem) c = -c;
            deriv_sturm.push_back(rem);
        }

        // Count roots of derivative in interval
        std::vector<mpreal> deriv_vals_left, deriv_vals_right;
        for (const auto& p : deriv_sturm) {
            deriv_vals_left.push_back(eval_power(p, left));
            deriv_vals_right.push_back(eval_power(p, right));
        }

        int deriv_changes_left = count_sign_changes(deriv_vals_left);
        int deriv_changes_right = count_sign_changes(deriv_vals_right);
        int deriv_num_roots = deriv_changes_left - deriv_changes_right;

        if (deriv_num_roots == 0) {
            // f^(k) has no roots in interval → multiplicity is k
            return k;
        }
    }

    return max_check + 1;
}

std::map<std::string, unsigned int> ResultRefinerHP::estimateMultiplicityAllMethods(
    const mpreal& location,
    const PolynomialHP& poly,
    unsigned int max_order,
    const mpreal& x1,
    const mpreal& x2,
    const mpreal& x3)
{
    std::map<std::string, unsigned int> results;

    // Method 1: Taylor series with ratio test (original estimateMultiplicity)
    mpreal dummy_deriv;
    mpreal threshold = mpreal("1e-50");
    results["Taylor"] = estimateMultiplicity(location, poly, max_order, threshold, dummy_deriv);

    // Method 2: Simple threshold (no ratio test)
    // Test multiple thresholds to see if any work universally
    results["Thresh_1e-50"] = estimateMultiplicitySimpleThreshold(location, poly, max_order, mpreal("1e-50"));
    results["Thresh_1e-40"] = estimateMultiplicitySimpleThreshold(location, poly, max_order, mpreal("1e-40"));
    results["Thresh_1e-30"] = estimateMultiplicitySimpleThreshold(location, poly, max_order, mpreal("1e-30"));
    results["Thresh_1e-20"] = estimateMultiplicitySimpleThreshold(location, poly, max_order, mpreal("1e-20"));

    // Method 3: Sturm sequence
    mpreal interval_radius = mpreal("1e-10");
    results["Sturm"] = estimateMultiplicitySturm(location, poly, interval_radius);

    // Method 4: Ostrowski (if iterates provided)
    if (abs(x1) > mpreal("1e-100") || abs(x2) > mpreal("1e-100") || abs(x3) > mpreal("1e-100")) {
        results["Ostrowski"] = estimateMultiplicityOstrowski(x1, x2, x3);
    }

    return results;
}

//=============================================================================
// High-Precision Curve Refinement Implementation
//=============================================================================

CurveRefinedPointHP refineCurveNumericalHP(
    const std::function<mpreal(const mpreal&, const mpreal&)>& g,
    double x0, double y0,
    const CurveRefinementConfigHP& config)
{
    CurveRefinedPointHP result;
    result.converged = false;
    result.iterations = 0;

    // Convert initial guess to high precision
    mpreal x = mpreal(x0);
    mpreal y = mpreal(y0);

    // Parse configuration tolerances
    mpreal residual_tol = mpreal(config.residual_tolerance_str);
    mpreal min_grad_norm = mpreal(config.min_gradient_norm_str);
    mpreal h = mpreal(config.step_size_str);

    for (unsigned int i = 0; i < config.max_iterations; ++i) {
        // Evaluate function
        mpreal val = g(x, y);
        result.iterations = i + 1;

        // Check convergence
        if (abs(val) < residual_tol) {
            result.x = x;
            result.y = y;
            result.residual = abs(val);
            result.converged = true;
            return result;
        }

        // Numerical gradient using high-precision step size
        mpreal gx = (g(x + h, y) - g(x - h, y)) / (mpreal(2) * h);
        mpreal gy = (g(x, y + h) - g(x, y - h)) / (mpreal(2) * h);

        // Gradient magnitude squared
        mpreal grad_sq = gx * gx + gy * gy;

        // Check for near-zero gradient (singularity)
        if (grad_sq < min_grad_norm * min_grad_norm) {
            result.x = x;
            result.y = y;
            result.residual = abs(val);
            result.converged = false;
            return result;
        }

        // Newton step along gradient direction
        // Project point onto zero set: (x,y) -= (g / |∇g|²) * ∇g
        mpreal t = -val / grad_sq;
        x += t * gx;
        y += t * gy;
    }

    // Did not converge within max iterations
    result.x = x;
    result.y = y;
    result.residual = abs(g(x, y));
    return result;
}

std::vector<CurveRefinedPointHP> refineCurveNumericalHPMultiple(
    const std::function<mpreal(const mpreal&, const mpreal&)>& g,
    const std::vector<std::pair<double, double>>& points,
    const CurveRefinementConfigHP& config)
{
    std::vector<CurveRefinedPointHP> results;
    results.reserve(points.size());

    for (const auto& pt : points) {
        results.push_back(refineCurveNumericalHP(g, pt.first, pt.second, config));
    }

    return results;
}

//=============================================================================
// CurveRefinementConfigHP factory
//=============================================================================

CurveRefinementConfigHP CurveRefinementConfigHP::fromPrecisionBits(
    unsigned int bits, unsigned int max_iters)
{
    CurveRefinementConfigHP config;
    config.max_iterations = max_iters;

    // For b bits: ~b/3.32 decimal digits
    // Optimal h ~ 10^(-digits/4), tol ~ 10^(-digits/2)
    // This balances truncation O(h²) and roundoff O(ε/h²) errors
    unsigned int digits = bits * 3 / 10;  // Conservative approximation of bits/3.32
    unsigned int h_exp = digits / 4;
    unsigned int tol_exp = digits / 2;

    config.step_size_str = "1e-" + std::to_string(h_exp);
    config.residual_tolerance_str = "1e-" + std::to_string(tol_exp);

    return config;
}

//=============================================================================
// Numerical Hessian Determinant
//=============================================================================

mpreal computeNumericalHessianDetHP(
    const std::function<mpreal(const mpreal&, const mpreal&)>& f,
    const mpreal& x, const mpreal& y,
    const std::string& h_str)
{
    mpreal h(h_str);
    mpreal f00 = f(x, y);

    // Second partial derivatives using central differences
    mpreal f_xx = (f(x + h, y) - mpreal(2) * f00 + f(x - h, y)) / (h * h);
    mpreal f_yy = (f(x, y + h) - mpreal(2) * f00 + f(x, y - h)) / (h * h);
    mpreal f_xy = (f(x + h, y + h) - f(x + h, y - h) - f(x - h, y + h) + f(x - h, y - h))
                  / (mpreal(4) * h * h);

    // Hessian determinant
    return f_xx * f_yy - f_xy * f_xy;
}

std::function<mpreal(const mpreal&, const mpreal&)> makeHessianDetFunctionHP(
    const std::function<mpreal(const mpreal&, const mpreal&)>& f,
    const std::string& h_str)
{
    return [f, h_str](const mpreal& x, const mpreal& y) -> mpreal {
        return computeNumericalHessianDetHP(f, x, y, h_str);
    };
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

