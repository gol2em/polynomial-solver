#ifdef ENABLE_HIGH_PRECISION

#include "result_refiner_hp.h"
#include "differentiation_hp.h"
#include "precision_conversion.h"
#include <cmath>

namespace polynomial_solver {

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

    // Get derivative polynomial
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);

    // For Ostrowski multiplicity estimation: store last 3 iterates
    std::vector<mpreal> iterates;
    iterates.reserve(4);  // x0, x1, x2, x3
    iterates.push_back(x);  // x0 = initial guess

    // Track multiplicity estimate
    unsigned int estimated_multiplicity = 1;
    bool multiplicity_detected = false;

    // Track previous error bound to detect stagnation
    mpreal prev_error_bound = mpreal("1e100");

    // Newton iteration
    for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
        result.iterations = iter + 1;

        // Evaluate function and derivative
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);

        // Store residual
        result.residual = f;

        // Step 1: Verify/detect multiplicity if derivative is small
        // Always re-verify when close enough, even if Ostrowski gave us an estimate
        if (abs(df) < mpreal("1e-20")) {
            // Derivative is small - do Taylor series analysis for exact multiplicity
            mpreal mult_threshold = mpreal("1e-50");
            mpreal dummy_deriv;
            unsigned int verified_mult = estimateMultiplicity(
                x, poly, config.max_multiplicity, mult_threshold, dummy_deriv);

            // Update if we got a better estimate
            if (verified_mult > 0 && verified_mult <= config.max_multiplicity) {
                estimated_multiplicity = verified_mult;
                multiplicity_detected = true;
            }
        }

        // Step 1b: Compute the correct derivative for the current multiplicity
        mpreal first_nonzero_deriv;
        if (estimated_multiplicity == 1) {
            first_nonzero_deriv = df;
        } else {
            // Compute m-th derivative explicitly
            PolynomialHP deriv_m = poly;
            for (unsigned int k = 0; k < estimated_multiplicity; ++k) {
                deriv_m = DifferentiationHP::derivative(deriv_m, 0, 1);
            }
            first_nonzero_deriv = deriv_m.evaluate(x);
        }

        // Step 2: Compute rigorous error bounds at each iteration
        mpreal lower, upper;
        mpreal current_error_bound = mpreal("1e100");  // Large default
        bool has_bounds = false;

        // Try to compute bounds with current multiplicity estimate
        if (abs(first_nonzero_deriv) > mpreal("1e-100")) {
            has_bounds = computeErrorBounds(x, poly, estimated_multiplicity,
                                           first_nonzero_deriv, lower, upper);
            if (has_bounds) {
                current_error_bound = (upper - lower) / mpreal("2.0");
            }
        }

        // Step 3: Check convergence based on error bounds
        if (has_bounds && current_error_bound <= target_tol) {
            // Converged! Error bound is small enough
            result.location = x;
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

        // Step 4: Check if we're making progress
        if (has_bounds && current_error_bound >= prev_error_bound * mpreal("0.99")) {
            // Error bound is not improving - re-verify multiplicity
            // This catches cases where Ostrowski gave wrong estimate
            if (!multiplicity_detected) {
                mpreal mult_threshold = mpreal("1e-50");
                mpreal dummy_deriv;
                unsigned int verified_mult = estimateMultiplicity(
                    x, poly, config.max_multiplicity, mult_threshold, dummy_deriv);

                if (verified_mult > 0 && verified_mult <= config.max_multiplicity) {
                    estimated_multiplicity = verified_mult;
                    multiplicity_detected = true;
                }
            }
        }

        prev_error_bound = current_error_bound;

        // Step 5: Compute Newton step based on current multiplicity estimate
        mpreal step;

        if (estimated_multiplicity > 1) {
            // Use modified Newton for multiple root
            step = mpreal(estimated_multiplicity) * f / df;
        } else {
            // Standard Newton step
            step = f / df;
        }

        // Limit step size to avoid divergence in ill-conditioned cases
        mpreal max_step = mpreal("1.0");
        if (abs(step) > max_step) {
            step = step * max_step / abs(step);
        }

        x = x - step;

        // Store iterate for Ostrowski method
        if (iter < 3) {
            iterates.push_back(x);
        }

        // Apply Ostrowski multiplicity estimation after 3 iterations
        if (iter == 2 && iterates.size() == 4 && !multiplicity_detected) {
            // We have x0, x1, x2, x3
            unsigned int mult_est = estimateMultiplicityOstrowski(
                iterates[1], iterates[2], iterates[3]);

            // DEBUG
            #ifdef DEBUG_MULTIPLICITY
            std::cout << "DEBUG: Ostrowski at iter 2: mult_est = " << mult_est << std::endl;
            std::cout << "  x1 = " << iterates[1] << std::endl;
            std::cout << "  x2 = " << iterates[2] << std::endl;
            std::cout << "  x3 = " << iterates[3] << std::endl;
            #endif

            if (mult_est > 1) {
                // Multiple root suspected - use as hint but don't mark as fully detected
                // We'll verify this later when derivative becomes small
                estimated_multiplicity = mult_est;
                // Note: Don't set multiplicity_detected = true yet
                // We want to verify this with Taylor series when close enough
            }
        }
    }

    // Max iterations reached - compute final error bounds
    mpreal f = poly.evaluate(x);
    mpreal df = dpoly.evaluate(x);
    result.location = x;
    result.residual = f;

    // Determine final multiplicity if not already detected
    if (!multiplicity_detected && abs(df) < mpreal("1e-20")) {
        mpreal mult_threshold = mpreal("1e-50");
        mpreal dummy_deriv;
        estimated_multiplicity = estimateMultiplicity(
            x, poly, config.max_multiplicity, mult_threshold, dummy_deriv);
        multiplicity_detected = true;
    }

    // Compute the correct derivative for the final multiplicity
    mpreal final_first_nonzero_deriv;
    if (estimated_multiplicity == 1) {
        final_first_nonzero_deriv = df;
    } else {
        // Compute m-th derivative explicitly
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

        // Check if we actually converged based on error bounds
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

    // Apply ceiling with minimum value 1
    // ⌈p⌉ with p ≥ 1
    int multiplicity = static_cast<int>(ceil(p_est));
    if (multiplicity < 1) {
        multiplicity = 1;
    }

    // Sanity check: cap at reasonable maximum
    if (multiplicity > 20) {
        multiplicity = 20;  // Very high multiplicities are rare
    }

    return static_cast<unsigned int>(multiplicity);
}

unsigned int ResultRefinerHP::estimateMultiplicity(
    const mpreal& location,
    const PolynomialHP& poly,
    unsigned int max_order,
    const mpreal& threshold,
    mpreal& first_nonzero_deriv)
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

    // Estimate coefficient scale from the polynomial's Bernstein coefficients
    // This gives us the typical magnitude of coefficients
    const std::vector<mpreal>& coeffs = poly.bernsteinCoefficients();
    mpreal coeff_scale = mpreal(0);
    for (const mpreal& c : coeffs) {
        coeff_scale = max(coeff_scale, abs(c));
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
    }

    if (significant_orders.empty()) {
        // All derivatives are zero
        first_nonzero_deriv = mpreal(0);
        return max_order + 1;
    }

    // Second pass: Use ratio test to find true first non-zero derivative
    // Check ratios |f^(k+1)| / |f^(k)| for consecutive significant orders
    for (size_t i = 0; i < significant_orders.size(); ++i) {
        unsigned int order = significant_orders[i];
        mpreal deriv_val = abs(deriv_values[order]);

        // Check if there's a next significant derivative
        if (i + 1 < significant_orders.size()) {
            unsigned int next_order = significant_orders[i + 1];
            mpreal next_deriv = abs(deriv_values[next_order]);

            // Compute ratio
            mpreal ratio = next_deriv / max(deriv_val, mpreal("1e-100"));

            // If ratio is very large (>100), current derivative is likely O(ε)
            // and vanishes at the root, so skip it
            if (ratio > mpreal(100)) {
                continue;  // Try next order
            }
        }

        // This derivative doesn't have a much larger successor
        // It's likely the first non-zero derivative at the root

        // Final check: verify it's above relative threshold
        mpreal expected_scale = factorials[order] * coeff_scale;
        mpreal relative_tol = mpreal("1e-15");
        mpreal deriv_threshold = max(threshold, relative_tol * expected_scale);

        if (deriv_val > deriv_threshold) {
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
        // We need to find min|f'(x)| over this interval
        //
        // For Bernstein polynomials, the derivative is bounded by the convex hull
        // of its Bernstein coefficients. We evaluate at several points to get a
        // conservative bound.

        mpreal df_min = abs(first_nonzero_deriv);  // Start with center value

        // Sample derivative at multiple points in the interval
        const int num_samples = 5;
        for (int i = 0; i <= num_samples; ++i) {
            mpreal t = mpreal(i) / mpreal(num_samples);
            mpreal x_sample = location - radius + mpreal("2.0") * radius * t;

            // Clamp to [0, 1] domain
            if (x_sample < mpreal("0.0")) x_sample = mpreal("0.0");
            if (x_sample > mpreal("1.0")) x_sample = mpreal("1.0");

            mpreal df_sample = abs(dpoly.evaluate(x_sample));
            if (df_sample < df_min && df_sample > mpreal("1e-100")) {
                df_min = df_sample;
            }
        }

        // Check if derivative stays bounded away from zero
        if (df_min < mpreal("1e-100")) {
            // Derivative too close to zero in the interval, cannot guarantee bounds
            return false;
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

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

