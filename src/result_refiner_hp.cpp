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

    // Newton iteration
    for (unsigned int iter = 0; iter < config.max_newton_iters; ++iter) {
        result.iterations = iter + 1;

        // Evaluate function and derivative
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);

        // Check if derivative is too small (multiple root or near-singular)
        // Do this BEFORE condition-aware check, because for multiple roots
        // the condition number formula breaks down (f' ≈ 0)
        mpreal df_threshold = mpreal("1e-50");
        if (abs(df) < df_threshold) {
            // Try to estimate multiplicity and use modified Newton
            mpreal first_nonzero_deriv = mpreal(0);
            mpreal mult_threshold = mpreal("1e-50");
            unsigned int mult = estimateMultiplicity(
                x, poly, config.max_multiplicity, mult_threshold, first_nonzero_deriv);

            result.multiplicity = mult;
            result.first_nonzero_derivative = first_nonzero_deriv;

            if (mult > 1) {
                // For multiple roots, check convergence based on residual only
                // Condition-aware check doesn't work because f' ≈ 0
                if (abs(f) < residual_tol) {
                    result.location = x;
                    result.residual = f;
                    result.condition_estimate = estimateConditionNumber1D(x, poly, first_nonzero_deriv);
                    result.converged = true;
                    return result;
                }

                // Use modified Newton for multiple root
                // For f(x) = (x-r)^m * g(x), we have f^(m)(x) ≈ m! * g(x)
                // So we can approximate: x_new = x - f(x) / (f^(m)(x) / m!)

                if (abs(first_nonzero_deriv) > df_threshold) {
                    // Compute factorial
                    mpreal factorial = mpreal(1);
                    for (unsigned int i = 2; i <= mult && i <= 10; ++i) {
                        factorial *= mpreal(i);
                    }

                    // Modified step using m-th derivative
                    mpreal step = f * factorial / first_nonzero_deriv;
                    x = x - step;
                    continue;
                }
            }

            // Can't make progress
            result.location = x;
            result.residual = f;
            result.converged = false;
            result.error_message = "Derivative too small, cannot make progress";
            return result;
        }

        // Check for convergence with condition-aware criterion (for simple roots)
        if (abs(f) < residual_tol) {
            // Residual is small - but is the error also small?
            // For ill-conditioned problems, small residual doesn't guarantee small error
            // Check condition number to estimate actual error

            // Compute second derivative for condition estimation
            PolynomialHP ddpoly = DifferentiationHP::derivative(dpoly, 0, 1);
            mpreal ddf = ddpoly.evaluate(x);

            // Estimate condition number: κ ≈ |f''| / |f'|²
            // This measures sensitivity of root to perturbations
            mpreal kappa = abs(ddf) / (abs(df) * abs(df) + mpreal("1e-100"));

            // Estimate actual error: |error| ≈ κ × |residual| / |f'|
            mpreal estimated_error = kappa * abs(f) / max(abs(df), mpreal("1e-50"));

            // Store results
            result.multiplicity = 1;  // Simple root (df is not small)
            result.first_nonzero_derivative = df;
            result.condition_estimate = kappa;

            // Accept root only if estimated error is within target tolerance
            if (estimated_error <= target_tol) {
                result.location = x;
                result.residual = f;
                result.converged = true;
                return result;
            }

            // Residual is small but estimated error is too large
            // This indicates an ill-conditioned problem
            // Continue iterating (may still improve with high precision)
        }

        // Standard Newton step
        mpreal x_new = x - f / df;
        x = x_new;
    }

    // Max iterations reached - check final residual
    mpreal f = poly.evaluate(x);
    result.location = x;
    result.residual = f;

    if (abs(f) < residual_tol) {
        // Check condition-aware convergence one more time
        mpreal df = dpoly.evaluate(x);
        PolynomialHP ddpoly = DifferentiationHP::derivative(dpoly, 0, 1);
        mpreal ddf = ddpoly.evaluate(x);

        mpreal kappa = abs(ddf) / (abs(df) * abs(df) + mpreal("1e-100"));
        mpreal estimated_error = kappa * abs(f) / max(abs(df), mpreal("1e-50"));

        result.multiplicity = 1;
        result.first_nonzero_derivative = df;
        result.condition_estimate = kappa;

        if (estimated_error <= target_tol) {
            result.converged = true;
            return result;
        }

        result.converged = false;
        result.error_message = "Residual small but estimated error too large (ill-conditioned)";
        return result;
    }

    result.converged = false;
    result.error_message = "Maximum iterations reached without convergence";
    return result;
}

unsigned int ResultRefinerHP::estimateMultiplicity(
    const mpreal& location,
    const PolynomialHP& poly,
    unsigned int max_order,
    const mpreal& threshold,
    mpreal& first_nonzero_deriv)
{
    first_nonzero_deriv = mpreal(0);

    // Check derivatives from order 1 to max_order
    for (unsigned int order = 1; order <= max_order; ++order) {
        PolynomialHP deriv = DifferentiationHP::derivative(poly, 0, order);
        mpreal deriv_val = deriv.evaluate(location);

        if (abs(deriv_val) > threshold) {
            // Found first non-zero derivative at order 'order'
            // This means multiplicity = order
            first_nonzero_deriv = deriv_val;
            return order;
        }
    }

    // All derivatives zero up to max_order
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

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

