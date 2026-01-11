#ifndef RESULT_REFINER_HP_H
#define RESULT_REFINER_HP_H

/**
 * @file result_refiner_hp.h
 * @brief High-precision result refinement for Tier 2 (fixed backend)
 *
 * This module provides high-precision Newton refinement for roots that are
 * ill-conditioned and cannot be accurately refined with double precision.
 *
 * ## Usage
 *
 * After solving with double precision and identifying roots with
 * needs_higher_precision=true, use this module to refine them:
 *
 * @code{.cpp}
 * #ifdef ENABLE_HIGH_PRECISION
 * #include "hp/result_refiner_hp.h"
 * using namespace polynomial_solver;
 *
 * // Set precision (256 bits = ~77 decimal digits)
 * PrecisionContext ctx(256);
 *
 * // Convert polynomial to high precision
 * PolynomialHP poly_hp(poly);
 *
 * // Refine root with high precision
 * RefinementConfigHP config_hp;
 * config_hp.target_tolerance_str = "1e-50";
 * 
 * RefinedRootHP result = ResultRefinerHP::refineRoot1D(
 *     initial_guess, poly_hp, config_hp);
 *
 * if (result.converged) {
 *     std::cout << "HP root: " << toString(result.location, 50) << std::endl;
 * }
 * #endif
 * @endcode
 *
 * Only available when ENABLE_HIGH_PRECISION is defined.
 */

#ifdef ENABLE_HIGH_PRECISION

#include "hp/polynomial_hp.h"
#include "hp/high_precision_types.h"
#include <vector>
#include <string>
#include <map>
#include <functional>

namespace polynomial_solver {

/**
 * @brief Multiplicity estimation method
 */
enum class MultiplicityMethod {
    NONE,           ///< No multiplicity detection (assume m=1)
    HINT,           ///< Use provided hint
    TAYLOR,         ///< Taylor series ratio test
    OSTROWSKI,      ///< Ostrowski's method (3 Newton iterates)
    TAYLOR_THEN_OSTROWSKI  ///< Try Taylor first, fall back to Ostrowski if needed
};

/**
 * @brief Iteration method for root refinement
 */
enum class IterationMethod {
    NEWTON,         ///< Standard Newton's method
    MODIFIED_NEWTON,///< Modified Newton for multiple roots: x_new = x - m*f/f'
    HALLEY,         ///< Halley's method (third-order)
    SCHRODER        ///< Schröder's method (third-order, good for multiple roots)
};

/**
 * @brief When to estimate multiplicity during iteration
 */
enum class MultiplicityTiming {
    ONCE_AT_START,  ///< Estimate once at the beginning
    EVERY_ITERATION,///< Re-estimate at every iteration
    WHEN_STAGNANT   ///< Re-estimate when convergence stagnates
};

/**
 * @brief Configuration for high-precision refinement
 */
struct RefinementConfigHP {
    std::string target_tolerance_str;    ///< Target error tolerance as string (default: "1e-50")
    std::string residual_tolerance_str;  ///< Residual tolerance as string (default: "1e-50")
    unsigned int max_newton_iters;       ///< Maximum Newton iterations (default: 100)
    unsigned int max_multiplicity;       ///< Maximum multiplicity to check (default: 10)
    double taylor_ratio_threshold;       ///< Ratio threshold for Taylor method (default: 10.0)
    unsigned int multiplicity_hint;      ///< Hint for multiplicity (0 = auto-detect)

    // Modular workflow control
    MultiplicityMethod multiplicity_method;  ///< How to estimate multiplicity
    IterationMethod iteration_method;        ///< Which iteration method to use
    MultiplicityTiming multiplicity_timing;  ///< When to estimate multiplicity

    RefinementConfigHP()
        : target_tolerance_str("1e-50"),
          residual_tolerance_str("1e-50"),
          max_newton_iters(100u),
          max_multiplicity(10u),
          taylor_ratio_threshold(10.0),
          multiplicity_hint(0u),
          multiplicity_method(MultiplicityMethod::TAYLOR),
          iteration_method(IterationMethod::MODIFIED_NEWTON),
          multiplicity_timing(MultiplicityTiming::ONCE_AT_START)
    {}
};

/**
 * @brief A refined root with high-precision results and guaranteed error bounds
 */
struct RefinedRootHP {
    mpreal location;                     ///< Root location (high precision)
    mpreal residual;                     ///< Residual f(x) at root
    unsigned int multiplicity;           ///< Estimated multiplicity (1 = simple root)
    mpreal first_nonzero_derivative;     ///< Value of first non-zero derivative
    mpreal condition_estimate;           ///< Estimated condition number
    unsigned int iterations;             ///< Newton iterations performed
    bool converged;                      ///< True if converged to target tolerance
    std::string error_message;           ///< Error message if not converged

    // Error bounds (interval containing the true root)
    mpreal max_error;                    ///< Maximum error bound (radius of interval)
    mpreal interval_lower;               ///< Lower bound of interval containing root
    mpreal interval_upper;               ///< Upper bound of interval containing root
    bool has_guaranteed_bounds;          ///< True if interval bounds are rigorous

    RefinedRootHP()
        : location(0), residual(0), multiplicity(1),
          first_nonzero_derivative(0), condition_estimate(1),
          iterations(0), converged(false),
          max_error(0), interval_lower(0), interval_upper(0),
          has_guaranteed_bounds(false)
    {}
};

/**
 * @class ResultRefinerHP
 * @brief High-precision root refinement using Newton's method
 *
 * This class provides high-precision versions of the refinement algorithms
 * from ResultRefiner. All arithmetic is performed in high precision (mpreal).
 *
 * This is the Tier 2 (non-template) version that works with PolynomialHP.
 */
class ResultRefinerHP {
public:
    /**
     * @brief Refine a 1D root to high precision using Newton's method
     *
     * Uses condition-aware convergence criterion that checks both residual
     * and estimated error. Handles multiple roots with modified Newton.
     *
     * @param initial_guess Initial guess (double precision is fine)
     * @param poly High-precision polynomial
     * @param config Refinement configuration
     * @return Refined root with convergence information
     */
    static RefinedRootHP refineRoot1D(
        double initial_guess,
        const PolynomialHP& poly,
        const RefinementConfigHP& config = RefinementConfigHP());

    /**
     * @brief Refine a 1D root using Schröder's method (third-order convergence)
     *
     * Schröder's method provides third-order convergence and is particularly
     * effective for multiple roots. It automatically handles multiplicity
     * without explicit detection.
     *
     * Formula: x_{n+1} = x_n - (f * f') / (f'^2 - f * f'')
     *
     * @param initial_guess Initial guess (double precision is fine)
     * @param poly High-precision polynomial
     * @param config Refinement configuration
     * @return Refined root with convergence information
     */
    static RefinedRootHP refineRoot1DSchroder(
        double initial_guess,
        const PolynomialHP& poly,
        const RefinementConfigHP& config = RefinementConfigHP());

    // ========================================================================
    // Modular Components for Flexible Workflows
    // ========================================================================

    /**
     * @brief Estimate multiplicity using the configured method
     *
     * @param x Current iterate
     * @param poly High-precision polynomial
     * @param config Refinement configuration
     * @param iterates Previous iterates (for Ostrowski method)
     * @return Estimated multiplicity (≥ 1)
     */
    static unsigned int estimateMultiplicityModular(
        const mpreal& x,
        const PolynomialHP& poly,
        const RefinementConfigHP& config,
        const std::vector<mpreal>& iterates = std::vector<mpreal>());

    /**
     * @brief Perform one iteration step using the configured method
     *
     * @param x Current iterate (will be updated)
     * @param poly High-precision polynomial
     * @param dpoly First derivative polynomial
     * @param ddpoly Second derivative polynomial (for Halley/Schröder)
     * @param multiplicity Current multiplicity estimate
     * @param method Iteration method to use
     * @return Step size taken
     */
    static mpreal performIterationStep(
        mpreal& x,
        const PolynomialHP& poly,
        const PolynomialHP& dpoly,
        const PolynomialHP& ddpoly,
        unsigned int multiplicity,
        IterationMethod method);

    // ========================================================================
    // Multiplicity Detection Methods
    // ========================================================================

    /**
     * @brief Estimate multiplicity using Ostrowski's method (1973)
     *
     * Uses 3 consecutive Newton iterates to estimate multiplicity:
     * p = ⌈1/2 + (x₁ - x₂)/(x₃ - 2x₂ + x₁)⌉
     *
     * This method works during convergence and doesn't require high-order derivatives.
     * Based on the asymptotic error ratio of Newton's method for multiple roots.
     *
     * @param x1 First Newton iterate
     * @param x2 Second Newton iterate
     * @param x3 Third Newton iterate
     * @return Estimated multiplicity (≥ 1)
     */
    static unsigned int estimateMultiplicityOstrowski(
        const mpreal& x1,
        const mpreal& x2,
        const mpreal& x3);

    /**
     * @brief Estimate multiplicity using Ostrowski's method from a starting point
     *
     * Performs 3 REGULAR Newton steps (not modified) from the starting point
     * and applies Ostrowski's formula to the resulting iterates.
     *
     * This is the correct way to use Ostrowski's method - it requires regular
     * Newton steps, not modified Newton steps.
     *
     * @param start Starting point (e.g., current iterate or interval midpoint)
     * @param poly High-precision polynomial
     * @return Estimated multiplicity (≥ 1)
     */
    static unsigned int estimateMultiplicityOstrowskiFromPoint(
        const mpreal& start,
        const PolynomialHP& poly);

    /**
     * @brief Estimate multiplicity of a root from derivatives (Taylor series method)
     *
     * Uses ratio test on consecutive derivatives to find first non-zero derivative.
     * For f(x) = c_m * (x-r)^m + ..., the ratio |f^(k+1)| / |f^(k)| reveals
     * whether f^(k) vanishes at the root.
     *
     * @param location Point to check (high precision)
     * @param poly High-precision polynomial
     * @param max_order Maximum order to check (default: 10)
     * @param threshold Threshold for "non-zero" (default: 1e-50)
     * @param first_nonzero_deriv Output: value of first non-zero derivative
     * @param ratio_threshold Ratio threshold for detecting vanishing derivatives (default: 10.0)
     * @return Estimated multiplicity
     */
    static unsigned int estimateMultiplicity(
        const mpreal& location,
        const PolynomialHP& poly,
        unsigned int max_order,
        const mpreal& threshold,
        mpreal& first_nonzero_deriv,
        double ratio_threshold = 10.0);

    /**
     * @brief Estimate multiplicity using simple threshold method
     *
     * Finds the first derivative |f^(k)| > threshold.
     * This is the simplest approach - just checks absolute values.
     *
     * @param location Point to check (high precision)
     * @param poly High-precision polynomial
     * @param max_order Maximum order to check
     * @param threshold Absolute threshold for detecting non-zero derivative
     * @return Estimated multiplicity (≥ 1)
     */
    static unsigned int estimateMultiplicitySimpleThreshold(
        const mpreal& location,
        const PolynomialHP& poly,
        unsigned int max_order,
        const mpreal& threshold);

    /**
     * @brief Estimate multiplicity using Sturm sequence
     *
     * Uses Sturm's theorem to count distinct roots in [location - radius, location + radius].
     * For a root of multiplicity m, f/gcd(f,f') has a simple root at that location.
     * We compute the Sturm sequence for f and f' to determine multiplicity.
     *
     * Key insight: If f has a root of multiplicity m at x=r, then:
     * - f has m roots at r (counted with multiplicity)
     * - f' has m-1 roots at r
     * - gcd(f, f') has m-1 roots at r
     * - f/gcd(f,f') has 1 simple root at r
     *
     * @param location Point to check (high precision)
     * @param poly High-precision polynomial
     * @param interval_radius Radius of interval to check
     * @return Estimated multiplicity (≥ 1)
     */
    static unsigned int estimateMultiplicitySturm(
        const mpreal& location,
        const PolynomialHP& poly,
        const mpreal& interval_radius);

    /**
     * @brief Comprehensive multiplicity detection using all methods
     *
     * Runs all available multiplicity detection methods and returns
     * detailed results for comparison and analysis.
     *
     * @param location Point to check (high precision)
     * @param poly High-precision polynomial
     * @param max_order Maximum multiplicity to check
     * @param x1 First Newton iterate (for Ostrowski, optional)
     * @param x2 Second Newton iterate (for Ostrowski, optional)
     * @param x3 Third Newton iterate (for Ostrowski, optional)
     * @return Map of method name to estimated multiplicity
     */
    static std::map<std::string, unsigned int> estimateMultiplicityAllMethods(
        const mpreal& location,
        const PolynomialHP& poly,
        unsigned int max_order,
        const mpreal& x1 = mpreal(0),
        const mpreal& x2 = mpreal(0),
        const mpreal& x3 = mpreal(0));

    /**
     * @brief Estimate condition number for root-finding problem
     *
     * Uses ratio of second derivative to first derivative squared.
     * Higher condition number indicates more ill-conditioned problem.
     *
     * @param location Point to evaluate (high precision)
     * @param poly High-precision polynomial
     * @param derivative_value First derivative value at location
     * @return Estimated condition number (>= 1.0)
     */
    static mpreal estimateConditionNumber1D(
        const mpreal& location,
        const PolynomialHP& poly,
        const mpreal& derivative_value);

    /**
     * @brief Compute rigorous error bounds for a refined root
     *
     * Uses interval Newton method or a posteriori error estimation to compute
     * a guaranteed interval containing the true root. For simple roots, uses
     * the Kantorovich theorem. For multiple roots, uses derivative-based bounds.
     *
     * @param location Refined root location
     * @param poly High-precision polynomial
     * @param multiplicity Estimated multiplicity
     * @param first_nonzero_deriv First non-zero derivative value
     * @param lower Output: lower bound of interval
     * @param upper Output: upper bound of interval
     * @return True if rigorous bounds were computed successfully
     */
    static bool computeErrorBounds(
        const mpreal& location,
        const PolynomialHP& poly,
        unsigned int multiplicity,
        const mpreal& first_nonzero_deriv,
        mpreal& lower,
        mpreal& upper);
};

//=============================================================================
// High-Precision Curve Refinement
//=============================================================================

/**
 * @brief Result of high-precision curve refinement
 */
struct CurveRefinedPointHP {
    mpreal x;                ///< Refined x coordinate
    mpreal y;                ///< Refined y coordinate
    mpreal residual;         ///< |g(x, y)| at refined point
    unsigned int iterations; ///< Number of iterations
    bool converged;          ///< True if converged within tolerance
};

/**
 * @brief Configuration for high-precision curve refinement
 */
struct CurveRefinementConfigHP {
    std::string residual_tolerance_str = "1e-50";  ///< Convergence criterion |g(x,y)| < tol
    unsigned int max_iterations = 100;              ///< Maximum Newton iterations
    std::string min_gradient_norm_str = "1e-100";   ///< Minimum |∇g| to avoid singularity
    std::string step_size_str = "1e-20";            ///< Step size h for numerical gradient

    /**
     * @brief Create config with optimal parameters for given precision bits
     *
     * Computes step size and tolerance to balance truncation O(h²) and roundoff O(ε/h²) errors.
     * For b bits (~b/3.32 decimal digits): h ≈ 10^(-digits/4), tol ≈ 10^(-digits/2)
     *
     * @param bits Precision bits (e.g., 128, 256)
     * @param max_iters Maximum iterations (default: 100)
     * @return Configured CurveRefinementConfigHP
     */
    static CurveRefinementConfigHP fromPrecisionBits(unsigned int bits, unsigned int max_iters = 100);
};

/**
 * @brief Refine a point onto a curve g(x,y) = 0 using high-precision arithmetic
 *
 * Uses gradient projection method with numerical gradient computed in high precision.
 * The smaller step size enabled by high precision gives much more accurate gradients.
 *
 * @param g Function g(x,y) using high-precision types
 * @param x0 Initial x coordinate (double precision, will be converted)
 * @param y0 Initial y coordinate (double precision, will be converted)
 * @param config Refinement configuration
 * @return Refinement result with converged point
 *
 * Usage:
 * @code
 * // Compute Hessian determinant with high precision
 * auto hessian_det_hp = [&f_hp](const mpreal& x, const mpreal& y) -> mpreal {
 *     mpreal h("1e-20");  // Much smaller step size in HP
 *     mpreal f00 = f_hp(x, y);
 *     mpreal f_xx = (f_hp(x+h, y) - 2*f00 + f_hp(x-h, y)) / (h*h);
 *     mpreal f_yy = (f_hp(x, y+h) - 2*f00 + f_hp(x, y-h)) / (h*h);
 *     mpreal f_xy = (f_hp(x+h, y+h) - f_hp(x+h, y-h) - f_hp(x-h, y+h) + f_hp(x-h, y-h)) / (4*h*h);
 *     return f_xx * f_yy - f_xy * f_xy;
 * };
 *
 * CurveRefinedPointHP result = refineCurveNumericalHP(hessian_det_hp, 0.5, 0.5);
 * @endcode
 */
CurveRefinedPointHP refineCurveNumericalHP(
    const std::function<mpreal(const mpreal&, const mpreal&)>& g,
    double x0, double y0,
    const CurveRefinementConfigHP& config = CurveRefinementConfigHP());

/**
 * @brief Refine multiple points onto a curve using high-precision arithmetic
 *
 * @param g Function g(x,y) using high-precision types
 * @param points Vector of (x, y) pairs to refine (double precision)
 * @param config Refinement configuration
 * @return Vector of high-precision refinement results
 */
std::vector<CurveRefinedPointHP> refineCurveNumericalHPMultiple(
    const std::function<mpreal(const mpreal&, const mpreal&)>& g,
    const std::vector<std::pair<double, double>>& points,
    const CurveRefinementConfigHP& config = CurveRefinementConfigHP());

/**
 * @brief Compute numerical Hessian determinant of a function at a point
 *
 * Uses central differences to compute second derivatives:
 *   f_xx = (f(x+h,y) - 2f(x,y) + f(x-h,y)) / h²
 *   f_yy = (f(x,y+h) - 2f(x,y) + f(x,y-h)) / h²
 *   f_xy = (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)) / (4h²)
 *   det(H) = f_xx * f_yy - f_xy²
 *
 * @param f Function f(x,y) to compute Hessian of
 * @param x X coordinate
 * @param y Y coordinate
 * @param h Step size (as string, e.g., "1e-20")
 * @return Hessian determinant det(H) at (x,y)
 */
mpreal computeNumericalHessianDetHP(
    const std::function<mpreal(const mpreal&, const mpreal&)>& f,
    const mpreal& x, const mpreal& y,
    const std::string& h_str);

/**
 * @brief Create a Hessian determinant function from a scalar function
 *
 * Returns a lambda that computes det(H_f) at any point using the given step size.
 * This is useful when you need to pass a curve function to refineCurveNumericalHP.
 *
 * @param f Function f(x,y) to compute Hessian of
 * @param h_str Step size for finite differences (e.g., "1e-20")
 * @return Function that computes det(H_f)(x,y)
 *
 * Usage:
 * @code
 * auto f_hp = [](const mpreal& x, const mpreal& y) { return exp(-(x*x + y*y)); };
 * auto hess_det = makeHessianDetFunctionHP(f_hp, "1e-20");
 * auto result = refineCurveNumericalHP(hess_det, x0, y0, config);
 * @endcode
 */
std::function<mpreal(const mpreal&, const mpreal&)> makeHessianDetFunctionHP(
    const std::function<mpreal(const mpreal&, const mpreal&)>& f,
    const std::string& h_str);

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // RESULT_REFINER_HP_H

