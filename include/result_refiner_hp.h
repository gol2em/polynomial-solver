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
 * #include "result_refiner_hp.h"
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

#include "polynomial_hp.h"
#include "high_precision_types.h"
#include <vector>
#include <string>

namespace polynomial_solver {

/**
 * @brief Configuration for high-precision refinement
 */
struct RefinementConfigHP {
    std::string target_tolerance_str;    ///< Target error tolerance as string (default: "1e-50")
    std::string residual_tolerance_str;  ///< Residual tolerance as string (default: "1e-50")
    unsigned int max_newton_iters;       ///< Maximum Newton iterations (default: 100)
    unsigned int max_multiplicity;       ///< Maximum multiplicity to check (default: 10)

    RefinementConfigHP()
        : target_tolerance_str("1e-50"),
          residual_tolerance_str("1e-50"),
          max_newton_iters(100u),
          max_multiplicity(10u)
    {}
};

/**
 * @brief A refined root with high-precision results
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

    RefinedRootHP()
        : location(0), residual(0), multiplicity(1),
          first_nonzero_derivative(0), condition_estimate(1),
          iterations(0), converged(false)
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
     * @brief Estimate multiplicity of a root from derivatives
     *
     * Checks derivatives from order 1 to max_order to find the first
     * non-zero derivative. The order of the first non-zero derivative
     * is the multiplicity.
     *
     * @param location Point to check (high precision)
     * @param poly High-precision polynomial
     * @param max_order Maximum order to check (default: 10)
     * @param threshold Threshold for "non-zero" (default: 1e-50)
     * @param first_nonzero_deriv Output: value of first non-zero derivative
     * @return Estimated multiplicity
     */
    static unsigned int estimateMultiplicity(
        const mpreal& location,
        const PolynomialHP& poly,
        unsigned int max_order,
        const mpreal& threshold,
        mpreal& first_nonzero_deriv);

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
};

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // RESULT_REFINER_HP_H

