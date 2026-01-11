#ifndef HIGH_PRECISION_REFINER_H
#define HIGH_PRECISION_REFINER_H

/**
 * @file high_precision_refiner.h
 * @brief High-precision refinement for ill-conditioned roots
 *
 * This module provides high-precision arithmetic support for refining roots
 * that are ill-conditioned and cannot be accurately computed with double precision.
 *
 * ## Usage
 *
 * 1. Build with ENABLE_HIGH_PRECISION=ON
 * 2. Solve with double precision to locate roots
 * 3. For roots with needs_higher_precision=true, use this module
 *
 * ## Example
 *
 * @code{.cpp}
 * #ifdef ENABLE_HIGH_PRECISION
 * #include "hp/high_precision_refiner.h"
 * using namespace polynomial_solver::high_precision;
 *
 * // After double-precision solving and refinement
 * for (const auto& root : refined.roots) {
 *     if (root.needs_higher_precision) {
 *         // Convert to high precision (256 bits = ~77 decimal digits)
 *         auto poly_hp = convertToHighPrecision(poly, 256);
 *         
 *         // Refine with high precision
 *         HighPrecisionConfig hp_config;
 *         hp_config.precision_bits = 256;
 *         hp_config.target_tolerance = 1e-50;
 *         
 *         auto root_hp = refineRootHighPrecision(
 *             root.location, {poly_hp}, hp_config);
 *         
 *         std::cout << "High-precision root: " << root_hp.location << std::endl;
 *     }
 * }
 * #endif
 * @endcode
 *
 * ## Requirements
 *
 * - Boost.Multiprecision library (header-only)
 * - Or: MPFR library (requires linking)
 *
 * ## Precision Levels
 *
 * - 128 bits: ~38 decimal digits
 * - 256 bits: ~77 decimal digits (recommended)
 * - 512 bits: ~154 decimal digits
 * - 1024 bits: ~308 decimal digits
 */

#ifdef ENABLE_HIGH_PRECISION

#include <vector>
#include <string>
#include "core/polynomial.h"

// Use Boost.Multiprecision with MPFR backend
// This is header-only if Boost is configured that way
#include <boost/multiprecision/mpfr.hpp>

namespace polynomial_solver {
namespace high_precision {

// Type alias for high-precision floating point
using mpreal = boost::multiprecision::mpfr_float;

/**
 * @brief Configuration for high-precision refinement
 */
struct HighPrecisionConfig {
    unsigned int precision_bits;      ///< Target precision in bits (default: 256 = ~77 digits)
    unsigned int max_newton_iters;    ///< Maximum Newton iterations (default: 100)
    std::string target_tolerance;     ///< Target relative error as string (default: "1e-50")
    std::string residual_tolerance;   ///< Target residual as string (default: "1e-50")
    
    HighPrecisionConfig()
        : precision_bits(256),
          max_newton_iters(100),
          target_tolerance("1e-50"),
          residual_tolerance("1e-50")
    {}
};

/**
 * @brief A refined root with high precision
 */
struct RefinedRootHP {
    std::string location;             ///< Root as decimal string (arbitrary precision)
    std::string residual;             ///< Residual as decimal string
    unsigned int precision_bits;      ///< Actual precision used
    unsigned int decimal_digits;      ///< Approximate decimal digits of precision
    unsigned int iterations;          ///< Newton iterations performed
    bool verified;                    ///< True if converged to target tolerance
    std::string error_message;        ///< Error message if not verified
};

/**
 * @brief Core high-precision operations
 *
 * This class provides high-precision versions of the core polynomial operations:
 * - Coefficient conversion (double → mpreal)
 * - Evaluation (De Casteljau algorithm)
 * - Differentiation (Bernstein derivative formula)
 */
class HighPrecisionOps {
public:
    /**
     * @brief Convert double-precision coefficients to high-precision
     *
     * @param coeffs Double-precision Bernstein coefficients
     * @param precision_bits Target precision in bits
     * @return High-precision coefficients
     */
    static std::vector<mpreal> convertCoefficients(
        const std::vector<double>& coeffs,
        unsigned int precision_bits = 256);

    /**
     * @brief Evaluate univariate Bernstein polynomial using De Casteljau
     *
     * High-precision version of DeCasteljau::evaluate1D
     *
     * @param coeffs High-precision Bernstein coefficients
     * @param t Parameter value
     * @return Function value at t
     */
    static mpreal evaluate1D(
        const std::vector<mpreal>& coeffs,
        const mpreal& t);

    /**
     * @brief Evaluate multivariate Bernstein polynomial using De Casteljau
     *
     * High-precision version of DeCasteljau::evaluateTensorProduct
     *
     * @param coeffs High-precision Bernstein coefficients (tensor-product layout)
     * @param degrees Polynomial degrees per variable
     * @param point Parameter values for each variable
     * @return Function value at point
     */
    static mpreal evaluateTensorProduct(
        const std::vector<mpreal>& coeffs,
        const std::vector<unsigned int>& degrees,
        const std::vector<mpreal>& point);

    /**
     * @brief Differentiate Bernstein polynomial along one axis
     *
     * High-precision version of Differentiation::differentiateAxis
     * Uses Bernstein derivative formula: d_i = n * (b_{i+1} - b_i)
     *
     * @param coeffs High-precision Bernstein coefficients
     * @param degrees Polynomial degrees per variable
     * @param axis Variable to differentiate with respect to
     * @return Coefficients of derivative polynomial (degree reduced by 1 along axis)
     */
    static std::vector<mpreal> differentiate(
        const std::vector<mpreal>& coeffs,
        const std::vector<unsigned int>& degrees,
        std::size_t axis);

    /**
     * @brief Compute gradient (all first partial derivatives)
     *
     * @param coeffs High-precision Bernstein coefficients
     * @param degrees Polynomial degrees per variable
     * @return Vector of derivative coefficients, one per dimension
     */
    static std::vector<std::vector<mpreal>> gradient(
        const std::vector<mpreal>& coeffs,
        const std::vector<unsigned int>& degrees);
};

/**
 * @brief High-precision refinement using Newton's method
 *
 * This class provides high-precision refinement for roots that are
 * ill-conditioned and cannot be accurately computed with double precision.
 */
class HighPrecisionRefiner {
public:
    /**
     * @brief Refine a univariate root to high precision
     *
     * Takes a double-precision polynomial and initial guess, converts to
     * high precision, and refines using Newton's method.
     *
     * @param initial_guess Initial guess from double-precision solver
     * @param poly Double-precision polynomial
     * @param config High-precision configuration
     * @return Refined root with high precision
     */
    static RefinedRootHP refineRoot1D(
        double initial_guess,
        const Polynomial& poly,
        const HighPrecisionConfig& config);

    /**
     * @brief Refine a multivariate root to high precision
     *
     * Takes a double-precision polynomial system and initial guess, converts to
     * high precision, and refines using Newton's method.
     *
     * @param initial_guess Initial guess from double-precision solver
     * @param system Double-precision polynomial system
     * @param config High-precision configuration
     * @return Refined root with high precision
     */
    static RefinedRootHP refineRootMultivariate(
        const std::vector<double>& initial_guess,
        const PolynomialSystem& system,
        const HighPrecisionConfig& config);

private:
    /**
     * @brief Newton iteration for univariate case (internal)
     *
     * @param x Initial guess (high precision)
     * @param coeffs Polynomial coefficients (high precision)
     * @param degrees Polynomial degrees
     * @param config Configuration
     * @param result Output result
     * @return True if converged
     */
    static bool newtonIteration1D(
        mpreal& x,
        const std::vector<mpreal>& coeffs,
        const std::vector<unsigned int>& degrees,
        const HighPrecisionConfig& config,
        RefinedRootHP& result);

    /**
     * @brief Newton iteration for multivariate case (internal)
     *
     * Uses Jacobian matrix and linear solve for Newton step
     *
     * @param x Initial guess (high precision)
     * @param system_coeffs Polynomial coefficients for each equation
     * @param system_degrees Polynomial degrees for each equation
     * @param config Configuration
     * @param result Output result
     * @return True if converged
     */
    static bool newtonIterationMultivariate(
        std::vector<mpreal>& x,
        const std::vector<std::vector<mpreal>>& system_coeffs,
        const std::vector<std::vector<unsigned int>>& system_degrees,
        const HighPrecisionConfig& config,
        RefinedRootHP& result);
};

/**
 * @brief Helper: Convert mpreal to decimal string with specified digits
 */
std::string mprealToString(const mpreal& value, unsigned int digits = 50);

/**
 * @brief Helper: Convert decimal string to mpreal
 */
mpreal stringToMpreal(const std::string& str, unsigned int precision_bits = 256);

} // namespace high_precision
} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // HIGH_PRECISION_REFINER_H

