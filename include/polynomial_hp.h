#ifndef POLYNOMIAL_HP_H
#define POLYNOMIAL_HP_H

/**
 * @file polynomial_hp.h
 * @brief High-precision polynomial class for Tier 2 (fixed backend)
 *
 * This module provides high-precision polynomial evaluation for root refinement.
 * It stores Bernstein coefficients in high-precision format and evaluates using
 * De Casteljau algorithm with high-precision arithmetic.
 *
 * Only available when ENABLE_HIGH_PRECISION is defined.
 */

#ifdef ENABLE_HIGH_PRECISION

#include "polynomial.h"
#include "high_precision_types.h"
#include <vector>
#include <cstddef>

namespace polynomial_solver {

/**
 * @class PolynomialHP
 * @brief High-precision polynomial for root refinement
 *
 * Stores Bernstein coefficients in high-precision format (mpreal).
 * Provides evaluation at high-precision points using De Casteljau algorithm.
 *
 * This class is designed for Tier 2 (fixed backend) - uses single mpreal type.
 * For Tier 3 (template-based), use Polynomial<T> instead.
 */
class PolynomialHP {
public:
    /**
     * @brief Default constructor creating an empty polynomial
     */
    PolynomialHP();

    /**
     * @brief Construct from high-precision Bernstein coefficients
     *
     * @param degrees  Degrees per variable (size = dimension)
     * @param bernstein_coeffs Coefficients in tensor-product Bernstein basis (HP)
     */
    PolynomialHP(const std::vector<unsigned int>& degrees,
                 const std::vector<mpreal>& bernstein_coeffs);

    /**
     * @brief Construct from double-precision polynomial
     *
     * Converts double-precision Bernstein coefficients to high-precision.
     * This is the primary way to create PolynomialHP from existing Polynomial.
     *
     * @param poly Double-precision polynomial to convert
     */
    explicit PolynomialHP(const Polynomial& poly);

    /**
     * @brief Destructor
     */
    ~PolynomialHP();

    /**
     * @brief Number of variables (polynomial dimension)
     */
    std::size_t dimension() const;

    /**
     * @brief Degrees per variable
     */
    const std::vector<unsigned int>& degrees() const;

    /**
     * @brief Total number of Bernstein coefficients
     */
    std::size_t coefficientCount() const;

    /**
     * @brief Access underlying high-precision Bernstein coefficients
     */
    const std::vector<mpreal>& bernsteinCoefficients() const;

    /**
     * @brief Evaluate polynomial at high-precision point using De Casteljau
     *
     * @param parameters Parameter values for each variable (size = dimension())
     * @return Value of polynomial at the given point (high-precision)
     */
    mpreal evaluate(const std::vector<mpreal>& parameters) const;

    /**
     * @brief Evaluate univariate polynomial at t using De Casteljau
     *
     * @param t Parameter value (must be 1D polynomial)
     * @return Value of polynomial at t (high-precision)
     */
    mpreal evaluate(const mpreal& t) const;

    /**
     * @brief Check if polynomial is empty
     */
    bool empty() const;

private:
    std::vector<unsigned int> degrees_;
    std::vector<mpreal> bernstein_coeffs_;
};

/**
 * @brief Convert double-precision polynomial to high-precision
 *
 * Convenience function for creating PolynomialHP from Polynomial.
 *
 * @param poly Double-precision polynomial
 * @return High-precision polynomial with converted coefficients
 */
PolynomialHP convertToHighPrecision(const Polynomial& poly);

/**
 * @brief Create PolynomialHP from power basis coefficients (high-precision)
 *
 * Converts power basis coefficients to Bernstein basis in high precision.
 * This avoids double precision limitations when creating test polynomials.
 *
 * @param degrees Degrees per variable (size = dimension)
 * @param power_coeffs Power basis coefficients in high precision
 * @return High-precision polynomial in Bernstein basis
 *
 * Example for 1D polynomial f(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ:
 * power_coeffs = [a₀, a₁, a₂, ..., aₙ]
 */
PolynomialHP fromPowerHP(const std::vector<unsigned int>& degrees,
                         const std::vector<mpreal>& power_coeffs);

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // POLYNOMIAL_HP_H

