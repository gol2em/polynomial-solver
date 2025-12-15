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
 * @brief High-precision polynomial with dual-representation support
 *
 * Stores both power and Bernstein coefficients in high-precision format (mpreal).
 * Follows the same dual-representation pattern as the double-precision Polynomial class.
 * Provides lazy conversion between representations with caching.
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
     * Converts double-precision coefficients to high-precision.
     * Uses high-precision conversion for power-to-Bernstein if needed.
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
     * @brief Total number of coefficients
     */
    std::size_t coefficientCount() const;

    /**
     * @brief Access high-precision Bernstein coefficients (lazy conversion)
     */
    const std::vector<mpreal>& bernsteinCoefficients() const;

    /**
     * @brief Access high-precision power coefficients (lazy conversion)
     */
    const std::vector<mpreal>& powerCoefficients() const;

    /**
     * @brief Get primary representation
     */
    PolynomialRepresentation primaryRepresentation() const;

    /**
     * @brief Check if power coefficients are available
     */
    bool hasPowerCoefficients() const;

    /**
     * @brief Check if Bernstein coefficients are available
     */
    bool hasBernsteinCoefficients() const;

    /**
     * @brief Ensure Bernstein is primary representation
     */
    void ensureBernsteinPrimary();

    /**
     * @brief Ensure power is primary representation
     */
    void ensurePowerPrimary();

    /**
     * @brief Explicitly convert power to Bernstein representation
     *
     * This is a public method to allow explicit conversion when needed.
     * Use this instead of relying on implicit conversion via bernsteinCoefficients().
     */
    void convertPowerToBernstein() const;

    /**
     * @brief Explicitly convert Bernstein to power representation
     *
     * This is a public method to allow explicit conversion when needed.
     * Use this instead of relying on implicit conversion via powerCoefficients().
     */
    void convertBernsteinToPower() const;

    /**
     * @brief Evaluate polynomial at high-precision point
     *
     * Uses Horner's method for power basis, De Casteljau for Bernstein.
     * Chooses based on primary representation for best accuracy.
     *
     * @param parameters Parameter values for each variable (size = dimension())
     * @return Value of polynomial at the given point (high-precision)
     */
    mpreal evaluate(const std::vector<mpreal>& parameters) const;

    /**
     * @brief Evaluate univariate polynomial at t
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

    /// Coefficients in tensor-product Bernstein basis (HP)
    mutable std::vector<mpreal> bernstein_coeffs_;

    /// Coefficients in tensor-product power basis (HP)
    mutable std::vector<mpreal> power_coeffs_;

    /// Which representation is primary (original/accurate)?
    PolynomialRepresentation primary_rep_;

    /// Is Bernstein representation currently valid/computed?
    mutable bool bernstein_valid_;

    /// Is power representation currently valid/computed?
    mutable bool power_valid_;

    // Friend function to allow fromPowerHP to access private members
    friend PolynomialHP fromPowerHP(const std::vector<unsigned int>&,
                                    const std::vector<mpreal>&);
};

/**
 * @brief Convert double-precision polynomial to high-precision
 *
 * Convenience function for creating PolynomialHP from Polynomial.
 * Preserves the primary representation from the double-precision polynomial.
 *
 * @param poly Double-precision polynomial
 * @return High-precision polynomial with converted coefficients
 */
PolynomialHP convertToHighPrecision(const Polynomial& poly);

/**
 * @brief Create PolynomialHP from power basis coefficients (high-precision)
 *
 * Stores power coefficients as primary, Bernstein computed lazily.
 * This avoids double precision limitations when creating test polynomials.
 *
 * @param degrees Degrees per variable (size = dimension)
 * @param power_coeffs Power basis coefficients in high precision
 * @return High-precision polynomial with power as primary representation
 *
 * Example for 1D polynomial f(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ:
 * power_coeffs = [a₀, a₁, a₂, ..., aₙ]
 */
PolynomialHP fromPowerHP(const std::vector<unsigned int>& degrees,
                         const std::vector<mpreal>& power_coeffs);

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // POLYNOMIAL_HP_H

