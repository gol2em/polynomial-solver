#ifndef DIFFERENTIATION_HP_H
#define DIFFERENTIATION_HP_H

/**
 * @file differentiation_hp.h
 * @brief High-precision differentiation utilities for Bernstein polynomials (Tier 2)
 *
 * This module provides high-precision versions of differentiation operations
 * for use in Tier 2 (fixed high-precision, no templates).
 *
 * The derivative of a Bernstein polynomial of degree n with coefficients
 * b_0, ..., b_n is a polynomial of degree n-1 with coefficients:
 * d_i = n * (b_{i+1} - b_i) for i = 0, ..., n-1.
 *
 * All arithmetic is performed in high precision (mpreal).
 */

#ifdef ENABLE_HIGH_PRECISION

#include "hp/polynomial_hp.h"
#include "core/polynomial.h"
#include <vector>

namespace polynomial_solver {

/**
 * @class DifferentiationHP
 * @brief Static utility class for computing high-precision derivatives.
 *
 * This is the Tier 2 (non-template) version that works with PolynomialHP.
 * All coefficient arithmetic is performed in high precision.
 */
class DifferentiationHP {
public:
    /**
     * @brief Compute the derivative of a high-precision polynomial.
     *
     * For univariate (1D): returns d^order f / dx^order
     * For multivariate: returns ∂^order f / ∂x_axis^order
     *
     * @param p The high-precision polynomial to differentiate
     * @param axis The variable index (0 <= axis < dimension)
     * @param order The derivative order (1, 2, 3, ...)
     * @return The derivative polynomial (degree reduced by 'order' along axis)
     */
    static PolynomialHP derivative(const PolynomialHP& p, std::size_t axis, unsigned int order = 1u);

    /**
     * @brief Compute the gradient of a high-precision polynomial.
     *
     * Returns a vector of first-order partial derivatives [∂f/∂x_0, ..., ∂f/∂x_{n-1}].
     *
     * @param p The high-precision polynomial
     * @return Vector of partial derivatives, one per dimension
     */
    static std::vector<PolynomialHP> gradient(const PolynomialHP& p);

    /**
     * @brief Core differentiation function: compute first derivative along one axis.
     *
     * This is the fundamental operation. All higher-order derivatives
     * are computed by iteratively applying this function.
     *
     * Uses the Bernstein derivative formula with high-precision arithmetic:
     * d_i = n * (b_{i+1} - b_i)
     *
     * @param p The high-precision polynomial to differentiate
     * @param axis The variable index
     * @return The first derivative along the specified axis
     */
    static PolynomialHP differentiateAxis(const PolynomialHP& p, std::size_t axis);

    /**
     * @brief Power-basis differentiation along one axis.
     *
     * Uses the power basis derivative formula with high-precision arithmetic:
     * d/dx(a_i * x^i) = i * a_i * x^(i-1)
     *
     * This avoids Bernstein conversion and preserves exact rational coefficients.
     * Use this when working with polynomials created via fromPowerHP().
     *
     * @param p The high-precision polynomial to differentiate (must have valid power coefficients)
     * @param axis The variable index
     * @return The first derivative along the specified axis (in power basis)
     */
    static PolynomialHP differentiateAxisPower(const PolynomialHP& p, std::size_t axis);

    /**
     * @brief Compute derivative from double-precision polynomial, return HP result.
     *
     * Convenience function that converts input to HP, then differentiates.
     * Useful when you have a double-precision polynomial but need HP derivative.
     *
     * @param p The double-precision polynomial
     * @param axis The variable index
     * @param order The derivative order
     * @return High-precision derivative polynomial
     */
    static PolynomialHP derivativeFromDouble(const Polynomial& p, std::size_t axis, unsigned int order = 1u);
};

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // DIFFERENTIATION_HP_H

