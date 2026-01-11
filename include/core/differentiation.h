#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H

/**
 * @file differentiation.h
 * @brief Differentiation utilities for Bernstein polynomials
 *
 * This module provides functionality for computing derivatives of polynomials
 * in Bernstein basis. The derivative of a Bernstein polynomial of degree n
 * with coefficients b_0, ..., b_n is a polynomial of degree n-1 with
 * coefficients: d_i = n * (b_{i+1} - b_i) for i = 0, ..., n-1.
 *
 * For multivariate tensor-product polynomials, partial derivatives are
 * computed dimension-by-dimension.
 */

#include "core/polynomial.h"
#include <map>
#include <vector>

namespace polynomial_solver {

/**
 * @class Differentiation
 * @brief Static utility class for computing derivatives of Bernstein polynomials.
 *
 * Provides stateless functions for one-off differentiation operations.
 * For repeated derivative computations, use DerivativeCache instead.
 */
class Differentiation {
public:
    /**
     * @brief Compute the derivative of a polynomial with respect to a variable.
     *
     * For univariate (1D): returns d^order f / dx^order
     * For multivariate: returns ∂^order f / ∂x_axis^order
     *
     * @param p The polynomial to differentiate
     * @param axis The variable index (0 <= axis < dimension)
     * @param order The derivative order (1, 2, 3, ...)
     * @return The derivative polynomial (degree reduced by 'order' along axis)
     */
    static Polynomial derivative(const Polynomial& p, std::size_t axis, unsigned int order = 1u);

    /**
     * @brief Compute the gradient of a polynomial.
     *
     * Returns a vector of first-order partial derivatives [∂f/∂x_0, ..., ∂f/∂x_{n-1}].
     *
     * @param p The polynomial
     * @return Vector of partial derivatives, one per dimension
     */
    static std::vector<Polynomial> gradient(const Polynomial& p);

    /**
     * @brief Compute the Hessian matrix of a polynomial.
     *
     * Returns a matrix H where H[i][j] = ∂²f/∂x_i∂x_j.
     * The Hessian is symmetric: H[i][j] = H[j][i].
     *
     * @param p The polynomial
     * @return 2D vector representing the Hessian matrix (dimension × dimension)
     */
    static std::vector<std::vector<Polynomial>> hessian(const Polynomial& p);

    /**
     * @brief Core differentiation function: compute first derivative along one axis.
     *
     * This is the fundamental operation. All higher-order and mixed derivatives
     * are computed by iteratively applying this function.
     *
     * @param p The polynomial to differentiate
     * @param axis The variable index
     * @return The first derivative along the specified axis
     */
    static Polynomial differentiateAxis(const Polynomial& p, std::size_t axis);

    /**
     * @brief Differentiate in power basis (more efficient for Newton methods).
     *
     * If the polynomial has power coefficients available, this method differentiates
     * directly in power basis without converting to Bernstein. The result has power
     * as primary representation.
     *
     * Power basis differentiation: d/dx(a_0 + a_1*x + a_2*x^2 + ...) = a_1 + 2*a_2*x + ...
     *
     * @param p The polynomial to differentiate
     * @param axis The variable index
     * @return The first derivative along the specified axis (power primary)
     */
    static Polynomial differentiateAxisPower(const Polynomial& p, std::size_t axis);
};

/**
 * @class DerivativeCache
 * @brief Caches computed derivatives of a polynomial for efficient reuse.
 *
 * This class stores a polynomial and all its computed derivatives in a map,
 * allowing efficient computation of higher-order and mixed partial derivatives
 * through iterative application of first-order derivatives.
 *
 * Derivatives are indexed by a multi-index (k_0, k_1, ..., k_{n-1}) where
 * k_i is the derivative order with respect to variable i.
 *
 * Example:
 *   - {0, 0, 0} = f (original polynomial)
 *   - {1, 0, 0} = ∂f/∂x_0
 *   - {0, 1, 0} = ∂f/∂x_1
 *   - {2, 0, 0} = ∂²f/∂x_0²
 *   - {1, 1, 0} = ∂²f/∂x_0∂x_1
 */
class DerivativeCache {
public:
    /**
     * @brief Construct a cache for the given polynomial.
     *
     * @param p The polynomial to cache derivatives for
     */
    explicit DerivativeCache(const Polynomial& p);

    /**
     * @brief Get a derivative by multi-index.
     *
     * The multi-index specifies the derivative order for each variable.
     * If the derivative hasn't been computed yet, it will be computed
     * and cached automatically.
     *
     * @param orders Multi-index of derivative orders (size = dimension)
     * @return Reference to the cached derivative polynomial
     */
    const Polynomial& get(const std::vector<unsigned int>& orders);

    /**
     * @brief Get a partial derivative with respect to a single variable.
     *
     * Convenience function for ∂^order f / ∂x_axis^order.
     *
     * @param axis The variable index
     * @param order The derivative order (default: 1)
     * @return Reference to the cached derivative polynomial
     */
    const Polynomial& getPartial(std::size_t axis, unsigned int order = 1u);

    /**
     * @brief Precompute all derivatives up to a given total order.
     *
     * This computes all derivatives where sum(orders) <= maxOrder.
     * Useful for precomputing all derivatives needed for multiplicity detection.
     *
     * @param maxOrder Maximum total derivative order
     */
    void precomputeUpToOrder(unsigned int maxOrder);

    /**
     * @brief Get the dimension of the cached polynomial.
     */
    std::size_t dimension() const;

private:
    /// The original polynomial (stored as derivative order {0, 0, ..., 0})
    std::size_t dimension_;

    /// Cache of computed derivatives, indexed by multi-index
    std::map<std::vector<unsigned int>, Polynomial> cache_;

    /**
     * @brief Compute a derivative if not already cached.
     *
     * Uses iterative differentiation: to compute ∂^k f / ∂x_i^k,
     * first get ∂^{k-1} f / ∂x_i^{k-1}, then differentiate once more.
     */
    const Polynomial& computeIfNeeded(const std::vector<unsigned int>& orders);
};

} // namespace polynomial_solver

#endif // DIFFERENTIATION_H

