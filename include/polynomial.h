#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

/**
 * @file polynomial.h
 * @brief Polynomial class for representing and manipulating polynomials
 *
 * This module provides functionality for polynomial operations including
 * evaluation, arithmetic operations, and root finding.
 */

#include <cstddef>
#include <vector>

namespace polynomial_solver {

/**
 * @class Polynomial
 * @brief Represents a multivariate polynomial with real coefficients
 *        stored in Bernstein form.
 *
 * Each polynomial is defined over a fixed number of variables
 * (the dimension). For each variable i, a non-negative degree d_i is
 * stored. The coefficients are stored in tensor-product Bernstein basis
 * order as a flat array.
 *
 * A polynomial can be constructed either from power basis coefficients
 * or directly from Bernstein coefficients. Internally, coefficients
 * are always stored in Bernstein form.
 */
class Polynomial {
public:
    /**
     * @brief Default constructor creating an empty polynomial.
     */
    Polynomial();

    /**
     * @brief Construct from Bernstein coefficients.
     *
     * @param degrees  Degrees per variable (size = dimension).
     * @param bernstein_coeffs Coefficients in tensor-product Bernstein basis.
     *
     * The number of coefficients should be
     *   prod_i (degrees[i] + 1).
     */
    Polynomial(const std::vector<unsigned int>& degrees,
               const std::vector<double>& bernstein_coeffs);

    /**
     * @brief Destructor.
     */
    ~Polynomial();

    /**
     * @brief Factory: construct from Bernstein coefficients.
     */
    static Polynomial fromBernstein(const std::vector<unsigned int>& degrees,
                                    const std::vector<double>& bernstein_coeffs);

    /**
     * @brief Factory: construct from power basis coefficients.
     *
     * The input coefficients are given in tensor-product power basis.
     * Internally they will be converted to Bernstein basis.
     */
    static Polynomial fromPower(const std::vector<unsigned int>& degrees,
                                const std::vector<double>& power_coeffs);

    /**
     * @brief Number of variables (polynomial dimension).
     */
    std::size_t dimension() const;

    /**
     * @brief Degrees per variable.
     */
    const std::vector<unsigned int>& degrees() const;

    /**
     * @brief Total number of Bernstein coefficients.
     */
    std::size_t coefficientCount() const;

    /**
     * @brief Access underlying Bernstein coefficients.
     */
    const std::vector<double>& bernsteinCoefficients() const;

    /**
     * @brief Evaluate the polynomial at a parameter point using De Casteljau.
     *
     * @param parameters Parameter values for each variable (size = dimension()).
     * @return Value of the polynomial at the given point.
     */
    double evaluate(const std::vector<double>& parameters) const;

    /**
     * @brief Evaluate a univariate polynomial at t using De Casteljau.
     *
     * This is a convenience overload for the case dimension() == 1.
     */
    double evaluate(double t) const;

    /**
     * @brief Build control points for the graph hypersurface (t, f(t)).
     *
     * The polynomial is assumed to be defined on the normalized domain
     * [0,1]^dimension(). For each Bernstein coefficient with multi-index
     * (i_0, ..., i_{n-1}), the corresponding control point is
     *
     *   ( i_0 / d_0, ..., i_{n-1} / d_{n-1}, c_{i_0,...,i_{n-1}} ),
     *
     * where d_k = degrees()[k] and c_* is the Bernstein coefficient.
     *
     * The output is a flat array of length coefficientCount() * (dimension() + 1).
     * For each coefficient index q in [0, coefficientCount()), the block
     *   control_points[q * (dimension()+1) + j]   (0 <= j < dimension())
     * stores the j-th coordinate in [0,1], and
     *   control_points[q * (dimension()+1) + dimension()]
     * stores the function value.
     */
    void graphControlPoints(std::vector<double>& control_points) const;

    /**
     * @brief Compute graph corners for geometric operations.
     *
     * For a polynomial f: [0,1]^n -> R, computes the corners of the graph
     * in R^{n+1}. The graph is the set {(x, f(x)) : x in [0,1]^n}.
     *
     * This evaluates f at all 2^n corners of the unit hypercube and returns
     * points in R^{n+1} of the form (x_1, ..., x_n, f(x_1,...,x_n)).
     *
     * This is different from graphControlPoints(), which returns only the
     * Bernstein control points. For polynomials with degree 0 in some dimensions,
     * the Bernstein control points don't span the full domain, so we need to
     * evaluate at all corners for correct geometric operations.
     *
     * @param corners Output vector of corner points, each of size dimension()+1
     */
    void graphCorners(std::vector<std::vector<double>>& corners) const;

    /**
     * @brief Restrict this polynomial to [a,b] along the given axis.
     *
     * For the specified axis (0 <= axis < dimension()), this constructs a new
     * polynomial q with the same degrees and dimension such that
     *
     *   q(u_0, ..., u_axis, ..., u_{n-1}) =
     *       p(u_0, ..., a + (b - a) * u_axis, ..., u_{n-1}),
     *
     * i.e. the original polynomial restricted to the parameter interval [a,b]
     * along that axis and reparameterized back to [0,1]. If the inputs do not
     * satisfy 0 <= a < b <= 1, the original polynomial is returned unchanged.
     */
    Polynomial restrictedToInterval(std::size_t axis, double a, double b) const;

private:
    /// Number of variables.
    std::size_t dimension_;

    /// Degree per variable.
    std::vector<unsigned int> degrees_;

    /// Coefficients in tensor-product Bernstein basis.
    std::vector<double> bernstein_coeffs_;

    /**
     * @brief Compute flattened coefficient index from a multi-index.
     *
     * The multi-index has one entry per variable, each in
     * [0, degrees_[i]].
     */
    std::size_t flattenIndex(const std::vector<unsigned int>& multi_index) const;
};

} // namespace polynomial_solver

#endif // POLYNOMIAL_H
