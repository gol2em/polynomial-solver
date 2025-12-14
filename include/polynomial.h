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
 * @enum PolynomialRepresentation
 * @brief Indicates which coefficient representation is primary (original/accurate)
 */
enum class PolynomialRepresentation {
    POWER,      ///< Power basis is primary (from user input or Newton methods)
    BERNSTEIN   ///< Bernstein basis is primary (from CAD or subdivision solver)
};

/**
 * @class Polynomial
 * @brief Represents a multivariate polynomial with dual coefficient storage.
 *
 * This class stores polynomials in BOTH power and Bernstein basis representations,
 * with lazy conversion between them. One representation is marked as "primary"
 * (the original, accurate one), and the other is computed on-demand when needed.
 *
 * Design rationale:
 * - CAD systems provide Bernstein coefficients (exact) → use Bernstein as primary
 * - User/test inputs provide power coefficients (exact) → use power as primary
 * - Subdivision solver REQUIRES Bernstein basis (for PP method)
 * - Newton refinement is MORE EFFICIENT with power basis (Horner's method)
 * - Differentiation is SIMPLER in power basis
 *
 * The dual storage allows each algorithm to use its preferred representation
 * without precision loss from unnecessary conversions.
 *
 * Memory overhead: Only stores primary representation initially. Secondary
 * representation is computed lazily on first access and cached.
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
     * @brief Access underlying Bernstein coefficients (with lazy conversion).
     *
     * If Bernstein coefficients are not yet computed, converts from power basis
     * using high-precision intermediate step (64-bit) for better accuracy.
     * The converted coefficients are cached for future access.
     *
     * @return Reference to Bernstein coefficients (always valid)
     */
    const std::vector<double>& bernsteinCoefficients() const;

    /**
     * @brief Access underlying power coefficients (with lazy conversion).
     *
     * If power coefficients are not yet computed, converts from Bernstein basis.
     * The converted coefficients are cached for future access.
     *
     * Note: Bernstein→Power conversion may lose some precision compared to
     * the original Bernstein coefficients.
     *
     * @return Reference to power coefficients (always valid)
     */
    const std::vector<double>& powerCoefficients() const;

    /**
     * @brief Query which representation is primary (original/accurate).
     *
     * @return POWER if created from power basis, BERNSTEIN if from Bernstein basis
     */
    PolynomialRepresentation primaryRepresentation() const;

    /**
     * @brief Check if power coefficients are currently available (computed).
     *
     * @return true if power coefficients are cached, false if they need conversion
     */
    bool hasPowerCoefficients() const;

    /**
     * @brief Check if Bernstein coefficients are currently available (computed).
     *
     * @return true if Bernstein coefficients are cached, false if they need conversion
     */
    bool hasBernsteinCoefficients() const;

    /**
     * @brief Switch primary representation to Bernstein.
     *
     * This method is called by the solver before using PP method (which requires Bernstein).
     * After this call, Bernstein becomes the authoritative representation.
     * If power was primary, it's converted to Bernstein and power becomes secondary.
     * Power coefficients remain cached and valid unless invalidated by operations.
     */
    void ensureBernsteinPrimary();

    /**
     * @brief Switch primary representation to Power.
     *
     * This method can be called by the refiner to use power basis for Newton iteration.
     * After this call, Power becomes the authoritative representation.
     * If Bernstein was primary, it's converted to power and Bernstein becomes secondary.
     * Bernstein coefficients remain cached and valid unless invalidated by operations.
     */
    void ensurePowerPrimary();

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
    mutable std::vector<double> bernstein_coeffs_;

    /// Coefficients in tensor-product power basis.
    mutable std::vector<double> power_coeffs_;

    /// Which representation is primary (original/accurate)?
    PolynomialRepresentation primary_rep_;

    /// Is Bernstein representation currently valid/computed?
    mutable bool bernstein_valid_;

    /// Is power representation currently valid/computed?
    mutable bool power_valid_;

    /**
     * @brief Compute flattened coefficient index from a multi-index.
     *
     * The multi-index has one entry per variable, each in
     * [0, degrees_[i]].
     */
    std::size_t flattenIndex(const std::vector<unsigned int>& multi_index) const;

    /**
     * @brief Convert power coefficients to Bernstein (via HP for precision).
     *
     * Uses 64-bit high-precision intermediate step to minimize conversion errors.
     * Updates bernstein_coeffs_ and sets bernstein_valid_ = true.
     */
    void convertPowerToBernstein() const;

    /**
     * @brief Convert Bernstein coefficients to power basis.
     *
     * Updates power_coeffs_ and sets power_valid_ = true.
     */
    void convertBernsteinToPower() const;
};

} // namespace polynomial_solver

#endif // POLYNOMIAL_H
