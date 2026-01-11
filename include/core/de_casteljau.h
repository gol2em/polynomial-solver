#ifndef DE_CASTELJAU_H
#define DE_CASTELJAU_H

/**
 * @file de_casteljau.h
 * @brief De Casteljau subdivision algorithm for polynomial curves
 *
 * This module implements the De Casteljau algorithm for subdividing
 * polynomial curves, which is useful for root isolation and refinement.
 */

#include <vector>

#include "core/polynomial.h"
#include "core/geometry.h"

namespace polynomial_solver {

/**
 * @class DeCasteljau
 * @brief Implements De Casteljau subdivision algorithm.
 *
 * The class provides static utility functions for evaluating polynomials
 * represented in Bernstein basis, both in the univariate case and in the
 * tensor-product multivariate case.
 */
class DeCasteljau {
public:
    /**
     * @brief Default constructor
     */
    DeCasteljau();

    /**
     * @brief Destructor
     */
    ~DeCasteljau();

    /**
     * @brief Evaluate a univariate Bernstein polynomial at parameter t.
     *
     * @param bernstein_coeffs Coefficients in Bernstein basis.
     * @param t Parameter value in [0, 1].
     * @return Value of the polynomial at t.
     */
    static double evaluate1D(const std::vector<double>& bernstein_coeffs,
                             double t);

    /**
     * @brief Subdivide a univariate Bernstein polynomial at parameter t.
     *
     * Given coefficients b[0..n] in Bernstein basis on [0,1], this computes
     * two new coefficient sets 'left' and 'right' (size n+1 each) that
     * represent the curve segments on [0,t] and [t,1], each reparameterized
     * back to [0,1].
     *
     * @param bernstein_coeffs Input coefficients in Bernstein basis.
     * @param t Subdivision parameter in [0, 1].
     * @param left Output coefficients for the [0,t] segment.
     * @param right Output coefficients for the [t,1] segment.
     */
    static void subdivide1D(const std::vector<double>& bernstein_coeffs,
                            double t,
                            std::vector<double>& left,
                            std::vector<double>& right);

    /**
     * @brief Evaluate a tensor-product Bernstein polynomial.
     *
     * @param degrees Degrees per variable.
     * @param bernstein_coeffs Coefficients in tensor-product Bernstein basis.
     * @param parameters Parameter values for each variable.
     * @return Value of the polynomial at the given point.
     */
    static double evaluateTensorProduct(const std::vector<unsigned int>& degrees,
                                        const std::vector<double>& bernstein_coeffs,
                                        const std::vector<double>& parameters);

    // TODO: Add further De Casteljau operations
    // - Subdivision at a parameter value
    // - Control point manipulation
    // - Bounding box computation

private:
    // TODO: Add member variables if needed
};

} // namespace polynomial_solver

#endif // DE_CASTELJAU_H
