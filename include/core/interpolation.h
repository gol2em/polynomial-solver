/**
 * @file interpolation.h
 * @brief Polynomial interpolation from function samples.
 *
 * This module provides polynomial interpolation for functions given as
 * function pointers (std::function). Supports 1D and multivariate cases
 * with configurable abscissae (uniform, Chebyshev, etc.).
 *
 * The interpolation works by sampling the function at specified abscissae
 * and fitting a polynomial in power basis using least-squares.
 */

#ifndef POLYNOMIAL_SOLVER_INTERPOLATION_H
#define POLYNOMIAL_SOLVER_INTERPOLATION_H

#include "core/polynomial.h"
#include "config.h"
#include <functional>
#include <vector>

#ifdef ENABLE_HIGH_PRECISION
#include "hp/high_precision_types.h"
#endif

namespace polynomial_solver {

/**
 * @brief Types of abscissae (sample points) for interpolation.
 */
enum class AbscissaeType {
    UNIFORM,    ///< Equally spaced points: x_i = i / n
    CHEBYSHEV,  ///< Chebyshev nodes: x_i = (1 + cos(π(2i+1)/(2n+2))) / 2
    CHEBYSHEV_EXTREMA  ///< Chebyshev extrema (Clenshaw-Curtis): x_i = (1 + cos(πi/n)) / 2
};

/**
 * @brief Polynomial interpolation utilities.
 *
 * Provides static methods for interpolating polynomials from function samples.
 * The function domain is always normalized to [0, 1]^d, with optional
 * transformation to a custom domain.
 */
class Interpolation {
public:
    // ========== 1D Interpolation ==========

    /**
     * @brief Interpolate a univariate function.
     *
     * @param func Function to interpolate, taking one double argument
     * @param degree Polynomial degree
     * @param x_min Lower bound of the domain
     * @param x_max Upper bound of the domain
     * @param abscissae Type of sample points to use
     * @return Interpolating polynomial on [0, 1]
     *
     * The returned polynomial is defined on [0, 1]. To evaluate at a point x
     * in [x_min, x_max], first transform: s = (x - x_min) / (x_max - x_min).
     */
    static Polynomial interpolate1D(
        const std::function<double(double)>& func,
        unsigned int degree,
        double x_min = 0.0,
        double x_max = 1.0,
        AbscissaeType abscissae = AbscissaeType::CHEBYSHEV);

    // ========== 2D Interpolation ==========

    /**
     * @brief Interpolate a bivariate function.
     *
     * @param func Function to interpolate, taking two double arguments
     * @param degree_x Polynomial degree in x direction
     * @param degree_y Polynomial degree in y direction
     * @param x_min Lower bound of x domain
     * @param x_max Upper bound of x domain
     * @param y_min Lower bound of y domain
     * @param y_max Upper bound of y domain
     * @param abscissae Type of sample points to use
     * @return Interpolating polynomial on [0, 1]^2
     */
    static Polynomial interpolate2D(
        const std::function<double(double, double)>& func,
        unsigned int degree_x,
        unsigned int degree_y,
        double x_min = 0.0,
        double x_max = 1.0,
        double y_min = 0.0,
        double y_max = 1.0,
        AbscissaeType abscissae = AbscissaeType::CHEBYSHEV);

    // ========== General Multivariate Interpolation ==========

    /**
     * @brief Interpolate a multivariate function.
     *
     * @param func Function to interpolate, taking a vector of doubles
     * @param degrees Polynomial degrees per variable
     * @param lower_bounds Lower bounds of the domain per variable
     * @param upper_bounds Upper bounds of the domain per variable
     * @param abscissae Type of sample points to use
     * @return Interpolating polynomial on [0, 1]^d
     */
    static Polynomial interpolateND(
        const std::function<double(const std::vector<double>&)>& func,
        const std::vector<unsigned int>& degrees,
        const std::vector<double>& lower_bounds,
        const std::vector<double>& upper_bounds,
        AbscissaeType abscissae = AbscissaeType::CHEBYSHEV);

    // ========== Utility Methods ==========

    /**
     * @brief Generate abscissae points in [0, 1].
     *
     * @param n Number of points (typically degree + 1)
     * @param type Type of abscissae
     * @return Vector of n points in [0, 1]
     */
    static std::vector<double> generateAbscissae(
        unsigned int n,
        AbscissaeType type);

    /**
     * @brief Transform a point from [0, 1] to [a, b].
     */
    static double transformFromUnit(double s, double a, double b) {
        return a + (b - a) * s;
    }

    /**
     * @brief Transform a point from [a, b] to [0, 1].
     */
    static double transformToUnit(double x, double a, double b) {
        return (x - a) / (b - a);
    }

private:
    /**
     * @brief Solve least-squares system A^T A x = A^T b using Gaussian elimination.
     *
     * @param A Vandermonde-like matrix (m x n)
     * @param b Right-hand side vector (m)
     * @return Solution vector (n)
     */
    static std::vector<double> solveLeastSquares(
        const std::vector<std::vector<double>>& A,
        const std::vector<double>& b);
};

//=============================================================================
// Domain Utilities
//=============================================================================

/**
 * @brief 2D domain for coordinate transformations between user space and [0,1]²
 *
 * The solver operates on [0,1]², but user functions are typically defined on
 * arbitrary rectangular domains. This struct handles the coordinate mapping.
 *
 * Usage:
 * @code
 * Domain2D domain(-1.5, 1.5, -1.5, 1.5);  // User domain [-1.5, 1.5]²
 *
 * // Transform user coordinates to unit domain
 * double u, v;
 * domain.toUnit(x, y, u, v);
 *
 * // Transform solver results back to user coordinates
 * double x_result, y_result;
 * domain.fromUnit(u_result, v_result, x_result, y_result);
 *
 * // Wrap user function for interpolation
 * auto f_unit = domain.wrapFunction(f_user);
 * Polynomial p = Interpolation::interpolate2D(f_unit, degree, degree,
 *                                              0.0, 1.0, 0.0, 1.0, AbscissaeType::CHEBYSHEV);
 * @endcode
 */
struct Domain2D {
    double x_min, x_max;  ///< X bounds of user domain
    double y_min, y_max;  ///< Y bounds of user domain

    Domain2D() : x_min(0), x_max(1), y_min(0), y_max(1) {}

    Domain2D(double xmin, double xmax, double ymin, double ymax)
        : x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax) {}

    /// Create symmetric domain [-hw, hw]²
    static Domain2D symmetric(double half_width) {
        return Domain2D(-half_width, half_width, -half_width, half_width);
    }

    /// Transform from user coordinates (x,y) to unit coordinates (u,v)
    void toUnit(double x, double y, double& u, double& v) const {
        u = (x - x_min) / (x_max - x_min);
        v = (y - y_min) / (y_max - y_min);
    }

    /// Transform from unit coordinates (u,v) to user coordinates (x,y)
    void fromUnit(double u, double v, double& x, double& y) const {
        x = x_min + (x_max - x_min) * u;
        y = y_min + (y_max - y_min) * v;
    }

    /// Wrap a user-space function for use with interpolation on [0,1]²
    std::function<double(double, double)> wrapFunction(
        const std::function<double(double, double)>& f_user) const
    {
        double xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max;
        return [=](double u, double v) -> double {
            double x = xmin + (xmax - xmin) * u;
            double y = ymin + (ymax - ymin) * v;
            return f_user(x, y);
        };
    }

#ifdef ENABLE_HIGH_PRECISION
    /// Wrap a high-precision user-space function for use on [0,1]²
    std::function<mpreal(const mpreal&, const mpreal&)> wrapFunctionHP(
        const std::function<mpreal(const mpreal&, const mpreal&)>& f_user) const
    {
        mpreal xmin(x_min), xmax(x_max), ymin(y_min), ymax(y_max);
        return [=](const mpreal& u, const mpreal& v) -> mpreal {
            mpreal x = xmin + (xmax - xmin) * u;
            mpreal y = ymin + (ymax - ymin) * v;
            return f_user(x, y);
        };
    }
#endif
};

} // namespace polynomial_solver

#endif // POLYNOMIAL_SOLVER_INTERPOLATION_H

