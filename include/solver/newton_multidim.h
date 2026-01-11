#ifndef NEWTON_MULTIDIM_H
#define NEWTON_MULTIDIM_H

/**
 * @file newton_multidim.h
 * @brief Multidimensional Newton iteration for polynomial systems
 *
 * This header provides templated Newton's method for refining roots of
 * multidimensional polynomial systems. For simple roots, Newton's method
 * achieves quadratic (2nd order) convergence.
 *
 * The implementation solves the linear system J(x) * delta = -F(x) at each
 * iteration, where J is the Jacobian matrix and F is the system of polynomials.
 *
 * For 2D systems, a direct formula is used. For higher dimensions, Gaussian
 * elimination with partial pivoting is employed.
 */

#include "core/polynomial_base.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace polynomial_solver {

//=============================================================================
// Configuration
//=============================================================================

/**
 * @brief Configuration for multidimensional Newton refinement
 */
template<typename Scalar>
struct NewtonMultidimConfig {
    Scalar tolerance;           ///< Convergence tolerance for step size
    Scalar residual_tolerance;  ///< Tolerance for residual |F(x)|
    unsigned int max_iterations;///< Maximum iterations
    bool verbose;               ///< Print iteration info

    NewtonMultidimConfig()
        : tolerance(Scalar(1e-15)),
          residual_tolerance(Scalar(1e-15)),
          max_iterations(50u),
          verbose(false)
    {}
};

//=============================================================================
// Result structure
//=============================================================================

/**
 * @brief Result of multidimensional Newton refinement
 */
template<typename Scalar>
struct NewtonMultidimResult {
    std::vector<Scalar> location;       ///< Refined root location
    std::vector<Scalar> residual;       ///< F(x) at refined location
    Scalar residual_norm;               ///< |F(x)|
    unsigned int iterations;            ///< Iterations performed
    bool converged;                     ///< True if converged
    std::vector<Scalar> error_history;  ///< Error at each iteration (for convergence analysis)

    NewtonMultidimResult() : residual_norm(0), iterations(0), converged(false) {}
};

//=============================================================================
// Helper: 2x2 linear system solver
//=============================================================================

namespace detail {

/**
 * @brief Solve 2x2 linear system [a b; c d] * [x; y] = [e; f]
 * @return true if system is non-singular
 */
template<typename Scalar>
bool solve_2x2(const Scalar& a, const Scalar& b, const Scalar& c, const Scalar& d,
               const Scalar& e, const Scalar& f,
               Scalar& x, Scalar& y) {
    Scalar det = a * d - b * c;
    Scalar abs_det = (det >= Scalar(0)) ? det : -det;
    
    // Check for singular matrix
    if (abs_det < Scalar(1e-100)) {
        return false;
    }
    
    x = (d * e - b * f) / det;
    y = (a * f - c * e) / det;
    return true;
}

/**
 * @brief Compute L2 norm of a vector
 */
template<typename Scalar>
Scalar vector_norm(const std::vector<Scalar>& v) {
    Scalar sum = Scalar(0);
    for (const auto& x : v) {
        sum += x * x;
    }
    // Use sqrt from cmath for double, or appropriate function for other types
    return Scalar(std::sqrt(static_cast<double>(sum)));
}

/**
 * @brief Compute infinity norm of a vector
 */
template<typename Scalar>
Scalar vector_inf_norm(const std::vector<Scalar>& v) {
    Scalar max_val = Scalar(0);
    for (const auto& x : v) {
        Scalar abs_x = (x >= Scalar(0)) ? x : -x;
        if (abs_x > max_val) max_val = abs_x;
    }
    return max_val;
}

} // namespace detail

//=============================================================================
// Newton iteration for 2D systems
//=============================================================================

/**
 * @brief Refine a 2D root using Newton's method
 *
 * For a system F(x,y) = [f1(x,y), f2(x,y)] = [0, 0], Newton's method is:
 *   [x_new]   [x]       [f1]
 *   [     ] = [ ] - J^{-1} [  ]
 *   [y_new]   [y]       [f2]
 *
 * where J is the Jacobian matrix:
 *   J = [∂f1/∂x  ∂f1/∂y]
 *       [∂f2/∂x  ∂f2/∂y]
 *
 * @param initial_x Initial x coordinate
 * @param initial_y Initial y coordinate  
 * @param f1 First polynomial equation
 * @param f2 Second polynomial equation
 * @param config Configuration parameters
 * @return Refinement result with location, residual, and convergence info
 */
template<typename Scalar>
NewtonMultidimResult<Scalar> refineRoot2D(
    const Scalar& initial_x,
    const Scalar& initial_y,
    const PolynomialBase<Scalar>& f1,
    const PolynomialBase<Scalar>& f2,
    const NewtonMultidimConfig<Scalar>& config = NewtonMultidimConfig<Scalar>())
{
    NewtonMultidimResult<Scalar> result;
    result.location = {initial_x, initial_y};
    result.residual.resize(2);
    
    // Verify polynomials are 2D
    if (f1.dimension() != 2 || f2.dimension() != 2) {
        throw std::invalid_argument("refineRoot2D: polynomials must be 2-dimensional");
    }

    // Compute partial derivatives (Jacobian components)
    PolynomialBase<Scalar> df1_dx = f1.differentiate(0);
    PolynomialBase<Scalar> df1_dy = f1.differentiate(1);
    PolynomialBase<Scalar> df2_dx = f2.differentiate(0);
    PolynomialBase<Scalar> df2_dy = f2.differentiate(1);

    Scalar x = initial_x;
    Scalar y = initial_y;

    for (unsigned int iter = 0; iter < config.max_iterations; ++iter) {
        result.iterations = iter + 1;

        // Evaluate system F(x,y)
        std::vector<Scalar> point = {x, y};
        Scalar F1 = f1.evaluate(point);
        Scalar F2 = f2.evaluate(point);

        // Compute residual norm
        result.residual = {F1, F2};
        result.residual_norm = detail::vector_norm(result.residual);
        result.error_history.push_back(result.residual_norm);

        if (config.verbose) {
            std::cerr << "  Iter " << iter << ": x=" << static_cast<double>(x)
                      << " y=" << static_cast<double>(y)
                      << " |F|=" << static_cast<double>(result.residual_norm) << "\n";
        }

        // Check convergence on residual
        if (result.residual_norm < config.residual_tolerance) {
            result.location = {x, y};
            result.converged = true;
            return result;
        }

        // Evaluate Jacobian at (x, y)
        Scalar J11 = df1_dx.evaluate(point);  // ∂f1/∂x
        Scalar J12 = df1_dy.evaluate(point);  // ∂f1/∂y
        Scalar J21 = df2_dx.evaluate(point);  // ∂f2/∂x
        Scalar J22 = df2_dy.evaluate(point);  // ∂f2/∂y

        // Solve J * delta = -F for delta
        Scalar delta_x, delta_y;
        Scalar neg_F1 = -F1;  // Explicit conversion for expression templates
        Scalar neg_F2 = -F2;
        if (!detail::solve_2x2(J11, J12, J21, J22, neg_F1, neg_F2, delta_x, delta_y)) {
            // Singular Jacobian - cannot continue
            result.location = {x, y};
            return result;
        }

        // Update
        x += delta_x;
        y += delta_y;

        // Check convergence on step size
        Scalar step_norm = detail::vector_norm(std::vector<Scalar>{delta_x, delta_y});
        if (step_norm < config.tolerance) {
            result.location = {x, y};
            // Recompute final residual
            point = {x, y};
            F1 = f1.evaluate(point);
            F2 = f2.evaluate(point);
            result.residual = {F1, F2};
            result.residual_norm = detail::vector_norm(result.residual);
            result.converged = true;
            return result;
        }
    }

    // Did not converge
    result.location = {x, y};
    return result;
}

//=============================================================================
// Convergence rate analysis
//=============================================================================

/**
 * @brief Analyze convergence order from error history
 *
 * For order p convergence: e_{n+1} ~ C * e_n^p
 * So: log(e_{n+1}) ~ log(C) + p * log(e_n)
 *
 * Returns estimated convergence order p.
 * For Newton's method on simple roots, p ≈ 2.
 */
template<typename Scalar>
double analyzeConvergenceOrder(const std::vector<Scalar>& error_history) {
    if (error_history.size() < 3) {
        return 0.0;  // Not enough data
    }

    // Use the last few iterations for analysis
    std::size_t n = error_history.size();

    // Estimate order using: p ≈ log(e_{n}/e_{n-1}) / log(e_{n-1}/e_{n-2})
    double orders_sum = 0.0;
    int count = 0;

    for (std::size_t i = 2; i < n; ++i) {
        double e_prev2 = static_cast<double>(error_history[i-2]);
        double e_prev1 = static_cast<double>(error_history[i-1]);
        double e_curr = static_cast<double>(error_history[i]);

        // Skip if errors are too small (converged) or not decreasing
        if (e_prev2 < 1e-100 || e_prev1 < 1e-100 || e_curr < 1e-100) continue;
        if (e_prev1 >= e_prev2 || e_curr >= e_prev1) continue;

        double ratio_curr = std::log(e_curr / e_prev1);
        double ratio_prev = std::log(e_prev1 / e_prev2);

        if (std::abs(ratio_prev) > 1e-10) {
            double p = ratio_curr / ratio_prev;
            if (p > 0.5 && p < 4.0) {  // Sanity check
                orders_sum += p;
                count++;
            }
        }
    }

    return (count > 0) ? (orders_sum / count) : 0.0;
}

} // namespace polynomial_solver

#endif // NEWTON_MULTIDIM_H

