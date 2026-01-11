/**
 * @file test_newton_multidim.cpp
 * @brief Test multidimensional Newton refinement with convergence analysis
 *
 * Tests the 2D Newton method on the circle-ellipse intersection problem:
 *   f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)
 *   f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)
 *
 * The intersection in [0,1]^2 is approximately (0.8944271909999159, 0.4472135954999579)
 * which is (2/√5, 1/√5).
 *
 * We verify:
 * 1. Newton's method converges to the correct root
 * 2. The convergence is quadratic (order ≈ 2) for this simple root
 */

#include "solver/newton_multidim.h"
#include "core/polynomial_base.h"
#include "hp/high_precision_types.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace polynomial_solver;

// Test with double precision
void test_circle_ellipse_double() {
    std::cout << "=== Circle-Ellipse Newton Refinement (double) ===\n\n";
    
    // f1(x,y) = x^2 + y^2 - 1
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> coeffs_f1 = {
        -1.0, 0.0, 1.0,   // constant, y, y^2
         0.0, 0.0, 0.0,   // x, xy, xy^2
         1.0, 0.0, 0.0    // x^2, x^2*y, x^2*y^2
    };
    auto f1 = PolynomialBase<double>::fromPower(degrees, coeffs_f1);
    
    // f2(x,y) = x^2/4 + 4*y^2 - 1
    std::vector<double> coeffs_f2 = {
        -1.0,  0.0, 4.0,   // constant, y, y^2
         0.0,  0.0, 0.0,   // x, xy, xy^2
         0.25, 0.0, 0.0    // x^2, x^2*y, x^2*y^2
    };
    auto f2 = PolynomialBase<double>::fromPower(degrees, coeffs_f2);
    
    // Exact root: (2/√5, 1/√5)
    double exact_x = 2.0 / std::sqrt(5.0);
    double exact_y = 1.0 / std::sqrt(5.0);
    
    std::cout << "Exact root: (" << std::setprecision(16) << exact_x 
              << ", " << exact_y << ")\n\n";
    
    // Initial guess from subdivision solver (somewhat close)
    double x0 = 0.89;
    double y0 = 0.45;
    
    std::cout << "Initial guess: (" << x0 << ", " << y0 << ")\n";
    std::cout << "Initial error: " << std::sqrt(std::pow(x0-exact_x,2) + std::pow(y0-exact_y,2)) << "\n\n";
    
    // Configure Newton
    NewtonMultidimConfig<double> config;
    config.tolerance = 1e-15;
    config.residual_tolerance = 1e-15;
    config.max_iterations = 20;
    config.verbose = true;
    
    std::cout << "Newton iterations:\n";
    auto result = refineRoot2D(x0, y0, f1, f2, config);
    
    std::cout << "\nResult:\n";
    std::cout << "  Converged: " << (result.converged ? "yes" : "no") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Location: (" << result.location[0] << ", " << result.location[1] << ")\n";
    std::cout << "  Residual norm: " << result.residual_norm << "\n";
    
    double final_error = std::sqrt(
        std::pow(result.location[0] - exact_x, 2) + 
        std::pow(result.location[1] - exact_y, 2));
    std::cout << "  Final error: " << final_error << "\n";
    
    // Print error history for analysis
    std::cout << "\nError history and convergence ratios:\n";
    std::cout << std::setw(6) << "Iter" << std::setw(15) << "|F|"
              << std::setw(15) << "ratio" << std::setw(15) << "order\n";
    std::cout << std::string(51, '-') << "\n";

    for (std::size_t i = 0; i < result.error_history.size(); ++i) {
        std::cout << std::setw(6) << i
                  << std::setw(15) << std::scientific << std::setprecision(3)
                  << result.error_history[i];

        if (i >= 1) {
            double ratio = result.error_history[i] / result.error_history[i-1];
            std::cout << std::setw(15) << ratio;
        } else {
            std::cout << std::setw(15) << "-";
        }

        if (i >= 2) {
            double e_prev2 = result.error_history[i-2];
            double e_prev1 = result.error_history[i-1];
            double e_curr = result.error_history[i];
            if (e_prev1 < e_prev2 && e_curr < e_prev1 && e_prev1 > 1e-100) {
                double p = std::log(e_curr / e_prev1) / std::log(e_prev1 / e_prev2);
                std::cout << std::setw(15) << std::fixed << std::setprecision(2) << p;
            } else {
                std::cout << std::setw(15) << "-";
            }
        } else {
            std::cout << std::setw(15) << "-";
        }
        std::cout << "\n";
    }

    // Manual convergence order calculation using middle iterations
    double order = 0;
    if (result.error_history.size() >= 3) {
        double e0 = result.error_history[0];
        double e1 = result.error_history[1];
        double e2 = result.error_history[2];
        if (e1 < e0 && e2 < e1 && e1 > 1e-100) {
            order = std::log(e2 / e1) / std::log(e1 / e0);
        }
    }
    std::cout << "\nConvergence order (iterations 0-1-2): " << std::setprecision(3) << order << "\n";
    std::cout << "(Expected ~2.0 for Newton's method on simple roots)\n";

    assert(result.converged);
    assert(final_error < 1e-14);
    assert(order > 1.8 && order < 2.5);  // Should be close to 2
    
    std::cout << "\nPASSED\n\n";
}

#ifdef ENABLE_HIGH_PRECISION
// Test with high precision to verify quadratic convergence more precisely
void test_circle_ellipse_high_precision() {
    std::cout << "=== Circle-Ellipse Newton Refinement (High Precision) ===\n\n";
    
    using HP = mpreal;
    
    // Set high precision (256 bits ≈ 77 decimal digits)
    setPrecision(256);
    
    std::cout << "Working precision: " << getPrecision() << " bits\n\n";
    
    // f1(x,y) = x^2 + y^2 - 1
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<HP> coeffs_f1 = {
        HP("-1"), HP("0"), HP("1"),
        HP("0"),  HP("0"), HP("0"),
        HP("1"),  HP("0"), HP("0")
    };
    auto f1 = PolynomialBase<HP>::fromPower(degrees, coeffs_f1);
    
    // f2(x,y) = x^2/4 + 4*y^2 - 1
    std::vector<HP> coeffs_f2 = {
        HP("-1"),   HP("0"), HP("4"),
        HP("0"),    HP("0"), HP("0"),
        HP("0.25"), HP("0"), HP("0")
    };
    auto f2 = PolynomialBase<HP>::fromPower(degrees, coeffs_f2);
    
    // Exact root: (2/√5, 1/√5)
    HP exact_x = HP("2") / sqrt(HP("5"));
    HP exact_y = HP("1") / sqrt(HP("5"));
    
    std::cout << "Exact root (high precision):\n";
    std::cout << "  x = " << exact_x.str(50) << "\n";
    std::cout << "  y = " << exact_y.str(50) << "\n\n";
    
    // Initial guess
    HP x0("0.89");
    HP y0("0.45");
    
    // Configure Newton with tight tolerance
    NewtonMultidimConfig<HP> config;
    config.tolerance = HP("1e-70");
    config.residual_tolerance = HP("1e-70");
    config.max_iterations = 30;
    config.verbose = false;

    auto result = refineRoot2D(x0, y0, f1, f2, config);

    std::cout << "Result:\n";
    std::cout << "  Converged: " << (result.converged ? "yes" : "no") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Final x = " << result.location[0].str(50) << "\n";
    std::cout << "  Final y = " << result.location[1].str(50) << "\n";

    HP error_x = abs(result.location[0] - exact_x);
    HP error_y = abs(result.location[1] - exact_y);
    HP final_error = sqrt(error_x*error_x + error_y*error_y);

    std::cout << "  Final error: " << final_error.str(10) << "\n\n";

    // Detailed convergence analysis
    std::cout << "Convergence analysis (quadratic verification):\n";
    std::cout << "For quadratic convergence: e_{n+1} ≈ C * e_n^2\n";
    std::cout << "So: log(e_{n+1}) / log(e_n) ≈ 2\n\n";

    std::cout << std::setw(6) << "Iter"
              << std::setw(20) << "|F(x)|"
              << std::setw(20) << "log ratio"
              << std::setw(15) << "est. order\n";
    std::cout << std::string(61, '-') << "\n";

    for (std::size_t i = 0; i < result.error_history.size(); ++i) {
        double e = static_cast<double>(result.error_history[i]);
        std::cout << std::setw(6) << i
                  << std::setw(20) << std::scientific << std::setprecision(6) << e;

        if (i >= 2) {
            double e_prev2 = static_cast<double>(result.error_history[i-2]);
            double e_prev1 = static_cast<double>(result.error_history[i-1]);
            double e_curr = static_cast<double>(result.error_history[i]);

            if (e_prev1 > 1e-100 && e_curr > 1e-100 && e_prev2 > 1e-100 &&
                e_prev1 < e_prev2 && e_curr < e_prev1) {
                double log_ratio = std::log(e_curr / e_prev1) / std::log(e_prev1 / e_prev2);
                std::cout << std::setw(20) << std::fixed << std::setprecision(6) << log_ratio
                          << std::setw(15) << std::fixed << std::setprecision(3) << log_ratio;
            }
        }
        std::cout << "\n";
    }

    double order = analyzeConvergenceOrder(result.error_history);
    std::cout << "\nOverall estimated order: " << std::setprecision(4) << order << "\n";

    assert(result.converged);
    assert(static_cast<double>(final_error) < 1e-60);
    assert(order > 1.8 && order < 2.2);

    std::cout << "\nPASSED\n\n";
}
#endif

int main() {
    std::cout << "======================================\n";
    std::cout << "Multidimensional Newton Refinement Test\n";
    std::cout << "======================================\n\n";

    test_circle_ellipse_double();

#ifdef ENABLE_HIGH_PRECISION
    test_circle_ellipse_high_precision();
#else
    std::cout << "High precision test skipped (ENABLE_HIGH_PRECISION not defined)\n";
#endif

    std::cout << "======================================\n";
    std::cout << "All tests PASSED\n";
    std::cout << "======================================\n";

    return 0;
}

