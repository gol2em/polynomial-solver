#include <iostream>
#include <iomanip>
#include <cmath>
#include "solver/solver.h"

using namespace polynomial_solver;

/**
 * @brief Test case for 2D ellipse intersection with geometry dump
 * 
 * System:
 *   x^2 + y^2 = 1         (circle)
 *   x^2/4 + 4*y^2 = 1     (ellipse)
 * 
 * In [0,1]^2, we solve:
 *   x^2 + y^2 - 1 = 0
 *   x^2/4 + 4*y^2 - 1 = 0
 */

int main() {
    std::cout << "=== 2D Ellipse Intersection Test with Geometry Dump ===" << std::endl;
    std::cout << std::endl;

    // System:
    //   f1(x,y) = x^2 + y^2 - 1 = 0
    //   f2(x,y) = x^2/4 + 4*y^2 - 1 = 0
    //
    // Convert to Bernstein form on [0,1]^2

    // f1(x,y) = x^2 + y^2 - 1
    // Degrees: (2, 2) where first is x, second is y
    // Storage order with y varying fastest: 1, y, y^2, x, xy, xy^2, x^2, x^2y, x^2y^2
    // f1 = -1 + 0*y + 1*y^2 + 0*x + 0*xy + 0*xy^2 + 1*x^2 + 0*x^2y + 0*x^2y^2
    std::vector<unsigned int> degrees1{2u, 2u};
    std::vector<double> power_coeffs1{
        -1.0,  0.0,  1.0,   // 1, y, y^2
         0.0,  0.0,  0.0,   // x, xy, xy^2
         1.0,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p1 = Polynomial::fromPower(degrees1, power_coeffs1);

    // f2(x,y) = x^2/4 + 4*y^2 - 1
    // f2 = -1 + 0*y + 4*y^2 + 0*x + 0*xy + 0*xy^2 + 0.25*x^2 + 0*x^2y + 0*x^2y^2
    std::vector<unsigned int> degrees2{2u, 2u};
    std::vector<double> power_coeffs2{
        -1.0,   0.0,  4.0,   // 1, y, y^2
         0.0,   0.0,  0.0,   // x, xy, xy^2
         0.25,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p2 = Polynomial::fromPower(degrees2, power_coeffs2);

    std::cout << "Equation 1: x^2 + y^2 - 1 = 0" << std::endl;
    std::cout << "Equation 2: x^2/4 + 4*y^2 - 1 = 0" << std::endl;
    std::cout << std::endl;

    // Create system
    PolynomialSystem system({p1, p2});

    // Configure solver with geometry dump
    SubdivisionConfig config;
    config.tolerance = 1e-6;
    config.max_depth = 50;

#ifdef ENABLE_GEOMETRY_DUMP
    config.dump_geometry = true;
    config.dump_prefix = "dumps/ellipse_test";
#endif

    std::cout << "Configuration:" << std::endl;
    std::cout << "  Tolerance: " << config.tolerance << std::endl;
    std::cout << "  Max depth: " << config.max_depth << std::endl;
#ifdef ENABLE_GEOMETRY_DUMP
    std::cout << "  Dump geometry: " << (config.dump_geometry ? "enabled" : "disabled") << std::endl;
    std::cout << "  Dump prefix: " << config.dump_prefix << std::endl;
#else
    std::cout << "  Dump geometry: disabled (compiled out in release mode)" << std::endl;
#endif
    std::cout << std::endl;

    // Solve with ProjectedPolyhedral method
    std::cout << "Solving with ProjectedPolyhedral method..." << std::endl;
    Solver solver;
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);

    std::cout << std::endl;
    std::cout << "=== Results ===" << std::endl;
    std::cout << "Resolved roots: " << result.num_resolved << std::endl;
    std::cout << "Unresolved boxes: " << (result.boxes.size() - result.num_resolved) << std::endl;
    std::cout << "Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    // Display resolved roots
    if (result.num_resolved > 0) {
        std::cout << "Resolved roots:" << std::endl;
        for (std::size_t i = 0; i < result.num_resolved; ++i) {
            const SubdivisionBoxResult& box = result.boxes[i];
            std::cout << "  Root " << (i + 1) << ":" << std::endl;
            std::cout << "    Center: (";
            for (std::size_t j = 0; j < box.center.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << std::setprecision(10) << box.center[j];
            }
            std::cout << ")" << std::endl;
            std::cout << "    Max error: (";
            for (std::size_t j = 0; j < box.max_error.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << std::scientific << std::setprecision(3) << box.max_error[j];
            }
            std::cout << ")" << std::endl;
            std::cout << "    Depth: " << box.depth << std::endl;

            // Evaluate at center
            std::vector<double> values;
            system.evaluate(box.center, values);
            std::cout << "    Evaluation: (";
            for (std::size_t j = 0; j < values.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << std::scientific << std::setprecision(3) << values[j];
            }
            std::cout << ")" << std::endl;
        }
    }

    std::cout << std::endl;
#ifdef ENABLE_GEOMETRY_DUMP
    std::cout << "Geometry dump written to: " << config.dump_prefix << "_geometry.txt" << std::endl;
#else
    std::cout << "Geometry dump: disabled (compiled out in release mode)" << std::endl;
#endif
    std::cout << std::endl;

    // Verify we found roots
    if (result.num_resolved == 0) {
        std::cerr << "ERROR: No roots found!" << std::endl;
        return 1;
    }

    std::cout << "Test completed successfully!" << std::endl;
    return 0;
}

