#include "core/polynomial.h"
#include "solver/solver.h"
#include <iostream>
#include <iomanip>

using namespace polynomial_solver;

/**
 * Test all three subdivision strategies on the same problem:
 * - ContractFirst (default)
 * - SubdivideFirst
 * - Simultaneous
 * 
 * Problem: Intersection of two ellipses in 2D
 * f1(x,y) = x^2 + y^2 - 1
 * f2(x,y) = x^2/4 + 4*y^2 - 1
 * 
 * Expected root in [0,1]^2: approximately (0.894, 0.447)
 */
int main() {
    std::cout << "Testing subdivision strategies on ellipse intersection problem\n";
    std::cout << "=============================================================\n\n";

    // f1(x,y) = x^2 + y^2 - 1
    // Power basis: -1 + 0*y + 1*y^2 + 0*x + 0*xy + 0*xy^2 + 1*x^2 + 0*x^2y + 0*x^2y^2
    std::vector<unsigned int> degrees_f1{2u, 2u};
    std::vector<double> power_coeffs1{
        -1.0,  0.0,  1.0,   // 1, y, y^2
         0.0,  0.0,  0.0,   // x, xy, xy^2
         1.0,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs1);

    // f2(x,y) = x^2/4 + 4*y^2 - 1
    // Power basis: -1 + 0*y + 4*y^2 + 0*x + 0*xy + 0*xy^2 + 0.25*x^2 + 0*x^2y + 0*x^2y^2
    std::vector<unsigned int> degrees_f2{2u, 2u};
    std::vector<double> power_coeffs2{
        -1.0,   0.0,  4.0,   // 1, y, y^2
         0.0,   0.0,  0.0,   // x, xy, xy^2
         0.25,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs2);

    PolynomialSystem system({p1, p2});

    // Test each strategy
    const char* strategy_names[] = {"ContractFirst", "SubdivideFirst", "Simultaneous"};
    SubdivisionStrategy strategies[] = {
        SubdivisionStrategy::ContractFirst,
        SubdivisionStrategy::SubdivideFirst,
        SubdivisionStrategy::Simultaneous
    };

    for (int s = 0; s < 3; ++s) {
        std::cout << "\n" << strategy_names[s] << " Strategy:\n";
        std::cout << std::string(50, '-') << "\n";

        Solver solver;
        SubdivisionConfig config;
        config.tolerance = 1e-6;
        config.max_depth = 100;
        config.contraction_threshold = 0.9;  // If new_width/old_width > 0.9, consider it insufficient
        config.strategy = strategies[s];

#ifdef ENABLE_GEOMETRY_DUMP
        config.dump_geometry = true;
        config.dump_prefix = std::string("dumps/strategy_") + strategy_names[s];
#endif

        SubdivisionSolverResult result = solver.subdivisionSolve(
            system, config, RootBoundingMethod::ProjectedPolyhedral);

        std::cout << "  Found " << result.boxes.size() << " box(es), "
                  << result.num_resolved << " resolved\n";
        std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";

        if (result.num_resolved > 0) {
            const SubdivisionBoxResult& box = result.boxes[0];
            std::cout << "  Root 1:\n";
            std::cout << "    Center: (" << std::setprecision(10) << box.center[0]
                      << ", " << box.center[1] << ")\n";
            std::cout << "    Max error: (" << std::scientific << std::setprecision(3)
                      << box.max_error[0] << ", " << box.max_error[1] << ")\n";
            std::cout << "    Depth: " << box.depth << "\n";

            // Evaluate at center
            std::vector<double> eval;
            system.evaluate(box.center, eval);
            std::cout << "    Evaluation: (" << eval[0] << ", " << eval[1] << ")\n";
        }

#ifdef ENABLE_GEOMETRY_DUMP
        std::cout << "  Geometry dump written to: " << config.dump_prefix << "_geometry.txt\n";

        // Suggest visualization command
        std::cout << "  Visualize: .venv/bin/python visualize_ellipse_dump.py "
                  << config.dump_prefix << "_geometry.txt --max-steps 5 --output-dir visualizations/viz_" << strategy_names[s] << "\n";
#else
        std::cout << "  Geometry dump: disabled (compiled out in release mode)\n";
#endif
    }

    std::cout << "\n=============================================================\n";
    std::cout << "All strategies tested successfully!\n";
#ifdef ENABLE_GEOMETRY_DUMP
    std::cout << "\nNext steps:\n";
    std::cout << "1. Visualize each strategy (first 5 steps):\n";
    std::cout << "   .venv/bin/python visualize_ellipse_dump.py dumps/strategy_ContractFirst_geometry.txt --max-steps 5 --output-dir visualizations/viz_ContractFirst\n";
    std::cout << "   .venv/bin/python visualize_ellipse_dump.py dumps/strategy_SubdivideFirst_geometry.txt --max-steps 5 --output-dir visualizations/viz_SubdivideFirst\n";
    std::cout << "   .venv/bin/python visualize_ellipse_dump.py dumps/strategy_Simultaneous_geometry.txt --max-steps 5 --output-dir visualizations/viz_Simultaneous\n";
    std::cout << "2. Compare the visualizations to understand the differences\n";
#else
    std::cout << "\nNote: Geometry dump is disabled (compiled out in release mode)\n";
    std::cout << "Rebuild in Debug mode to enable geometry visualization.\n";
#endif

    return 0;
}

