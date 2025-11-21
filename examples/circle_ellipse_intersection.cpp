/**
 * Example: Circle-Ellipse Intersection
 * 
 * This example demonstrates solving a system of polynomial equations
 * representing the intersection of a circle and an ellipse.
 * 
 * Problem:
 *   f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)
 *   f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)
 * 
 * Expected root in [0,1]^2: approximately (0.894, 0.447)
 * 
 * This example tests all three subdivision strategies:
 * - ContractFirst: Contract bounds before subdividing
 * - SubdivideFirst: Subdivide before contracting
 * - Simultaneous: Balance both approaches
 */

#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>

using namespace polynomial_solver;

void print_result(const std::string& strategy_name,
                  const SubdivisionSolverResult& result,
                  const PolynomialSystem& system) {
    std::cout << "\n" << strategy_name << " Strategy:\n";
    std::cout << std::string(60, '=') << "\n";

    std::cout << "  Total boxes: " << result.boxes.size() << "\n";
    std::cout << "  Resolved: " << result.num_resolved << "\n";
    std::cout << "  Unresolved: " << (result.boxes.size() - result.num_resolved) << "\n";
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";

    if (result.num_resolved > 0) {
        std::cout << "\n  Resolved roots:\n";
        for (size_t i = 0; i < result.num_resolved; ++i) {
            const SubdivisionBoxResult& box = result.boxes[i];
            std::cout << "    Root " << (i+1) << ":\n";
            std::cout << "      Center: (" << std::setprecision(10) << box.center[0]
                      << ", " << box.center[1] << ")\n";
            std::cout << "      Max error: (" << std::scientific << std::setprecision(3)
                      << box.max_error[0] << ", " << box.max_error[1] << ")\n";
            std::cout << "      Depth: " << box.depth << "\n";

            // Evaluate at center to verify
            std::vector<double> eval;
            system.evaluate(box.center, eval);
            std::cout << "      Residual: (" << eval[0] << ", " << eval[1] << ")\n";

            // Calculate actual distance from expected root
            double expected_x = 2.0 / std::sqrt(5.0);  // ≈ 0.894427191
            double expected_y = 1.0 / std::sqrt(5.0);  // ≈ 0.447213595
            double error_x = std::abs(box.center[0] - expected_x);
            double error_y = std::abs(box.center[1] - expected_y);
            std::cout << "      Error from expected: (" << error_x << ", " << error_y << ")\n";
        }
    }

    if (result.boxes.size() > result.num_resolved) {
        std::cout << "\n  Unresolved boxes:\n";
        for (size_t i = result.num_resolved; i < result.boxes.size(); ++i) {
            const SubdivisionBoxResult& box = result.boxes[i];
            std::cout << "    Box " << (i - result.num_resolved + 1) << ":\n";
            std::cout << "      Center: (" << std::setprecision(10) << box.center[0]
                      << ", " << box.center[1] << ")\n";
            std::cout << "      Box size: (" << std::scientific << std::setprecision(3)
                      << (box.upper[0] - box.lower[0]) << ", "
                      << (box.upper[1] - box.lower[1]) << ")\n";
            std::cout << "      Depth: " << box.depth << "\n";
        }
    }
}

void dump_result(const SubdivisionSolverResult& result, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }

    out << "# Solver Result Dump\n";
    out << "# Dimension: 2\n";
    out << "# Total boxes: " << result.boxes.size() << "\n";
    out << "# Resolved: " << result.num_resolved << "\n";
    out << "# Unresolved: " << (result.boxes.size() - result.num_resolved) << "\n";
    out << "# Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";
    out << "\n";

    // Dump resolved boxes
    out << "## Resolved Boxes\n";
    for (size_t i = 0; i < result.num_resolved; ++i) {
        const SubdivisionBoxResult& box = result.boxes[i];
        out << "Box " << i << "\n";
        out << "  Lower: " << std::setprecision(17) << box.lower[0] << " " << box.lower[1] << "\n";
        out << "  Upper: " << box.upper[0] << " " << box.upper[1] << "\n";
        out << "  Center: " << box.center[0] << " " << box.center[1] << "\n";
        out << "  Depth: " << box.depth << "\n";
        out << "  Converged: " << (box.converged ? "yes" : "no") << "\n";
        out << "\n";
    }

    // Dump unresolved boxes
    if (result.boxes.size() > result.num_resolved) {
        out << "## Unresolved Boxes\n";
        for (size_t i = result.num_resolved; i < result.boxes.size(); ++i) {
            const SubdivisionBoxResult& box = result.boxes[i];
            out << "Box " << (i - result.num_resolved) << "\n";
            out << "  Lower: " << std::setprecision(17) << box.lower[0] << " " << box.lower[1] << "\n";
            out << "  Upper: " << box.upper[0] << " " << box.upper[1] << "\n";
            out << "  Center: " << box.center[0] << " " << box.center[1] << "\n";
            out << "  Depth: " << box.depth << "\n";
            out << "  Converged: " << (box.converged ? "yes" : "no") << "\n";
            out << "\n";
        }
    }

    out.close();
    std::cout << "  Result dump saved to: " << filename << "\n";
}

int main(int argc, char* argv[]) {
    // Default parameters
    bool dump_geometry = false;
    double tolerance = 1e-8;
    unsigned int max_depth = 100;
    double degeneracy_multiplier = 5.0;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--dump-geometry") == 0) {
            dump_geometry = true;
        } else if (std::strcmp(argv[i], "--tolerance") == 0 || std::strcmp(argv[i], "-t") == 0) {
            if (i + 1 < argc) {
                tolerance = std::atof(argv[++i]);
            }
        } else if (std::strcmp(argv[i], "--max-depth") == 0 || std::strcmp(argv[i], "-d") == 0) {
            if (i + 1 < argc) {
                max_depth = std::atoi(argv[++i]);
            }
        } else if (std::strcmp(argv[i], "--degeneracy-multiplier") == 0 || std::strcmp(argv[i], "-m") == 0) {
            if (i + 1 < argc) {
                degeneracy_multiplier = std::atof(argv[++i]);
            }
        } else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            std::cout << "Usage: " << argv[0] << " [options]\n";
            std::cout << "Options:\n";
            std::cout << "  --dump-geometry              Enable geometry dump for visualization (default: off)\n";
            std::cout << "  --tolerance, -t <value>      Box size tolerance for convergence (default: 1e-6)\n";
            std::cout << "  --max-depth, -d <value>      Maximum subdivision depth (default: 100)\n";
            std::cout << "  --degeneracy-multiplier, -m <value>  Multiplier for degeneracy detection (default: 2.0)\n";
            std::cout << "  --help, -h                   Show this help message\n";
            return 0;
        }
    }

    std::cout << "Circle-Ellipse Intersection Example\n";
    std::cout << std::string(60, '=') << "\n\n";

    std::cout << "Problem:\n";
    std::cout << "  f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)\n";
    std::cout << "  f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)\n";
    std::cout << "  Domain: [0, 1] × [0, 1]\n\n";

    std::cout << "Parameters:\n";
    std::cout << "  Tolerance:             " << std::scientific << std::setprecision(3) << tolerance << "\n";
    std::cout << "  Max depth:             " << max_depth << "\n";
    std::cout << "  Degeneracy multiplier: " << std::fixed << std::setprecision(1) << degeneracy_multiplier << "\n\n";
    
    // Define f1(x,y) = x^2 + y^2 - 1
    std::vector<unsigned int> degrees_f1{2u, 2u};
    std::vector<double> power_coeffs1{
        -1.0,  0.0,  1.0,   // 1, y, y^2
         0.0,  0.0,  0.0,   // x, xy, xy^2
         1.0,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs1);
    
    // Define f2(x,y) = x^2/4 + 4*y^2 - 1
    std::vector<unsigned int> degrees_f2{2u, 2u};
    std::vector<double> power_coeffs2{
        -1.0,   0.0,  4.0,   // 1, y, y^2
         0.0,   0.0,  0.0,   // x, xy, xy^2
         0.25,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs2);
    
    PolynomialSystem system({p1, p2});
    
    // Test all three strategies
    const char* strategy_names[] = {"ContractFirst", "SubdivideFirst", "Simultaneous"};
    SubdivisionStrategy strategies[] = {
        SubdivisionStrategy::ContractFirst,
        SubdivisionStrategy::SubdivideFirst,
        SubdivisionStrategy::Simultaneous
    };
    
    Solver solver;
    
    for (int s = 0; s < 3; ++s) {
        SubdivisionConfig config;
        config.tolerance = tolerance;
        config.max_depth = max_depth;
        config.degeneracy_multiplier = degeneracy_multiplier;
        config.contraction_threshold = 0.9;
        config.strategy = strategies[s];
#ifdef ENABLE_GEOMETRY_DUMP
        config.dump_geometry = dump_geometry;
        config.dump_prefix = std::string("dumps/example_") + strategy_names[s];
#endif

        SubdivisionSolverResult result = solver.subdivisionSolve(
            system, config, RootBoundingMethod::ProjectedPolyhedral);

        print_result(strategy_names[s], result, system);

        // Dump result to file
        std::string result_file = std::string("dumps/example_") + strategy_names[s] + "_result.txt";
        dump_result(result, result_file);
    }

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Example completed successfully!\n";
#ifdef ENABLE_GEOMETRY_DUMP
    if (dump_geometry) {
        std::cout << "\nGeometry dumps saved to dumps/example_*_geometry.txt\n";
    }
#endif
    std::cout << "Result dumps saved to dumps/example_*_result.txt\n";
    std::cout << "Visualize with: tools/visualize_solver.py\n";

    return 0;
}

