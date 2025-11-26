/**
 * @file circle_ellipse_intersection.cpp
 * @brief Example demonstrating 2D polynomial system solving
 *
 * This example shows how to solve a system of polynomial equations
 * representing the intersection of a circle and an ellipse.
 *
 * Problem:
 *   f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)
 *   f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)
 *
 * Expected root in [0,1]^2: approximately (0.894, 0.447)
 *
 * Note: 2D refinement is not yet implemented, so this example only
 * demonstrates the solver step.
 */

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace polynomial_solver;

// Helper function to parse command-line arguments
struct Config {
    // Solver parameters
    double tolerance = 1e-8;
    unsigned int max_depth = 100;
    double degeneracy_multiplier = 5.0;
    bool dump_geometry = false;
    bool dump_result = true;

    bool show_help = false;
};

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tolerance") == 0) {
            if (i + 1 < argc) config.tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--max-depth") == 0) {
            if (i + 1 < argc) config.max_depth = std::atoi(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--degeneracy-multiplier") == 0) {
            if (i + 1 < argc) config.degeneracy_multiplier = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--dump-geometry") == 0) {
            config.dump_geometry = true;
        } else if (strcmp(argv[i], "--no-dump") == 0) {
            config.dump_result = false;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            config.show_help = true;
        }
    }
    return config;
}

void print_help() {
    std::cout << "Usage: circle_ellipse_intersection [OPTIONS]\n\n";
    std::cout << "Solve 2D system: circle x^2+y^2=1 and ellipse x^2/4+4y^2=1\n\n";
    std::cout << "Solver Options:\n";
    std::cout << "  -t, --tolerance <value>           Box size tolerance (default: 1e-8)\n";
    std::cout << "  -d, --max-depth <value>           Maximum subdivision depth (default: 100)\n";
    std::cout << "  -m, --degeneracy-multiplier <val> Degeneracy detection multiplier (default: 5.0)\n";
    std::cout << "  --dump-geometry                   Enable geometry dump for visualization\n";
    std::cout << "  --no-dump                         Disable result dump file\n\n";
    std::cout << "Other Options:\n";
    std::cout << "  -h, --help                        Show this help message\n\n";
    std::cout << "Note: 2D refinement is not yet implemented.\n";
    std::cout << "See docs/PARAMETERS.md for detailed parameter documentation.\n";
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
    // Parse command-line arguments
    Config config = parse_args(argc, argv);

    if (config.show_help) {
        print_help();
        return 0;
    }

    std::cout << "========================================\n";
    std::cout << "Circle-Ellipse Intersection (2D)\n";
    std::cout << "========================================\n\n";

    std::cout << "Problem:\n";
    std::cout << "  f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)\n";
    std::cout << "  f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)\n";
    std::cout << "  Domain: [0, 1] × [0, 1]\n";
    std::cout << "  Expected root: (0.894, 0.447)\n\n";

    // Print configuration
    std::cout << "Solver Configuration:\n";
    std::cout << "  Tolerance: " << std::scientific << config.tolerance << "\n";
    std::cout << "  Max depth: " << config.max_depth << "\n";
    std::cout << "  Degeneracy multiplier: " << std::fixed << std::setprecision(1)
              << config.degeneracy_multiplier << "\n";
    std::cout << "  Geometry dump: " << (config.dump_geometry ? "enabled" : "disabled") << "\n\n";
    
    // Define f1(x,y) = x^2 + y^2 - 1
    std::vector<unsigned int> degrees_f1 = {2, 2};
    std::vector<double> power_coeffs1 = {
        -1.0,  0.0,  1.0,   // 1, y, y^2
         0.0,  0.0,  0.0,   // x, xy, xy^2
         1.0,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs1);

    // Define f2(x,y) = x^2/4 + 4*y^2 - 1
    std::vector<unsigned int> degrees_f2 = {2, 2};
    std::vector<double> power_coeffs2 = {
        -1.0,   0.0,  4.0,   // 1, y, y^2
         0.0,   0.0,  0.0,   // x, xy, xy^2
         0.25,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs2);

    PolynomialSystem system(std::vector<Polynomial>{p1, p2});

    // Configure solver
    SubdivisionConfig solver_config;
    solver_config.tolerance = config.tolerance;
    solver_config.max_depth = config.max_depth;
    solver_config.degeneracy_multiplier = config.degeneracy_multiplier;

#ifdef ENABLE_GEOMETRY_DUMP
    if (config.dump_geometry) {
        solver_config.dump_geometry = true;
        solver_config.dump_prefix = "dumps/circle_ellipse";
    }
#endif

    // ============================================================
    // LINE 1: SOLVE (fast, double precision)
    // ============================================================
    Solver solver;
    auto result = solver.subdivisionSolve(system, solver_config, RootBoundingMethod::ProjectedPolyhedral);

    std::cout << "Step 1: Solve (fast, double precision)\n";
    std::cout << "  Found " << result.num_resolved << " root(s)\n";
    std::cout << "  Unresolved boxes: " << (result.boxes.size() - result.num_resolved) << "\n";
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n\n";

    // Dump solver result to file
    if (config.dump_result) {
        dump_result(result, "dumps/circle_ellipse_result.txt");
    }

    // ============================================================
    // LINE 2: REFINE (not implemented for 2D yet)
    // ============================================================
    ResultRefiner refiner;
    RefinementConfig refine_config = defaultRefinementConfig();
    auto refined = refiner.refine(result, system, refine_config);

    std::cout << "Step 2: Refine (not implemented for 2D yet)\n";
    std::cout << "  2D refinement is not yet implemented.\n";
    std::cout << "  Using solver results directly.\n\n";

    // ============================================================
    // RESULTS
    // ============================================================
    std::cout << "========================================\n";
    std::cout << "Results\n";
    std::cout << "========================================\n\n";

    if (result.num_resolved > 0) {
        std::cout << "Resolved Roots:\n";
        for (std::size_t i = 0; i < result.num_resolved; ++i) {
            const auto& box = result.boxes[i];
            std::cout << "  Root " << (i + 1) << ":\n";
            std::cout << "    Center: (" << std::setprecision(10) << std::fixed
                      << box.center[0] << ", " << box.center[1] << ")\n";
            std::cout << "    Box size: (" << std::scientific << std::setprecision(3)
                      << (box.upper[0] - box.lower[0]) << ", "
                      << (box.upper[1] - box.lower[1]) << ")\n";

            // Evaluate residual
            std::vector<double> eval;
            system.evaluate(box.center, eval);
            std::cout << "    Residual: (" << std::setprecision(4) << eval[0]
                      << ", " << eval[1] << ")\n";

            // Calculate error from expected root
            double expected_x = 2.0 / std::sqrt(5.0);  // ≈ 0.894427191
            double expected_y = 1.0 / std::sqrt(5.0);  // ≈ 0.447213595
            double error_x = std::abs(box.center[0] - expected_x);
            double error_y = std::abs(box.center[1] - expected_y);
            std::cout << "    Error from expected: (" << error_x << ", " << error_y << ")\n\n";
        }
    }

    std::cout << "========================================\n";
    std::cout << "Summary\n";
    std::cout << "========================================\n\n";
    std::cout << "✅ Found " << result.num_resolved << " root(s) in [0,1]×[0,1]\n";
    std::cout << "Note: 2D refinement is not yet implemented.\n";
    std::cout << "For high-precision 2D roots, consider using multi-precision arithmetic.\n\n";

    return 0;
}

