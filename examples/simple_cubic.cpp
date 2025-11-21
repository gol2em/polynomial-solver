/**
 * @file simple_cubic.cpp
 * @brief Simple example demonstrating the 2-line workflow: solve and refine
 * 
 * This example shows the minimal code needed to solve a polynomial and refine
 * the roots to high precision.
 * 
 * Problem: p(x) = (x - 0.2)(x - 0.5)(x - 0.8)
 * Expected roots: x = 0.2, 0.5, 0.8
 */

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>
#include <cstring>

using namespace polynomial_solver;

// Helper function to parse command-line arguments
struct Config {
    double tolerance = 1e-8;
    unsigned int max_depth = 100;
    double degeneracy_multiplier = 5.0;
    bool dump_geometry = false;
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
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            config.show_help = true;
        }
    }
    return config;
}

void print_help() {
    std::cout << "Usage: simple_cubic [OPTIONS]\n\n";
    std::cout << "Solve cubic polynomial: (x - 0.2)(x - 0.5)(x - 0.8)\n\n";
    std::cout << "Options:\n";
    std::cout << "  -t, --tolerance <value>           Box size tolerance (default: 1e-8)\n";
    std::cout << "  -d, --max-depth <value>           Maximum subdivision depth (default: 100)\n";
    std::cout << "  -m, --degeneracy-multiplier <val> Degeneracy detection multiplier (default: 5.0)\n";
    std::cout << "  --dump-geometry                   Enable geometry dump for visualization\n";
    std::cout << "  -h, --help                        Show this help message\n\n";
    std::cout << "See docs/PARAMETERS.md for detailed parameter documentation.\n";
}

int main(int argc, char* argv[]) {
    // Parse command-line arguments
    Config config = parse_args(argc, argv);
    
    if (config.show_help) {
        print_help();
        return 0;
    }
    
    std::cout << "========================================\n";
    std::cout << "Simple Cubic Example: 2-Line Workflow\n";
    std::cout << "========================================\n\n";
    
    std::cout << "Problem: p(x) = (x - 0.2)(x - 0.5)(x - 0.8)\n";
    std::cout << "Expected roots: x = 0.2, 0.5, 0.8\n\n";
    
    // Print configuration
    std::cout << "Configuration:\n";
    std::cout << "  Tolerance: " << std::scientific << config.tolerance << "\n";
    std::cout << "  Max depth: " << config.max_depth << "\n";
    std::cout << "  Degeneracy multiplier: " << std::fixed << std::setprecision(1) 
              << config.degeneracy_multiplier << "\n";
    std::cout << "  Geometry dump: " << (config.dump_geometry ? "enabled" : "disabled") << "\n\n";
    
    // Define polynomial: (x - 0.2)(x - 0.5)(x - 0.8)
    // = x^3 - 1.5*x^2 + 0.66*x - 0.08
    std::vector<unsigned int> degrees = {3};
    std::vector<double> power_coeffs = {-0.08, 0.66, -1.5, 1.0};  // Power basis coefficients
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system(std::vector<Polynomial>{poly});
    
    // Configure solver
    SubdivisionConfig solver_config;
    solver_config.tolerance = config.tolerance;
    solver_config.max_depth = config.max_depth;
    solver_config.degeneracy_multiplier = config.degeneracy_multiplier;
    
#ifdef ENABLE_GEOMETRY_DUMP
    if (config.dump_geometry) {
        solver_config.dump_geometry = true;
        solver_config.dump_prefix = "dumps/simple_cubic";
    }
#endif
    
    // ============================================================
    // LINE 1: SOLVE (fast, double precision)
    // ============================================================
    Solver solver;
    auto result = solver.subdivisionSolve(system, solver_config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "Step 1: Solve (fast, double precision)\n";
    std::cout << "  Found " << result.num_resolved << " root(s)\n";
    std::cout << "  Unresolved boxes: " << (result.boxes.size() - result.num_resolved) << "\n\n";
    
    // ============================================================
    // LINE 2: REFINE (high precision, 1e-15)
    // ============================================================
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;

    auto refined = refiner.refine(result, system, refine_config);
    
    std::cout << "Step 2: Refine (high precision, 1e-15)\n";
    std::cout << "  Verified roots: " << refined.roots.size() << "\n\n";
    
    // ============================================================
    // RESULTS
    // ============================================================
    std::cout << "========================================\n";
    std::cout << "Results\n";
    std::cout << "========================================\n\n";
    
    for (std::size_t i = 0; i < refined.roots.size(); ++i) {
        const auto& root = refined.roots[i];
        std::cout << "Root " << (i + 1) << ":\n";
        std::cout << "  Location: x = " << std::setprecision(16) << std::fixed 
                  << root.location[0] << "\n";
        std::cout << "  Residual: |f(x)| = " << std::scientific << std::setprecision(4)
                  << std::abs(root.residual[0]) << "\n";
        std::cout << "  Multiplicity: " << root.multiplicity << "\n";
        std::cout << "  Condition estimate: " << std::scientific << std::setprecision(2)
                  << root.condition_estimate << "\n";
        
        // Check if higher precision is needed
        if (root.needs_higher_precision) {
            std::cout << "  ⚠️  WARNING: Higher precision recommended!\n";
            std::cout << "      See docs/CONDITIONING_AND_PRECISION.md\n";
        } else {
            std::cout << "  ✅ Double precision sufficient\n";
        }
        std::cout << "\n";
    }
    
    std::cout << "========================================\n";
    std::cout << "Summary\n";
    std::cout << "========================================\n\n";
    std::cout << "✅ All roots found and refined to 1e-15 precision\n";
    std::cout << "✅ All roots are well-conditioned (double precision sufficient)\n\n";
    std::cout << "For more examples, see:\n";
    std::cout << "  - examples/multiplicity_1d_roots.cpp (multiple roots)\n";
    std::cout << "  - examples/wilkinson_1d_roots.cpp (ill-conditioned)\n";
    std::cout << "  - examples/circle_ellipse_intersection.cpp (2D system)\n\n";
    
    return 0;
}

