/**
 * @file circle_ellipse_templated.cpp
 * @brief Templated solver + 2D Newton refinement workflow
 * 
 * Problem:
 *   f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)
 *   f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)
 * 
 * Expected root in [0,1]^2: (2/√5, 1/√5) ≈ (0.8944, 0.4472)
 */

#include "core/polynomial_base.h"
#include "solver/solver_base.h"
#include "solver/newton_multidim.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <cassert>

using namespace polynomial_solver;

// Type aliases
using Poly = PolynomialBase<double>;
using System = PolynomialSystemBase<double>;
using Solver = SolverBase<double>;
using SolverConfig = SubdivisionConfigBase<double>;
using NewtonConfig = NewtonMultidimConfig<double>;

struct Config {
    double tolerance = 1e-8;
    unsigned int max_depth = 100;
    double refine_tolerance = 1e-15;
    bool test_mode = false;
    bool show_help = false;
};

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tolerance") == 0) {
            if (i + 1 < argc) config.tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--max-depth") == 0) {
            if (i + 1 < argc) config.max_depth = std::atoi(argv[++i]);
        } else if (strcmp(argv[i], "--refine-tolerance") == 0) {
            if (i + 1 < argc) config.refine_tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--test") == 0) {
            config.test_mode = true;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            config.show_help = true;
        }
    }
    return config;
}

void print_help() {
    std::cout << "Usage: circle_ellipse_templated [OPTIONS]\n\n";
    std::cout << "Solve 2D circle-ellipse intersection using templated solver\n\n";
    std::cout << "Options:\n";
    std::cout << "  -t, --tolerance <value>       Box size tolerance (default: 1e-8)\n";
    std::cout << "  -d, --max-depth <value>       Maximum subdivision depth (default: 100)\n";
    std::cout << "  --refine-tolerance <value>    Newton refinement tolerance (default: 1e-15)\n";
    std::cout << "  --test                        Run in test mode (assertions enabled)\n";
    std::cout << "  -h, --help                    Show this help\n";
}

int main(int argc, char* argv[]) {
    Config config = parse_args(argc, argv);
    if (config.show_help) { print_help(); return 0; }
    
    std::cout << "=== Circle-Ellipse Intersection (Templated Workflow) ===\n\n";
    std::cout << "f1(x,y) = x^2 + y^2 - 1 = 0\n";
    std::cout << "f2(x,y) = x^2/4 + 4*y^2 - 1 = 0\n";
    std::cout << "Expected root: (2/√5, 1/√5) ≈ (0.8944, 0.4472)\n\n";
    
    // Define f1(x,y) = x^2 + y^2 - 1
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> coeffs_f1 = {
        -1.0, 0.0, 1.0,   // 1, y, y^2
         0.0, 0.0, 0.0,   // x, xy, xy^2
         1.0, 0.0, 0.0    // x^2, x^2y, x^2y^2
    };
    Poly f1 = Poly::fromPower(degrees, coeffs_f1);
    
    // Define f2(x,y) = x^2/4 + 4*y^2 - 1
    std::vector<double> coeffs_f2 = {
        -1.0,  0.0, 4.0,   // 1, y, y^2
         0.0,  0.0, 0.0,   // x, xy, xy^2
         0.25, 0.0, 0.0    // x^2, x^2y, x^2y^2
    };
    Poly f2 = Poly::fromPower(degrees, coeffs_f2);
    
    System system({f1, f2});
    
    // Step 1: Solve
    SolverConfig solver_config;
    solver_config.tolerance = config.tolerance;
    solver_config.max_depth = config.max_depth;
    
    Solver solver;
    auto result = solver.subdivisionSolve(system, solver_config, 
                                          RootBoundingMethodBase::ProjectedPolyhedral);
    
    std::cout << "Step 1: Solve\n";
    std::cout << "  Found " << result.num_resolved << " root box(es)\n\n";
    
    // Step 2: Refine using 2D Newton
    NewtonConfig newton_config;
    newton_config.tolerance = config.refine_tolerance;
    newton_config.residual_tolerance = config.refine_tolerance;
    newton_config.max_iterations = 20;
    
    double exact_x = 2.0 / std::sqrt(5.0);
    double exact_y = 1.0 / std::sqrt(5.0);
    
    std::cout << "Step 2: Refine (2D Newton)\n";
    for (std::size_t i = 0; i < result.num_resolved; ++i) {
        double x0 = result.boxes[i].center[0];
        double y0 = result.boxes[i].center[1];
        
        auto refined = refineRoot2D(x0, y0, f1, f2, newton_config);
        
        double error = std::sqrt(std::pow(refined.location[0] - exact_x, 2) + 
                                 std::pow(refined.location[1] - exact_y, 2));
        
        std::cout << "  Root " << (i+1) << ":\n";
        std::cout << "    Location: (" << std::setprecision(16) << refined.location[0]
                  << ", " << refined.location[1] << ")\n";
        std::cout << "    Residual: " << std::scientific << std::setprecision(3) 
                  << refined.residual_norm << "\n";
        std::cout << "    Iterations: " << refined.iterations << "\n";
        std::cout << "    Error from exact: " << error << "\n";
        std::cout << "    " << (refined.converged ? "✓ Converged" : "✗ Not converged") << "\n\n";
        
        if (config.test_mode) {
            assert(refined.converged && "Should converge");
            assert(error < 1e-14 && "Error should be small");
            assert(refined.residual_norm < 1e-14 && "Residual should be small");
        }
    }
    
    if (config.test_mode) {
        assert(result.num_resolved == 1 && "Should find exactly 1 root in [0,1]^2");
        std::cout << "✓ All assertions passed\n";
    }
    
    std::cout << "=== Complete ===\n";
    return 0;
}

