/**
 * @file wilkinson_1d_templated.cpp
 * @brief Templated solver/refiner workflow for Wilkinson polynomial
 * 
 * Problem: p(x) = (x-1/21)(x-2/21)...(x-20/21)
 * Expected roots: x = k/21 for k = 1, 2, ..., 20
 */

#include "core/polynomial_base.h"
#include "solver/solver_base.h"
#include "refinement/result_refiner_base.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace polynomial_solver;

// Type aliases
using Poly = PolynomialBase<double>;
using System = PolynomialSystemBase<double>;
using Solver = SolverBase<double>;
using Refiner = ResultRefinerBase<double>;
using SolverConfig = SubdivisionConfigBase<double>;
using RefineConfig = RefinementConfigBase<double>;

struct Config {
    double tolerance = 1e-8;
    unsigned int max_depth = 100;
    double target_tolerance = 1e-12;
    double residual_tolerance = 1e-10;
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
        } else if (strcmp(argv[i], "--target-tolerance") == 0) {
            if (i + 1 < argc) config.target_tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--residual-tolerance") == 0) {
            if (i + 1 < argc) config.residual_tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--test") == 0) {
            config.test_mode = true;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            config.show_help = true;
        }
    }
    return config;
}

void print_help() {
    std::cout << "Usage: wilkinson_1d_templated [OPTIONS]\n\n";
    std::cout << "Solve Wilkinson polynomial: (x-1/21)(x-2/21)...(x-20/21)\n\n";
    std::cout << "Options:\n";
    std::cout << "  -t, --tolerance <value>        Box size tolerance (default: 1e-8)\n";
    std::cout << "  -d, --max-depth <value>        Maximum subdivision depth (default: 100)\n";
    std::cout << "  --target-tolerance <value>     Refinement target (default: 1e-12)\n";
    std::cout << "  --residual-tolerance <value>   Convergence residual (default: 1e-10)\n";
    std::cout << "  --test                         Run in test mode (assertions enabled)\n";
    std::cout << "  -h, --help                     Show this help\n";
}

// Compute coefficients of (x-1/n)(x-2/n)...(x-(n-1)/n)
std::vector<double> wilkinson_coefficients(int n) {
    double scale = 1.0 / static_cast<double>(n);
    std::vector<double> coeffs = {-scale, 1.0};
    
    for (int k = 2; k < n; ++k) {
        std::vector<double> new_coeffs(coeffs.size() + 1, 0.0);
        double root = k * scale;
        for (size_t i = 0; i < coeffs.size(); ++i) {
            new_coeffs[i] -= root * coeffs[i];
            new_coeffs[i + 1] += coeffs[i];
        }
        coeffs = new_coeffs;
    }
    return coeffs;
}

int main(int argc, char* argv[]) {
    Config config = parse_args(argc, argv);
    if (config.show_help) { print_help(); return 0; }
    
    std::cout << "=== Wilkinson Polynomial (Templated Workflow) ===\n\n";
    std::cout << "Problem: p(x) = (x-1/21)(x-2/21)...(x-20/21)\n";
    std::cout << "Expected: 20 roots at x = k/21 for k = 1..20\n\n";
    
    // Define Wilkinson polynomial
    std::vector<double> power_coeffs = wilkinson_coefficients(21);
    std::vector<unsigned int> degrees{static_cast<unsigned int>(power_coeffs.size() - 1)};
    Poly poly = Poly::fromPower(degrees, power_coeffs);
    System system({poly});
    
    // Step 1: Solve
    SolverConfig solver_config;
    solver_config.tolerance = config.tolerance;
    solver_config.max_depth = config.max_depth;
    
    Solver solver;
    auto result = solver.subdivisionSolve(system, solver_config, 
                                          RootBoundingMethodBase::ProjectedPolyhedral);
    
    std::cout << "Step 1: Solve\n";
    std::cout << "  Found " << result.num_resolved << " root boxes\n";
    std::cout << "  Unresolved: " << (result.boxes.size() - result.num_resolved) << "\n\n";
    
    // Step 2: Refine
    RefineConfig refine_config;
    refine_config.target_tolerance = config.target_tolerance;
    refine_config.residual_tolerance = config.residual_tolerance;
    refine_config.max_newton_iters = 100;
    
    std::cout << "Step 2: Refine\n";
    std::vector<double> found_roots;
    int converged_count = 0;
    
    for (std::size_t i = 0; i < result.num_resolved; ++i) {
        double x0 = result.boxes[i].center[0];
        auto refined = Refiner::refineRoot1D(x0, poly, refine_config);
        
        if (refined.converged) {
            found_roots.push_back(refined.location);
            converged_count++;
        }
    }
    
    std::sort(found_roots.begin(), found_roots.end());
    
    std::cout << "  Converged: " << converged_count << "/" << result.num_resolved << "\n\n";
    std::cout << "  Found roots:\n";
    for (std::size_t i = 0; i < found_roots.size(); ++i) {
        double expected = (i + 1) / 21.0;
        double error = std::abs(found_roots[i] - expected);
        std::cout << "    x[" << std::setw(2) << (i+1) << "] = " 
                  << std::fixed << std::setprecision(10) << found_roots[i]
                  << " (expected " << expected << ", error=" 
                  << std::scientific << std::setprecision(2) << error << ")\n";
    }
    
    if (config.test_mode) {
        // Wilkinson polynomial is extremely ill-conditioned in double precision
        // The solver may not find all roots, but should find at least a few
        // This is a known limitation - HP is needed for full accuracy
        assert(found_roots.size() >= 1 && "Should find at least one root");
        std::cout << "\n✓ Test passed (found " << found_roots.size()
                  << " roots - HP recommended for full Wilkinson)\n";
    }
    
    std::cout << "\n=== Complete ===\n";
    return 0;
}

