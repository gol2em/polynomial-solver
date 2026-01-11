/**
 * @file simple_cubic_templated.cpp
 * @brief Templated solver/refiner workflow for simple cubic
 *
 * Demonstrates the batch refinement workflow with automatic merging
 * of nearby roots and degeneracy detection.
 *
 * Problem: p(x) = (x - 0.2)(x - 0.5)(x - 0.8)
 * Expected roots: x = 0.2, 0.5, 0.8
 */

#include "core/polynomial_base.h"
#include "solver/solver_base.h"
#include "refinement/result_refiner_base.h"
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
using Refiner = ResultRefinerBase<double>;
using SolverConfig = SubdivisionConfigBase<double>;
using RefineConfig = RefinementConfigBase<double>;
using BatchResult = RefinementResultBase<double>;

struct Config {
    double tolerance = 1e-8;
    unsigned int max_depth = 100;
    double target_tolerance = 1e-15;
    double residual_tolerance = 1e-15;
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
    std::cout << "Usage: simple_cubic_templated [OPTIONS]\n\n";
    std::cout << "Solve cubic polynomial using templated solver: (x - 0.2)(x - 0.5)(x - 0.8)\n\n";
    std::cout << "Options:\n";
    std::cout << "  -t, --tolerance <value>        Box size tolerance (default: 1e-8)\n";
    std::cout << "  -d, --max-depth <value>        Maximum subdivision depth (default: 100)\n";
    std::cout << "  --target-tolerance <value>     Refinement target (default: 1e-15)\n";
    std::cout << "  --residual-tolerance <value>   Convergence residual (default: 1e-15)\n";
    std::cout << "  --test                         Run in test mode (assertions enabled)\n";
    std::cout << "  -h, --help                     Show this help\n";
}

int main(int argc, char* argv[]) {
    Config config = parse_args(argc, argv);
    if (config.show_help) { print_help(); return 0; }
    
    std::cout << "=== Simple Cubic (Templated Workflow) ===\n\n";
    std::cout << "Problem: p(x) = (x - 0.2)(x - 0.5)(x - 0.8)\n";
    std::cout << "Expected roots: 0.2, 0.5, 0.8\n\n";
    
    // Define polynomial: (x - 0.2)(x - 0.5)(x - 0.8) = x^3 - 1.5*x^2 + 0.66*x - 0.08
    std::vector<unsigned int> degrees = {3};
    std::vector<double> power_coeffs = {-0.08, 0.66, -1.5, 1.0};
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
    if (result.degeneracy_detected) {
        std::cout << "  Warning: Degeneracy detected during solving\n";
    }
    std::cout << "\n";

    // Step 2: Batch refine with automatic merging
    RefineConfig refine_config;
    refine_config.target_tolerance = config.target_tolerance;
    refine_config.residual_tolerance = config.residual_tolerance;

    std::cout << "Step 2: Batch Refine (with automatic merging)\n";
    BatchResult batch = Refiner::refine(result, poly, refine_config);

    std::cout << "  Unique roots found: " << batch.roots.size() << "\n";
    std::cout << "  Cancelled (merged) boxes: " << batch.cancelled_boxes.size() << "\n";
    std::cout << "  Unverified boxes: " << batch.unverified_boxes.size() << "\n";
    std::cout << "  Problematic regions: " << batch.problematic_regions.size() << "\n";
    if (batch.any_needs_higher_precision) {
        std::cout << "  Warning: Some roots may need higher precision\n";
    }
    std::cout << "\n";

    std::vector<double> expected_roots = {0.2, 0.5, 0.8};
    std::vector<double> found_roots;

    std::cout << "  Refined roots:\n";
    for (std::size_t i = 0; i < batch.roots.size(); ++i) {
        const auto& root = batch.roots[i];
        std::cout << "    Root " << (i+1) << ": x = " << std::setprecision(16) << root.location;
        std::cout << " |f(x)| = " << std::scientific << std::setprecision(2)
                  << std::abs(root.residual);
        std::cout << " mult=" << root.multiplicity;
        std::cout << " (merged " << root.source_boxes.size() << " box(es))";
        std::cout << (root.converged ? " ✓" : " ✗");
        if (root.needs_higher_precision) std::cout << " [HP]";
        std::cout << "\n";

        found_roots.push_back(root.location);

        if (config.test_mode) {
            assert(root.converged && "Root should converge");
            assert(std::abs(root.residual) < 1e-14 && "Residual should be small");
        }
    }

    // Verify all expected roots found
    if (config.test_mode) {
        std::sort(found_roots.begin(), found_roots.end());

        assert(found_roots.size() == 3 && "Should find 3 unique roots");
        for (std::size_t i = 0; i < 3; ++i) {
            assert(std::abs(found_roots[i] - expected_roots[i]) < 1e-10
                   && "Root location should match expected");
        }
        std::cout << "\n✓ All assertions passed\n";
    }

    std::cout << "\n=== Complete ===\n";
    return 0;
}

