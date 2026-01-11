/**
 * @file test_refiner_wilkinson.cpp
 * @brief Test result refiner on Wilkinson polynomial
 */

#include "refinement/result_refiner.h"
#include "core/polynomial.h"
#include "solver/solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

int main(int argc, char** argv) {
    double tolerance = 1e-8;
    if (argc > 1) {
        tolerance = std::atof(argv[1]);
    }
    
    std::cout << "Wilkinson Polynomial (5 roots) with Result Refiner" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << std::endl;
    
    // p(x) = (x-0.1)(x-0.3)(x-0.5)(x-0.7)(x-0.9)
    std::vector<unsigned int> degrees{5u};
    std::vector<double> power_coeffs{-0.00945, 0.1689, -0.95, 2.3, -2.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    std::cout << "Problem: p(x) = (x-0.1)(x-0.3)(x-0.5)(x-0.7)(x-0.9)" << std::endl;
    std::cout << "Expected: 5 simple roots" << std::endl;
    std::cout << std::endl;
    
    // Solve
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = tolerance;
    config.max_depth = 100;
    
    std::cout << "Step 1: Subdivision solver" << std::endl;
    std::cout << "  Tolerance: " << std::scientific << config.tolerance << std::endl;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "  Found " << result.num_resolved << " resolved boxes" << std::endl;
    std::cout << "  Found " << (result.boxes.size() - result.num_resolved) 
              << " unresolved boxes" << std::endl;
    
    // Show solver box errors
    if (result.num_resolved > 0) {
        double max_solver_error = 0.0;
        for (std::size_t i = 0; i < result.num_resolved; ++i) {
            if (result.boxes[i].max_error[0] > max_solver_error) {
                max_solver_error = result.boxes[i].max_error[0];
            }
        }
        std::cout << "  Max solver box error: " << std::scientific << max_solver_error << std::endl;
    }
    
    // Refine
    std::cout << std::endl;
    std::cout << "Step 2: Newton refinement" << std::endl;
    std::cout << "  Target tolerance: 1e-15" << std::endl;
    std::cout << "  Residual tolerance: 1e-12" << std::endl;
    
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;
    refine_config.max_newton_iters = 50;
    
    RefinementResult refined = refiner.refine(result, system, refine_config);
    
    std::cout << "  Refined to " << refined.roots.size() << " unique roots" << std::endl;
    std::cout << "  Cancelled " << refined.cancelled_boxes.size() << " boxes" << std::endl;
    std::cout << "  Unverified " << refined.unverified_boxes.size() << " boxes" << std::endl;
    
    // Print refined roots
    std::cout << std::endl;
    std::cout << "Refined roots:" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    double expected_roots[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    double max_residual = 0.0;
    double max_error = 0.0;
    
    for (std::size_t i = 0; i < refined.roots.size(); ++i) {
        const auto& root = refined.roots[i];
        double res = std::abs(root.residual[0]);
        if (res > max_residual) max_residual = res;
        
        // Find closest expected root
        double min_dist = 1.0;
        double expected = 0.0;
        for (int j = 0; j < 5; ++j) {
            double dist = std::abs(root.location[0] - expected_roots[j]);
            if (dist < min_dist) {
                min_dist = dist;
                expected = expected_roots[j];
            }
        }
        if (min_dist > max_error) max_error = min_dist;
        
        std::cout << "  Root " << (i+1) << ": x=" << std::fixed << std::setprecision(16) 
                  << root.location[0]
                  << " (expected " << expected << ")" << std::endl;
        std::cout << "          error=" << std::scientific << std::setprecision(2) << min_dist
                  << ", residual=" << res
                  << ", mult=" << root.multiplicity << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Summary:" << std::endl;
    std::cout << "  Roots found: " << refined.roots.size() << " / 5 expected" << std::endl;
    std::cout << "  Max error: " << std::scientific << max_error << std::endl;
    std::cout << "  Max residual: " << max_residual << std::endl;
    std::cout << "  Status: " << (refined.roots.size() == 5 && max_residual < 1e-12 ? "PASS ✓" : "FAIL ✗") 
              << std::endl;
    std::cout << std::endl;
    
    return (refined.roots.size() == 5 && max_residual < 1e-12) ? 0 : 1;
}

