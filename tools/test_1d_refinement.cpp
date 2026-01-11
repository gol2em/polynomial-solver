/**
 * @file test_1d_refinement.cpp
 * @brief Test systematic 1D root refinement on all examples
 */

#include "refinement/result_refiner.h"
#include "core/polynomial.h"
#include "solver/solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

void test_polynomial(const std::string& name,
                     const Polynomial& poly,
                     double solver_tolerance = 1e-8) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Test: " << name << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    PolynomialSystem system({poly});
    
    // Step 1: Solve with moderate precision
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = solver_tolerance;
    config.max_depth = 100;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "Solver (tolerance=" << solver_tolerance << "):" << std::endl;
    std::cout << "  Resolved boxes: " << result.num_resolved << std::endl;
    std::cout << "  Unresolved boxes: " << (result.boxes.size() - result.num_resolved) << std::endl;
    
    // Print resolved boxes
    for (std::size_t i = 0; i < result.num_resolved; ++i) {
        const auto& box = result.boxes[i];
        std::cout << "  Box " << i << ": x ∈ [" << box.lower[0] << ", " << box.upper[0] 
                  << "], center=" << box.center[0]
                  << ", max_error=" << box.max_error[0]
                  << ", depth=" << box.depth << std::endl;
    }
    
    // Step 2: Refine to high precision
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;
    refine_config.max_multiplicity = 10;
    refine_config.exclusion_multiplier = 3.0;
    
    RefinementResult refined = refiner.refine(result, system, refine_config);
    
    std::cout << "\nRefinement Results:" << std::endl;
    std::cout << "  Verified roots: " << refined.roots.size() << std::endl;
    std::cout << "  Cancelled boxes: " << refined.cancelled_boxes.size() << std::endl;
    std::cout << "  Unverified boxes: " << refined.unverified_boxes.size() << std::endl;
    
    // Print verified roots with details
    for (std::size_t i = 0; i < refined.roots.size(); ++i) {
        const auto& root = refined.roots[i];
        std::cout << "\n  Root " << i << ":" << std::endl;
        std::cout << "    Location: x = " << std::setprecision(16) << root.location[0] << std::endl;
        std::cout << "    Residual: |f(x)| = " << std::scientific << root.residual[0] << std::endl;
        std::cout << "    Multiplicity: " << root.multiplicity << std::endl;
        std::cout << "    First non-zero derivative: " << root.first_nonzero_derivative << std::endl;
        std::cout << "    Max error: " << root.max_error[0] << std::endl;
        std::cout << "    Source boxes: " << root.source_boxes.size() << std::endl;
        std::cout << "    Depth: " << root.depth << std::endl;
        std::cout << "    Exclusion radius: " << root.exclusion_radius << std::endl;
    }
    
    // Print unverified boxes
    if (!refined.unverified_boxes.empty()) {
        std::cout << "\n  Unverified boxes:" << std::endl;
        for (std::size_t idx : refined.unverified_boxes) {
            const auto& box = result.boxes[idx];
            std::cout << "    Box " << idx << ": x ∈ [" << box.lower[0] << ", " << box.upper[0] 
                      << "], center=" << box.center[0] << std::endl;
        }
    }
}

int main() {
    std::cout << "Systematic 1D Root Refinement Test" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Test 1: Cubic polynomial with 3 simple roots
    // p(x) = (x - 0.2)(x - 0.5)(x - 0.8)
    {
        std::vector<unsigned int> degrees{3u};
        std::vector<double> power_coeffs{-0.08, 0.66, -1.5, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        test_polynomial("Cubic: (x-0.2)(x-0.5)(x-0.8)", p);
    }
    
    // Test 2: Multiplicity polynomial
    // p(x) = (x-0.2)(x-0.6)^6
    {
        std::vector<unsigned int> degrees{7u};
        std::vector<double> power_coeffs{
            -0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0
        };
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        test_polynomial("Multiplicity: (x-0.2)(x-0.6)^6", p);
    }
    
    // Test 3: Wilkinson-like polynomial with 5 roots
    {
        std::vector<unsigned int> degrees{5u};
        std::vector<double> power_coeffs{0.00378, -0.0945, 0.8925, -3.95, 7.5, -5.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        test_polynomial("Wilkinson-5: roots at 0.1, 0.3, 0.5, 0.7, 0.9", p);
    }
    
    // Test 4: Double root
    // p(x) = (x-0.5)^2
    {
        std::vector<unsigned int> degrees{2u};
        std::vector<double> power_coeffs{0.25, -1.0, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        test_polynomial("Double root: (x-0.5)^2", p);
    }
    
    // Test 5: Triple root
    // p(x) = (x-0.5)^3
    {
        std::vector<unsigned int> degrees{3u};
        std::vector<double> power_coeffs{-0.125, 0.75, -1.5, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        test_polynomial("Triple root: (x-0.5)^3", p);
    }
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "All tests completed!" << std::endl;
    
    return 0;
}

