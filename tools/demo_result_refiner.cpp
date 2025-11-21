/**
 * @file demo_result_refiner.cpp
 * @brief Demonstration of the result refiner tool
 *
 * Shows how the result refiner consolidates solver results at high precision (1e-15).
 */

#include "result_refiner.h"
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

void print_separator() {
    std::cout << std::string(70, '=') << std::endl;
}

void demo_cubic() {
    std::cout << "\nDemo 1: Cubic Polynomial with Simple Roots" << std::endl;
    print_separator();
    
    // p(x) = (x-0.2)(x-0.5)(x-0.8)
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{-0.08, 0.66, -1.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    std::cout << "Problem: p(x) = (x-0.2)(x-0.5)(x-0.8)" << std::endl;
    std::cout << "Expected: 3 simple roots at x = 0.2, 0.5, 0.8" << std::endl;
    std::cout << std::endl;
    
    // Solve with high precision
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-14;  // Need tight tolerance for 1e-15 verification
    config.max_depth = 100;
    
    std::cout << "Step 1: Solve with tolerance = " << config.tolerance << std::endl;
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "  Solver found " << result.num_resolved << " resolved boxes" << std::endl;
    
    // Show raw solver results
    std::cout << "\n  Raw solver boxes:" << std::endl;
    for (std::size_t i = 0; i < result.num_resolved; ++i) {
        const auto& box = result.boxes[i];
        std::cout << "    Box " << i << ": x=" << std::setprecision(16) << box.center[0]
                  << ", error=" << std::scientific << box.max_error[0]
                  << ", depth=" << box.depth << std::endl;
    }
    
    // Refine results
    std::cout << "\nStep 2: Refine with verification_tolerance = 1e-15" << std::endl;
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.verification_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;
    refine_config.exclusion_multiplier = 3.0;
    
    RefinementResult refined = refiner.refine(result, system, refine_config);
    
    std::cout << "  Refined to " << refined.roots.size() << " unique roots" << std::endl;
    std::cout << "  Cancelled " << refined.cancelled_boxes.size() << " duplicate boxes" << std::endl;
    
    // Show refined results
    std::cout << "\n  Verified roots:" << std::endl;
    for (std::size_t i = 0; i < refined.roots.size(); ++i) {
        const auto& root = refined.roots[i];
        std::cout << "    Root " << i << ": x=" << std::setprecision(16) << root.location[0]
                  << ", multiplicity=" << root.multiplicity
                  << ", residual=" << std::scientific << std::setprecision(2) << root.residual[0]
                  << ", error=" << root.max_error[0]
                  << std::endl;
    }
    
    std::cout << "\n  Result: Successfully consolidated " << result.num_resolved 
              << " boxes → " << refined.roots.size() << " unique roots" << std::endl;
}

void demo_multiplicity() {
    std::cout << "\n\nDemo 2: Polynomial with Multiple Root" << std::endl;
    print_separator();
    
    // p(x) = (x-0.5)^3
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{-0.125, 0.75, -1.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    std::cout << "Problem: p(x) = (x-0.5)^3" << std::endl;
    std::cout << "Expected: 1 triple root at x = 0.5" << std::endl;
    std::cout << std::endl;
    
    // Solve with default precision
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;
    
    std::cout << "Step 1: Solve with tolerance = " << config.tolerance << std::endl;
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "  Solver found " << result.num_resolved << " resolved boxes" << std::endl;
    std::cout << "  Solver found " << (result.boxes.size() - result.num_resolved) 
              << " unresolved boxes (degeneracy)" << std::endl;
    
    // Refine results
    std::cout << "\nStep 2: Refine with verification_tolerance = 1e-15" << std::endl;
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.verification_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;
    refine_config.max_multiplicity = 10;
    
    RefinementResult refined = refiner.refine(result, system, refine_config);
    
    std::cout << "  Refined to " << refined.roots.size() << " verified roots" << std::endl;
    std::cout << "  Unverified boxes: " << refined.unverified_boxes.size() 
              << " (did not meet 1e-15 precision)" << std::endl;
    
    if (refined.roots.size() > 0) {
        std::cout << "\n  Verified roots:" << std::endl;
        for (const auto& root : refined.roots) {
            std::cout << "    x=" << std::setprecision(16) << root.location[0]
                      << ", multiplicity=" << root.multiplicity << std::endl;
        }
    }
    
    std::cout << "\n  Result: Multiple roots are difficult to resolve to 1e-15 precision" << std::endl;
    std::cout << "          Solver stops at degeneracy threshold, not precision threshold" << std::endl;
}

int main() {
    std::cout << "\nResult Refiner Demonstration" << std::endl;
    print_separator();
    std::cout << "This tool demonstrates the result refiner's ability to:" << std::endl;
    std::cout << "  1. Verify roots at high precision (1e-15)" << std::endl;
    std::cout << "  2. Consolidate duplicate boxes" << std::endl;
    std::cout << "  3. Estimate root multiplicity" << std::endl;
    std::cout << "  4. Cancel nearby boxes within exclusion radius" << std::endl;
    
    demo_cubic();
    demo_multiplicity();
    
    std::cout << "\n";
    print_separator();
    std::cout << "Summary:" << std::endl;
    std::cout << "- For simple roots: Solver with tolerance 1e-14 produces boxes" << std::endl;
    std::cout << "  that pass 1e-15 verification" << std::endl;
    std::cout << "- For multiple roots: Degeneracy detection prevents achieving" << std::endl;
    std::cout << "  1e-15 precision (expected behavior)" << std::endl;
    std::cout << std::endl;
    
    return 0;
}

