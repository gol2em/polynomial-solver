#include "result_refiner.h"
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

// Test 1: Simple cubic with 3 roots
int test_cubic_refinement() {
    std::cout << "Test 1: Cubic polynomial refinement" << std::endl;
    
    // p(x) = (x-0.2)(x-0.5)(x-0.8) = x^3 - 1.5*x^2 + 0.66*x - 0.08
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{-0.08, 0.66, -1.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    // Solve with high precision
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-10;
    config.max_depth = 100;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "  Solver found " << result.num_resolved << " boxes" << std::endl;
    
    // Refine results
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;
    
    RefinementResult refined = refiner.refine(result, system, refine_config);
    
    std::cout << "  Refined to " << refined.roots.size() << " unique roots" << std::endl;
    std::cout << "  Cancelled " << refined.cancelled_boxes.size() << " boxes" << std::endl;
    
    // Check results
    if (refined.roots.size() != 3) {
        std::cerr << "  FAIL: Expected 3 unique roots, got " << refined.roots.size() << std::endl;
        return 1;
    }
    
    // Check that all roots are simple (multiplicity 1)
    for (const auto& root : refined.roots) {
        if (root.multiplicity != 1) {
            std::cerr << "  FAIL: Expected multiplicity 1, got " << root.multiplicity << std::endl;
            return 1;
        }
        
        // Check residual
        double max_residual = 0.0;
        for (double r : root.residual) {
            max_residual = std::max(max_residual, std::abs(r));
        }
        
        if (max_residual >= refine_config.residual_tolerance) {
            std::cerr << "  FAIL: Residual too large: " << max_residual << std::endl;
            return 1;
        }
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 2: Multiplicity root
int test_multiplicity_refinement() {
    std::cout << "\nTest 2: Multiplicity polynomial refinement" << std::endl;
    
    // p(x) = (x-0.2)(x-0.6)^6
    std::vector<unsigned int> degrees{7u};
    std::vector<double> power_coeffs{
        -0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0
    };
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    // Solve with moderate precision
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "  Solver found " << result.num_resolved << " resolved boxes" << std::endl;
    std::cout << "  Solver found " << (result.boxes.size() - result.num_resolved) 
              << " unresolved boxes" << std::endl;
    
    // Refine results
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;
    refine_config.max_multiplicity = 10;

    RefinementResult refined = refiner.refine(result, system, refine_config);

    std::cout << "  Refined to " << refined.roots.size() << " unique roots" << std::endl;
    std::cout << "  Cancelled " << refined.cancelled_boxes.size() << " boxes" << std::endl;
    std::cout << "  Unverified boxes: " << refined.unverified_boxes.size() << std::endl;

    // Print details of resolved boxes that didn't pass verification
    for (std::size_t idx : refined.unverified_boxes) {
        if (idx < result.num_resolved) {
            const auto& box = result.boxes[idx];
            std::cout << "  Unverified resolved box at x=" << box.center[0]
                      << ", max_error=" << box.max_error[0]
                      << ", depth=" << box.depth << std::endl;
        }
    }

    // Print verified roots
    for (const auto& root : refined.roots) {
        std::cout << "  Verified root at x=" << root.location[0]
                  << " with multiplicity=" << root.multiplicity
                  << ", residual=" << root.residual[0] << std::endl;
    }

    // Check that we found at least the simple root at 0.2
    bool found_simple = false;
    for (const auto& root : refined.roots) {
        double dist_to_02 = std::abs(root.location[0] - 0.2);
        if (dist_to_02 < 1e-6) {
            found_simple = true;
            if (root.multiplicity != 1) {
                std::cerr << "  FAIL: Root at 0.2 should have multiplicity 1, got "
                          << root.multiplicity << std::endl;
                return 1;
            }
        }
    }

    if (!found_simple) {
        std::cerr << "  FAIL: Did not find simple root at 0.2 with precision 1e-15" << std::endl;
        std::cerr << "  This indicates the solver needs higher precision for simple roots in multiplicity polynomials" << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 3: Multiplicity detection
int test_multiplicity_detection() {
    std::cout << "\nTest 3: Multiplicity detection from derivatives" << std::endl;

    ResultRefiner refiner;

    // Test case 1: Simple root at x=0.5
    // p(x) = (x - 0.5) = x - 0.5
    {
        std::vector<unsigned int> degrees{1u};
        std::vector<double> power_coeffs{-0.5, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        PolynomialSystem system({p});

        std::vector<double> point{0.5};
        double first_deriv = 0.0;
        unsigned int mult = refiner.estimateMultiplicity(point, system, 10, 1e-10, first_deriv);

        std::cout << "  Simple root (x-0.5): multiplicity = " << mult;
        if (mult == 1) {
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗ (expected 1)" << std::endl;
            return 1;
        }
    }

    // Test case 2: Double root at x=0.5
    // p(x) = (x - 0.5)^2 = x^2 - x + 0.25
    {
        std::vector<unsigned int> degrees{2u};
        std::vector<double> power_coeffs{0.25, -1.0, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        PolynomialSystem system({p});

        std::vector<double> point{0.5};
        double first_deriv = 0.0;
        unsigned int mult = refiner.estimateMultiplicity(point, system, 10, 1e-10, first_deriv);

        std::cout << "  Double root (x-0.5)^2: multiplicity = " << mult;
        if (mult == 2) {
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗ (expected 2)" << std::endl;
            return 1;
        }
    }

    // Test case 3: Triple root at x=0.5
    // p(x) = (x - 0.5)^3 = x^3 - 1.5x^2 + 0.75x - 0.125
    {
        std::vector<unsigned int> degrees{3u};
        std::vector<double> power_coeffs{-0.125, 0.75, -1.5, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        PolynomialSystem system({p});

        std::vector<double> point{0.5};
        double first_deriv = 0.0;
        unsigned int mult = refiner.estimateMultiplicity(point, system, 10, 1e-10, first_deriv);

        std::cout << "  Triple root (x-0.5)^3: multiplicity = " << mult;
        if (mult == 3) {
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗ (expected 3)" << std::endl;
            return 1;
        }
    }

    // Test case 4: Quadruple root at x=0.3
    // p(x) = (x - 0.3)^4
    {
        std::vector<unsigned int> degrees{4u};
        std::vector<double> power_coeffs{0.0081, -0.108, 0.54, -1.2, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        PolynomialSystem system({p});

        std::vector<double> point{0.3};
        double first_deriv = 0.0;
        unsigned int mult = refiner.estimateMultiplicity(point, system, 10, 1e-10, first_deriv);

        std::cout << "  Quadruple root (x-0.3)^4: multiplicity = " << mult;
        if (mult == 4) {
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗ (expected 4)" << std::endl;
            return 1;
        }
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

int main() {
    std::cout << "Result Refiner Tests" << std::endl;
    std::cout << "====================" << std::endl;

    int failures = 0;
    failures += test_cubic_refinement();
    failures += test_multiplicity_refinement();
    failures += test_multiplicity_detection();

    std::cout << "\n====================" << std::endl;
    if (failures == 0) {
        std::cout << "All tests passed!" << std::endl;
    } else {
        std::cout << failures << " test(s) failed!" << std::endl;
    }

    return failures;
}

