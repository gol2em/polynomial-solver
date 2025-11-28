/**
 * @file high_precision_example.cpp
 * @brief Example demonstrating high-precision refinement for ill-conditioned roots
 * 
 * This example shows the complete workflow:
 * 1. Solve with double precision (fast)
 * 2. Refine with double precision
 * 3. Detect ill-conditioned roots
 * 4. Refine ill-conditioned roots with high precision
 * 
 * Problem: Wilkinson-like polynomial with closely-spaced roots
 * Expected: Double precision fails, high precision succeeds
 * 
 * Build:
 *   cmake -DENABLE_HIGH_PRECISION=ON ..
 *   make example_high_precision
 */

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>

#ifdef ENABLE_HIGH_PRECISION
#include "high_precision_refiner.h"
#endif

using namespace polynomial_solver;

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "High-Precision Refinement Example" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;

    // Create an ill-conditioned polynomial: (x-0.5)^5
    // This has a root of multiplicity 5 at x=0.5
    // Condition number is very high near this root
    std::vector<unsigned int> degrees{5};
    std::vector<double> power_coeffs{
        -0.03125,   // (0.5)^5
         0.15625,   // -5*(0.5)^4
        -0.3125,    // 10*(0.5)^3
         0.3125,    // -10*(0.5)^2
        -0.15625,   // 5*(0.5)
         0.03125    // -1
    };
    
    // Expand to get actual coefficients
    // (x-0.5)^5 = x^5 - 2.5x^4 + 2.5x^3 - 1.25x^2 + 0.3125x - 0.03125
    power_coeffs[0] = -0.03125;
    power_coeffs[1] = 0.3125;
    power_coeffs[2] = -1.25;
    power_coeffs[3] = 2.5;
    power_coeffs[4] = -2.5;
    power_coeffs[5] = 1.0;
    
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({poly});

    std::cout << "Polynomial: (x - 0.5)^5" << std::endl;
    std::cout << "Expected root: x = 0.5 (multiplicity 5)" << std::endl;
    std::cout << "Condition number: Very high (ill-conditioned)" << std::endl;
    std::cout << std::endl;

    // ========================================
    // STEP 1: Solve with double precision
    // ========================================
    std::cout << "=== Step 1: Solve with double precision ===" << std::endl;
    
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;
    
    auto result = solver.subdivisionSolve(system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "Resolved boxes: " << result.num_resolved << std::endl;
    std::cout << "Unresolved boxes: " << (result.boxes.size() - result.num_resolved) << std::endl;
    std::cout << std::endl;

    // ========================================
    // STEP 2: Refine with double precision
    // ========================================
    std::cout << "=== Step 2: Refine with double precision ===" << std::endl;
    
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-15;
    
    auto refined = refiner.refine(result, system, refine_config);
    
    std::cout << "Refined roots: " << refined.roots.size() << std::endl;
    
    for (size_t i = 0; i < refined.roots.size(); ++i) {
        const auto& root = refined.roots[i];
        std::cout << "  Root " << (i+1) << ":" << std::endl;
        std::cout << "    Location: " << std::fixed << std::setprecision(15) 
                  << root.location[0] << std::endl;
        std::cout << "    Residual: " << std::scientific << std::setprecision(3)
                  << root.residual[0] << std::endl;
        std::cout << "    Multiplicity: " << root.multiplicity << std::endl;
        std::cout << "    Condition estimate: " << root.condition_estimate << std::endl;
        std::cout << "    Needs higher precision: " 
                  << (root.needs_higher_precision ? "YES" : "NO") << std::endl;
    }
    std::cout << std::endl;

    // ========================================
    // STEP 3: High-precision refinement
    // ========================================
#ifdef ENABLE_HIGH_PRECISION
    std::cout << "=== Step 3: High-precision refinement ===" << std::endl;
    std::cout << "High-precision support: ENABLED" << std::endl;
    std::cout << std::endl;
    
    using namespace high_precision;
    
    for (size_t i = 0; i < refined.roots.size(); ++i) {
        const auto& root = refined.roots[i];
        
        if (root.needs_higher_precision) {
            std::cout << "Refining root " << (i+1) << " with high precision..." << std::endl;
            
            // Convert polynomial to high precision
            auto coeffs_hp = convertCoefficientsToHighPrecision(poly, 256);
            
            // Configure high-precision refinement
            HighPrecisionConfig hp_config;
            hp_config.precision_bits = 256;  // ~77 decimal digits
            hp_config.target_tolerance = "1e-50";
            hp_config.residual_tolerance = "1e-50";
            hp_config.max_newton_iters = 100;
            
            // Refine with high precision
            auto root_hp = refineRootHighPrecision(
                root.location,
                {coeffs_hp},
                {poly.degrees()},
                hp_config);
            
            std::cout << "  High-precision result:" << std::endl;
            std::cout << "    Location: " << root_hp.location << std::endl;
            std::cout << "    Residual: " << root_hp.residual << std::endl;
            std::cout << "    Precision: " << root_hp.precision_bits << " bits ("
                      << root_hp.decimal_digits << " digits)" << std::endl;
            std::cout << "    Iterations: " << root_hp.iterations << std::endl;
            std::cout << "    Verified: " << (root_hp.verified ? "YES" : "NO") << std::endl;
            
            if (!root_hp.verified) {
                std::cout << "    Error: " << root_hp.error_message << std::endl;
            }
            std::cout << std::endl;
        }
    }
#else
    std::cout << "=== Step 3: High-precision refinement ===" << std::endl;
    std::cout << "High-precision support: DISABLED" << std::endl;
    std::cout << "Rebuild with -DENABLE_HIGH_PRECISION=ON to enable" << std::endl;
    std::cout << std::endl;
#endif

    std::cout << "========================================" << std::endl;
    std::cout << "Example completed!" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}

