/**
 * @file test_auto_precision_lift.cpp
 * @brief Test automatic precision lifting workflow for isolated multiple roots
 *
 * Workflow:
 * 1. Detect multiplicity in double precision using Ostrowski method
 * 2. Estimate condition number to determine if HP is needed
 * 3. Automatically lift to high precision based on condition number
 * 4. Use modified Newton with detected multiplicity
 * 5. Compute rigorous error bounds
 * 6. Return when error < double precision threshold (1e-15)
 */

#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "core/polynomial.h"
#include "hp/polynomial_hp.h"
#include "solver/solver.h"
#include "refinement/result_refiner.h"
#include "hp/result_refiner_hp.h"
#include "core/differentiation.h"
#include "hp/differentiation_hp.h"
#include "hp/high_precision_types.h"
#include "hp/precision_context.h"
#include "hp/precision_conversion.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

// Helper to create (x - r)^m polynomial in power form
Polynomial createMultipleRootPolynomial(double root, unsigned int multiplicity) {
    // Create (x - r)^m in power form
    std::vector<double> power_coeffs = {-root, 1.0};
    
    // Multiply (m-1) times to get (x - r)^m
    for (unsigned int i = 1; i < multiplicity; ++i) {
        std::vector<double> new_coeffs(power_coeffs.size() + 1, 0.0);
        
        // Multiply by (x - r)
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-root);
            new_coeffs[j+1] += power_coeffs[j];
        }
        
        power_coeffs = new_coeffs;
    }
    
    // Create polynomial in power form
    std::vector<unsigned int> degrees = {multiplicity};
    return Polynomial::fromPower(degrees, power_coeffs);
}

// Test automatic precision lifting for a single isolated multiple root
void testAutoPrecisionLift(unsigned int multiplicity, double root, double initial_guess) {
    std::cout << "\n========================================\n";
    std::cout << "Test: (x - " << root << ")^" << multiplicity << "\n";
    std::cout << "Initial guess: " << initial_guess << "\n";
    std::cout << "========================================\n\n";
    
    // Step 1: Create polynomial in double precision (power form)
    Polynomial poly = createMultipleRootPolynomial(root, multiplicity);
    
    // Step 2: Detect multiplicity in double precision using Ostrowski
    std::cout << "Step 1: Detecting multiplicity in double precision...\n";
    
    ResultRefiner refiner;
    unsigned int detected_mult = refiner.estimateMultiplicityOstrowskiFromPoint(initial_guess, poly);
    
    std::cout << "  True multiplicity: " << multiplicity << "\n";
    std::cout << "  Detected multiplicity: " << detected_mult << "\n";
    
    if (detected_mult != multiplicity) {
        std::cout << "  ⚠ Warning: Multiplicity detection mismatch!\n";
    } else {
        std::cout << "  ✓ Multiplicity correctly detected\n";
    }
    
    // Step 3: Refine in double precision to get initial approximation
    std::cout << "\nStep 2: Refining in double precision...\n";
    
    RefinementConfig config;
    config.max_newton_iters = 50;
    config.target_tolerance = 1e-15;
    config.residual_tolerance = 1e-15;
    
    double refined_x;
    double residual;
    bool converged = refiner.refineRoot1D_fromPoint(initial_guess, poly, config, refined_x, residual);
    
    if (!converged) {
        std::cout << "  ✗ Failed to converge in double precision\n";
        return;
    }
    
    std::cout << "  Refined location: " << std::setprecision(16) << refined_x << "\n";
    std::cout << "  Residual: " << std::scientific << std::setprecision(2) << residual << "\n";
    
    // Step 4: Estimate condition number
    std::cout << "\nStep 3: Estimating condition number...\n";
    
    double first_nonzero_deriv = 0.0;
    std::vector<double> point = {refined_x};
    PolynomialSystem system({poly});

    unsigned int mult_check = refiner.estimateMultiplicity(
        point, system, 10, 1e-10, first_nonzero_deriv);
    
    double condition = refiner.estimateConditionNumber1D(refined_x, poly, first_nonzero_deriv);
    double est_error = condition * std::abs(residual) / std::max(std::abs(first_nonzero_deriv), 1e-14);
    
    std::cout << "  Condition number: " << std::scientific << std::setprecision(2) << condition << "\n";
    std::cout << "  Estimated error: " << est_error << "\n";
    std::cout << "  First nonzero derivative: " << first_nonzero_deriv << "\n";
    
    // Step 5: Decide if high precision is needed
    bool needs_hp = (est_error > 1e-10 || condition > 1e5);
    
    std::cout << "\nStep 4: Precision decision...\n";
    if (!needs_hp) {
        std::cout << "  ✓ Double precision is sufficient (error < 1e-10, condition < 1e5)\n";
        std::cout << "  Final root: " << std::fixed << std::setprecision(16) << refined_x << "\n";
        std::cout << "  Error from true root: " << std::scientific << std::abs(refined_x - root) << "\n";
        return;
    }
    
    std::cout << "  → High precision needed (error > 1e-10 or condition > 1e5)\n";

    // Step 6: Select precision based on condition number
    unsigned int precision_bits;
    if (condition < 1e5) {
        precision_bits = 256;  // ~77 decimal digits
    } else if (condition < 1e10) {
        precision_bits = 512;  // ~154 decimal digits
    } else {
        precision_bits = 1024; // ~308 decimal digits
    }

    std::cout << "  Selected precision: " << precision_bits << " bits\n";

    // Step 7: Lift to high precision and refine with modified Newton
    std::cout << "\nStep 5: Refining in high precision...\n";

    // Set precision (use global function from high_precision_types.h)
    setPrecision(precision_bits);

    // Convert polynomial to high precision
    PolynomialHP poly_hp(poly);

    // Configure HP refinement
    RefinementConfigHP config_hp;
    config_hp.max_newton_iters = 100;

    // Set tolerance based on precision (target: double precision threshold)
    config_hp.target_tolerance_str = "1e-15";
    config_hp.residual_tolerance_str = "1e-15";
    config_hp.max_multiplicity = 10;

    // Refine with high precision using modified Newton
    RefinedRootHP result_hp = ResultRefinerHP::refineRoot1D(refined_x, poly_hp, config_hp);

    if (!result_hp.converged) {
        std::cout << "  ✗ Failed to converge in high precision\n";
        std::cout << "  Error: " << result_hp.error_message << "\n";
        return;
    }

    std::cout << "  ✓ Converged in high precision\n";
    std::cout << "  Iterations: " << result_hp.iterations << "\n";
    std::cout << "  Detected multiplicity: " << result_hp.multiplicity << "\n";

    // Step 8: Check error bounds
    std::cout << "\nStep 6: Verifying error bounds...\n";

    if (result_hp.has_guaranteed_bounds) {
        std::cout << "  ✓ Rigorous error bounds computed\n";
        std::cout << "  Interval: [" << result_hp.interval_lower << ", " << result_hp.interval_upper << "]\n";
        std::cout << "  Max error: " << result_hp.max_error << "\n";

        // Check if error is below double precision threshold
        mpreal error_threshold = mpreal("1e-15");
        if (result_hp.max_error <= error_threshold) {
            std::cout << "  ✓ Error < 1e-15 (double precision threshold)\n";
        } else {
            std::cout << "  ⚠ Error > 1e-15 (may need more iterations)\n";
        }
    } else {
        std::cout << "  ⚠ Could not compute rigorous error bounds\n";
    }

    // Step 9: Convert back to double and report
    std::cout << "\nStep 7: Final result...\n";

    double final_root = toDouble(result_hp.location);
    double true_error = std::abs(final_root - root);

    std::cout << "  Final root: " << std::fixed << std::setprecision(16) << final_root << "\n";
    std::cout << "  True error: " << std::scientific << std::setprecision(2) << true_error << "\n";
    std::cout << "  Condition number: " << result_hp.condition_estimate << "\n";

    // Verify error is below double precision threshold
    if (true_error < 1e-15) {
        std::cout << "  ✓ SUCCESS: Error < 1e-15 (double precision achieved)\n";
    } else {
        std::cout << "  ⚠ Warning: Error > 1e-15\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "Automatic Precision Lifting Test\n";
    std::cout << "========================================\n";
    std::cout << "\nWorkflow:\n";
    std::cout << "1. Detect multiplicity in double precision (Ostrowski)\n";
    std::cout << "2. Estimate condition number\n";
    std::cout << "3. Automatically lift to high precision if needed\n";
    std::cout << "4. Use modified Newton with detected multiplicity\n";
    std::cout << "5. Compute rigorous error bounds\n";
    std::cout << "6. Return when error < 1e-15 (double precision)\n";

    // Test cases: isolated multiple roots with varying multiplicity
    testAutoPrecisionLift(1, 0.5, 0.48);  // Simple root (should not need HP)
    testAutoPrecisionLift(2, 0.5, 0.48);  // Double root
    testAutoPrecisionLift(3, 0.5, 0.48);  // Triple root
    testAutoPrecisionLift(5, 0.5, 0.48);  // Multiplicity 5
    testAutoPrecisionLift(8, 0.5, 0.48);  // Multiplicity 8 (challenging)

    std::cout << "\n========================================\n";
    std::cout << "All tests completed\n";
    std::cout << "========================================\n";

    return 0;
}

#else
#include <iostream>
int main() {
    std::cout << "This test requires ENABLE_HIGH_PRECISION to be enabled.\n";
    return 0;
}
#endif

