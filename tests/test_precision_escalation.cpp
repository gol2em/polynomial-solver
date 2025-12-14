/**
 * @file test_precision_escalation.cpp
 * @brief Test automatic precision escalation for ill-conditioned roots
 */

#ifdef ENABLE_HIGH_PRECISION

#include "result_refiner.h"
#include "polynomial.h"
#include "precision_context.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace polynomial_solver;

Polynomial createMultipleRootPolynomial(double root, unsigned int multiplicity) {
    std::vector<double> power_coeffs = {-root, 1.0};
    for (unsigned int i = 1; i < multiplicity; ++i) {
        std::vector<double> new_coeffs(power_coeffs.size() + 1, 0.0);
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-root);
            new_coeffs[j+1] += power_coeffs[j];
        }
        power_coeffs = new_coeffs;
    }
    std::vector<unsigned int> degrees = {multiplicity};
    return Polynomial::fromPower(degrees, power_coeffs);
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Precision Escalation Test\n";
    std::cout << "========================================\n";
    
    ResultRefiner refiner;
    RefinementConfig config;
    config.target_tolerance = 1e-15;
    config.residual_tolerance = 1e-15;
    config.max_newton_iters = 50;
    config.max_multiplicity = 15;
    
    std::cout << "\nTesting roots with different multiplicities:\n";
    std::cout << "m | Double OK | HP Used | Final Error | Multiplicity\n";
    std::cout << "--+-----------+---------+-------------+-------------\n";
    
    for (unsigned int m = 2; m <= 8; ++m) {
        Polynomial poly = createMultipleRootPolynomial(0.5, m);
        double initial_guess = 0.48;

        std::cout << "\n--- Testing m=" << m << " ---\n";

        RefinedRoot result;
        bool success = refiner.refineRoot1DWithPrecisionEscalation(
            initial_guess, poly, config, result);

        std::cout << std::setw(2) << m << " | ";
        std::cout << std::setw(9) << (result.needs_higher_precision ? "No" : "Yes") << " | ";
        std::cout << std::setw(7) << (result.needs_higher_precision ? "Yes" : "No") << " | ";

        if (success && !result.max_error.empty()) {
            std::cout << std::scientific << std::setprecision(2) << result.max_error[0] << " | ";
        } else {
            std::cout << "N/A         | ";
        }

        std::cout << std::setw(12) << result.multiplicity << "\n";
    }
    
    // Test a very ill-conditioned case
    std::cout << "\n=== Testing very high multiplicity (m=10) ===\n";
    Polynomial poly10 = createMultipleRootPolynomial(0.5, 10);
    RefinedRoot result10;
    
    bool success = refiner.refineRoot1DWithPrecisionEscalation(
        0.48, poly10, config, result10);
    
    if (success) {
        std::cout << "Success!\n";
        std::cout << "  Root location: " << std::setprecision(15) << result10.location[0] << "\n";
        std::cout << "  Multiplicity: " << result10.multiplicity << "\n";
        std::cout << "  Condition: " << std::scientific << result10.condition_estimate << "\n";
        if (!result10.max_error.empty()) {
            std::cout << "  Error bound: " << result10.max_error[0] << "\n";
        }
    } else {
        std::cout << "Failed to refine even with high precision\n";
    }
    
    // Test cluster of roots (ill-conditioned)
    std::cout << "\n=== Testing close roots (ill-conditioned) ===\n";
    // (x-0.5)*(x-0.50001) - two very close simple roots
    std::vector<double> close_roots_coeffs = {0.5 * 0.50001, -(0.5 + 0.50001), 1.0};
    std::vector<unsigned int> degrees = {2};
    Polynomial close_poly = Polynomial::fromPower(degrees, close_roots_coeffs);
    
    RefinedRoot result_close;
    success = refiner.refineRoot1DWithPrecisionEscalation(
        0.499, close_poly, config, result_close);
    
    if (success) {
        std::cout << "Success!\n";
        std::cout << "  Root location: " << std::setprecision(15) << result_close.location[0] << "\n";
        std::cout << "  Multiplicity: " << result_close.multiplicity << "\n";
        std::cout << "  Condition: " << std::scientific << result_close.condition_estimate << "\n";
        std::cout << "  HP used: " << (result_close.needs_higher_precision ? "No (double OK)" : "Yes") << "\n";
    } else {
        std::cout << "Failed\n";
    }
    
    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    std::cout << "Rebuild with -DENABLE_HIGH_PRECISION to test precision escalation\n";
    return 1;
}
#endif

