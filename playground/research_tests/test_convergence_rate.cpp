/**
 * @file test_convergence_rate.cpp
 * @brief Test convergence rate of modified Newton method for multiple roots
 * 
 * Verifies that modified Newton achieves quadratic convergence when:
 * 1. Multiplicity is correctly detected
 * 2. Polynomial coefficients are in high precision
 * 3. Sufficient precision is available
 */

#ifdef ENABLE_HIGH_PRECISION

#include "result_refiner_hp.h"
#include "polynomial_hp.h"
#include "polynomial.h"
#include "differentiation_hp.h"
#include "precision_context.h"
#include "precision_conversion.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

using namespace polynomial_solver;

// Create polynomial in HIGH PRECISION
PolynomialHP createMultipleRootPolynomialHP(const mpreal& root, unsigned int multiplicity) {
    std::vector<mpreal> power_coeffs = {-root, mpreal(1)};
    for (unsigned int i = 1; i < multiplicity; ++i) {
        std::vector<mpreal> new_coeffs(power_coeffs.size() + 1, mpreal(0));
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-root);
            new_coeffs[j+1] += power_coeffs[j];
        }
        power_coeffs = new_coeffs;
    }
    std::vector<unsigned int> degrees = {multiplicity};
    return fromPowerHP(degrees, power_coeffs);
}

void testConvergenceRate(unsigned int multiplicity, unsigned int precision_bits) {
    PrecisionContext ctx(precision_bits);
    
    mpreal true_root = mpreal("0.5");
    PolynomialHP poly = createMultipleRootPolynomialHP(true_root, multiplicity);
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    
    // Start from initial guess (further away to see convergence progression)
    mpreal x = mpreal("0.3");
    
    std::cout << "\n========================================\n";
    std::cout << "Multiplicity: " << multiplicity << ", Precision: " << precision_bits << " bits\n";
    std::cout << "========================================\n";
    std::cout << "Iter | Error          | Error Ratio    | Convergence\n";
    std::cout << "-----+----------------+----------------+------------\n";
    
    std::vector<mpreal> errors;
    bool quadratic_convergence = true;
    std::string failure_reason;
    
    for (unsigned int iter = 0; iter < 20; ++iter) {
        mpreal error = abs(x - true_root);
        errors.push_back(error);
        
        std::cout << std::setw(4) << iter << " | ";
        std::cout << std::scientific << std::setprecision(6) << toDouble(error) << " | ";
        
        // Check convergence rate
        if (iter >= 2) {
            // For quadratic convergence: e_{n+1} ≈ C * e_n^2
            // So: e_{n+1} / e_n^2 ≈ C (constant)
            // And: log(e_{n+1}) ≈ 2*log(e_n) + log(C)
            // Convergence order ≈ log(e_{n+1}/e_n) / log(e_n/e_{n-1})
            
            mpreal ratio = errors[iter] / errors[iter-1];
            mpreal prev_ratio = errors[iter-1] / errors[iter-2];
            
            std::cout << std::setw(14) << toDouble(ratio) << " | ";
            
            if (abs(prev_ratio) > mpreal("1e-100") && abs(errors[iter-1]) > mpreal("1e-100")) {
                mpreal order = log(abs(ratio)) / log(abs(prev_ratio));
                std::cout << std::setw(10) << std::setprecision(2) << toDouble(order);
                
                // Check if convergence order is close to 2 (quadratic)
                if (iter >= 3 && abs(errors[iter]) > mpreal("1e-100")) {
                    if (toDouble(order) < 1.5) {
                        quadratic_convergence = false;
                        if (failure_reason.empty()) {
                            failure_reason = "Convergence order < 1.5 at iteration " + std::to_string(iter);
                        }
                    }
                }
            } else {
                std::cout << "    N/A    ";
            }
        } else {
            std::cout << "      N/A      |     N/A    ";
        }
        std::cout << "\n";
        
        // Check if converged
        if (error < mpreal("1e-100")) {
            std::cout << "\nConverged to machine precision!\n";
            break;
        }
        
        // Modified Newton step
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        
        if (abs(df) < mpreal("1e-200")) {
            std::cout << "\nDerivative too small, stopping.\n";
            quadratic_convergence = false;
            failure_reason = "Derivative vanished at iteration " + std::to_string(iter);
            break;
        }
        
        mpreal step = mpreal(multiplicity) * f / df;
        x = x - step;
    }
    
    if (quadratic_convergence) {
        std::cout << "\n✓ QUADRATIC CONVERGENCE VERIFIED\n";
    } else {
        std::cout << "\n✗ QUADRATIC CONVERGENCE FAILED\n";
        std::cout << "Reason: " << failure_reason << "\n";
        std::cout << "*** MARKED FOR FURTHER ANALYSIS ***\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Convergence Rate Test\n";
    std::cout << "  Modified Newton for Multiple Roots\n";
    std::cout << "========================================\n";
    
    // Test different multiplicities with appropriate precision (rule: 128*m bits)
    std::vector<std::pair<unsigned int, unsigned int>> test_cases = {
        {2, 256},    // m=2, 256 bits (128*2)
        {3, 384},    // m=3, 384 bits (128*3)
        {4, 512},    // m=4, 512 bits (128*4)
        {5, 640},    // m=5, 640 bits (128*5)
        {6, 768},    // m=6, 768 bits (128*6)
        {7, 896},    // m=7, 896 bits (128*7)
        {8, 1024},   // m=8, 1024 bits (128*8)
        {10, 1280}   // m=10, 1280 bits (128*10)
    };
    
    for (const auto& test : test_cases) {
        testConvergenceRate(test.first, test.second);
    }
    
    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

