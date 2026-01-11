/**
 * @file test_convergence_comparison.cpp
 * @brief Compare convergence rates using Taylor vs Ostrowski multiplicity estimates
 *
 * This test runs modified Newton with multiplicity estimates from different methods
 * to see which gives faster convergence in practice.
 */

#ifdef ENABLE_HIGH_PRECISION

#include "hp/result_refiner_hp.h"
#include "hp/polynomial_hp.h"
#include "core/polynomial.h"
#include "hp/differentiation_hp.h"
#include "hp/high_precision_types.h"
#include "hp/precision_context.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

// Helper to create (x - r)^m polynomial
PolynomialHP createMultipleRootPolynomial(double root, unsigned int multiplicity) {
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
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    return PolynomialHP(poly);
}

// Test convergence with a specific multiplicity estimation method
void testConvergenceWithMethod(const std::string& method_name,
                                const PolynomialHP& poly,
                                double initial_guess,
                                unsigned int true_multiplicity,
                                double true_root,
                                bool use_taylor)
{
    mpreal x = mpreal(initial_guess);
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    
    unsigned int max_iters = 50;
    
    std::cout << "\n" << method_name << ":\n";
    std::cout << "Iter | Error          | m_est | Convergence\n";
    std::cout << "-----+----------------+-------+------------\n";
    
    mpreal prev_error = abs(x - mpreal(true_root));
    
    for (unsigned int iter = 0; iter < max_iters; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        mpreal error = abs(x - mpreal(true_root));
        
        // Estimate multiplicity
        unsigned int m_est;
        if (use_taylor) {
            mpreal dummy;
            m_est = ResultRefinerHP::estimateMultiplicity(x, poly, 10, mpreal("1e-50"), dummy);
        } else {
            m_est = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);
        }
        
        if (m_est == 0) m_est = 1;
        
        // Compute convergence rate
        std::string conv_rate = "-";
        if (iter > 0 && prev_error > mpreal("1e-100")) {
            mpreal rate = error / prev_error;
            if (rate < mpreal("1e-10")) {
                conv_rate = "superlinear";
            } else if (rate < mpreal("0.1")) {
                conv_rate = "fast";
            } else if (rate < mpreal("0.9")) {
                conv_rate = "linear";
            } else {
                conv_rate = "slow";
            }
        }
        
        std::cout << std::setw(4) << iter << " | ";
        std::cout << std::scientific << std::setprecision(6) << error << " | ";
        std::cout << std::setw(5) << m_est << " | ";
        std::cout << conv_rate << "\n";
        
        if (error < mpreal("1e-50")) {
            std::cout << "Converged in " << iter + 1 << " iterations!\n";
            break;
        }
        
        // Modified Newton step
        mpreal step = mpreal(m_est) * f / df;
        prev_error = error;
        x = x - step;
    }
}

// Compare methods side by side
void compareConvergence(const std::string& test_name,
                        unsigned int true_multiplicity,
                        double true_root,
                        double initial_guess)
{
    std::cout << "\n========================================\n";
    std::cout << "Test: " << test_name << "\n";
    std::cout << "True multiplicity: " << true_multiplicity << "\n";
    std::cout << "True root: " << true_root << "\n";
    std::cout << "Initial guess: " << initial_guess << "\n";
    std::cout << "========================================\n";
    
    PolynomialHP poly = createMultipleRootPolynomial(true_root, true_multiplicity);
    
    testConvergenceWithMethod("Taylor Method", poly, initial_guess, true_multiplicity, true_root, true);
    testConvergenceWithMethod("Ostrowski Method", poly, initial_guess, true_multiplicity, true_root, false);
}

int main() {
    PrecisionContext ctx(256);
    
    std::cout << "========================================\n";
    std::cout << "  Convergence Rate Comparison\n";
    std::cout << "  Taylor vs Ostrowski\n";
    std::cout << "========================================\n";
    
    // Test different multiplicities
    compareConvergence("Double root: (x - 0.5)^2", 2, 0.5, 0.48);
    compareConvergence("Triple root: (x - 0.5)^3", 3, 0.5, 0.48);
    compareConvergence("Quadruple root: (x - 0.5)^4", 4, 0.5, 0.48);
    compareConvergence("Quintuple root: (x - 0.5)^5", 5, 0.5, 0.48);
    
    std::cout << "\n========================================\n";
    std::cout << "  All tests completed!\n";
    std::cout << "========================================\n";
    
    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

