/**
 * @file test_multiplicity_methods.cpp
 * @brief Comprehensive test of all multiplicity detection methods
 *
 * Tests all multiplicity detection methods at each Newton iteration
 * to compare robustness and accuracy.
 */

#ifdef ENABLE_HIGH_PRECISION

#include "result_refiner_hp.h"
#include "polynomial_hp.h"
#include "polynomial.h"
#include "differentiation_hp.h"
#include "high_precision_types.h"
#include "precision_context.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>

using namespace polynomial_solver;

// Helper to create (x - r)^m polynomial in power form, then convert to Bernstein
PolynomialHP createMultipleRootPolynomial(double root, unsigned int multiplicity) {
    // Create (x - r)^m in power form
    // Start with coefficients of (x - r)
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

    // Create double-precision polynomial
    std::vector<unsigned int> degrees = {multiplicity};
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);

    // Convert to high precision
    return PolynomialHP(poly);
}

// Test a single polynomial with all methods at each iteration
void testMultiplicityDetection(const std::string& test_name,
                               const PolynomialHP& poly,
                               double initial_guess,
                               unsigned int true_multiplicity,
                               double true_root)
{
    std::cout << "\n========================================\n";
    std::cout << "Test: " << test_name << "\n";
    std::cout << "True multiplicity: " << true_multiplicity << "\n";
    std::cout << "True root: " << true_root << "\n";
    std::cout << "========================================\n\n";
    
    // Perform Newton iterations manually and test at each step
    mpreal x = mpreal(initial_guess);
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    
    std::vector<mpreal> iterates;
    iterates.push_back(x);
    
    unsigned int max_iters = 50;
    
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "Iter | Error      | Taylor | T_1e-50 | T_1e-40 | T_1e-30 | T_1e-20 | Sturm | Ostr\n";
    std::cout << "-----+------------+--------+---------+---------+---------+---------+-------+------\n";

    for (unsigned int iter = 0; iter < max_iters; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);

        // Compute error from true root
        mpreal error = abs(x - mpreal(true_root));

        // Test all methods
        std::map<std::string, unsigned int> estimates;

        // Run all methods
        estimates = ResultRefinerHP::estimateMultiplicityAllMethods(x, poly, 10);

        // Add corrected Ostrowski method (performs 3 regular Newton steps from current point)
        estimates["Ostrowski_corrected"] = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);

        // Print results
        std::cout << std::setw(4) << iter << " | ";
        std::cout << std::scientific << std::setprecision(2) << error << " | ";
        std::cout << std::setw(6) << estimates["Taylor"] << " | ";
        std::cout << std::setw(7) << estimates["Thresh_1e-50"] << " | ";
        std::cout << std::setw(7) << estimates["Thresh_1e-40"] << " | ";
        std::cout << std::setw(7) << estimates["Thresh_1e-30"] << " | ";
        std::cout << std::setw(7) << estimates["Thresh_1e-20"] << " | ";
        std::cout << std::setw(5) << estimates["Sturm"] << " | ";
        std::cout << std::setw(4) << estimates["Ostrowski_corrected"];
        std::cout << "\n";
        
        // Check convergence
        if (error < mpreal("1e-50")) {
            std::cout << "\nConverged to root!\n";
            break;
        }
        
        // Newton step with estimated multiplicity (use Taylor estimate)
        unsigned int m = estimates["Taylor"];
        if (m == 0) m = 1;
        
        mpreal step = mpreal(m) * f / df;
        x = x - step;
        iterates.push_back(x);
    }
}

int main() {
    // Set precision
    PrecisionContext ctx(256);

    std::cout << "========================================\n";
    std::cout << "  Multiplicity Detection Methods Test\n";
    std::cout << "========================================\n";

    // Test 1: Simple root (x - 0.5)
    {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, 1);
        testMultiplicityDetection("Simple root: (x - 0.5)", poly, 0.48, 1, 0.5);
    }

    // Test 2: Double root (x - 0.5)^2
    {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, 2);
        testMultiplicityDetection("Double root: (x - 0.5)^2", poly, 0.48, 2, 0.5);
    }

    // Test 3: Triple root (x - 0.5)^3
    {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, 3);
        testMultiplicityDetection("Triple root: (x - 0.5)^3", poly, 0.48, 3, 0.5);
    }

    // Test 4: Quadruple root (x - 0.5)^4
    {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, 4);
        testMultiplicityDetection("Quadruple root: (x - 0.5)^4", poly, 0.48, 4, 0.5);
    }

    // Test 5: Quintuple root (x - 0.5)^5
    {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, 5);
        testMultiplicityDetection("Quintuple root: (x - 0.5)^5", poly, 0.48, 5, 0.5);
    }

    std::cout << "\n========================================\n";
    std::cout << "  All tests completed!\n";
    std::cout << "========================================\n";

    return 0;
}

#else

#include <iostream>
int main() {
    std::cout << "High precision support not enabled. Skipping test.\n";
    return 0;
}

#endif

