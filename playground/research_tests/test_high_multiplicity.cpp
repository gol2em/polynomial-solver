/**
 * @file test_high_multiplicity.cpp
 * @brief Test limits of multiplicity detection methods
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

using namespace polynomial_solver;

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

void testOstrowskiDetail(unsigned int m) {
    std::cout << "\n========== Ostrowski Detail for m=" << m << " ==========\n";
    
    PolynomialHP poly = createMultipleRootPolynomial(0.5, m);
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    
    mpreal x = mpreal(0.48);
    mpreal true_root = mpreal(0.5);
    
    std::vector<mpreal> iterates;
    iterates.push_back(x);
    
    for (int i = 0; i < 5; ++i) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        if (abs(df) < mpreal("1e-100")) break;
        x = x - f / df;
        iterates.push_back(x);
    }
    
    if (iterates.size() >= 4) {
        mpreal x1 = iterates[1];
        mpreal x2 = iterates[2];
        mpreal x3 = iterates[3];
        
        mpreal num = x1 - x2;
        mpreal denom = x3 - mpreal(2)*x2 + x1;
        mpreal p = mpreal("0.5") + num / denom;
        
        std::cout << "x1 = " << x1 << "\n";
        std::cout << "x2 = " << x2 << "\n";
        std::cout << "x3 = " << x3 << "\n";
        std::cout << "numerator = " << num << "\n";
        std::cout << "denominator = " << denom << "\n";
        std::cout << "p = " << p << "\n";
        std::cout << "floor(p) = " << static_cast<int>(floor(p)) << "\n";
        
        // Check error ratios
        mpreal e1 = abs(x1 - true_root);
        mpreal e2 = abs(x2 - true_root);
        mpreal e3 = abs(x3 - true_root);
        std::cout << "e1 = " << e1 << "\n";
        std::cout << "e2 = " << e2 << "\n";
        std::cout << "e3 = " << e3 << "\n";
        std::cout << "e2/e1 = " << (e2/e1) << " (should be " << (m-1.0)/m << ")\n";
        std::cout << "e3/e2 = " << (e3/e2) << " (should be " << (m-1.0)/m << ")\n";
    }
}

int main() {
    PrecisionContext ctx(256);
    
    std::cout << "========================================\n";
    std::cout << "  High Multiplicity Detection Test\n";
    std::cout << "========================================\n";
    
    // Test with different precisions
    std::cout << "\n=== Testing with 256-bit precision ===\n";
    std::cout << "True m | Taylor(max=15) | Taylor(max=20) | Ostrowski\n";
    std::cout << "-------+----------------+----------------+-----------\n";
    
    for (unsigned int m = 2; m <= 15; ++m) {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, m);
        mpreal x = mpreal(0.48);
        
        mpreal first_nonzero;
        unsigned int taylor15 = ResultRefinerHP::estimateMultiplicity(x, poly, 15, mpreal("1e-50"), first_nonzero);
        unsigned int taylor20 = ResultRefinerHP::estimateMultiplicity(x, poly, 20, mpreal("1e-50"), first_nonzero);
        unsigned int ostro = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);
        
        std::cout << std::setw(6) << m << " | ";
        std::cout << std::setw(14) << taylor15 << " | ";
        std::cout << std::setw(14) << taylor20 << " | ";
        std::cout << std::setw(9) << ostro << "\n";
    }
    
    // Investigate Ostrowski failures
    testOstrowskiDetail(9);
    testOstrowskiDetail(10);
    
    // Test with higher precision
    std::cout << "\n=== Testing with 512-bit precision ===\n";
    PrecisionContext ctx512(512);

    std::cout << "True m | Taylor(max=20) | Ostrowski\n";
    std::cout << "-------+----------------+-----------\n";

    for (unsigned int m = 2; m <= 15; ++m) {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, m);
        mpreal x = mpreal(0.48);

        mpreal first_nonzero;
        unsigned int taylor = ResultRefinerHP::estimateMultiplicity(x, poly, 20, mpreal("1e-50"), first_nonzero);
        unsigned int ostro = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);

        std::cout << std::setw(6) << m << " | ";
        std::cout << std::setw(14) << taylor << " | ";
        std::cout << std::setw(9) << ostro << "\n";
    }

    // Test configurable ratio threshold for extreme multiplicities
    std::cout << "\n=== Testing configurable ratio threshold (256-bit) ===\n";
    PrecisionContext ctx256_2(256);

    std::cout << "True m | Thresh=10 | Thresh=50 | Thresh=100\n";
    std::cout << "-------+-----------+-----------+-----------\n";

    for (unsigned int m = 8; m <= 15; ++m) {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, m);
        mpreal x = mpreal(0.48);

        mpreal dummy;
        unsigned int t10 = ResultRefinerHP::estimateMultiplicity(x, poly, 20, mpreal("1e-50"), dummy, 10.0);
        unsigned int t50 = ResultRefinerHP::estimateMultiplicity(x, poly, 20, mpreal("1e-50"), dummy, 50.0);
        unsigned int t100 = ResultRefinerHP::estimateMultiplicity(x, poly, 20, mpreal("1e-50"), dummy, 100.0);

        std::cout << std::setw(6) << m << " | ";
        std::cout << std::setw(9) << t10 << " | ";
        std::cout << std::setw(9) << t50 << " | ";
        std::cout << std::setw(9) << t100 << "\n";
    }

    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

