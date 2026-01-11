/**
 * @file test_hp_conversion_benefit.cpp
 * @brief Test if power→HP→Bernstein→double is better than power→Bernstein(double)
 * 
 * Compare two approaches:
 * 1. Direct: power coeffs → Polynomial::fromPower() → Bernstein in double
 * 2. Via HP: power coeffs → HP → Bernstein in HP → convert to double
 * 
 * This shows how much precision is lost in double-precision Bernstein conversion.
 */

#ifdef ENABLE_HIGH_PRECISION

#include "core/polynomial.h"
#include "hp/polynomial_hp.h"
#include "hp/precision_conversion.h"
#include "hp/precision_context.h"
#include "refinement/result_refiner.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

// Create (x - r)^m directly in double
Polynomial createMultipleRootDirect(double root, unsigned int multiplicity) {
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

// Create (x - r)^m via HP then convert to double
Polynomial createMultipleRootViaHP(double root, unsigned int multiplicity, unsigned int precision_bits) {
    PrecisionContext ctx(precision_bits);
    
    // Create in HP
    std::vector<mpreal> power_coeffs = {-mpreal(root), mpreal(1)};
    for (unsigned int i = 1; i < multiplicity; ++i) {
        std::vector<mpreal> new_coeffs(power_coeffs.size() + 1, mpreal(0));
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-mpreal(root));
            new_coeffs[j+1] += power_coeffs[j];
        }
        power_coeffs = new_coeffs;
    }
    
    // Convert to Bernstein in HP
    std::vector<unsigned int> degrees = {multiplicity};
    PolynomialHP poly_hp = fromPowerHP(degrees, power_coeffs);

    // Convert HP Bernstein coefficients to double
    const auto& hp_coeffs = poly_hp.bernsteinCoefficients();
    std::vector<double> double_coeffs = toDouble(hp_coeffs);

    return Polynomial(degrees, double_coeffs);
}

void testMultiplicityDetection() {
    std::cout << "\n=== Multiplicity Detection Comparison ===\n";
    std::cout << "Polynomial: (x - 0.5)^m, test at x = 0.48\n\n";
    
    ResultRefiner refiner;
    
    std::cout << "m  | Direct | Via 64-bit HP | Via 128-bit HP | Via 256-bit HP\n";
    std::cout << "---+--------+---------------+----------------+----------------\n";
    
    for (unsigned int m = 1; m <= 10; ++m) {
        Polynomial poly_direct = createMultipleRootDirect(0.5, m);
        Polynomial poly_hp64 = createMultipleRootViaHP(0.5, m, 64);
        Polynomial poly_hp128 = createMultipleRootViaHP(0.5, m, 128);
        Polynomial poly_hp256 = createMultipleRootViaHP(0.5, m, 256);
        
        double x0 = 0.48;
        unsigned int detected_direct = refiner.estimateMultiplicityOstrowskiFromPoint(x0, poly_direct);
        unsigned int detected_hp64 = refiner.estimateMultiplicityOstrowskiFromPoint(x0, poly_hp64);
        unsigned int detected_hp128 = refiner.estimateMultiplicityOstrowskiFromPoint(x0, poly_hp128);
        unsigned int detected_hp256 = refiner.estimateMultiplicityOstrowskiFromPoint(x0, poly_hp256);
        
        std::cout << std::setw(2) << m << " | ";
        
        // Direct
        if (detected_direct == m) {
            std::cout << "   ✓   | ";
        } else {
            std::cout << "   " << detected_direct << "   | ";
        }
        
        // Via 64-bit HP
        if (detected_hp64 == m) {
            std::cout << "      ✓       | ";
        } else {
            std::cout << "      " << detected_hp64 << "       | ";
        }
        
        // Via 128-bit HP
        if (detected_hp128 == m) {
            std::cout << "       ✓        | ";
        } else {
            std::cout << "       " << detected_hp128 << "        | ";
        }
        
        // Via 256-bit HP
        if (detected_hp256 == m) {
            std::cout << "       ✓\n";
        } else {
            std::cout << "       " << detected_hp256 << "\n";
        }
    }
}

void testCoefficientDifference() {
    std::cout << "\n=== Bernstein Coefficient Difference ===\n";
    std::cout << "Polynomial: (x - 0.5)^m\n";
    std::cout << "Compare Bernstein coefficients from different conversion paths\n\n";

    for (unsigned int m : {5, 8, 9, 10}) {
        std::cout << "\n--- m = " << m << " ---\n";

        Polynomial poly_direct = createMultipleRootDirect(0.5, m);
        Polynomial poly_hp64 = createMultipleRootViaHP(0.5, m, 64);
        Polynomial poly_hp128 = createMultipleRootViaHP(0.5, m, 128);

        const auto& coeffs_direct = poly_direct.bernsteinCoefficients();
        const auto& coeffs_hp64 = poly_hp64.bernsteinCoefficients();
        const auto& coeffs_hp128 = poly_hp128.bernsteinCoefficients();

        // Compute max difference
        double max_diff_direct = 0.0;
        double max_diff_hp64 = 0.0;

        for (size_t i = 0; i < coeffs_direct.size(); ++i) {
            double diff_direct = std::abs(coeffs_direct[i] - coeffs_hp128[i]);
            double diff_hp64 = std::abs(coeffs_hp64[i] - coeffs_hp128[i]);
            max_diff_direct = std::max(max_diff_direct, diff_direct);
            max_diff_hp64 = std::max(max_diff_hp64, diff_hp64);
        }

        std::cout << "Max coefficient difference from 128-bit HP:\n";
        std::cout << "  Direct (double):  " << std::scientific << std::setprecision(2) << max_diff_direct << "\n";
        std::cout << "  Via 64-bit HP:    " << max_diff_hp64 << "\n";

        if (max_diff_hp64 < max_diff_direct * 0.1) {
            std::cout << "  ✓ HP conversion gives " << (max_diff_direct / max_diff_hp64) << "x better precision\n";
        } else if (max_diff_hp64 < max_diff_direct) {
            std::cout << "  ⚠ HP conversion slightly better\n";
        } else {
            std::cout << "  ✗ No benefit from HP conversion\n";
        }
    }
}

void testMultipleRootsDetection() {
    std::cout << "\n=== Multiple Roots Detection Comparison ===\n";
    std::cout << "Polynomial: (x - 0.2) * (x - 0.5)^3 * (x - 0.8)^2\n\n";

    // Pre-expanded coefficients
    std::vector<double> power_coeffs = {
        0.008,      // constant
        -0.13568,   // x
        0.8732,     // x^2
        -2.717,     // x^3
        4.29,       // x^4
        -3.3,       // x^5
        1.0         // x^6
    };
    std::vector<unsigned int> degrees = {6};

    // Create via direct and HP conversion
    Polynomial poly_direct = Polynomial::fromPower(degrees, power_coeffs);

    PrecisionContext ctx(128);
    std::vector<mpreal> power_coeffs_hp = toHighPrecision(power_coeffs);
    PolynomialHP poly_hp = fromPowerHP(degrees, power_coeffs_hp);
    std::vector<double> hp_coeffs_double = toDouble(poly_hp.bernsteinCoefficients());
    Polynomial poly_via_hp(degrees, hp_coeffs_double);

    ResultRefiner refiner;

    struct TestCase {
        double root;
        unsigned int mult;
        double test_x;
    };

    std::vector<TestCase> tests = {
        {0.2, 1, 0.19},
        {0.5, 3, 0.49},
        {0.8, 2, 0.79},
    };

    std::cout << "Root | True m | Test x | Direct | Via 128-bit HP\n";
    std::cout << "-----+--------+--------+--------+----------------\n";

    for (const auto& test : tests) {
        unsigned int detected_direct = refiner.estimateMultiplicityOstrowskiFromPoint(test.test_x, poly_direct);
        unsigned int detected_hp = refiner.estimateMultiplicityOstrowskiFromPoint(test.test_x, poly_via_hp);

        std::cout << std::fixed << std::setprecision(1) << test.root << " | ";
        std::cout << std::setw(6) << test.mult << " | ";
        std::cout << std::setprecision(2) << test.test_x << " | ";

        if (detected_direct == test.mult) {
            std::cout << "   ✓   | ";
        } else {
            std::cout << "   " << detected_direct << "   | ";
        }

        if (detected_hp == test.mult) {
            std::cout << "       ✓\n";
        } else {
            std::cout << "       " << detected_hp << "\n";
        }
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  HP Conversion Benefit Test\n";
    std::cout << "  Power→HP→Bernstein→Double vs\n";
    std::cout << "  Power→Bernstein(Double)\n";
    std::cout << "========================================\n";

    testMultiplicityDetection();
    testCoefficientDifference();
    testMultipleRootsDetection();

    std::cout << "\n========================================\n";
    std::cout << "  Tests completed!\n";
    std::cout << "========================================\n";

    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

