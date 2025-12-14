/**
 * @file test_detection_limits.cpp
 * @brief Test limits of multiplicity detection with 64-bit precision
 * 
 * Based on findings that Taylor/Ostrowski work at 64-bit precision,
 * this test explores:
 * 1. Polynomials with multiple roots (some simple, some multiple)
 * 2. How far from the root can we detect correctly
 * 3. Maximum multiplicity detectable at 64-bit precision
 */

#ifdef ENABLE_HIGH_PRECISION

#include "result_refiner_hp.h"
#include "polynomial_hp.h"
#include "differentiation_hp.h"
#include "precision_context.h"
#include "precision_conversion.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace polynomial_solver;

// Create (x - r)^m in native high precision
PolynomialHP createMultipleRootHP(const mpreal& root, unsigned int multiplicity) {
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

// Multiply two polynomials in power form
std::vector<mpreal> multiplyPowerPolynomials(const std::vector<mpreal>& p1, const std::vector<mpreal>& p2) {
    std::vector<mpreal> result(p1.size() + p2.size() - 1, mpreal(0));
    for (size_t i = 0; i < p1.size(); ++i) {
        for (size_t j = 0; j < p2.size(); ++j) {
            result[i + j] += p1[i] * p2[j];
        }
    }
    return result;
}

// Test 1: Polynomial with multiple roots
void testMultipleRoots() {
    PrecisionContext ctx(64);
    
    std::cout << "\n=== Test 1: Polynomial with Multiple Roots (64-bit) ===\n";
    std::cout << "Polynomial: (x - 0.2) * (x - 0.5)^3 * (x - 0.8)^2\n\n";
    
    // Create each factor in power form
    std::vector<mpreal> p1 = {-mpreal("0.2"), mpreal(1)};  // (x - 0.2)
    
    std::vector<mpreal> p2 = {-mpreal("0.5"), mpreal(1)};  // (x - 0.5)
    for (int i = 1; i < 3; ++i) {
        std::vector<mpreal> new_coeffs(p2.size() + 1, mpreal(0));
        for (size_t j = 0; j < p2.size(); ++j) {
            new_coeffs[j] += p2[j] * (-mpreal("0.5"));
            new_coeffs[j+1] += p2[j];
        }
        p2 = new_coeffs;
    }
    
    std::vector<mpreal> p3 = {-mpreal("0.8"), mpreal(1)};  // (x - 0.8)
    for (int i = 1; i < 2; ++i) {
        std::vector<mpreal> new_coeffs(p3.size() + 1, mpreal(0));
        for (size_t j = 0; j < p3.size(); ++j) {
            new_coeffs[j] += p3[j] * (-mpreal("0.8"));
            new_coeffs[j+1] += p3[j];
        }
        p3 = new_coeffs;
    }
    
    // Multiply: p1 * p2 * p3
    std::vector<mpreal> product = multiplyPowerPolynomials(p1, p2);
    product = multiplyPowerPolynomials(product, p3);
    
    std::vector<unsigned int> degrees = {6};  // degree 1+3+2 = 6
    PolynomialHP poly = fromPowerHP(degrees, product);
    
    std::cout << "Root | True m | Test x  | Distance | Taylor | Ostrowski | Status\n";
    std::cout << "-----+--------+---------+----------+--------+-----------+--------\n";
    
    struct TestCase {
        double root;
        unsigned int mult;
        double test_x;
    };
    
    std::vector<TestCase> tests = {
        {0.2, 1, 0.19},
        {0.2, 1, 0.18},
        {0.5, 3, 0.49},
        {0.5, 3, 0.48},
        {0.5, 3, 0.45},
        {0.8, 2, 0.79},
        {0.8, 2, 0.78},
    };
    
    for (const auto& test : tests) {
        mpreal x = mpreal(test.test_x);
        mpreal dummy;
        unsigned int taylor = ResultRefinerHP::estimateMultiplicity(
            x, poly, 15, mpreal("1e-50"), dummy, 10.0);
        unsigned int ostro = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);
        
        double dist = std::abs(test.test_x - test.root);
        
        std::cout << std::fixed << std::setprecision(1) << test.root << " | ";
        std::cout << std::setw(6) << test.mult << " | ";
        std::cout << std::setprecision(2) << test.test_x << " | ";
        std::cout << std::scientific << std::setprecision(2) << dist << " | ";
        std::cout << std::setw(6) << taylor << " | ";
        std::cout << std::setw(9) << ostro << " | ";
        
        bool correct = (taylor == test.mult && ostro == test.mult);
        if (correct) {
            std::cout << "✓";
        } else {
            std::cout << "✗ (T=" << taylor << ", O=" << ostro << ")";
        }
        std::cout << "\n";
    }
}

// Test 2: Distance from root
void testDistanceFromRoot() {
    PrecisionContext ctx(64);
    
    std::cout << "\n=== Test 2: Detection vs Distance from Root (64-bit) ===\n";
    std::cout << "Polynomial: (x - 0.5)^5\n\n";
    
    unsigned int true_m = 5;
    PolynomialHP poly = createMultipleRootHP(mpreal("0.5"), true_m);
    
    std::vector<double> distances = {0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
    
    std::cout << "Distance | x      | Taylor | Ostrowski | Status\n";
    std::cout << "---------+--------+--------+-----------+--------\n";
    
    for (double dist : distances) {
        mpreal x = mpreal("0.5") - mpreal(dist);
        mpreal dummy;
        unsigned int taylor = ResultRefinerHP::estimateMultiplicity(
            x, poly, 15, mpreal("1e-50"), dummy, 10.0);
        unsigned int ostro = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);
        
        std::cout << std::setw(8) << std::fixed << std::setprecision(3) << dist << " | ";
        std::cout << std::setw(6) << std::fixed << std::setprecision(3) << toDouble(x) << " | ";
        std::cout << std::setw(6) << taylor << " | ";
        std::cout << std::setw(9) << ostro << " | ";
        
        bool correct = (taylor == true_m && ostro == true_m);
        if (correct) {
            std::cout << "✓ Both correct";
        } else if (taylor == true_m) {
            std::cout << "⚠ Taylor OK, Ostro=" << ostro;
        } else if (ostro == true_m) {
            std::cout << "⚠ Ostro OK, Taylor=" << taylor;
        } else {
            std::cout << "✗ Both wrong";
        }
        std::cout << "\n";
    }
}

// Test 3: Maximum multiplicity at 64-bit
void testMaxMultiplicity() {
    PrecisionContext ctx(64);

    std::cout << "\n=== Test 3: Maximum Multiplicity (64-bit) ===\n";
    std::cout << "Polynomial: (x - 0.5)^m, testing at x = 0.48\n\n";
    std::cout << "m  | Taylor | Ostrowski | Status\n";
    std::cout << "---+--------+-----------+--------\n";

    for (unsigned int m = 1; m <= 20; ++m) {
        PolynomialHP poly = createMultipleRootHP(mpreal("0.5"), m);
        mpreal x = mpreal("0.48");

        mpreal dummy;
        unsigned int taylor = ResultRefinerHP::estimateMultiplicity(
            x, poly, 15, mpreal("1e-50"), dummy, 10.0);
        unsigned int ostro = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);

        std::cout << std::setw(2) << m << " | ";
        std::cout << std::setw(6) << taylor << " | ";
        std::cout << std::setw(9) << ostro << " | ";

        bool taylor_correct = (taylor == m);
        bool ostro_correct = (ostro == m);

        if (taylor_correct && ostro_correct) {
            std::cout << "✓ Both correct";
        } else if (taylor_correct) {
            std::cout << "⚠ Taylor OK, Ostro=" << ostro;
        } else if (ostro_correct) {
            std::cout << "⚠ Ostro OK, Taylor=" << taylor;
        } else {
            std::cout << "✗ Both wrong (T=" << taylor << ", O=" << ostro << ")";
        }
        std::cout << "\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Detection Limits Test (64-bit HP)\n";
    std::cout << "========================================\n";

    testMultipleRoots();
    testDistanceFromRoot();
    testMaxMultiplicity();

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

