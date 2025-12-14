/**
 * @file test_double_ostrowski.cpp
 * @brief Test Ostrowski-based multiplicity detection in double precision
 *
 * Focus: Test DETECTION only, not convergence
 * Convergence in double precision is expected to fail for multiple roots
 * We just want to verify that Ostrowski can detect multiplicity correctly
 */

#include "result_refiner.h"
#include "polynomial.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

// Create (x - r)^m in power form, then convert to Bernstein
Polynomial createMultipleRoot(double root, unsigned int multiplicity) {
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

// Test 1: Multiplicity detection for single root
void testMultiplicityDetection() {
    std::cout << "\n=== Test 1: Multiplicity Detection (Double Precision) ===\n";
    std::cout << "Polynomial: (x - 0.5)^m\n";
    std::cout << "Test point: x = 0.48 (distance = 0.02 from root)\n\n";

    ResultRefiner refiner;

    std::cout << "m  | Detected m | Status\n";
    std::cout << "---+------------+--------\n";

    for (unsigned int m = 1; m <= 10; ++m) {
        Polynomial poly = createMultipleRoot(0.5, m);

        double x0 = 0.48;  // Start near root
        unsigned int detected_m = refiner.estimateMultiplicityOstrowskiFromPoint(x0, poly);

        std::cout << std::setw(2) << m << " | ";
        std::cout << std::setw(10) << detected_m << " | ";

        if (detected_m == m) {
            std::cout << "✓ Correct";
        } else {
            std::cout << "✗ Wrong (detected " << detected_m << ", expected " << m << ")";
        }
        std::cout << "\n";
    }
}

// Test 2: Distance tolerance
void testDistanceTolerance() {
    std::cout << "\n=== Test 2: Distance Tolerance ===\n";
    std::cout << "Polynomial: (x - 0.5)^5\n\n";

    Polynomial poly = createMultipleRoot(0.5, 5);
    ResultRefiner refiner;

    std::vector<double> distances = {0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2};

    std::cout << "Distance | Test x | Detected m | Status\n";
    std::cout << "---------+--------+------------+--------\n";

    for (double dist : distances) {
        double x = 0.5 - dist;
        unsigned int detected_m = refiner.estimateMultiplicityOstrowskiFromPoint(x, poly);

        std::cout << std::setw(8) << std::fixed << std::setprecision(3) << dist << " | ";
        std::cout << std::setw(6) << std::setprecision(3) << x << " | ";
        std::cout << std::setw(10) << detected_m << " | ";

        if (detected_m == 5) {
            std::cout << "✓ Correct";
        } else {
            std::cout << "✗ Wrong (detected " << detected_m << ", expected 5)";
        }
        std::cout << "\n";
    }
}

// Test 3: Multiple roots in same polynomial
void testMultipleRootsDetection() {
    std::cout << "\n=== Test 3: Multiple Roots in Same Polynomial ===\n";
    std::cout << "Polynomial: (x - 0.2) * (x - 0.5)^3 * (x - 0.8)^2\n\n";

    // Pre-expanded coefficients to avoid precision loss from multiplication
    // (x - 0.2)(x - 0.5)^3(x - 0.8)^2
    // Expanded form: x^6 - 3.3x^5 + 4.29x^4 - 2.717x^3 + 0.8732x^2 - 0.13568x + 0.008
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
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);

    ResultRefiner refiner;

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

    std::cout << "Root | True m | Test x | Detected m | Status\n";
    std::cout << "-----+--------+--------+------------+--------\n";

    for (const auto& test : tests) {
        unsigned int detected_m = refiner.estimateMultiplicityOstrowskiFromPoint(test.test_x, poly);

        std::cout << std::fixed << std::setprecision(1) << test.root << " | ";
        std::cout << std::setw(6) << test.mult << " | ";
        std::cout << std::setprecision(2) << test.test_x << " | ";
        std::cout << std::setw(10) << detected_m << " | ";

        if (detected_m == test.mult) {
            std::cout << "✓ Correct";
        } else {
            std::cout << "✗ Wrong (detected " << detected_m << ", expected " << test.mult << ")";
        }
        std::cout << "\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Double Precision Ostrowski Test\n";
    std::cout << "  (Detection Only - Not Convergence)\n";
    std::cout << "========================================\n";

    testMultiplicityDetection();
    testDistanceTolerance();
    testMultipleRootsDetection();

    std::cout << "\n========================================\n";
    std::cout << "  Tests completed!\n";
    std::cout << "========================================\n";

    return 0;
}

