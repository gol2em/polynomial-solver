/**
 * @file test_double_precision_limits.cpp
 * @brief Test limits of multiplicity detection in DOUBLE precision
 * 
 * This test explores:
 * 1. Maximum multiplicity detectable in double precision
 * 2. How far from the root can we detect multiplicity correctly
 * 3. Polynomials with multiple roots (some simple, some multiple)
 */

#include "polynomial.h"
#include "differentiation.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

// Create (x - r)^m in double precision
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

// Simple Taylor ratio test in double precision
unsigned int detectMultiplicityTaylor(double x, const Polynomial& poly, unsigned int max_order = 15) {
    std::vector<double> deriv_values;
    
    // Evaluate derivatives
    Polynomial p = poly;
    for (unsigned int k = 0; k <= max_order; ++k) {
        double val = std::abs(p.evaluate(x));
        deriv_values.push_back(val);
        
        if (k < max_order) {
            p = Differentiation::derivative(p, 0, 1);
        }
    }
    
    // Find first non-zero derivative using ratio test
    double ratio_threshold = 10.0;
    for (unsigned int k = 0; k < max_order; ++k) {
        if (deriv_values[k] < 1e-100) continue;  // Skip zeros
        
        // Check ratio with next derivative
        if (k + 1 < deriv_values.size() && deriv_values[k+1] > 1e-100) {
            double ratio = deriv_values[k] / deriv_values[k+1];
            if (ratio > ratio_threshold) {
                return k + 1;  // f^(k+1) is first non-zero
            }
        }
    }
    
    return 1;  // Default to simple root
}

// Test 1: Maximum multiplicity in double precision
void testMaxMultiplicity() {
    std::cout << "\n=== Test 1: Maximum Multiplicity in Double Precision ===\n";
    std::cout << "Polynomial: (x - 0.5)^m\n";
    std::cout << "Strategy: Do a few Newton iterations first, then detect multiplicity\n\n";
    std::cout << "m  | After iters | Distance | Detected | Status\n";
    std::cout << "---+-------------+----------+----------+--------\n";

    for (unsigned int m = 1; m <= 20; ++m) {
        Polynomial poly = createMultipleRoot(0.5, m);
        double x = 0.3;  // Start far away

        // Do 3 standard Newton iterations to get closer
        for (int iter = 0; iter < 3; ++iter) {
            double f = poly.evaluate(x);
            Polynomial dpoly = Differentiation::derivative(poly, 0, 1);
            double df = dpoly.evaluate(x);
            if (std::abs(df) > 1e-100) {
                x = x - f / df;
            }
        }

        double distance = std::abs(x - 0.5);
        unsigned int detected = detectMultiplicityTaylor(x, poly);

        std::cout << std::setw(2) << m << " | ";
        std::cout << std::setw(11) << "3 Newton" << " | ";
        std::cout << std::scientific << std::setprecision(2) << distance << " | ";
        std::cout << std::setw(8) << detected << " | ";

        if (detected == m) {
            std::cout << "✓ Correct";
        } else {
            std::cout << "✗ Wrong (detected " << detected << ")";
        }
        std::cout << "\n";
    }
}

// Test 2: Distance from root
void testDistanceFromRoot() {
    std::cout << "\n=== Test 2: Detection vs Distance from Root ===\n";
    std::cout << "Polynomial: (x - 0.5)^5\n\n";
    
    unsigned int true_m = 5;
    Polynomial poly = createMultipleRoot(0.5, true_m);
    
    std::vector<double> distances = {0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.49};
    
    std::cout << "Distance | x      | Detected | f(x)        | Status\n";
    std::cout << "---------+--------+----------+-------------+--------\n";
    
    for (double dist : distances) {
        double x = 0.5 - dist;
        unsigned int detected = detectMultiplicityTaylor(x, poly);
        double f = poly.evaluate(x);
        
        std::cout << std::setw(8) << std::fixed << std::setprecision(3) << dist << " | ";
        std::cout << std::setw(6) << std::fixed << std::setprecision(3) << x << " | ";
        std::cout << std::setw(8) << detected << " | ";
        std::cout << std::scientific << std::setprecision(2) << f << " | ";
        
        if (detected == true_m) {
            std::cout << "✓";
        } else {
            std::cout << "✗ (detected " << detected << ")";
        }
        std::cout << "\n";
    }
}

// Test 3: Polynomial with multiple roots
void testMultipleRoots() {
    std::cout << "\n=== Test 3: Polynomial with Multiple Roots ===\n";
    std::cout << "Polynomial: (x - 0.3)^3 * (x - 0.7)^2\n\n";
    
    // Create (x - 0.3)^3
    Polynomial p1 = createMultipleRoot(0.3, 3);
    // Create (x - 0.7)^2
    Polynomial p2 = createMultipleRoot(0.7, 2);
    
    // Multiply them
    // Get Bernstein coefficients
    auto coeffs1 = p1.bernsteinCoefficients();
    auto coeffs2 = p2.bernsteinCoefficients();
    
    // For multiplication, we need to convert to power, multiply, then back
    // This is complex, so let's create it directly in power form
    std::vector<double> power1 = {-0.3, 1.0};
    for (int i = 1; i < 3; ++i) {
        std::vector<double> new_coeffs(power1.size() + 1, 0.0);
        for (size_t j = 0; j < power1.size(); ++j) {
            new_coeffs[j] += power1[j] * (-0.3);
            new_coeffs[j+1] += power1[j];
        }
        power1 = new_coeffs;
    }
    
    std::vector<double> power2 = {-0.7, 1.0};
    for (int i = 1; i < 2; ++i) {
        std::vector<double> new_coeffs(power2.size() + 1, 0.0);
        for (size_t j = 0; j < power2.size(); ++j) {
            new_coeffs[j] += power2[j] * (-0.7);
            new_coeffs[j+1] += power2[j];
        }
        power2 = new_coeffs;
    }
    
    // Multiply power1 and power2
    std::vector<double> power_product(power1.size() + power2.size() - 1, 0.0);
    for (size_t i = 0; i < power1.size(); ++i) {
        for (size_t j = 0; j < power2.size(); ++j) {
            power_product[i + j] += power1[i] * power2[j];
        }
    }
    
    std::vector<unsigned int> degrees = {5};  // degree 3 + 2 = 5
    Polynomial poly = Polynomial::fromPower(degrees, power_product);
    
    std::cout << "Testing near root at x=0.3 (multiplicity 3):\n";
    std::cout << "Distance | Detected | Status\n";
    std::cout << "---------+----------+--------\n";
    
    std::vector<double> test_points = {0.28, 0.29, 0.295, 0.298, 0.299};
    for (double x : test_points) {
        unsigned int detected = detectMultiplicityTaylor(x, poly);
        double dist = std::abs(x - 0.3);
        
        std::cout << std::setw(8) << std::fixed << std::setprecision(3) << dist << " | ";
        std::cout << std::setw(8) << detected << " | ";
        
        if (detected == 3) {
            std::cout << "✓";
        } else {
            std::cout << "✗ (expected 3)";
        }
        std::cout << "\n";
    }
    
    std::cout << "\nTesting near root at x=0.7 (multiplicity 2):\n";
    std::cout << "Distance | Detected | Status\n";
    std::cout << "---------+----------+--------\n";
    
    test_points = {0.68, 0.69, 0.695, 0.698, 0.699};
    for (double x : test_points) {
        unsigned int detected = detectMultiplicityTaylor(x, poly);
        double dist = std::abs(x - 0.7);
        
        std::cout << std::setw(8) << std::fixed << std::setprecision(3) << dist << " | ";
        std::cout << std::setw(8) << detected << " | ";
        
        if (detected == 2) {
            std::cout << "✓";
        } else {
            std::cout << "✗ (expected 2)";
        }
        std::cout << "\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Double Precision Limits Test\n";
    std::cout << "========================================\n";
    
    testMaxMultiplicity();
    testDistanceFromRoot();
    testMultipleRoots();
    
    std::cout << "\n========================================\n";
    std::cout << "  Tests completed!\n";
    std::cout << "========================================\n";
    
    return 0;
}

