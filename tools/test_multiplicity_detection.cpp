/**
 * @file test_multiplicity_detection.cpp
 * @brief Test multiplicity detection on various polynomial roots
 */

#include "refinement/result_refiner.h"
#include "core/polynomial.h"
#include "solver/solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

void test_multiplicity(const std::string& name, 
                       const std::vector<double>& power_coeffs,
                       double root_location,
                       unsigned int expected_mult) {
    std::vector<unsigned int> degrees{static_cast<unsigned int>(power_coeffs.size() - 1)};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    ResultRefiner refiner;
    std::vector<double> point{root_location};
    double first_deriv = 0.0;
    unsigned int mult = refiner.estimateMultiplicity(point, system, 10, 1e-10, first_deriv);

    std::cout << "  " << std::setw(30) << std::left << name
              << ": mult=" << mult 
              << " (expected " << expected_mult << ")";
    
    if (mult == expected_mult) {
        std::cout << " ✓" << std::endl;
    } else {
        std::cout << " ✗" << std::endl;
    }
}

int main() {
    std::cout << "Multiplicity Detection Test" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << std::endl;
    
    // Simple roots
    std::cout << "Simple roots (multiplicity = 1):" << std::endl;
    test_multiplicity("(x-0.5)", {-0.5, 1.0}, 0.5, 1);
    test_multiplicity("(x-0.3)", {-0.3, 1.0}, 0.3, 1);
    test_multiplicity("(x-0.7)", {-0.7, 1.0}, 0.7, 1);
    
    std::cout << std::endl;
    std::cout << "Double roots (multiplicity = 2):" << std::endl;
    test_multiplicity("(x-0.5)^2", {0.25, -1.0, 1.0}, 0.5, 2);
    test_multiplicity("(x-0.3)^2", {0.09, -0.6, 1.0}, 0.3, 2);
    
    std::cout << std::endl;
    std::cout << "Triple roots (multiplicity = 3):" << std::endl;
    test_multiplicity("(x-0.5)^3", {-0.125, 0.75, -1.5, 1.0}, 0.5, 3);
    test_multiplicity("(x-0.4)^3", {-0.064, 0.48, -1.2, 1.0}, 0.4, 3);
    
    std::cout << std::endl;
    std::cout << "Quadruple roots (multiplicity = 4):" << std::endl;
    test_multiplicity("(x-0.5)^4", {0.0625, -0.5, 1.5, -2.0, 1.0}, 0.5, 4);
    test_multiplicity("(x-0.3)^4", {0.0081, -0.108, 0.54, -1.2, 1.0}, 0.3, 4);
    
    std::cout << std::endl;
    std::cout << "Quintuple roots (multiplicity = 5):" << std::endl;
    test_multiplicity("(x-0.5)^5", {-0.03125, 0.3125, -1.25, 2.5, -2.5, 1.0}, 0.5, 5);
    
    std::cout << std::endl;
    std::cout << "Sextuple roots (multiplicity = 6):" << std::endl;
    test_multiplicity("(x-0.6)^6", {0.046656, -0.46656, 1.9440, -4.32, 5.4, -3.6, 1.0}, 0.6, 6);
    
    std::cout << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "All multiplicity detection tests completed!" << std::endl;
    
    return 0;
}

