/**
 * @file test_sturm_debug.cpp
 * @brief Debug Sturm sequence multiplicity detection
 */

#ifdef ENABLE_HIGH_PRECISION

#include "result_refiner_hp.h"
#include "polynomial_hp.h"
#include "differentiation_hp.h"
#include "precision_context.h"
#include "precision_conversion.h"
#include <iostream>
#include <iomanip>

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

int main() {
    PrecisionContext ctx(256);
    
    std::cout << "========================================\n";
    std::cout << "  Sturm Sequence Debug\n";
    std::cout << "========================================\n\n";
    
    // Test case: (x - 0.5)^3
    unsigned int m = 3;
    mpreal root = mpreal("0.5");
    PolynomialHP poly = createMultipleRootPolynomialHP(root, m);
    
    std::cout << "Polynomial: (x - 0.5)^" << m << "\n";
    std::cout << "True root: " << root << "\n";
    std::cout << "True multiplicity: " << m << "\n\n";
    
    // Test at the root
    mpreal x = mpreal("0.5");
    std::cout << "Testing at x = " << x << "\n";
    std::cout << "f(x) = " << poly.evaluate(x) << "\n\n";
    
    // Test derivatives
    std::cout << "Derivatives at x = " << x << ":\n";
    for (unsigned int k = 0; k <= 5; ++k) {
        PolynomialHP deriv = DifferentiationHP::derivative(poly, 0, k);
        mpreal val = deriv.evaluate(x);
        std::cout << "  f^(" << k << ")(x) = " << val << "\n";
    }
    std::cout << "\n";
    
    // Test Sturm with different radii
    std::cout << "Sturm method with different interval radii:\n";
    std::cout << "Radius      | Detected m\n";
    std::cout << "------------+-----------\n";
    
    std::vector<std::string> radii = {"1e-5", "1e-6", "1e-8", "1e-10", "1e-12", "1e-15"};
    for (const auto& r_str : radii) {
        mpreal radius = mpreal(r_str);
        unsigned int detected = ResultRefinerHP::estimateMultiplicitySturm(x, poly, radius);
        std::cout << std::setw(11) << r_str << " | " << detected << "\n";
    }
    std::cout << "\n";
    
    // Test at nearby point
    mpreal x_near = mpreal("0.48");
    std::cout << "Testing at nearby point x = " << x_near << "\n";
    std::cout << "f(x) = " << poly.evaluate(x_near) << "\n\n";
    
    std::cout << "Sturm method with different interval radii:\n";
    std::cout << "Radius      | Detected m\n";
    std::cout << "------------+-----------\n";
    
    for (const auto& r_str : radii) {
        mpreal radius = mpreal(r_str);
        unsigned int detected = ResultRefinerHP::estimateMultiplicitySturm(x_near, poly, radius);
        std::cout << std::setw(11) << r_str << " | " << detected << "\n";
    }
    std::cout << "\n";
    
    // Compare with Taylor method
    std::cout << "Comparison with Taylor method:\n";
    mpreal dummy;
    unsigned int taylor = ResultRefinerHP::estimateMultiplicity(
        x_near, poly, 15, mpreal("1e-50"), dummy, 10.0);
    std::cout << "Taylor: " << taylor << "\n";
    
    unsigned int ostro = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x_near, poly);
    std::cout << "Ostrowski: " << ostro << "\n";
    
    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

