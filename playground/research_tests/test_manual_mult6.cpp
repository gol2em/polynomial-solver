#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "hp/polynomial_hp.h"
#include "core/polynomial.h"
#include "hp/differentiation_hp.h"
#include "hp/result_refiner_hp.h"
#include "hp/high_precision_types.h"
#include "hp/precision_context.h"
#include <iostream>
#include <iomanip>

using namespace polynomial_solver;

// Helper to create (x - r)^m polynomial in power form
PolynomialHP createMultipleRootPolynomial(double root, unsigned int multiplicity) {
    // Create (x - r)^m in power form by repeated multiplication
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

    // Create polynomial in power form (double precision first)
    std::vector<unsigned int> degrees = {multiplicity};
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);

    // Convert to high precision
    return PolynomialHP(poly);
}

int main() {
    PrecisionContext ctx(512);  // 512 bits precision

    std::cout << "========================================\n";
    std::cout << "Manual Test: Modified Newton with m=6\n";
    std::cout << "========================================\n\n";
    
    // Create (x - 0.6)^6
    PolynomialHP poly = createMultipleRootPolynomial(0.6, 6);
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    
    std::cout << "Polynomial: (x - 0.6)^6\n";
    std::cout << "True root: 0.6\n";
    std::cout << "True multiplicity: 6\n\n";
    
    // Test 1: Modified Newton with m=6 from x=0.5999
    std::cout << "Test 1: Modified Newton with m=6 from x=0.5999\n";
    std::cout << "------------------------------------------------\n";
    mpreal x = mpreal("0.5999");
    for (int iter = 0; iter < 10; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        
        std::cout << "Iter " << iter << ": x = " << x 
                  << ", f = " << f << ", df = " << df << std::endl;
        
        if (abs(f) < mpreal("1e-100")) {
            std::cout << "Converged!\n";
            break;
        }
        
        // Modified Newton: x_new = x - m * f / f'
        mpreal step = mpreal(6) * f / df;
        x = x - step;
    }
    std::cout << "Final x = " << x << std::endl;
    std::cout << "Error = " << abs(x - mpreal("0.6")) << "\n\n";
    
    // Test 2: Modified Newton with m=5 from x=0.5999
    std::cout << "Test 2: Modified Newton with m=5 from x=0.5999 (WRONG multiplicity)\n";
    std::cout << "----------------------------------------------------------------------\n";
    x = mpreal("0.5999");
    for (int iter = 0; iter < 10; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        
        std::cout << "Iter " << iter << ": x = " << x 
                  << ", f = " << f << ", df = " << df << std::endl;
        
        if (abs(f) < mpreal("1e-100")) {
            std::cout << "Converged!\n";
            break;
        }
        
        // Modified Newton with WRONG multiplicity: x_new = x - 5 * f / f'
        mpreal step = mpreal(5) * f / df;
        x = x - step;
    }
    std::cout << "Final x = " << x << std::endl;
    std::cout << "Error = " << abs(x - mpreal("0.6")) << "\n\n";
    
    // Test 3: Check derivatives at x=0.6
    std::cout << "Test 3: Derivatives at x=0.6\n";
    std::cout << "-----------------------------\n";
    for (unsigned int k = 1; k <= 8; ++k) {
        PolynomialHP deriv = DifferentiationHP::derivative(poly, 0, k);
        mpreal val = abs(deriv.evaluate(mpreal("0.6")));
        std::cout << "f^(" << k << ")(0.6) = " << val << std::endl;
    }
    
    return 0;
}

#else
int main() {
    std::cout << "High precision support not enabled. Rebuild with -DENABLE_HIGH_PRECISION=ON\n";
    return 1;
}
#endif

