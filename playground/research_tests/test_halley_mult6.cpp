#ifdef ENABLE_HIGH_PRECISION

#include "config.h"
#include "polynomial_hp.h"
#include "polynomial.h"
#include "differentiation_hp.h"
#include "precision_context.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

// Create (x - 0.6)^6 in power form by repeated multiplication
// Using high-precision arithmetic throughout to avoid double->HP conversion errors
PolynomialHP createMult6Polynomial() {
    mpreal root = mpreal("0.6");  // Use exact high-precision value
    unsigned int mult = 6;

    std::vector<mpreal> power_coeffs = {-root, mpreal("1.0")};

    for (unsigned int i = 1; i < mult; ++i) {
        std::vector<mpreal> new_coeffs(power_coeffs.size() + 1, mpreal("0.0"));
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-root);
            new_coeffs[j+1] += power_coeffs[j];
        }
        power_coeffs = new_coeffs;
    }

    // Create polynomial directly in high precision using fromPowerHP
    std::vector<unsigned int> degrees = {mult};
    return fromPowerHP(degrees, power_coeffs);
}

int main() {
    // Set precision to 512 bits
    PrecisionContext ctx(512);

    std::cout << std::setprecision(16);
    std::cout << "========================================" << std::endl;
    std::cout << "Halley's Method vs Modified Newton" << std::endl;
    std::cout << "Test polynomial: (x - 0.6)^6" << std::endl;
    std::cout << "Initial guess: x0 = 0.5999" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Create polynomial (x - 0.6)^6
    PolynomialHP poly = createMult6Polynomial();

    // Derivatives
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    PolynomialHP ddpoly = DifferentiationHP::derivative(poly, 0, 2);
    
    mpreal x0("0.5999");
    mpreal true_root("0.6");
    unsigned int max_iters = 10;
    
    // Test 1: Halley's method (standard, for simple roots)
    std::cout << "Test 1: Halley's Method (standard formula)" << std::endl;
    std::cout << "Formula: x_{n+1} = x_n - (2*f*f') / (2*(f')^2 - f*f'')" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    mpreal x = x0;
    for (unsigned int iter = 0; iter < max_iters; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        mpreal ddf = ddpoly.evaluate(x);
        
        mpreal error = abs(x - true_root);
        std::cout << "Iter " << iter << ": x = " << x 
                  << ", error = " << error << std::endl;
        
        if (error < mpreal("1e-50")) {
            std::cout << "  ✓ Converged!" << std::endl;
            break;
        }
        
        // Halley's formula
        mpreal numerator = mpreal("2.0") * f * df;
        mpreal denominator = mpreal("2.0") * df * df - f * ddf;
        
        if (abs(denominator) < mpreal("1e-100")) {
            std::cout << "  ✗ Denominator too small!" << std::endl;
            break;
        }
        
        x = x - numerator / denominator;
    }
    
    mpreal final_error1 = abs(x - true_root);
    std::cout << "\nFinal error: " << final_error1 << "\n" << std::endl;
    
    // Test 2: Modified Halley's Method with m=6
    std::cout << "\nTest 2: Modified Halley's Method with m=6" << std::endl;
    std::cout << "Formula: x_{n+1} = x_n - m*(2*f*f') / (2*(f')^2 - m*f*f'')" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    x = x0;
    unsigned int m = 6;
    for (unsigned int iter = 0; iter < max_iters; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        mpreal ddf = ddpoly.evaluate(x);

        mpreal error = abs(x - true_root);
        std::cout << "Iter " << iter << ": x = " << x
                  << ", error = " << error << std::endl;

        if (error < mpreal("1e-50")) {
            std::cout << "  ✓ Converged!" << std::endl;
            break;
        }

        // Modified Halley's formula
        mpreal numerator = mpreal(m) * mpreal("2.0") * f * df;
        mpreal denominator = mpreal("2.0") * df * df - mpreal(m) * f * ddf;

        if (abs(denominator) < mpreal("1e-100")) {
            std::cout << "  ✗ Denominator too small!" << std::endl;
            break;
        }

        x = x - numerator / denominator;
    }

    mpreal final_error_mod_halley = abs(x - true_root);
    std::cout << "\nFinal error: " << final_error_mod_halley << "\n" << std::endl;

    // Test 3: Modified Newton with m=6
    std::cout << "\nTest 3: Modified Newton with m=6" << std::endl;
    std::cout << "Formula: x_{n+1} = x_n - m * f / f'" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    x = x0;
    m = 6;
    for (unsigned int iter = 0; iter < max_iters; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        
        mpreal error = abs(x - true_root);
        std::cout << "Iter " << iter << ": x = " << x 
                  << ", error = " << error << std::endl;
        
        if (error < mpreal("1e-50")) {
            std::cout << "  ✓ Converged!" << std::endl;
            break;
        }
        
        if (abs(df) < mpreal("1e-100")) {
            std::cout << "  ✗ Derivative too small!" << std::endl;
            break;
        }
        
        x = x - mpreal(m) * f / df;
    }

    mpreal final_error_mod_newton = abs(x - true_root);
    std::cout << "\nFinal error: " << final_error_mod_newton << "\n" << std::endl;

    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Halley's method final error:          " << final_error1 << std::endl;
    std::cout << "Modified Halley (m=6) final error:    " << final_error_mod_halley << std::endl;
    std::cout << "Modified Newton (m=6) final error:    " << final_error_mod_newton << std::endl;
    
    return 0;
}

#else
int main() {
    std::cout << "High precision support not enabled. Rebuild with -DENABLE_HIGH_PRECISION=ON" << std::endl;
    return 1;
}
#endif

