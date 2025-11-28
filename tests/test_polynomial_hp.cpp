#ifdef ENABLE_HIGH_PRECISION

#include "polynomial_hp.h"
#include "polynomial.h"
#include "precision_conversion.h"
#include "precision_context.h"
#include <iostream>
#include <cmath>
#include <cassert>

using namespace polynomial_solver;

void test_polynomial_hp_1d() {
    std::cout << "Test: PolynomialHP 1D evaluation" << std::endl;
    
    // Set precision to 256 bits (~77 decimal digits)
    PrecisionContext ctx(256);
    
    // Create polynomial: f(x) = x^3 - 2x + 1
    // Power coefficients: [1, -2, 0, 1]
    std::vector<unsigned int> degrees = {3};
    std::vector<double> power_coeffs = {1.0, -2.0, 0.0, 1.0};
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    
    // Convert to high-precision
    PolynomialHP poly_hp(poly);
    
    // Test evaluation at x = 0.5
    mpreal x = toHighPrecision(0.5);
    mpreal result_hp = poly_hp.evaluate(x);
    double result_double = poly.evaluate(0.5);
    
    // Expected: 0.5^3 - 2*0.5 + 1 = 0.125 - 1.0 + 1.0 = 0.125
    double expected = 0.125;
    
    std::cout << "  x = 0.5" << std::endl;
    std::cout << "  f(x) double = " << result_double << std::endl;
    std::cout << "  f(x) HP     = " << toString(result_hp, 20) << std::endl;
    std::cout << "  Expected    = " << expected << std::endl;
    
    double error_double = std::abs(result_double - expected);
    double error_hp = std::abs(toDouble(result_hp) - expected);
    
    std::cout << "  Error (double) = " << error_double << std::endl;
    std::cout << "  Error (HP)     = " << error_hp << std::endl;
    
    assert(error_double < 1e-14);
    assert(error_hp < 1e-14);
    
    std::cout << "  PASSED" << std::endl << std::endl;
}

void test_polynomial_hp_precision() {
    std::cout << "Test: PolynomialHP high-precision evaluation" << std::endl;
    
    // Set precision to 512 bits (~154 decimal digits)
    PrecisionContext ctx(512);
    
    // Create polynomial with ill-conditioned root
    // f(x) = (x - 0.5)^5 = x^5 - 2.5x^4 + 2.5x^3 - 1.25x^2 + 0.3125x - 0.03125
    std::vector<unsigned int> degrees = {5};
    std::vector<double> power_coeffs = {-0.03125, 0.3125, -1.25, 2.5, -2.5, 1.0};
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    
    // Convert to high-precision
    PolynomialHP poly_hp(poly);
    
    // Evaluate at x = 0.5 (the root)
    mpreal x = toHighPrecision(0.5);
    mpreal result_hp = poly_hp.evaluate(x);
    double result_double = poly.evaluate(0.5);
    
    std::cout << "  Polynomial: (x - 0.5)^5" << std::endl;
    std::cout << "  x = 0.5 (exact root)" << std::endl;
    std::cout << "  f(x) double = " << result_double << std::endl;
    std::cout << "  f(x) HP     = " << toString(result_hp, 50) << std::endl;
    
    // Both should be very close to zero
    assert(std::abs(result_double) < 1e-14);
    assert(std::abs(toDouble(result_hp)) < 1e-14);
    
    // Evaluate near the root: x = 0.5 + 1e-10
    mpreal x_near = toHighPrecision(0.5) + toHighPrecision(1e-10);
    mpreal result_near_hp = poly_hp.evaluate(x_near);
    double result_near_double = poly.evaluate(0.5 + 1e-10);
    
    std::cout << "  x = 0.5 + 1e-10" << std::endl;
    std::cout << "  f(x) double = " << result_near_double << std::endl;
    std::cout << "  f(x) HP     = " << toString(result_near_hp, 30) << std::endl;
    
    // Expected: (1e-10)^5 = 1e-50
    double expected = 1e-50;
    std::cout << "  Expected    ≈ " << expected << std::endl;
    
    // HP should be more accurate
    double error_hp = std::abs(toDouble(result_near_hp) - expected);
    std::cout << "  Error (HP)  = " << error_hp << std::endl;
    
    std::cout << "  PASSED" << std::endl << std::endl;
}

void test_polynomial_hp_conversion() {
    std::cout << "Test: Polynomial to PolynomialHP conversion" << std::endl;
    
    PrecisionContext ctx(256);
    
    // Create polynomial
    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, -3.0, 2.0};  // 2x^2 - 3x + 1
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    
    // Convert using constructor
    PolynomialHP poly_hp1(poly);
    
    // Convert using utility function
    PolynomialHP poly_hp2 = convertToHighPrecision(poly);
    
    // Both should give same results
    mpreal x = toHighPrecision(0.7);
    mpreal result1 = poly_hp1.evaluate(x);
    mpreal result2 = poly_hp2.evaluate(x);
    
    std::cout << "  Polynomial: 2x^2 - 3x + 1" << std::endl;
    std::cout << "  x = 0.7" << std::endl;
    std::cout << "  Result (constructor) = " << toString(result1, 20) << std::endl;
    std::cout << "  Result (utility)     = " << toString(result2, 20) << std::endl;
    
    double diff = std::abs(toDouble(result1 - result2));
    std::cout << "  Difference = " << diff << std::endl;
    
    assert(diff < 1e-50);
    
    std::cout << "  PASSED" << std::endl << std::endl;
}

int main() {
    std::cout << "=== PolynomialHP Tests ===" << std::endl << std::endl;
    
    test_polynomial_hp_1d();
    test_polynomial_hp_precision();
    test_polynomial_hp_conversion();
    
    std::cout << "=== All PolynomialHP tests passed ===" << std::endl;
    return 0;
}

#else

#include <iostream>

int main() {
    std::cout << "PolynomialHP tests skipped (ENABLE_HIGH_PRECISION not defined)" << std::endl;
    return 0;
}

#endif // ENABLE_HIGH_PRECISION

