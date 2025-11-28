#ifdef ENABLE_HIGH_PRECISION

#include "differentiation_hp.h"
#include "differentiation.h"
#include "polynomial_hp.h"
#include "polynomial.h"
#include "precision_context.h"
#include "precision_conversion.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace polynomial_solver;

// Test 1: First derivative of 1D polynomial - verify true HP precision
bool test_first_derivative_1d() {
    std::cout << "\nTest: First derivative of 1D polynomial - verify true HP precision\n";

    // Set precision to 256 bits (~77 decimal digits)
    PrecisionContext ctx(256);

    // Create polynomial DIRECTLY in HP: f(x) = x^3 - 2x + 1
    std::vector<unsigned int> degrees = {3};
    std::vector<mpreal> power_coeffs_hp = {
        mpreal(1), mpreal(-2), mpreal(0), mpreal(1)
    };
    PolynomialHP poly_hp = fromPowerHP(degrees, power_coeffs_hp);

    // Differentiate in HP
    PolynomialHP dpoly_hp = DifferentiationHP::derivative(poly_hp, 0, 1);

    // Evaluate at x = 0.5 (exactly representable)
    mpreal x_hp = mpreal("0.5");
    mpreal df_hp = dpoly_hp.evaluate(x_hp);

    // Expected: f'(x) = 3x^2 - 2, so f'(0.5) = 3*0.25 - 2 = -1.25 (exactly representable)
    mpreal expected_hp = mpreal("-1.25");

    // Compute error in high precision
    mpreal error_hp = abs(df_hp - expected_hp);

    std::cout << "  f(x) = x^3 - 2x + 1\n";
    std::cout << "  f'(0.5) HP     = " << toString(df_hp, 80) << "\n";
    std::cout << "  Expected       = " << toString(expected_hp, 80) << "\n";
    std::cout << "  Error (HP)     = " << toString(error_hp, 10) << "\n";

    // For 256-bit precision, machine epsilon is ~2^-256 ≈ 1e-77
    // Since we constructed polynomial directly in HP, error should be near machine precision
    std::cout << "  HP machine epsilon: ~1e-77 for 256-bit precision\n";
    std::cout << "  Note: Polynomial constructed directly in HP (no double precision limitation)\n";

    // Error should be near HP machine precision
    bool pass = error_hp < mpreal(1e-70);  // Should be ~1e-77, allow slack
    std::cout << "  " << (pass ? "PASSED" : "FAILED") << " (HP error < 1e-70)\n";
    return pass;
}

// Test 2: Second derivative - verify HP precision
bool test_second_derivative_1d() {
    std::cout << "\nTest: Second derivative of 1D polynomial - verify HP precision\n";

    PrecisionContext ctx(256);

    // Polynomial: f(x) = x^4
    // Bernstein coefficients for degree 4: [0, 0, 0, 0, 1]
    std::vector<double> coeffs_double = {0.0, 0.0, 0.0, 0.0, 1.0};
    Polynomial poly(std::vector<unsigned int>{4}, coeffs_double);

    // Convert to HP and differentiate twice
    PolynomialHP poly_hp(poly);
    PolynomialHP dpoly_hp = DifferentiationHP::derivative(poly_hp, 0, 1);
    PolynomialHP ddpoly_hp = DifferentiationHP::derivative(dpoly_hp, 0, 1);

    // Evaluate at x = 0.5 (exactly representable)
    mpreal x_hp = toHighPrecision(0.5);
    mpreal ddf_hp = ddpoly_hp.evaluate(x_hp);

    // Expected: f''(x) = 12x^2, so f''(0.5) = 12*0.25 = 3.0 (exactly representable)
    mpreal expected_hp = mpreal(3.0);

    // Compute error in high precision
    mpreal error_hp = abs(ddf_hp - expected_hp);

    std::cout << "  f(x) = x^4\n";
    std::cout << "  f''(0.5) HP = " << toString(ddf_hp, 30) << "\n";
    std::cout << "  Expected    = " << toString(expected_hp, 30) << "\n";
    std::cout << "  Error (HP)  = " << toString(error_hp, 10) << "\n";

    // For 256-bit precision, machine epsilon is ~2^-256 ≈ 1e-77
    std::cout << "  Expected HP epsilon: ~1e-77 for 256-bit precision\n";

    // Error should be near HP machine precision
    bool pass = error_hp < mpreal(1e-70);  // Should be ~1e-77, allow slack
    std::cout << "  " << (pass ? "PASSED" : "FAILED") << " (HP error < 1e-70)\n";
    return pass;
}

// Test 3: Gradient of 2D polynomial
bool test_gradient_2d() {
    std::cout << "\nTest: Gradient of 2D polynomial\n";

    PrecisionContext ctx(256);

    // Use a simple polynomial that we can verify with double precision first
    // f(x,y) = x*y (degree [1,1])
    // Bernstein coefficients: [0, 0, 0, 1] (tensor product of [0,1] and [0,1])
    std::vector<double> coeffs_double = {0.0, 0.0, 0.0, 1.0};
    Polynomial poly(std::vector<unsigned int>{1, 1}, coeffs_double);

    // First verify with double precision
    std::vector<Polynomial> grad_double = Differentiation::gradient(poly);
    double df_dx_double = grad_double[0].evaluate(std::vector<double>{0.5, 0.5});
    double df_dy_double = grad_double[1].evaluate(std::vector<double>{0.5, 0.5});

    // Convert to HP and compute gradient
    PolynomialHP poly_hp(poly);
    std::vector<PolynomialHP> grad_hp = DifferentiationHP::gradient(poly_hp);

    // Evaluate at (0.5, 0.5)
    std::vector<mpreal> point_hp = {toHighPrecision(0.5), toHighPrecision(0.5)};
    mpreal df_dx_hp = grad_hp[0].evaluate(point_hp);
    mpreal df_dy_hp = grad_hp[1].evaluate(point_hp);

    // Expected: ∇f = [y, x], so ∇f(0.5, 0.5) = [0.5, 0.5]
    double expected_x = 0.5;
    double expected_y = 0.5;

    std::cout << "  f(x,y) = x*y\n";
    std::cout << "  ∂f/∂x(0.5, 0.5) double = " << df_dx_double << "\n";
    std::cout << "  ∂f/∂y(0.5, 0.5) double = " << df_dy_double << "\n";
    std::cout << "  ∂f/∂x(0.5, 0.5) HP     = " << toString(df_dx_hp, 20) << "\n";
    std::cout << "  ∂f/∂y(0.5, 0.5) HP     = " << toString(df_dy_hp, 20) << "\n";
    std::cout << "  Expected               = [" << expected_x << ", " << expected_y << "]\n";

    bool pass = std::abs(toDouble(df_dx_hp) - expected_x) < 1e-10 &&
                std::abs(toDouble(df_dy_hp) - expected_y) < 1e-10;
    std::cout << "  " << (pass ? "PASSED" : "FAILED") << "\n";
    return pass;
}

// Test 4: derivativeFromDouble convenience function
bool test_derivative_from_double() {
    std::cout << "\nTest: derivativeFromDouble convenience function\n";

    PrecisionContext ctx(256);

    // Polynomial: f(x) = x^2
    // Bernstein coefficients for x^2 on [0,1]: [0, 0, 1]
    // (Not [0, 0.5, 1] - that would be x, not x^2)
    std::vector<double> coeffs_double = {0.0, 0.0, 1.0};
    Polynomial poly(std::vector<unsigned int>{2}, coeffs_double);

    // First verify with double precision
    Polynomial dpoly_double = Differentiation::derivative(poly, 0, 1);
    double df_double = dpoly_double.evaluate(0.7);

    // Use convenience function
    PolynomialHP dpoly_hp = DifferentiationHP::derivativeFromDouble(poly, 0, 1);

    // Evaluate at x = 0.7
    mpreal x_hp = toHighPrecision(0.7);
    mpreal df_hp = dpoly_hp.evaluate(x_hp);

    // Expected: f'(x) = 2x, so f'(0.7) = 1.4
    double expected = 1.4;

    std::cout << "  f(x) = x^2\n";
    std::cout << "  f'(0.7) double = " << df_double << "\n";
    std::cout << "  f'(0.7) HP     = " << toString(df_hp, 20) << "\n";
    std::cout << "  Expected       = " << expected << "\n";
    std::cout << "  Error (double) = " << std::abs(df_double - expected) << "\n";
    std::cout << "  Error (HP)     = " << toDouble(abs(df_hp - toHighPrecision(expected))) << "\n";

    bool pass = std::abs(toDouble(df_hp) - expected) < 1e-10;
    std::cout << "  " << (pass ? "PASSED" : "FAILED") << "\n";
    return pass;
}

int main() {
    std::cout << "=== DifferentiationHP Tests ===\n";
    
    int passed = 0;
    int total = 4;
    
    if (test_first_derivative_1d()) ++passed;
    if (test_second_derivative_1d()) ++passed;
    if (test_gradient_2d()) ++passed;
    if (test_derivative_from_double()) ++passed;
    
    std::cout << "\n=== " << passed << "/" << total << " tests passed ===\n";
    
    return (passed == total) ? 0 : 1;
}

#else

int main() {
    std::cout << "High-precision support not enabled\n";
    return 0;
}

#endif // ENABLE_HIGH_PRECISION

