#include "core/differentiation.h"
#include "core/polynomial.h"
#include <iostream>
#include <cmath>
#include <vector>

/**
 * @file test_differentiation.cpp
 * @brief Test suite for differentiation of Bernstein polynomials
 */

using namespace polynomial_solver;

// Helper: Check if two values are approximately equal
bool approx_equal(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

// Test 1: Derivative of x^3 (univariate)
int test_derivative_cubic()
{
    std::cout << "Test 1: Derivative of x^3" << std::endl;
    
    // p(x) = x^3 in power basis
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{0.0, 0.0, 0.0, 1.0}; // x^3
    
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    
    // First derivative: 3x^2
    Polynomial dp = Differentiation::derivative(p, 0, 1);
    
    // Check at several points
    std::vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double x : test_points) {
        double expected = 3.0 * x * x;
        double actual = dp.evaluate(x);
        
        if (!approx_equal(actual, expected, 1e-9)) {
            std::cerr << "  FAIL at x=" << x << ": expected " << expected 
                      << ", got " << actual << std::endl;
            return 1;
        }
    }
    
    // Second derivative: 6x
    Polynomial d2p = Differentiation::derivative(p, 0, 2);
    for (double x : test_points) {
        double expected = 6.0 * x;
        double actual = d2p.evaluate(x);
        
        if (!approx_equal(actual, expected, 1e-9)) {
            std::cerr << "  FAIL (2nd deriv) at x=" << x << ": expected " << expected 
                      << ", got " << actual << std::endl;
            return 1;
        }
    }
    
    // Third derivative: 6 (constant)
    Polynomial d3p = Differentiation::derivative(p, 0, 3);
    for (double x : test_points) {
        double expected = 6.0;
        double actual = d3p.evaluate(x);
        
        if (!approx_equal(actual, expected, 1e-9)) {
            std::cerr << "  FAIL (3rd deriv) at x=" << x << ": expected " << expected 
                      << ", got " << actual << std::endl;
            return 1;
        }
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 2: Derivative of x^5 (higher degree)
int test_derivative_quintic()
{
    std::cout << "Test 2: Derivative of x^5" << std::endl;
    
    // p(x) = x^5
    std::vector<unsigned int> degrees{5u};
    std::vector<double> power_coeffs{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    
    // First derivative: 5x^4
    Polynomial dp = Differentiation::derivative(p, 0, 1);
    
    std::vector<double> test_points = {0.0, 0.3, 0.5, 0.7, 1.0};
    for (double x : test_points) {
        double x2 = x * x;
        double x4 = x2 * x2;
        double expected = 5.0 * x4;
        double actual = dp.evaluate(x);
        
        if (!approx_equal(actual, expected, 1e-8)) {
            std::cerr << "  FAIL at x=" << x << ": expected " << expected 
                      << ", got " << actual << std::endl;
            return 1;
        }
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 3: Gradient of 2D polynomial
int test_gradient_2d()
{
    std::cout << "Test 3: Gradient of f(x,y) = x^2 + y^2" << std::endl;

    // f(x,y) = x^2 + y^2
    // In Bernstein basis for degree (2,2):
    // We'll use power basis and convert
    std::vector<unsigned int> degrees{2u, 2u};

    // Power basis: coefficients for x^i * y^j
    // Layout: last dimension (y) varies fastest
    // Index order: (0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)
    std::vector<double> power_coeffs(9, 0.0);
    power_coeffs[6] = 1.0;  // x^2 * y^0 -> index (2,0) = 2*3 + 0 = 6
    power_coeffs[2] = 1.0;  // x^0 * y^2 -> index (0,2) = 0*3 + 2 = 2

    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    
    // Gradient: [2x, 2y]
    std::vector<Polynomial> grad = Differentiation::gradient(p);
    
    if (grad.size() != 2) {
        std::cerr << "  FAIL: gradient size should be 2, got " << grad.size() << std::endl;
        return 1;
    }
    
    // Test at several points
    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0}, {0.5, 0.5}, {0.3, 0.7}, {1.0, 0.0}, {0.0, 1.0}
    };
    
    for (const auto& pt : test_points) {
        double x = pt[0], y = pt[1];
        
        double df_dx_expected = 2.0 * x;
        double df_dy_expected = 2.0 * y;
        
        double df_dx_actual = grad[0].evaluate(pt);
        double df_dy_actual = grad[1].evaluate(pt);
        
        if (!approx_equal(df_dx_actual, df_dx_expected, 1e-9) ||
            !approx_equal(df_dy_actual, df_dy_expected, 1e-9)) {
            std::cerr << "  FAIL at (" << x << "," << y << "): expected (" 
                      << df_dx_expected << "," << df_dy_expected << "), got ("
                      << df_dx_actual << "," << df_dy_actual << ")" << std::endl;
            return 1;
        }
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 4: DerivativeCache with iterative computation
int test_derivative_cache()
{
    std::cout << "Test 4: DerivativeCache for x^4" << std::endl;

    // p(x) = x^4
    std::vector<unsigned int> degrees{4u};
    std::vector<double> power_coeffs{0.0, 0.0, 0.0, 0.0, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    DerivativeCache cache(p);

    // Get derivatives using cache
    const Polynomial& dp = cache.getPartial(0, 1);   // 4x^3
    const Polynomial& d2p = cache.getPartial(0, 2);  // 12x^2
    const Polynomial& d3p = cache.getPartial(0, 3);  // 24x
    const Polynomial& d4p = cache.getPartial(0, 4);  // 24

    std::vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};

    for (double x : test_points) {
        double x2 = x * x;
        double x3 = x2 * x;

        double expected_dp = 4.0 * x3;
        double expected_d2p = 12.0 * x2;
        double expected_d3p = 24.0 * x;
        double expected_d4p = 24.0;

        if (!approx_equal(dp.evaluate(x), expected_dp, 1e-9) ||
            !approx_equal(d2p.evaluate(x), expected_d2p, 1e-9) ||
            !approx_equal(d3p.evaluate(x), expected_d3p, 1e-9) ||
            !approx_equal(d4p.evaluate(x), expected_d4p, 1e-9)) {
            std::cerr << "  FAIL at x=" << x << std::endl;
            return 1;
        }
    }

    // Verify caching: getting the same derivative again should return same reference
    const Polynomial& dp_again = cache.getPartial(0, 1);
    if (&dp != &dp_again) {
        std::cerr << "  FAIL: cache not returning same reference" << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 5: Mixed partial derivatives (2D)
int test_mixed_partials()
{
    std::cout << "Test 5: Mixed partials of f(x,y) = x^2*y^2" << std::endl;

    // f(x,y) = x^2 * y^2
    std::vector<unsigned int> degrees{2u, 2u};
    std::vector<double> power_coeffs(9, 0.0);
    // x^2 * y^2 -> index (2,2) = 2*3 + 2 = 8
    power_coeffs[8] = 1.0;

    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    DerivativeCache cache(p);

    // ∂f/∂x = 2x*y^2
    const Polynomial& df_dx = cache.get({1, 0});

    // ∂f/∂y = 2x^2*y
    const Polynomial& df_dy = cache.get({0, 1});

    // ∂²f/∂x² = 2y^2
    const Polynomial& d2f_dx2 = cache.get({2, 0});

    // ∂²f/∂y² = 2x^2
    const Polynomial& d2f_dy2 = cache.get({0, 2});

    // ∂²f/∂x∂y = 4xy (mixed partial)
    const Polynomial& d2f_dxdy = cache.get({1, 1});

    std::vector<std::vector<double>> test_points = {
        {0.5, 0.5}, {0.3, 0.7}, {0.8, 0.2}
    };

    for (const auto& pt : test_points) {
        double x = pt[0], y = pt[1];
        double x2 = x * x, y2 = y * y;

        double expected_df_dx = 2.0 * x * y2;
        double expected_df_dy = 2.0 * x2 * y;
        double expected_d2f_dx2 = 2.0 * y2;
        double expected_d2f_dy2 = 2.0 * x2;
        double expected_d2f_dxdy = 4.0 * x * y;

        if (!approx_equal(df_dx.evaluate(pt), expected_df_dx, 1e-9) ||
            !approx_equal(df_dy.evaluate(pt), expected_df_dy, 1e-9) ||
            !approx_equal(d2f_dx2.evaluate(pt), expected_d2f_dx2, 1e-9) ||
            !approx_equal(d2f_dy2.evaluate(pt), expected_d2f_dy2, 1e-9) ||
            !approx_equal(d2f_dxdy.evaluate(pt), expected_d2f_dxdy, 1e-9)) {
            std::cerr << "  FAIL at (" << x << "," << y << ")" << std::endl;
            return 1;
        }
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 6: Hessian matrix
int test_hessian()
{
    std::cout << "Test 6: Hessian of f(x,y) = x^2 + y^2" << std::endl;

    // f(x,y) = x^2 + y^2
    std::vector<unsigned int> degrees{2u, 2u};
    std::vector<double> power_coeffs(9, 0.0);
    power_coeffs[6] = 1.0;  // x^2 * y^0 -> index (2,0) = 6
    power_coeffs[2] = 1.0;  // x^0 * y^2 -> index (0,2) = 2

    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    // Hessian should be [[2, 0], [0, 2]]
    std::vector<std::vector<Polynomial>> H = Differentiation::hessian(p);

    if (H.size() != 2 || H[0].size() != 2 || H[1].size() != 2) {
        std::cerr << "  FAIL: Hessian size should be 2x2" << std::endl;
        return 1;
    }

    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0}, {0.5, 0.5}, {0.3, 0.7}
    };

    for (const auto& pt : test_points) {
        // H[0][0] = ∂²f/∂x² = 2
        // H[0][1] = ∂²f/∂x∂y = 0
        // H[1][0] = ∂²f/∂y∂x = 0
        // H[1][1] = ∂²f/∂y² = 2

        if (!approx_equal(H[0][0].evaluate(pt), 2.0, 1e-9) ||
            !approx_equal(H[0][1].evaluate(pt), 0.0, 1e-9) ||
            !approx_equal(H[1][0].evaluate(pt), 0.0, 1e-9) ||
            !approx_equal(H[1][1].evaluate(pt), 2.0, 1e-9)) {
            std::cerr << "  FAIL at (" << pt[0] << "," << pt[1] << ")" << std::endl;
            return 1;
        }
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 7: Precompute up to order
int test_precompute()
{
    std::cout << "Test 7: Precompute derivatives up to order 3" << std::endl;

    // p(x) = x^5
    std::vector<unsigned int> degrees{5u};
    std::vector<double> power_coeffs{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    DerivativeCache cache(p);

    // Precompute all derivatives up to order 3
    cache.precomputeUpToOrder(3);

    // Now access them - should be already cached
    const Polynomial& d1 = cache.getPartial(0, 1);
    const Polynomial& d2 = cache.getPartial(0, 2);
    const Polynomial& d3 = cache.getPartial(0, 3);

    // Verify correctness at x=0.5
    double x = 0.5;
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x2 * x2;

    double expected_d1 = 5.0 * x4;
    double expected_d2 = 20.0 * x3;
    double expected_d3 = 60.0 * x2;

    if (!approx_equal(d1.evaluate(x), expected_d1, 1e-9) ||
        !approx_equal(d2.evaluate(x), expected_d2, 1e-9) ||
        !approx_equal(d3.evaluate(x), expected_d3, 1e-9)) {
        std::cerr << "  FAIL: precomputed derivatives incorrect" << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

int main()
{
    std::cout << "========================================" << std::endl;
    std::cout << "Differentiation Test Suite" << std::endl;
    std::cout << "========================================" << std::endl;

    int failures = 0;

    failures += test_derivative_cubic();
    failures += test_derivative_quintic();
    failures += test_gradient_2d();
    failures += test_derivative_cache();
    failures += test_mixed_partials();
    failures += test_hessian();
    failures += test_precompute();

    std::cout << "========================================" << std::endl;
    if (failures == 0) {
        std::cout << "All tests PASSED!" << std::endl;
    } else {
        std::cout << failures << " test(s) FAILED!" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    return failures;
}

