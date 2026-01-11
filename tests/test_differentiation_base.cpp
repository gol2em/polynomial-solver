/**
 * @file test_differentiation_base.cpp
 * @brief Tests for templated differentiation in PolynomialBase
 *
 * Tests the differentiate(), derivative(), gradient(), and hessian() methods
 * for both Bernstein and Power bases.
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

#include "core/polynomial_base.h"

using namespace polynomial_solver;
using PolyDouble = PolynomialBase<double>;
using DerivCacheDouble = DerivativeCacheBase<double>;

static bool approx_equal(double a, double b, double tol = 1e-10) {
    return std::fabs(a - b) < tol;
}

//=============================================================================
// Test 1: Derivative of x^3 (univariate, power basis)
//=============================================================================
void test_derivative_cubic_power() {
    std::cout << "  test_derivative_cubic_power... ";

    // p(x) = x^3 in power basis
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{0.0, 0.0, 0.0, 1.0}; // x^3

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    // First derivative: 3x^2
    PolyDouble dp = p.differentiate(0);

    // Check at several points
    std::vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double x : test_points) {
        double expected = 3.0 * x * x;
        double actual = dp.evaluate(x);
        assert(approx_equal(actual, expected));
    }

    // Second derivative: 6x
    PolyDouble d2p = p.derivative(0, 2);
    for (double x : test_points) {
        double expected = 6.0 * x;
        double actual = d2p.evaluate(x);
        assert(approx_equal(actual, expected));
    }

    // Third derivative: 6 (constant)
    PolyDouble d3p = p.derivative(0, 3);
    for (double x : test_points) {
        double expected = 6.0;
        double actual = d3p.evaluate(x);
        assert(approx_equal(actual, expected));
    }

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 2: Derivative of Bernstein polynomial
//=============================================================================
void test_derivative_bernstein() {
    std::cout << "  test_derivative_bernstein... ";

    // Create p(x) = x^2 = sum of Bernstein coefficients [0, 0, 1] for degree 2
    // Actually, x^2 in Bernstein basis of degree 2 is [0, 0, 1]
    // Let's use p(x) = x which is [0, 1] in degree 1 Bernstein
    std::vector<unsigned int> degrees{1u};
    std::vector<double> bern_coeffs{0.0, 1.0}; // B_0,1(x)=1-x, B_1,1(x)=x => p(x) = x

    PolyDouble p(degrees, bern_coeffs);

    // Derivative of x is 1
    PolyDouble dp = p.differentiate(0);

    for (double x = 0.0; x <= 1.0; x += 0.1) {
        double actual = dp.evaluate(x);
        assert(approx_equal(actual, 1.0));
    }

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 3: Gradient of multivariate polynomial
//=============================================================================
void test_gradient_2d() {
    std::cout << "  test_gradient_2d... ";

    // p(x,y) = x^2 + y^2 in power basis
    // Coefficients for degree (2,2):
    // c[i][j] = coefficient of x^i * y^j
    // We want c[2][0] = 1, c[0][2] = 1, rest = 0
    std::vector<unsigned int> degrees{2u, 2u};
    // Layout: [c00, c01, c02, c10, c11, c12, c20, c21, c22] (last dim fastest)
    std::vector<double> coeffs(9, 0.0);
    coeffs[0*3 + 0] = 0.0;  // constant
    coeffs[2*3 + 0] = 1.0;  // x^2
    coeffs[0*3 + 2] = 1.0;  // y^2

    PolyDouble p = PolyDouble::fromPower(degrees, coeffs);

    auto grad = p.gradient();
    assert(grad.size() == 2);

    // ∂p/∂x = 2x
    // ∂p/∂y = 2y
    std::vector<double> test_point{0.5, 0.3};
    double dp_dx = grad[0].evaluate(test_point);
    double dp_dy = grad[1].evaluate(test_point);

    assert(approx_equal(dp_dx, 2.0 * 0.5));  // 2x = 1.0
    assert(approx_equal(dp_dy, 2.0 * 0.3));  // 2y = 0.6

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 4: Derivative of constant is zero
//=============================================================================
void test_derivative_constant() {
    std::cout << "  test_derivative_constant... ";

    // p(x) = 5 (constant)
    std::vector<unsigned int> degrees{0u};
    std::vector<double> coeffs{5.0};

    PolyDouble p = PolyDouble::fromPower(degrees, coeffs);
    PolyDouble dp = p.differentiate(0);

    // Derivative of constant is 0
    assert(approx_equal(dp.evaluate(0.5), 0.0));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 5: Hessian of quadratic
//=============================================================================
void test_hessian_2d() {
    std::cout << "  test_hessian_2d... ";

    // p(x,y) = x^2 + 2xy + 3y^2
    // ∂p/∂x = 2x + 2y
    // ∂p/∂y = 2x + 6y
    // ∂²p/∂x² = 2
    // ∂²p/∂y² = 6
    // ∂²p/∂x∂y = 2
    std::vector<unsigned int> degrees{2u, 2u};
    std::vector<double> coeffs(9, 0.0);
    coeffs[2*3 + 0] = 1.0;  // x^2
    coeffs[1*3 + 1] = 2.0;  // 2xy
    coeffs[0*3 + 2] = 3.0;  // 3y^2

    PolyDouble p = PolyDouble::fromPower(degrees, coeffs);
    auto hess = p.hessian();

    assert(hess.size() == 2);
    assert(hess[0].size() == 2);

    std::vector<double> pt{0.5, 0.3};

    // Check Hessian values (should be constant for quadratic)
    assert(approx_equal(hess[0][0].evaluate(pt), 2.0));  // ∂²p/∂x²
    assert(approx_equal(hess[1][1].evaluate(pt), 6.0));  // ∂²p/∂y²
    assert(approx_equal(hess[0][1].evaluate(pt), 2.0));  // ∂²p/∂x∂y
    assert(approx_equal(hess[1][0].evaluate(pt), 2.0));  // ∂²p/∂y∂x (symmetric)

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 6: Region preservation after differentiation
//=============================================================================
void test_region_preservation() {
    std::cout << "  test_region_preservation... ";

    std::vector<unsigned int> degrees{2u};
    std::vector<double> coeffs{1.0, 2.0, 3.0};

    PolyDouble p(degrees, coeffs);  // Bernstein basis

    // Subdivide to create a non-trivial region
    PolyDouble sub = p.restrictedToInterval(0, 0.25, 0.75);

    auto [lower, upper] = sub.getOriginalBox();
    assert(approx_equal(lower[0], 0.25));
    assert(approx_equal(upper[0], 0.75));

    // Differentiate - region should be preserved
    PolyDouble dsub = sub.differentiate(0);

    auto [d_lower, d_upper] = dsub.getOriginalBox();
    assert(approx_equal(d_lower[0], 0.25));
    assert(approx_equal(d_upper[0], 0.75));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 7: Compare derivative evaluation
//=============================================================================
void test_compare_derivative_eval() {
    std::cout << "  test_compare_derivative_eval... ";

    // Create polynomial: 1 + 2x + 3x^2 + 4x^3
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{1.0, 2.0, 3.0, 4.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    // First derivative: 2 + 6x + 12x^2
    PolyDouble dp = p.differentiate(0);

    for (double x = 0.0; x <= 1.0; x += 0.1) {
        double expected = 2.0 + 6.0*x + 12.0*x*x;
        double actual = dp.evaluate(x);
        assert(approx_equal(actual, expected, 1e-9));
    }

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 8: Order-0 derivative returns copy
//=============================================================================
void test_order_zero() {
    std::cout << "  test_order_zero... ";

    std::vector<unsigned int> degrees{2u};
    std::vector<double> coeffs{1.0, 2.0, 3.0};

    PolyDouble p = PolyDouble::fromPower(degrees, coeffs);
    PolyDouble d0 = p.derivative(0, 0);  // Order 0 = identity

    for (double x = 0.0; x <= 1.0; x += 0.1) {
        assert(approx_equal(p.evaluate(x), d0.evaluate(x)));
    }

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 9: DerivativeCacheBase
//=============================================================================
void test_derivative_cache() {
    std::cout << "  test_derivative_cache... ";

    // p(x,y) = x^2 + xy + y^2
    std::vector<unsigned int> degrees{2u, 2u};
    std::vector<double> coeffs(9, 0.0);
    coeffs[2*3 + 0] = 1.0;  // x^2
    coeffs[1*3 + 1] = 1.0;  // xy
    coeffs[0*3 + 2] = 1.0;  // y^2

    PolyDouble p = PolyDouble::fromPower(degrees, coeffs);

    DerivativeCacheBase<double> cache(p);

    // Test getPartial
    const PolyDouble& dp_dx = cache.getPartial(0, 1);  // ∂p/∂x = 2x + y
    const PolyDouble& dp_dy = cache.getPartial(1, 1);  // ∂p/∂y = x + 2y

    std::vector<double> pt{0.5, 0.3};
    assert(approx_equal(dp_dx.evaluate(pt), 2.0*0.5 + 0.3));  // 1.3
    assert(approx_equal(dp_dy.evaluate(pt), 0.5 + 2.0*0.3));  // 1.1

    // Test get with multi-index
    std::vector<unsigned int> orders{1, 1};  // ∂²p/∂x∂y = 1
    const PolyDouble& d2p_dxdy = cache.get(orders);
    assert(approx_equal(d2p_dxdy.evaluate(pt), 1.0));

    // Test precompute
    cache.precomputeUpToOrder(3);

    // After precompute, third derivatives should be zero (quadratic polynomial)
    std::vector<unsigned int> third_x{3, 0};
    const PolyDouble& d3p_dx3 = cache.get(third_x);
    assert(approx_equal(d3p_dx3.evaluate(pt), 0.0));

    std::cout << "PASSED\n";
}

//=============================================================================
// Main
//=============================================================================
int main() {
    std::cout << "=== Testing PolynomialBase Differentiation ===" << std::endl;

    test_derivative_cubic_power();
    test_derivative_bernstein();
    test_gradient_2d();
    test_derivative_constant();
    test_hessian_2d();
    test_region_preservation();
    test_compare_derivative_eval();
    test_order_zero();
    test_derivative_cache();

    std::cout << "\n=== All Differentiation Tests PASSED ===" << std::endl;
    return 0;
}
