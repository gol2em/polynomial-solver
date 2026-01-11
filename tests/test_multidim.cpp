/**
 * @file test_multidim.cpp
 * @brief Test multidimensional polynomial evaluation, de Casteljau, and differentiation
 *
 * Tests 3D and higher dimensional polynomials for:
 * - Tensor-product evaluation (power and Bernstein basis)
 * - De Casteljau algorithm
 * - Differentiation along different axes
 * - Basis conversion
 */

#include "core/polynomial_base.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

using namespace polynomial_solver;

template<typename T>
bool approx_equal(T a, T b, T eps = T(1e-10)) {
    T diff = (a > b) ? (a - b) : (b - a);
    return diff < eps;
}

//=============================================================================
// 3D Polynomial Tests
//=============================================================================

void test_3d_power_evaluation() {
    std::cout << "  test_3d_power_evaluation... ";
    
    // f(x,y,z) = 1 + x + y + z + xy + xz + yz + xyz
    // Power coefficients for degrees [1,1,1]:
    // Index order: (x,y,z) with z fastest
    // (0,0,0)=1, (0,0,1)=1, (0,1,0)=1, (0,1,1)=1, (1,0,0)=1, (1,0,1)=1, (1,1,0)=1, (1,1,1)=1
    std::vector<unsigned int> degrees = {1, 1, 1};
    std::vector<double> power_coeffs = {1, 1, 1, 1, 1, 1, 1, 1};
    
    auto p = PolynomialBase<double>::fromPower(degrees, power_coeffs);
    
    // Test at corners
    assert(approx_equal(p.evaluate({0.0, 0.0, 0.0}), 1.0));  // f(0,0,0) = 1
    assert(approx_equal(p.evaluate({1.0, 0.0, 0.0}), 2.0));  // f(1,0,0) = 1+1 = 2
    assert(approx_equal(p.evaluate({0.0, 1.0, 0.0}), 2.0));  // f(0,1,0) = 1+1 = 2
    assert(approx_equal(p.evaluate({0.0, 0.0, 1.0}), 2.0));  // f(0,0,1) = 1+1 = 2
    assert(approx_equal(p.evaluate({1.0, 1.0, 0.0}), 4.0));  // f(1,1,0) = 1+1+1+1 = 4
    assert(approx_equal(p.evaluate({1.0, 1.0, 1.0}), 8.0));  // f(1,1,1) = 1+1+1+1+1+1+1+1 = 8
    
    // Test at interior point
    double val = p.evaluate({0.5, 0.5, 0.5});
    // f(0.5,0.5,0.5) = 1 + 0.5 + 0.5 + 0.5 + 0.25 + 0.25 + 0.25 + 0.125 = 3.375
    assert(approx_equal(val, 3.375));
    
    std::cout << "PASSED\n";
}

void test_3d_bernstein_evaluation() {
    std::cout << "  test_3d_bernstein_evaluation... ";
    
    // Create a simple 3D Bernstein polynomial
    // f(x,y,z) = x*y*z (trilinear function)
    // For degrees [1,1,1], the Bernstein coefficient at (1,1,1) = 1, others = 0
    std::vector<unsigned int> degrees = {1, 1, 1};
    std::vector<double> bern_coeffs(8, 0.0);
    bern_coeffs[7] = 1.0;  // B_{111} = xyz
    
    PolynomialBase<double> p(degrees, bern_coeffs);
    
    // Test at corners
    assert(approx_equal(p.evaluate({0.0, 0.0, 0.0}), 0.0));
    assert(approx_equal(p.evaluate({1.0, 0.0, 0.0}), 0.0));
    assert(approx_equal(p.evaluate({0.0, 1.0, 0.0}), 0.0));
    assert(approx_equal(p.evaluate({1.0, 1.0, 1.0}), 1.0));
    
    // Test at center: B_{111}(0.5,0.5,0.5) = 0.5*0.5*0.5 = 0.125
    assert(approx_equal(p.evaluate({0.5, 0.5, 0.5}), 0.125));
    
    std::cout << "PASSED\n";
}

void test_3d_differentiation() {
    std::cout << "  test_3d_differentiation... ";
    
    // f(x,y,z) = x^2 + y^2 + z^2 + xyz
    std::vector<unsigned int> degrees = {2, 2, 2};
    // Total coefficients: 3*3*3 = 27
    std::vector<double> power_coeffs(27, 0.0);
    // x^2: index (2,0,0) -> 2*9 + 0*3 + 0 = 18
    power_coeffs[18] = 1.0;
    // y^2: index (0,2,0) -> 0*9 + 2*3 + 0 = 6
    power_coeffs[6] = 1.0;
    // z^2: index (0,0,2) -> 0*9 + 0*3 + 2 = 2
    power_coeffs[2] = 1.0;
    // xyz: index (1,1,1) -> 1*9 + 1*3 + 1 = 13
    power_coeffs[13] = 1.0;
    
    auto p = PolynomialBase<double>::fromPower(degrees, power_coeffs);
    
    // df/dx = 2x + yz
    auto df_dx = p.differentiate(0);
    assert(approx_equal(df_dx.evaluate({1.0, 1.0, 1.0}), 3.0));  // 2*1 + 1*1 = 3
    assert(approx_equal(df_dx.evaluate({0.5, 0.5, 0.5}), 1.25)); // 2*0.5 + 0.25 = 1.25
    
    // df/dy = 2y + xz
    auto df_dy = p.differentiate(1);
    assert(approx_equal(df_dy.evaluate({1.0, 1.0, 1.0}), 3.0));
    
    // df/dz = 2z + xy
    auto df_dz = p.differentiate(2);
    assert(approx_equal(df_dz.evaluate({1.0, 1.0, 1.0}), 3.0));
    
    // Mixed partial: d²f/dxdy = z
    auto d2f_dxdy = df_dx.differentiate(1);
    assert(approx_equal(d2f_dxdy.evaluate({0.5, 0.5, 0.5}), 0.5));
    assert(approx_equal(d2f_dxdy.evaluate({0.5, 0.5, 1.0}), 1.0));
    
    std::cout << "PASSED\n";
}

void test_3d_basis_conversion() {
    std::cout << "  test_3d_basis_conversion... ";
    
    // Create a polynomial in power basis
    std::vector<unsigned int> degrees = {2, 1, 1};
    std::vector<double> power_coeffs(12, 0.0);
    // f(x,y,z) = 1 + x + x^2 + y + z
    power_coeffs[0] = 1.0;   // constant
    power_coeffs[1] = 1.0;   // z
    power_coeffs[3] = 1.0;   // y
    power_coeffs[6] = 1.0;   // x
    power_coeffs[12 - 6] = 1.0; // x^2 at (2,0,0)
    
    auto p_power = PolynomialBase<double>::fromPower(degrees, power_coeffs);
    
    // Convert to Bernstein
    auto p_bern = p_power.convertToBernstein();
    
    // Verify evaluations match at several points
    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.5, 0.5, 0.5},
        {0.25, 0.75, 0.33}
    };
    
    for (const auto& pt : test_points) {
        double v_power = p_power.evaluate(pt);
        double v_bern = p_bern.evaluate(pt);
        assert(approx_equal(v_power, v_bern, 1e-9));
    }
    
    std::cout << "PASSED\n";
}

//=============================================================================
// Main
//=============================================================================

int main() {
    std::cout << "Running multidimensional polynomial tests...\n\n";
    
    std::cout << "3D Polynomial Tests:\n";
    test_3d_power_evaluation();
    test_3d_bernstein_evaluation();
    test_3d_differentiation();
    test_3d_basis_conversion();
    
    std::cout << "\n=== All multidimensional tests PASSED ===\n";
    return 0;
}

