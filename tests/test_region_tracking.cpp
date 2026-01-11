/**
 * @file test_region_tracking.cpp
 * @brief Test region tracking for PolynomialBase
 *
 * Verifies that:
 * 1. Default region is [0,1]^n (backward compatibility)
 * 2. Subdivision updates original_box correctly
 * 3. getOriginalBox returns correct values
 * 4. Multiple subdivisions track regions correctly
 */

#include "core/polynomial_base.h"
#include <iostream>
#include <cmath>
#include <cassert>

using namespace polynomial_solver;

bool approx_equal(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

void test_default_region() {
    std::cout << "  test_default_region... ";

    // Create polynomial - should default to [0,1]^n
    std::vector<unsigned int> degrees = {2};
    std::vector<double> bern_coeffs = {1.0, 2.0, 3.0};

    PolynomialBase<double> p(degrees, bern_coeffs);

    auto [lower, upper] = p.getOriginalBox();

    assert(lower.size() == 1);
    assert(upper.size() == 1);
    assert(approx_equal(lower[0], 0.0));
    assert(approx_equal(upper[0], 1.0));

    std::cout << "PASSED\n";
}

void test_subdivision_updates_region() {
    std::cout << "  test_subdivision_updates_region... ";

    // Create polynomial on [0,1]
    std::vector<unsigned int> degrees = {2};
    std::vector<double> bern_coeffs = {1.0, 2.0, 3.0};

    PolynomialBase<double> p(degrees, bern_coeffs);

    // Subdivide to [0.3, 0.7]
    auto sub = p.restrictedToInterval(0, 0.3, 0.7);

    auto [lower, upper] = sub.getOriginalBox();

    assert(lower.size() == 1);
    assert(upper.size() == 1);
    assert(approx_equal(lower[0], 0.3));
    assert(approx_equal(upper[0], 0.7));

    std::cout << "PASSED\n";
}

void test_multiple_subdivisions() {
    std::cout << "  test_multiple_subdivisions... ";

    // Create 2D polynomial on [0,1]^2
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> bern_coeffs = {
        1.0, 2.0, 3.0,
        2.0, 3.0, 4.0,
        3.0, 4.0, 5.0
    };

    PolynomialBase<double> p(degrees, bern_coeffs);

    // Subdivide axis 0 to [0.0, 0.5]
    auto sub1 = p.restrictedToInterval(0, 0.0, 0.5);

    auto [lower1, upper1] = sub1.getOriginalBox();
    assert(approx_equal(lower1[0], 0.0));
    assert(approx_equal(upper1[0], 0.5));
    assert(approx_equal(lower1[1], 0.0));
    assert(approx_equal(upper1[1], 1.0));

    // Subdivide axis 1 to [0.5, 1.0]
    auto sub2 = sub1.restrictedToInterval(1, 0.5, 1.0);

    auto [lower2, upper2] = sub2.getOriginalBox();
    assert(approx_equal(lower2[0], 0.0));
    assert(approx_equal(upper2[0], 0.5));
    assert(approx_equal(lower2[1], 0.5));
    assert(approx_equal(upper2[1], 1.0));

    std::cout << "PASSED\n";
}

void test_nested_subdivisions() {
    std::cout << "  test_nested_subdivisions... ";

    // Create 1D polynomial
    std::vector<unsigned int> degrees = {3};
    std::vector<double> bern_coeffs = {0.0, 1.0/3.0, 2.0/3.0, 1.0};

    PolynomialBase<double> p(degrees, bern_coeffs);

    // Subdivide [0,1] -> [0.2, 0.8]
    auto sub1 = p.restrictedToInterval(0, 0.2, 0.8);
    auto [l1, u1] = sub1.getOriginalBox();
    assert(approx_equal(l1[0], 0.2));
    assert(approx_equal(u1[0], 0.8));

    // Subdivide [0.2, 0.8] -> [0.3, 0.7] in local coords
    // In global coords: 0.2 + 0.3*(0.8-0.2) = 0.38
    //                   0.2 + 0.7*(0.8-0.2) = 0.62
    auto sub2 = sub1.restrictedToInterval(0, 0.3, 0.7);
    auto [l2, u2] = sub2.getOriginalBox();
    assert(approx_equal(l2[0], 0.38));
    assert(approx_equal(u2[0], 0.62));

    std::cout << "PASSED\n";
}

void test_fromPower_has_region() {
    std::cout << "  test_fromPower_has_region... ";

    // Create from power basis
    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, 2.0, 3.0};

    auto p = PolynomialBase<double>::fromPower(degrees, power_coeffs);

    auto [lower, upper] = p.getOriginalBox();

    assert(lower.size() == 1);
    assert(upper.size() == 1);
    assert(approx_equal(lower[0], 0.0));
    assert(approx_equal(upper[0], 1.0));

    std::cout << "PASSED\n";
}

int main() {
    std::cout << "=== Testing Region Tracking ===\n\n";

    test_default_region();
    test_subdivision_updates_region();
    test_multiple_subdivisions();
    test_nested_subdivisions();
    test_fromPower_has_region();

    std::cout << "\n=== All Region Tracking Tests PASSED ===\n";

    return 0;
}

