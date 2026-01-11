/**
 * @file test_polynomial_base.cpp
 * @brief Unit tests for PolynomialBase<Scalar> template class (Tier 3)
 *
 * Tests the templated polynomial implementation to ensure it has equivalent
 * functionality to the existing Polynomial class and works with different
 * scalar types.
 */

#include "core/polynomial_base.h"
// polynomial.h is included by polynomial_base.h for PolynomialRepresentation enum

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace polynomial_solver;

// Type alias for testing
using PolyDouble = PolynomialBase<double>;

// Tolerance for floating point comparisons
constexpr double EPSILON = 1e-10;

bool approx_equal(double a, double b, double tol = EPSILON) {
    return std::fabs(a - b) < tol;
}

//=============================================================================
// Test: Construction and Basic Accessors
//=============================================================================

void test_default_constructor() {
    std::cout << "  test_default_constructor... ";

    PolyDouble p;
    assert(p.dimension() == 0);
    assert(p.empty());

    std::cout << "PASSED\n";
}

void test_bernstein_constructor() {
    std::cout << "  test_bernstein_constructor... ";

    // 1D quadratic polynomial with Bernstein coefficients [1, 2, 3]
    std::vector<unsigned int> degrees = {2};
    std::vector<double> coeffs = {1.0, 2.0, 3.0};

    PolyDouble p(degrees, coeffs);

    assert(p.dimension() == 1);
    assert(p.degrees().size() == 1);
    assert(p.degrees()[0] == 2);
    assert(p.coefficientCount() == 3);
    assert(!p.empty());
    assert(p.hasBernsteinCoefficients());
    assert(!p.hasPowerCoefficients());
    assert(p.primaryRepresentation() == PolynomialRepresentation::BERNSTEIN);

    std::cout << "PASSED\n";
}

void test_from_power_factory() {
    std::cout << "  test_from_power_factory... ";

    // 1D polynomial: f(x) = 1 + 2x + 3x^2 (power basis)
    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, 2.0, 3.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    assert(p.dimension() == 1);
    assert(p.degrees()[0] == 2);
    assert(p.hasPowerCoefficients());
    assert(!p.hasBernsteinCoefficients());
    assert(p.primaryRepresentation() == PolynomialRepresentation::POWER);

    std::cout << "PASSED\n";
}

void test_from_bernstein_factory() {
    std::cout << "  test_from_bernstein_factory... ";

    std::vector<unsigned int> degrees = {2};
    std::vector<double> bern_coeffs = {1.0, 2.0, 3.0};

    PolyDouble p = PolyDouble::fromBernstein(degrees, bern_coeffs);

    assert(p.dimension() == 1);
    assert(p.hasBernsteinCoefficients());
    assert(p.primaryRepresentation() == PolynomialRepresentation::BERNSTEIN);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test: Evaluation
//=============================================================================

void test_evaluate_power_1d() {
    std::cout << "  test_evaluate_power_1d... ";

    // f(x) = 1 + 2x + 3x^2
    // f(0) = 1, f(1) = 6, f(0.5) = 1 + 1 + 0.75 = 2.75
    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, 2.0, 3.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    assert(approx_equal(p.evaluate(0.0), 1.0));
    assert(approx_equal(p.evaluate(1.0), 6.0));
    assert(approx_equal(p.evaluate(0.5), 2.75));

    std::cout << "PASSED\n";
}

void test_evaluate_bernstein_1d() {
    std::cout << "  test_evaluate_bernstein_1d... ";

    // Bernstein polynomial with coefficients [0, 1, 0] and degree 2
    // B_0^2(t) = (1-t)^2, B_1^2(t) = 2t(1-t), B_2^2(t) = t^2
    // f(t) = 0*(1-t)^2 + 1*2t(1-t) + 0*t^2 = 2t(1-t)
    // f(0) = 0, f(1) = 0, f(0.5) = 0.5
    std::vector<unsigned int> degrees = {2};
    std::vector<double> bern_coeffs = {0.0, 1.0, 0.0};

    PolyDouble p = PolyDouble::fromBernstein(degrees, bern_coeffs);

    assert(approx_equal(p.evaluate(0.0), 0.0));
    assert(approx_equal(p.evaluate(1.0), 0.0));
    assert(approx_equal(p.evaluate(0.5), 0.5));

    std::cout << "PASSED\n";
}

void test_evaluate_2d() {
    std::cout << "  test_evaluate_2d... ";

    // f(x,y) = 1 + x + y (power basis)
    // degrees = [1, 1], coefficients arranged as:
    // [c_00, c_01, c_10, c_11] = [1, 1, 1, 0]
    // (x^0 y^0)=1, (x^0 y^1)=1, (x^1 y^0)=1, (x^1 y^1)=0
    std::vector<unsigned int> degrees = {1, 1};
    std::vector<double> power_coeffs = {1.0, 1.0, 1.0, 0.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    assert(approx_equal(p.evaluate({0.0, 0.0}), 1.0));
    assert(approx_equal(p.evaluate({1.0, 0.0}), 2.0));
    assert(approx_equal(p.evaluate({0.0, 1.0}), 2.0));
    assert(approx_equal(p.evaluate({1.0, 1.0}), 3.0));
    assert(approx_equal(p.evaluate({0.5, 0.5}), 2.0));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test: Basis Conversion
//=============================================================================

void test_power_to_bernstein_conversion() {
    std::cout << "  test_power_to_bernstein_conversion... ";

    // f(x) = 1 + 2x + 3x^2 (power basis)
    // Expected Bernstein coefficients for degree 2:
    // b_0 = f(0) = 1
    // b_1 = f(0) + f'(0)/2 = 1 + 2/2 = 2
    // b_2 = f(1) = 1 + 2 + 3 = 6
    // Actually: b_k = sum_{i=0}^k a_i * C(k,i) / C(n,i)
    // b_0 = a_0 = 1
    // b_1 = a_0 + a_1 * C(1,1)/C(2,1) = 1 + 2 * 1/2 = 2
    // b_2 = a_0 + a_1 + a_2 = 1 + 2 + 3 = 6

    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, 2.0, 3.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    // Trigger conversion
    p.convertPowerToBernstein();

    assert(p.hasBernsteinCoefficients());
    const std::vector<double>& bern = p.bernsteinCoefficients();
    assert(bern.size() == 3);
    assert(approx_equal(bern[0], 1.0));
    assert(approx_equal(bern[1], 2.0));
    assert(approx_equal(bern[2], 6.0));

    std::cout << "PASSED\n";
}

void test_ensure_bernstein_primary() {
    std::cout << "  test_ensure_bernstein_primary... ";

    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, 2.0, 3.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);
    assert(p.primaryRepresentation() == PolynomialRepresentation::POWER);

    // New design: conversion returns new polynomial, or modifies in-place
    p.ensureBernsteinPrimary();

    assert(p.primaryRepresentation() == PolynomialRepresentation::BERNSTEIN);
    assert(p.hasBernsteinCoefficients());
    // New design: single representation, no dual storage
    // assert(p.hasPowerCoefficients()); // Old dual-rep test, no longer valid

    std::cout << "PASSED\n";
}

//=============================================================================
// Test: Arithmetic Operations
//=============================================================================

void test_polynomial_addition() {
    std::cout << "  test_polynomial_addition... ";

    // p(x) = 1 + 2x
    // q(x) = 3 + 4x
    // (p + q)(x) = 4 + 6x
    std::vector<unsigned int> degrees = {1};

    PolyDouble p = PolyDouble::fromPower(degrees, {1.0, 2.0});
    PolyDouble q = PolyDouble::fromPower(degrees, {3.0, 4.0});

    PolyDouble sum = p + q;

    assert(approx_equal(sum.evaluate(0.0), 4.0));
    assert(approx_equal(sum.evaluate(1.0), 10.0));
    assert(approx_equal(sum.evaluate(0.5), 7.0));

    std::cout << "PASSED\n";
}

void test_polynomial_subtraction() {
    std::cout << "  test_polynomial_subtraction... ";

    // p(x) = 5 + 3x
    // q(x) = 2 + x
    // (p - q)(x) = 3 + 2x
    std::vector<unsigned int> degrees = {1};

    PolyDouble p = PolyDouble::fromPower(degrees, {5.0, 3.0});
    PolyDouble q = PolyDouble::fromPower(degrees, {2.0, 1.0});

    PolyDouble diff = p - q;

    assert(approx_equal(diff.evaluate(0.0), 3.0));
    assert(approx_equal(diff.evaluate(1.0), 5.0));

    std::cout << "PASSED\n";
}

void test_polynomial_multiplication() {
    std::cout << "  test_polynomial_multiplication... ";

    // p(x) = 1 + x
    // q(x) = 2 + 3x
    // (p * q)(x) = 2 + 5x + 3x^2
    std::vector<unsigned int> degrees = {1};

    PolyDouble p = PolyDouble::fromPower(degrees, {1.0, 1.0});
    PolyDouble q = PolyDouble::fromPower(degrees, {2.0, 3.0});

    PolyDouble prod = p * q;

    assert(prod.degrees()[0] == 2); // degree 1 + 1
    assert(approx_equal(prod.evaluate(0.0), 2.0));
    assert(approx_equal(prod.evaluate(1.0), 10.0)); // 2 + 5 + 3

    std::cout << "PASSED\n";
}

void test_scalar_multiplication() {
    std::cout << "  test_scalar_multiplication... ";

    // p(x) = 1 + 2x
    // 3 * p(x) = 3 + 6x
    std::vector<unsigned int> degrees = {1};

    PolyDouble p = PolyDouble::fromPower(degrees, {1.0, 2.0});
    PolyDouble scaled = p * 3.0;

    assert(approx_equal(scaled.evaluate(0.0), 3.0));
    assert(approx_equal(scaled.evaluate(1.0), 9.0));

    std::cout << "PASSED\n";
}

void test_unary_negation() {
    std::cout << "  test_unary_negation... ";

    std::vector<unsigned int> degrees = {1};
    PolyDouble p = PolyDouble::fromPower(degrees, {1.0, 2.0});
    PolyDouble neg = -p;

    assert(approx_equal(neg.evaluate(0.0), -1.0));
    assert(approx_equal(neg.evaluate(1.0), -3.0));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test: Subdivision / Restriction
//=============================================================================

void test_restricted_to_interval() {
    std::cout << "  test_restricted_to_interval... ";

    // p(x) = x (identity in Bernstein: coefficients [0, 1])
    std::vector<unsigned int> degrees = {1};
    std::vector<double> bern_coeffs = {0.0, 1.0};

    PolyDouble p = PolyDouble::fromBernstein(degrees, bern_coeffs);

    // Restrict to [0.25, 0.75]
    // q(u) = p(0.25 + 0.5*u) = 0.25 + 0.5*u
    // In Bernstein on [0,1]: q(0) = 0.25, q(1) = 0.75
    PolyDouble q = p.restrictedToInterval(0, 0.25, 0.75);

    assert(approx_equal(q.evaluate(0.0), 0.25, 1e-6));
    assert(approx_equal(q.evaluate(1.0), 0.75, 1e-6));
    assert(approx_equal(q.evaluate(0.5), 0.5, 1e-6));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test: Graph Control Points
//=============================================================================

void test_graph_control_points_1d() {
    std::cout << "  test_graph_control_points_1d... ";

    // Bernstein coefficients [1, 2, 3] for degree 2
    // Control points: (0/2, 1), (1/2, 2), (2/2, 3) = (0, 1), (0.5, 2), (1, 3)
    std::vector<unsigned int> degrees = {2};
    std::vector<double> bern_coeffs = {1.0, 2.0, 3.0};

    PolyDouble p = PolyDouble::fromBernstein(degrees, bern_coeffs);

    std::vector<double> cpts;
    p.graphControlPoints(cpts);

    // 3 control points, each with 2 coordinates (t, f)
    assert(cpts.size() == 6);

    // Point 0: (0, 1)
    assert(approx_equal(cpts[0], 0.0));
    assert(approx_equal(cpts[1], 1.0));

    // Point 1: (0.5, 2)
    assert(approx_equal(cpts[2], 0.5));
    assert(approx_equal(cpts[3], 2.0));

    // Point 2: (1, 3)
    assert(approx_equal(cpts[4], 1.0));
    assert(approx_equal(cpts[5], 3.0));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test: Graph Corners
//=============================================================================

void test_graph_corners_1d() {
    std::cout << "  test_graph_corners_1d... ";

    // p(x) = 2x (power basis)
    std::vector<unsigned int> degrees = {1};
    std::vector<double> power_coeffs = {0.0, 2.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    std::vector<std::vector<double>> corners;
    p.graphCorners(corners);

    // 2 corners: (0, f(0))=(0, 0), (1, f(1))=(1, 2)
    assert(corners.size() == 2);
    assert(corners[0].size() == 2);

    assert(approx_equal(corners[0][0], 0.0));
    assert(approx_equal(corners[0][1], 0.0));

    assert(approx_equal(corners[1][0], 1.0));
    assert(approx_equal(corners[1][1], 2.0));

    std::cout << "PASSED\n";
}

void test_graph_corners_2d() {
    std::cout << "  test_graph_corners_2d... ";

    // p(x,y) = x + y (power basis)
    std::vector<unsigned int> degrees = {1, 1};
    std::vector<double> power_coeffs = {0.0, 1.0, 1.0, 0.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    std::vector<std::vector<double>> corners;
    p.graphCorners(corners);

    // 4 corners
    assert(corners.size() == 4);
    assert(corners[0].size() == 3);

    // Corner (0,0): f=0
    // Corner (1,0): f=1
    // Corner (0,1): f=1
    // Corner (1,1): f=2

    // Note: order is based on bit pattern: (x,y) = (bit0, bit1)
    // i=0: (0,0), i=1: (1,0), i=2: (0,1), i=3: (1,1)
    assert(approx_equal(corners[0][2], 0.0)); // f(0,0)
    assert(approx_equal(corners[1][2], 1.0)); // f(1,0)
    assert(approx_equal(corners[2][2], 1.0)); // f(0,1)
    assert(approx_equal(corners[3][2], 2.0)); // f(1,1)

    std::cout << "PASSED\n";
}

//=============================================================================
// Main
//=============================================================================

int main() {
    std::cout << "=== Testing PolynomialBase<double> (Tier 3 Template) ===\n\n";

    std::cout << "Construction and Basic Accessors:\n";
    test_default_constructor();
    test_bernstein_constructor();
    test_from_power_factory();
    test_from_bernstein_factory();

    std::cout << "\nEvaluation:\n";
    test_evaluate_power_1d();
    test_evaluate_bernstein_1d();
    test_evaluate_2d();

    std::cout << "\nBasis Conversion:\n";
    test_power_to_bernstein_conversion();
    test_ensure_bernstein_primary();

    std::cout << "\nArithmetic Operations:\n";
    test_polynomial_addition();
    test_polynomial_subtraction();
    test_polynomial_multiplication();
    test_scalar_multiplication();
    test_unary_negation();

    std::cout << "\nSubdivision / Restriction:\n";
    test_restricted_to_interval();

    std::cout << "\nGraph Control Points:\n";
    test_graph_control_points_1d();

    std::cout << "\nGraph Corners:\n";
    test_graph_corners_1d();
    test_graph_corners_2d();

    std::cout << "\n=== All Tests PASSED ===\n";
    return 0;
}
