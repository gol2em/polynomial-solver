/**
 * @file test_result_refiner_base.cpp
 * @brief Tests for templated ResultRefinerBase
 *
 * Tests the refineRoot1D method with various iteration and multiplicity
 * detection methods.
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

#include "core/polynomial_base.h"
#include "refinement/result_refiner_base.h"

using namespace polynomial_solver;
using PolyDouble = PolynomialBase<double>;
using RefinerDouble = ResultRefinerBase<double>;
using ConfigDouble = RefinementConfigBase<double>;
using ResultDouble = RefinedRootBase<double>;

static bool approx_equal(double a, double b, double tol = 1e-10) {
    return std::fabs(a - b) < tol;
}

//=============================================================================
// Test 1: Simple root with standard Newton
//=============================================================================
void test_simple_root_newton() {
    std::cout << "  test_simple_root_newton... ";

    // p(x) = x^2 - 2 (roots at ±√2)
    std::vector<unsigned int> degrees{2u};
    std::vector<double> power_coeffs{-2.0, 0.0, 1.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    ConfigDouble config;
    config.target_tolerance = 1e-14;
    config.residual_tolerance = 1e-14;
    config.iteration_method = IterationMethodBase::NEWTON;
    config.multiplicity_method = MultiplicityMethodBase::NONE;

    ResultDouble result = RefinerDouble::refineRoot1D(1.5, p, config);

    assert(result.converged);
    assert(approx_equal(result.location, std::sqrt(2.0), 1e-12));
    assert(result.multiplicity == 1);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 2: Double root with modified Newton
//=============================================================================
void test_double_root_modified_newton() {
    std::cout << "  test_double_root_modified_newton... ";

    // p(x) = (x - 0.5)^2 = x^2 - x + 0.25
    std::vector<unsigned int> degrees{2u};
    std::vector<double> power_coeffs{0.25, -1.0, 1.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    ConfigDouble config;
    config.target_tolerance = 1e-12;
    config.residual_tolerance = 1e-14;
    config.iteration_method = IterationMethodBase::MODIFIED_NEWTON;
    config.multiplicity_method = MultiplicityMethodBase::TAYLOR;

    ResultDouble result = RefinerDouble::refineRoot1D(0.6, p, config);

    assert(result.converged);
    assert(approx_equal(result.location, 0.5, 1e-10));
    assert(result.multiplicity == 2);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 3: Triple root with Schröder's method
//=============================================================================
void test_triple_root_schroder() {
    std::cout << "  test_triple_root_schroder... ";

    // p(x) = (x - 0.3)^3 = x^3 - 0.9x^2 + 0.27x - 0.027
    // Expanded: x^3 - 0.9*x^2 + 0.27*x - 0.027
    double r = 0.3;
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{
        -r*r*r,                // constant: -0.027
        3.0*r*r,               // x^1: 0.27
        -3.0*r,                // x^2: -0.9
        1.0                    // x^3: 1.0
    };

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    ConfigDouble config;
    config.target_tolerance = 1e-10;
    config.residual_tolerance = 1e-14;
    config.iteration_method = IterationMethodBase::SCHRODER;
    config.multiplicity_method = MultiplicityMethodBase::TAYLOR;

    ResultDouble result = RefinerDouble::refineRoot1D(0.35, p, config);

    assert(result.converged);
    assert(approx_equal(result.location, 0.3, 1e-8));
    // Note: multiplicity detection may give 3 or higher
    assert(result.multiplicity >= 2);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 4: Ostrowski multiplicity detection
//=============================================================================
void test_ostrowski_multiplicity() {
    std::cout << "  test_ostrowski_multiplicity... ";

    // p(x) = (x - 0.5)^2 = x^2 - x + 0.25
    std::vector<unsigned int> degrees{2u};
    std::vector<double> power_coeffs{0.25, -1.0, 1.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    // Test Ostrowski from point
    unsigned int mult = RefinerDouble::estimateMultiplicityOstrowskiFromPoint(0.6, p);

    // Should detect multiplicity 2
    assert(mult == 2);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 5: Halley's method
//=============================================================================
void test_halley_method() {
    std::cout << "  test_halley_method... ";

    // p(x) = x^3 - 2 (root at 2^(1/3) ≈ 1.2599)
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{-2.0, 0.0, 0.0, 1.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    ConfigDouble config;
    config.target_tolerance = 1e-14;
    config.residual_tolerance = 1e-14;
    config.iteration_method = IterationMethodBase::HALLEY;
    config.multiplicity_method = MultiplicityMethodBase::NONE;

    ResultDouble result = RefinerDouble::refineRoot1D(1.5, p, config);

    assert(result.converged);
    assert(approx_equal(result.location, std::cbrt(2.0), 1e-12));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 6: Condition number estimation
//=============================================================================
void test_condition_number() {
    std::cout << "  test_condition_number... ";

    // Well-conditioned: p(x) = x - 0.5
    std::vector<unsigned int> degrees1{1u};
    std::vector<double> coeffs1{-0.5, 1.0};
    PolyDouble p1 = PolyDouble::fromPower(degrees1, coeffs1);

    double kappa1 = static_cast<double>(
        RefinerDouble::estimateConditionNumber1D(0.5, p1, 1.0));
    // Should be close to 1
    assert(kappa1 < 10.0);

    // Ill-conditioned: p(x) = (x - 0.5)^5
    std::vector<unsigned int> degrees2{5u};
    double r = 0.5;
    std::vector<double> coeffs2{
        -std::pow(r, 5),
        5*std::pow(r, 4),
        -10*std::pow(r, 3),
        10*std::pow(r, 2),
        -5*r,
        1.0
    };
    PolyDouble p2 = PolyDouble::fromPower(degrees2, coeffs2);

    // Near the root, first derivative is very small
    double df_at_root = p2.differentiate(0).evaluate(0.5);
    double kappa2 = static_cast<double>(
        RefinerDouble::estimateConditionNumber1D(0.5, p2, df_at_root));
    // Should be larger due to multiple root
    assert(kappa2 > 1.0);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 7: Error bounds computation
//=============================================================================
void test_error_bounds() {
    std::cout << "  test_error_bounds... ";

    // p(x) = x^2 - 2
    std::vector<unsigned int> degrees{2u};
    std::vector<double> power_coeffs{-2.0, 0.0, 1.0};
    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    double root = std::sqrt(2.0);
    double first_deriv = 2.0 * root;  // f'(x) = 2x

    double lower, upper;
    bool success = RefinerDouble::computeErrorBounds(
        root, p, 1, first_deriv, lower, upper);

    assert(success);
    assert(lower < root);
    assert(upper > root);
    assert(lower < upper);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 8: Wilkinson-like polynomial
//=============================================================================
void test_wilkinson_like() {
    std::cout << "  test_wilkinson_like... ";

    // p(x) = (x - 0.1)(x - 0.2)(x - 0.3) - clustered roots
    // = x^3 - 0.6x^2 + 0.11x - 0.006
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{-0.006, 0.11, -0.6, 1.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);

    ConfigDouble config;
    config.target_tolerance = 1e-12;
    config.residual_tolerance = 1e-14;
    config.iteration_method = IterationMethodBase::NEWTON;

    // Find root near 0.2
    ResultDouble result = RefinerDouble::refineRoot1D(0.19, p, config);

    assert(result.converged);
    assert(approx_equal(result.location, 0.2, 1e-10));
    assert(result.multiplicity == 1);

    std::cout << "PASSED\n";
}

//=============================================================================
// Main
//=============================================================================
int main() {
    std::cout << "=== Testing ResultRefinerBase ===" << std::endl;

    test_simple_root_newton();
    test_double_root_modified_newton();
    test_triple_root_schroder();
    test_ostrowski_multiplicity();
    test_halley_method();
    test_condition_number();
    test_error_bounds();
    test_wilkinson_like();

    std::cout << "\n=== All ResultRefinerBase Tests PASSED ===" << std::endl;
    return 0;
}

