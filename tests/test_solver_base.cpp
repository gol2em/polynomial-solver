/**
 * @file test_solver_base.cpp
 * @brief Tests for templated SolverBase<Scalar>
 *
 * Tests the templated solver implementation to ensure it produces
 * equivalent results to the existing Solver class.
 */

#include "solver/solver_base.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace polynomial_solver;

// Type aliases
using PolyDouble = PolynomialBase<double>;
using SystemDouble = PolynomialSystemBase<double>;
using SolverDouble = SolverBase<double>;
using ConfigDouble = SubdivisionConfigBase<double>;

// Tolerance for floating point comparisons
constexpr double EPSILON = 1e-6;

bool approx_equal(double a, double b, double tol = EPSILON) {
    return std::fabs(a - b) < tol;
}

//=============================================================================
// Test 1: Empty system
//=============================================================================

void test_empty_system() {
    std::cout << "  test_empty_system... ";

    SystemDouble system;
    SolverDouble solver;
    ConfigDouble config;

    auto result = solver.subdivisionSolve(system, config);
    
    assert(result.boxes.empty());
    assert(result.num_resolved == 0);
    assert(!result.degeneracy_detected);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 2: Simple 1D polynomial (x^2 - 0.25 = 0, roots at x=0.5, x=-0.5)
//=============================================================================

void test_simple_1d_quadratic() {
    std::cout << "  test_simple_1d_quadratic... ";

    // p(x) = x^2 - 0.25, root at x = 0.5 in [0,1]
    std::vector<unsigned int> degrees{2u};
    std::vector<double> power_coeffs{-0.25, 0.0, 1.0};

    PolyDouble p = PolyDouble::fromPower(degrees, power_coeffs);
    SystemDouble system({p});

    SolverDouble solver;
    ConfigDouble config;
    config.tolerance = 1e-8;

    auto result = solver.subdivisionSolve(system, config, 
                                          RootBoundingMethodBase::ProjectedPolyhedral);

    assert(result.num_resolved >= 1);
    
    // Check that we found a root near 0.5
    bool found_root = false;
    for (std::size_t i = 0; i < result.num_resolved; ++i) {
        double center = result.boxes[i].center[0];
        if (approx_equal(center, 0.5, 1e-6)) {
            found_root = true;
            break;
        }
    }
    assert(found_root);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 3: PolynomialSystemBase construction and evaluation
//=============================================================================

void test_polynomial_system_base() {
    std::cout << "  test_polynomial_system_base... ";

    // p1(x,y) = x - 0.3
    // p2(x,y) = y - 0.7
    // Root at (0.3, 0.7)
    //
    // Coefficient ordering for 2D (degrees [1,1]):
    // [c_00, c_01, c_10, c_11] where c_ij = coeff of x^i * y^j
    // So for p(x,y) = a + b*x + c*y + d*xy:
    //   coeffs = [a, c, b, d]
    //   (i.e., y-index varies fastest)

    std::vector<unsigned int> degrees{1u, 1u};

    // p1(x,y) = x - 0.3 = -0.3 + 1*x + 0*y + 0*xy
    // coeffs: [c_00, c_01, c_10, c_11] = [-0.3, 0, 1, 0]
    std::vector<double> power1 = {-0.3, 0.0, 1.0, 0.0};
    PolyDouble p1 = PolyDouble::fromPower(degrees, power1);

    // p2(x,y) = y - 0.7 = -0.7 + 0*x + 1*y + 0*xy
    // coeffs: [c_00, c_01, c_10, c_11] = [-0.7, 1, 0, 0]
    std::vector<double> power2 = {-0.7, 1.0, 0.0, 0.0};
    PolyDouble p2 = PolyDouble::fromPower(degrees, power2);

    SystemDouble system({p1, p2});

    assert(system.dimension() == 2);
    assert(system.equationCount() == 2);

    // Evaluate at root
    std::vector<double> root_point = {0.3, 0.7};
    std::vector<double> values;
    system.evaluate(root_point, values);

    assert(approx_equal(values[0], 0.0, 1e-10));
    assert(approx_equal(values[1], 0.0, 1e-10));
    assert(system.isApproximateRoot(root_point, 1e-6));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 4: 2D system solve
//=============================================================================

void test_2d_system_solve() {
    std::cout << "  test_2d_system_solve... ";

    // Same system as above - root at (0.3, 0.7)
    // p1(x,y) = x - 0.3, p2(x,y) = y - 0.7
    std::vector<unsigned int> degrees{1u, 1u};
    std::vector<double> power1 = {-0.3, 0.0, 1.0, 0.0};  // x - 0.3
    std::vector<double> power2 = {-0.7, 1.0, 0.0, 0.0};  // y - 0.7

    PolyDouble p1 = PolyDouble::fromPower(degrees, power1);
    PolyDouble p2 = PolyDouble::fromPower(degrees, power2);
    SystemDouble system({p1, p2});

    SolverDouble solver;
    ConfigDouble config;
    config.tolerance = 1e-8;

    auto result = solver.subdivisionSolve(system, config,
                                          RootBoundingMethodBase::ProjectedPolyhedral);

    assert(result.num_resolved >= 1);

    // Verify root is near (0.3, 0.7)
    bool found_root = false;
    for (std::size_t i = 0; i < result.num_resolved; ++i) {
        double cx = result.boxes[i].center[0];
        double cy = result.boxes[i].center[1];
        if (approx_equal(cx, 0.3, 1e-5) && approx_equal(cy, 0.7, 1e-5)) {
            found_root = true;
            break;
        }
    }
    assert(found_root);

    std::cout << "PASSED\n";
}

//=============================================================================
// Main
//=============================================================================

int main() {
    std::cout << "Running SolverBase<double> tests...\n";

    test_empty_system();
    test_simple_1d_quadratic();
    test_polynomial_system_base();
    test_2d_system_solve();

    std::cout << "\nAll SolverBase tests PASSED!\n";
    return 0;
}

