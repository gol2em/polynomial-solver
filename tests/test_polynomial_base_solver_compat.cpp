/**
 * @file test_polynomial_base_solver_compat.cpp
 * @brief Shadow tests to verify PolynomialBase<double> is compatible with Solver operations
 * 
 * This test simulates the exact operations that the Solver performs on Polynomial objects:
 * - dimension(), degrees(), coefficientCount()
 * - evaluate(point)
 * - graphControlPoints(control_points)
 * - ensureBernsteinPrimary()
 * - restrictedToInterval(axis, a, b)
 * 
 * Both Polynomial and PolynomialBase<double> are tested side-by-side and their results
 * are compared to ensure they produce identical outputs.
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <sstream>

#include "core/polynomial.h"
#include "core/polynomial_base.h"

using namespace polynomial_solver;

// Test utilities
static double max_diff = 0.0;  // Track maximum difference across all tests
static int test_count = 0;
static int pass_count = 0;

bool double_equal(double a, double b, double tol = 1e-14) {
    double diff = std::fabs(a - b);
    if (diff > max_diff) max_diff = diff;
    return diff < tol;
}

bool vectors_equal(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-14) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!double_equal(a[i], b[i], tol)) return false;
    }
    return true;
}

bool vectors_equal_ui(const std::vector<unsigned int>& a, const std::vector<unsigned int>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

#define TEST(name) \
    void name(); \
    struct name##_register { \
        name##_register() { \
            test_count++; \
            std::cout << "Running " << #name << "... "; \
            try { \
                name(); \
                pass_count++; \
                std::cout << "PASSED\n"; \
            } catch (const std::exception& e) { \
                std::cout << "FAILED: " << e.what() << "\n"; \
            } \
        } \
    } name##_instance; \
    void name()

#define ASSERT_TRUE(cond) \
    if (!(cond)) { \
        std::ostringstream oss; \
        oss << "Assertion failed: " << #cond << " at line " << __LINE__; \
        throw std::runtime_error(oss.str()); \
    }

#define ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::ostringstream oss; \
        oss << "Assertion failed: " << #a << " == " << #b << " (" << (a) << " != " << (b) << ") at line " << __LINE__; \
        throw std::runtime_error(oss.str()); \
    }

// ============================================================================
// Test Framework: Create parallel polynomials and compare operations
// ============================================================================

/**
 * @brief Create parallel Polynomial and PolynomialBase<double> from Bernstein coefficients
 */
std::pair<Polynomial, PolynomialBase<double>> create_parallel_bernstein(
    const std::vector<unsigned int>& degrees,
    const std::vector<double>& coeffs)
{
    Polynomial poly(degrees, coeffs);
    PolynomialBase<double> polyBase = PolynomialBase<double>::fromBernstein(degrees, coeffs);
    return {poly, polyBase};
}

/**
 * @brief Create parallel polynomials from power coefficients
 */
std::pair<Polynomial, PolynomialBase<double>> create_parallel_power(
    const std::vector<unsigned int>& degrees,
    const std::vector<double>& coeffs)
{
    Polynomial poly = Polynomial::fromPower(degrees, coeffs);
    PolynomialBase<double> polyBase = PolynomialBase<double>::fromPower(degrees, coeffs);
    return {poly, polyBase};
}

// ============================================================================
// PART 1: Basic property equivalence
// ============================================================================

TEST(test_dimension_equivalence_1d)
{
    std::vector<unsigned int> degrees = {3};
    std::vector<double> coeffs = {1.0, -2.0, 3.0, -1.0};
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    ASSERT_EQ(poly.dimension(), polyBase.dimension());
    ASSERT_EQ(poly.dimension(), 1u);
}

TEST(test_dimension_equivalence_2d)
{
    std::vector<unsigned int> degrees = {2, 3};
    std::vector<double> coeffs(12);
    for (int i = 0; i < 12; ++i) coeffs[i] = static_cast<double>(i + 1);
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    ASSERT_EQ(poly.dimension(), polyBase.dimension());
    ASSERT_EQ(poly.dimension(), 2u);
}

TEST(test_degrees_equivalence)
{
    std::vector<unsigned int> degrees = {4, 2, 3};
    std::vector<double> coeffs(60);
    for (int i = 0; i < 60; ++i) coeffs[i] = static_cast<double>(i);
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    ASSERT_TRUE(vectors_equal_ui(poly.degrees(), polyBase.degrees()));
}

TEST(test_coefficient_count_equivalence)
{
    std::vector<unsigned int> degrees = {3, 2};
    std::vector<double> coeffs(12);
    for (int i = 0; i < 12; ++i) coeffs[i] = static_cast<double>(i);
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    ASSERT_EQ(poly.coefficientCount(), polyBase.coefficientCount());
    ASSERT_EQ(poly.coefficientCount(), 12u);
}

// ============================================================================
// PART 2: Evaluation equivalence (critical for solver)
// ============================================================================

TEST(test_evaluate_1d_bernstein)
{
    // f(t) = t^2 in Bernstein form on [0,1]: [0, 0, 1]
    std::vector<unsigned int> degrees = {2};
    std::vector<double> coeffs = {0.0, 0.0, 1.0};
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    // Test at multiple points
    std::vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double t : test_points) {
        std::vector<double> params = {t};
        double val_poly = poly.evaluate(params);
        double val_base = polyBase.evaluate(params);
        ASSERT_TRUE(double_equal(val_poly, val_base));
    }
}

TEST(test_evaluate_2d_bernstein)
{
    // Bilinear patch: f(s,t) = s*t on [0,1]^2
    // Bernstein coeffs: [0,0,0,1] for degrees [1,1]
    std::vector<unsigned int> degrees = {1, 1};
    std::vector<double> coeffs = {0.0, 0.0, 0.0, 1.0};  // (0,0)=0, (1,0)=0, (0,1)=0, (1,1)=1
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    // Test at grid of points
    for (double s = 0.0; s <= 1.0; s += 0.25) {
        for (double t = 0.0; t <= 1.0; t += 0.25) {
            std::vector<double> params = {s, t};
            double val_poly = poly.evaluate(params);
            double val_base = polyBase.evaluate(params);
            ASSERT_TRUE(double_equal(val_poly, val_base));
        }
    }
}

TEST(test_evaluate_power_basis)
{
    // f(x) = x^3 - 2x + 1 in power basis
    std::vector<unsigned int> degrees = {3};
    std::vector<double> power_coeffs = {1.0, -2.0, 0.0, 1.0};  // 1 - 2x + x^3
    
    auto [poly, polyBase] = create_parallel_power(degrees, power_coeffs);
    
    std::vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double t : test_points) {
        std::vector<double> params = {t};
        double val_poly = poly.evaluate(params);
        double val_base = polyBase.evaluate(params);
        ASSERT_TRUE(double_equal(val_poly, val_base));
    }
}

// ============================================================================
// PART 3: graphControlPoints equivalence (critical for solver bounding)
// ============================================================================

TEST(test_graph_control_points_1d)
{
    std::vector<unsigned int> degrees = {3};
    std::vector<double> coeffs = {1.0, 0.5, -0.5, 2.0};
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    std::vector<double> cp_poly;
    std::vector<double> cp_base;
    
    poly.graphControlPoints(cp_poly);
    polyBase.graphControlPoints(cp_base);
    
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
}

TEST(test_graph_control_points_2d)
{
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> coeffs(9);
    for (int i = 0; i < 9; ++i) coeffs[i] = std::sin(static_cast<double>(i));
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    std::vector<double> cp_poly;
    std::vector<double> cp_base;
    
    poly.graphControlPoints(cp_poly);
    polyBase.graphControlPoints(cp_base);
    
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
}

TEST(test_graph_control_points_from_power)
{
    // Test that power -> Bernstein conversion gives same control points
    std::vector<unsigned int> degrees = {3};
    std::vector<double> power_coeffs = {1.0, -1.0, 2.0, -0.5};
    
    auto [poly, polyBase] = create_parallel_power(degrees, power_coeffs);
    
    // Ensure both have Bernstein representation
    poly.ensureBernsteinPrimary();
    polyBase.ensureBernsteinPrimary();
    
    std::vector<double> cp_poly;
    std::vector<double> cp_base;
    
    poly.graphControlPoints(cp_poly);
    polyBase.graphControlPoints(cp_base);
    
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
}

// ============================================================================
// PART 4: restrictedToInterval equivalence (critical for subdivision)
// ============================================================================

TEST(test_restricted_to_interval_1d)
{
    std::vector<unsigned int> degrees = {4};
    std::vector<double> coeffs = {1.0, 2.0, -1.0, 3.0, -2.0};
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    // Restrict to [0.25, 0.75]
    Polynomial restricted_poly = poly.restrictedToInterval(0, 0.25, 0.75);
    PolynomialBase<double> restricted_base = polyBase.restrictedToInterval(0, 0.25, 0.75);
    
    // Compare control points
    std::vector<double> cp_poly;
    std::vector<double> cp_base;
    restricted_poly.graphControlPoints(cp_poly);
    restricted_base.graphControlPoints(cp_base);
    
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Compare evaluation at multiple points
    for (double t = 0.0; t <= 1.0; t += 0.1) {
        std::vector<double> params = {t};
        double val_poly = restricted_poly.evaluate(params);
        double val_base = restricted_base.evaluate(params);
        ASSERT_TRUE(double_equal(val_poly, val_base));
    }
}

TEST(test_restricted_to_interval_2d_axis0)
{
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> coeffs(9);
    for (int i = 0; i < 9; ++i) coeffs[i] = static_cast<double>(i + 1);
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    // Restrict axis 0 to [0.0, 0.5]
    Polynomial restricted_poly = poly.restrictedToInterval(0, 0.0, 0.5);
    PolynomialBase<double> restricted_base = polyBase.restrictedToInterval(0, 0.0, 0.5);
    
    // Compare at grid points
    for (double s = 0.0; s <= 1.0; s += 0.25) {
        for (double t = 0.0; t <= 1.0; t += 0.25) {
            std::vector<double> params = {s, t};
            double val_poly = restricted_poly.evaluate(params);
            double val_base = restricted_base.evaluate(params);
            ASSERT_TRUE(double_equal(val_poly, val_base));
        }
    }
}

TEST(test_restricted_to_interval_2d_axis1)
{
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> coeffs(9);
    for (int i = 0; i < 9; ++i) coeffs[i] = static_cast<double>(i + 1);
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    // Restrict axis 1 to [0.5, 1.0]
    Polynomial restricted_poly = poly.restrictedToInterval(1, 0.5, 1.0);
    PolynomialBase<double> restricted_base = polyBase.restrictedToInterval(1, 0.5, 1.0);
    
    // Compare at grid points
    for (double s = 0.0; s <= 1.0; s += 0.25) {
        for (double t = 0.0; t <= 1.0; t += 0.25) {
            std::vector<double> params = {s, t};
            double val_poly = restricted_poly.evaluate(params);
            double val_base = restricted_base.evaluate(params);
            ASSERT_TRUE(double_equal(val_poly, val_base));
        }
    }
}

// ============================================================================
// PART 5: Simulated Solver Workflow
// Tests the exact sequence of operations the Solver performs
// ============================================================================

TEST(test_solver_workflow_1d)
{
    // Simulate solver workflow for 1D polynomial
    // f(x) = (x - 0.3)(x - 0.7) = x^2 - x + 0.21 in power basis
    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {0.21, -1.0, 1.0};  // c0 + c1*x + c2*x^2
    
    auto [poly, polyBase] = create_parallel_power(degrees, power_coeffs);
    
    // Step 1: Solver converts to Bernstein (once at the beginning)
    poly.ensureBernsteinPrimary();
    polyBase.ensureBernsteinPrimary();
    
    // Step 2: Get initial graph control points for bounding
    std::vector<double> cp_poly, cp_base;
    poly.graphControlPoints(cp_poly);
    polyBase.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Step 3: Subdivision - restrict to [0, 0.5] (left child)
    Polynomial left_poly = poly.restrictedToInterval(0, 0.0, 0.5);
    PolynomialBase<double> left_base = polyBase.restrictedToInterval(0, 0.0, 0.5);
    
    // Step 4: Get control points for left child
    left_poly.graphControlPoints(cp_poly);
    left_base.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Step 5: Subdivision - restrict to [0.5, 1.0] (right child)
    Polynomial right_poly = poly.restrictedToInterval(0, 0.5, 1.0);
    PolynomialBase<double> right_base = polyBase.restrictedToInterval(0, 0.5, 1.0);
    
    // Step 6: Get control points for right child
    right_poly.graphControlPoints(cp_poly);
    right_base.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Step 7: Further subdivision on right child [0.5, 0.75]
    Polynomial right_left_poly = right_poly.restrictedToInterval(0, 0.0, 0.5);
    PolynomialBase<double> right_left_base = right_base.restrictedToInterval(0, 0.0, 0.5);
    
    right_left_poly.graphControlPoints(cp_poly);
    right_left_base.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Step 8: Final evaluation check
    for (double t = 0.0; t <= 1.0; t += 0.2) {
        std::vector<double> params = {t};
        ASSERT_TRUE(double_equal(
            right_left_poly.evaluate(params),
            right_left_base.evaluate(params)
        ));
    }
}

TEST(test_solver_workflow_2d_circle_ellipse)
{
    // Simulate 2D solver workflow for system:
    // f1(x,y) = x^2 + y^2 - 1 (circle)
    // f2(x,y) = x^2/4 + y^2 - 1 (ellipse)
    
    // We'll test just f1 here (circle equation)
    // In power basis: -1 + x^2 + y^2
    // degrees [2,2], coeffs indexed as [i + j*3] for i,j in 0..2
    // c[0,0] = -1, c[2,0] = 1, c[0,2] = 1, others = 0
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> power_coeffs(9, 0.0);
    power_coeffs[0] = -1.0;  // constant
    power_coeffs[2] = 1.0;   // x^2 (index = 2 + 0*3 = 2)
    power_coeffs[6] = 1.0;   // y^2 (index = 0 + 2*3 = 6)
    
    auto [poly, polyBase] = create_parallel_power(degrees, power_coeffs);
    
    // Convert to Bernstein
    poly.ensureBernsteinPrimary();
    polyBase.ensureBernsteinPrimary();
    
    // Check control points
    std::vector<double> cp_poly, cp_base;
    poly.graphControlPoints(cp_poly);
    polyBase.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Subdivision: lower-left quadrant [0, 0.5] x [0, 0.5]
    Polynomial sub_poly = poly.restrictedToInterval(0, 0.0, 0.5);
    PolynomialBase<double> sub_base = polyBase.restrictedToInterval(0, 0.0, 0.5);
    
    sub_poly = sub_poly.restrictedToInterval(1, 0.0, 0.5);
    sub_base = sub_base.restrictedToInterval(1, 0.0, 0.5);
    
    sub_poly.graphControlPoints(cp_poly);
    sub_base.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Evaluate at center of quadrant
    std::vector<double> center = {0.5, 0.5};  // Center of [0,0.5]x[0,0.5] in local coords
    ASSERT_TRUE(double_equal(sub_poly.evaluate(center), sub_base.evaluate(center)));
}

TEST(test_deep_subdivision_chain)
{
    // Test many subdivision levels (simulates deep solver recursion)
    std::vector<unsigned int> degrees = {5};
    std::vector<double> coeffs = {1.0, -2.0, 3.0, -2.0, 1.0, 0.5};
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    Polynomial current_poly = poly;
    PolynomialBase<double> current_base = polyBase;
    
    // 10 levels of subdivision (alternating left/right)
    for (int level = 0; level < 10; ++level) {
        double a = (level % 2 == 0) ? 0.0 : 0.5;
        double b = (level % 2 == 0) ? 0.5 : 1.0;
        
        current_poly = current_poly.restrictedToInterval(0, a, b);
        current_base = current_base.restrictedToInterval(0, a, b);
        
        // Check equivalence at each level
        std::vector<double> cp_poly, cp_base;
        current_poly.graphControlPoints(cp_poly);
        current_base.graphControlPoints(cp_base);
        ASSERT_TRUE(vectors_equal(cp_poly, cp_base, 1e-12));  // Slightly relaxed tolerance for accumulated error
    }
    
    // Final evaluation check
    for (double t = 0.0; t <= 1.0; t += 0.25) {
        std::vector<double> params = {t};
        ASSERT_TRUE(double_equal(
            current_poly.evaluate(params),
            current_base.evaluate(params),
            1e-12
        ));
    }
}

// ============================================================================
// PART 6: Edge cases and numerical stability
// ============================================================================

TEST(test_degree_zero_polynomial)
{
    // Constant polynomial
    std::vector<unsigned int> degrees = {0};
    std::vector<double> coeffs = {3.14159};
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    ASSERT_EQ(poly.dimension(), polyBase.dimension());
    ASSERT_TRUE(vectors_equal_ui(poly.degrees(), polyBase.degrees()));
    
    std::vector<double> cp_poly, cp_base;
    poly.graphControlPoints(cp_poly);
    polyBase.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    std::vector<double> params = {0.5};
    ASSERT_TRUE(double_equal(poly.evaluate(params), polyBase.evaluate(params)));
}

TEST(test_high_degree_wilkinson)
{
    // Wilkinson-like polynomial: high degree, needs good numerical properties
    std::vector<unsigned int> degrees = {10};
    std::vector<double> coeffs(11);
    for (int i = 0; i <= 10; ++i) {
        coeffs[i] = std::pow(-1.0, i) * static_cast<double>(i + 1);
    }
    
    auto [poly, polyBase] = create_parallel_bernstein(degrees, coeffs);
    
    // Check control points
    std::vector<double> cp_poly, cp_base;
    poly.graphControlPoints(cp_poly);
    polyBase.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base));
    
    // Subdivision
    Polynomial sub_poly = poly.restrictedToInterval(0, 0.3, 0.7);
    PolynomialBase<double> sub_base = polyBase.restrictedToInterval(0, 0.3, 0.7);
    
    sub_poly.graphControlPoints(cp_poly);
    sub_base.graphControlPoints(cp_base);
    ASSERT_TRUE(vectors_equal(cp_poly, cp_base, 1e-12));
}

// NOTE: Arithmetic operations (+, -, *) are NOT used by the solver.
// The solver only uses: evaluate, graphControlPoints, restrictedToInterval,
// ensureBernsteinPrimary, dimension, degrees. So we skip arithmetic tests here.
// If arithmetic is needed for other modules, add a separate test file.

// ============================================================================
// Main
// ============================================================================

int main() {
    std::cout << "\n========================================\n";
    std::cout << "PolynomialBase<double> Solver Compatibility Tests\n";
    std::cout << "========================================\n\n";
    
    // Tests are auto-registered and run by static initializers
    
    std::cout << "\n========================================\n";
    std::cout << "Results: " << pass_count << "/" << test_count << " tests passed\n";
    std::cout << "Maximum numerical difference observed: " << max_diff << "\n";
    std::cout << "========================================\n";
    
    return (pass_count == test_count) ? 0 : 1;
}
