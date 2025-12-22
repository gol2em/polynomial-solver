/**
 * @file test_interval_arithmetic.cpp
 * @brief Test interval arithmetic for polynomial evaluation
 */

#ifdef ENABLE_HIGH_PRECISION

#include "interval_arithmetic.h"
#include "polynomial_hp.h"
#include "precision_context.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace polynomial_solver;

// Helper for comparing mpreal values
bool approxEqual(const mpreal& a, const mpreal& b, const mpreal& tol = mpreal("1e-50")) {
    return abs(a - b) <= tol;
}

// Test basic interval operations
void test_interval_addition() {
    std::cout << "\n=== Test: Interval Addition ===" << std::endl;
    PrecisionContext ctx(256);

    Interval a(mpreal("1.0"), mpreal("2.0"));
    Interval b(mpreal("3.0"), mpreal("5.0"));
    Interval c = a + b;

    assert(c.lower == mpreal("4.0"));
    assert(c.upper == mpreal("7.0"));
    std::cout << "[1,2] + [3,5] = [" << c.lower << ", " << c.upper << "]" << std::endl;
    std::cout << "PASSED" << std::endl;
}

void test_interval_subtraction() {
    std::cout << "\n=== Test: Interval Subtraction ===" << std::endl;
    PrecisionContext ctx(256);

    Interval a(mpreal("1.0"), mpreal("2.0"));
    Interval b(mpreal("3.0"), mpreal("5.0"));
    Interval c = a - b;
    // [1,2] - [3,5] = [1-5, 2-3] = [-4, -1]
    assert(c.lower == mpreal("-4.0"));
    assert(c.upper == mpreal("-1.0"));
    std::cout << "[1,2] - [3,5] = [" << c.lower << ", " << c.upper << "]" << std::endl;
    std::cout << "PASSED" << std::endl;
}

void test_interval_multiplication() {
    std::cout << "\n=== Test: Interval Multiplication ===" << std::endl;
    PrecisionContext ctx(256);

    // Positive * Positive
    Interval a(mpreal("2.0"), mpreal("3.0"));
    Interval b(mpreal("4.0"), mpreal("5.0"));
    Interval c = a * b;
    assert(c.lower == mpreal("8.0"));
    assert(c.upper == mpreal("15.0"));
    std::cout << "[2,3] * [4,5] = [" << c.lower << ", " << c.upper << "]" << std::endl;

    // Mixed signs
    Interval d(mpreal("-1.0"), mpreal("2.0"));
    Interval e(mpreal("3.0"), mpreal("4.0"));
    Interval f = d * e;
    // min(-3, -4, 6, 8) = -4, max = 8
    assert(f.lower == mpreal("-4.0"));
    assert(f.upper == mpreal("8.0"));
    std::cout << "[-1,2] * [3,4] = [" << f.lower << ", " << f.upper << "]" << std::endl;
    std::cout << "PASSED" << std::endl;
}

void test_interval_contains_zero() {
    std::cout << "\n=== Test: Interval Contains Zero ===" << std::endl;
    PrecisionContext ctx(256);

    Interval a(mpreal("-1.0"), mpreal("1.0"));
    assert(a.containsZero());
    std::cout << "[-1,1] contains zero: true" << std::endl;

    Interval b(mpreal("0.5"), mpreal("1.5"));
    assert(!b.containsZero());
    std::cout << "[0.5,1.5] contains zero: false" << std::endl;

    Interval c(mpreal("-2.0"), mpreal("-0.5"));
    assert(!c.containsZero());
    std::cout << "[-2,-0.5] contains zero: false" << std::endl;
    std::cout << "PASSED" << std::endl;
}

// Test interval power
void test_interval_pow_odd() {
    std::cout << "\n=== Test: Interval Power (Odd) ===" << std::endl;
    PrecisionContext ctx(256);

    // Odd power: monotonic
    Interval x(mpreal("2.0"), mpreal("3.0"));
    Interval x3 = intervalPow(x, 3);
    assert(x3.lower == mpreal("8.0"));   // 2^3 = 8
    assert(x3.upper == mpreal("27.0"));  // 3^3 = 27
    std::cout << "[2,3]^3 = [" << x3.lower << ", " << x3.upper << "]" << std::endl;

    // Negative interval, odd power
    Interval y(mpreal("-3.0"), mpreal("-2.0"));
    Interval y3 = intervalPow(y, 3);
    assert(y3.lower == mpreal("-27.0")); // (-3)^3 = -27
    assert(y3.upper == mpreal("-8.0"));  // (-2)^3 = -8
    std::cout << "[-3,-2]^3 = [" << y3.lower << ", " << y3.upper << "]" << std::endl;
    std::cout << "PASSED" << std::endl;
}

void test_interval_pow_even() {
    std::cout << "\n=== Test: Interval Power (Even) ===" << std::endl;
    PrecisionContext ctx(256);

    // Even power, positive interval: monotonic increasing
    Interval x(mpreal("2.0"), mpreal("3.0"));
    Interval x2 = intervalPow(x, 2);
    assert(x2.lower == mpreal("4.0"));   // 2^2 = 4
    assert(x2.upper == mpreal("9.0"));   // 3^2 = 9
    std::cout << "[2,3]^2 = [" << x2.lower << ", " << x2.upper << "]" << std::endl;

    // Even power, negative interval: monotonic decreasing
    Interval y(mpreal("-3.0"), mpreal("-2.0"));
    Interval y2 = intervalPow(y, 2);
    assert(y2.lower == mpreal("4.0"));   // (-2)^2 = 4
    assert(y2.upper == mpreal("9.0"));   // (-3)^2 = 9
    std::cout << "[-3,-2]^2 = [" << y2.lower << ", " << y2.upper << "]" << std::endl;

    // Even power, interval containing zero: minimum at 0
    Interval z(mpreal("-2.0"), mpreal("3.0"));
    Interval z2 = intervalPow(z, 2);
    assert(z2.lower == mpreal("0.0"));   // minimum at x=0
    assert(z2.upper == mpreal("9.0"));   // max(4, 9) = 9
    std::cout << "[-2,3]^2 = [" << z2.lower << ", " << z2.upper << "]" << std::endl;
    std::cout << "PASSED" << std::endl;
}

// Test polynomial evaluation on interval
void test_linear_polynomial() {
    std::cout << "\n=== Test: Linear Polynomial Evaluation ===" << std::endl;
    PrecisionContext ctx(256);

    // p(x) = 2x + 1
    std::vector<unsigned int> degrees = {1};
    std::vector<mpreal> coeffs = {mpreal("1.0"), mpreal("2.0")};
    PolynomialHP poly = fromPowerHP(degrees, coeffs);

    Interval x(mpreal("0.0"), mpreal("1.0"));
    Interval result = evaluateOnInterval(poly, x);

    // p([0,1]) = [2*0+1, 2*1+1] = [1, 3]
    assert(result.lower == mpreal("1.0"));
    assert(result.upper == mpreal("3.0"));
    std::cout << "p(x) = 2x + 1 on [0,1] = [" << result.lower << ", " << result.upper << "]" << std::endl;
    std::cout << "PASSED" << std::endl;
}

void test_quadratic_polynomial() {
    std::cout << "\n=== Test: Quadratic Polynomial Evaluation ===" << std::endl;
    PrecisionContext ctx(256);

    // p(x) = x^2 - 1
    std::vector<unsigned int> degrees = {2};
    std::vector<mpreal> coeffs = {mpreal("-1.0"), mpreal("0.0"), mpreal("1.0")};
    PolynomialHP poly = fromPowerHP(degrees, coeffs);

    // Interval [0.5, 1.5]: x^2 in [0.25, 2.25], so p in [-0.75, 1.25]
    Interval x1(mpreal("0.5"), mpreal("1.5"));
    Interval r1 = evaluateOnInterval(poly, x1);
    assert(approxEqual(r1.lower, mpreal("-0.75")));
    assert(approxEqual(r1.upper, mpreal("1.25")));
    std::cout << "p(x) = x^2 - 1 on [0.5,1.5] = [" << r1.lower << ", " << r1.upper << "]" << std::endl;

    // Interval [-1, 1]: x^2 in [0, 1], so p in [-1, 0]
    Interval x2(mpreal("-1.0"), mpreal("1.0"));
    Interval r2 = evaluateOnInterval(poly, x2);
    assert(approxEqual(r2.lower, mpreal("-1.0")));
    assert(approxEqual(r2.upper, mpreal("0.0")));
    std::cout << "p(x) = x^2 - 1 on [-1,1] = [" << r2.lower << ", " << r2.upper << "]" << std::endl;
    std::cout << "PASSED" << std::endl;
}

void test_min_abs_on_interval() {
    std::cout << "\n=== Test: Minimum Absolute Value on Interval ===" << std::endl;
    PrecisionContext ctx(256);

    // p(x) = x - 0.5 crosses zero at x=0.5
    std::vector<unsigned int> degrees = {1};
    std::vector<mpreal> coeffs = {mpreal("-0.5"), mpreal("1.0")};
    PolynomialHP poly = fromPowerHP(degrees, coeffs);

    Interval x(mpreal("0.0"), mpreal("1.0"));
    mpreal minAbs = minAbsOnInterval(poly, x);

    // p([0,1]) = [-0.5, 0.5] contains zero
    assert(minAbs == mpreal("0.0"));
    std::cout << "p(x) = x - 0.5 on [0,1]: min|p| = " << minAbs << std::endl;

    // p(x) = x - 0.5 on [0.6, 0.8] doesn't cross zero
    Interval x2(mpreal("0.6"), mpreal("0.8"));
    mpreal minAbs2 = minAbsOnInterval(poly, x2);
    // p([0.6, 0.8]) = [0.1, 0.3], so min|p| = 0.1
    assert(approxEqual(minAbs2, mpreal("0.1")));
    std::cout << "p(x) = x - 0.5 on [0.6,0.8]: min|p| = " << minAbs2 << std::endl;
    std::cout << "PASSED" << std::endl;
}

// Test interval Newton on simple cubic
void test_interval_newton_simple_cubic() {
    std::cout << "\n=== Test: Interval Newton on Simple Cubic ===" << std::endl;
    PrecisionContext ctx(256);

    // p(x) = (x - 1/5)(x - 1/2)(x - 4/5)
    // Using exact rational arithmetic for coefficients:
    // = x^3 - (1/5 + 1/2 + 4/5)*x^2 + (1/5*1/2 + 1/5*4/5 + 1/2*4/5)*x - 1/5*1/2*4/5
    // = x^3 - 3/2*x^2 + 33/50*x - 2/25
    //
    // c0 = -2/25 = -0.08 (exact)
    // c1 = 33/50 = 0.66 (exact)
    // c2 = -3/2 = -1.5 (exact)
    // c3 = 1 (exact)
    //
    // The roots 1/5, 1/2, 4/5 have exact representations:
    // 1/5 = 0.2 (but 0.2 has no exact binary representation!)
    // 1/2 = 0.5 (exact in binary)
    // 4/5 = 0.8 (no exact binary representation)

    // Use exact fractions for coefficients
    mpreal c0 = mpreal(-2) / mpreal(25);   // -0.08
    mpreal c1 = mpreal(33) / mpreal(50);   // 0.66
    mpreal c2 = mpreal(-3) / mpreal(2);    // -1.5
    mpreal c3 = mpreal(1);

    std::vector<mpreal> coeffs = {c0, c1, c2, c3};
    PolynomialHP poly = fromPowerHP({3}, coeffs);

    // Derivative: p'(x) = 3*x^2 - 3*x + 0.66
    mpreal d0 = mpreal(33) / mpreal(50);   // 0.66
    mpreal d1 = mpreal(-3);                 // -3
    mpreal d2 = mpreal(3);                  // 3

    std::vector<mpreal> dcoeffs = {d0, d1, d2};
    PolynomialHP dpoly = fromPowerHP({2}, dcoeffs);

    // Expected roots as exact fractions
    mpreal root1 = mpreal(1) / mpreal(5);   // 0.2
    mpreal root2 = mpreal(1) / mpreal(2);   // 0.5
    mpreal root3 = mpreal(4) / mpreal(5);   // 0.8

    std::cout << "Polynomial: p(x) = (x - 1/5)(x - 1/2)(x - 4/5)" << std::endl;
    std::cout << "Expected roots: 1/5 = " << root1 << ", 1/2 = " << root2 << ", 4/5 = " << root3 << std::endl;
    std::cout << std::endl;

    // Note: Interval Newton requires tight initial intervals to avoid the
    // "dependency problem" where interval arithmetic overestimates ranges.
    // The solver typically provides boxes of size ~1e-8, which work well.
    // Wide intervals like [0.1, 0.3] cause interval evaluation of f'(x) to
    // include zero even when f'(x) is bounded away from zero on that interval.

    // Test root at 1/5 with tight initial interval (simulating solver output)
    {
        // Solver would give us a box like [0.2 - 1e-8, 0.2 + 1e-8]
        Interval init(mpreal("0.19999999"), mpreal("0.20000001"));
        IntervalNewtonResult result = intervalNewton(poly, init, dpoly, 100, mpreal("1e-100"));

        std::cout << "Root near 1/5 = 0.2:" << std::endl;
        std::cout << "  Initial interval: [0.19999999, 0.20000001]" << std::endl;
        std::cout << std::setprecision(80);
        std::cout << "  Final enclosure: [" << result.enclosure.lower << ", " << result.enclosure.upper << "]" << std::endl;
        std::cout << "  Root estimate: " << result.root << std::endl;
        std::cout << "  Exact root 1/5: " << root1 << std::endl;
        std::cout << std::setprecision(6);
        std::cout << "  Enclosure width: " << result.enclosure.width() << std::endl;
        std::cout << "  Error bound: " << result.error_bound << std::endl;
        std::cout << "  Distance from lower: " << (root1 - result.enclosure.lower) << std::endl;
        std::cout << "  Distance from upper: " << (result.enclosure.upper - root1) << std::endl;
        std::cout << "  Unique root verified: " << (result.unique_root ? "yes" : "no") << std::endl;
        std::cout << "  Converged: " << (result.converged ? "yes" : "no") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Status: " << result.status << std::endl;

        mpreal error = abs(result.root - root1);
        std::cout << "  Actual error from exact 1/5: " << error << std::endl;

        // Verify we found a unique root and converged
        assert(result.unique_root);
        assert(result.converged);
        // The error should be very small (< 1e-70 for 256-bit precision)
        assert(error < mpreal("1e-70"));
        std::cout << "  PASSED" << std::endl;
    }

    // Test root at 1/2 with tight initial interval
    {
        Interval init(mpreal("0.49999999"), mpreal("0.50000001"));
        IntervalNewtonResult result = intervalNewton(poly, init, dpoly, 100, mpreal("1e-100"));

        std::cout << "\nRoot near 1/2 = 0.5:" << std::endl;
        std::cout << "  Initial interval: [0.49999999, 0.50000001]" << std::endl;
        std::cout << std::setprecision(80);
        std::cout << "  Final enclosure: [" << result.enclosure.lower << ", " << result.enclosure.upper << "]" << std::endl;
        std::cout << "  Root estimate: " << result.root << std::endl;
        std::cout << std::setprecision(6);
        std::cout << "  Enclosure width: " << result.enclosure.width() << std::endl;
        std::cout << "  Error bound: " << result.error_bound << std::endl;
        std::cout << "  Unique root verified: " << (result.unique_root ? "yes" : "no") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;

        mpreal error = abs(result.root - root2);
        std::cout << "  Actual error from exact 1/2: " << error << std::endl;

        // Verify we found a unique root and converged with small error
        assert(result.unique_root);
        assert(result.converged);
        assert(error < mpreal("1e-70"));
        std::cout << "  PASSED" << std::endl;
    }

    // Test root at 4/5 with tight initial interval
    {
        Interval init(mpreal("0.79999999"), mpreal("0.80000001"));
        IntervalNewtonResult result = intervalNewton(poly, init, dpoly, 100, mpreal("1e-100"));

        std::cout << "\nRoot near 4/5 = 0.8:" << std::endl;
        std::cout << "  Initial interval: [0.79999999, 0.80000001]" << std::endl;
        std::cout << std::setprecision(80);
        std::cout << "  Final enclosure: [" << result.enclosure.lower << ", " << result.enclosure.upper << "]" << std::endl;
        std::cout << "  Root estimate: " << result.root << std::endl;
        std::cout << std::setprecision(6);
        std::cout << "  Enclosure width: " << result.enclosure.width() << std::endl;
        std::cout << "  Error bound: " << result.error_bound << std::endl;
        std::cout << "  Unique root verified: " << (result.unique_root ? "yes" : "no") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;

        mpreal error = abs(result.root - root3);
        std::cout << "  Actual error from exact 4/5: " << error << std::endl;

        // Verify we found a unique root and converged with small error
        assert(result.unique_root);
        assert(result.converged);
        assert(error < mpreal("1e-70"));
        std::cout << "  PASSED" << std::endl;
    }

    std::cout << "\nAll interval Newton tests PASSED" << std::endl;
}

#endif // ENABLE_HIGH_PRECISION

int main() {
#ifdef ENABLE_HIGH_PRECISION
    std::cout << "========================================" << std::endl;
    std::cout << "Interval Arithmetic Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    test_interval_addition();
    test_interval_subtraction();
    test_interval_multiplication();
    test_interval_contains_zero();
    test_interval_pow_odd();
    test_interval_pow_even();
    test_linear_polynomial();
    test_quadratic_polynomial();
    test_min_abs_on_interval();
    test_interval_newton_simple_cubic();

    std::cout << "\n========================================" << std::endl;
    std::cout << "All interval arithmetic tests PASSED!" << std::endl;
    std::cout << "========================================" << std::endl;
#else
    std::cout << "High precision not enabled. Test skipped." << std::endl;
#endif
    return 0;
}

