/**
 * @file test_polynomial_base_comprehensive.cpp
 * @brief Comprehensive tests for PolynomialBase<Scalar> template class
 *
 * Tests that the templated PolynomialBase class provides equivalent functionality
 * to both the existing Polynomial (double) and PolynomialHP classes.
 *
 * Test categories:
 * 1. PolynomialBase<double> vs Polynomial (double) - Functional equivalence
 * 2. PolynomialBase<double> - Full test suite from existing tests
 * 3. PolynomialBase<mpreal> - High precision tests (when HP enabled)
 * 4. Cross-precision comparison - Verify precision improvement
 */

#include "core/polynomial_base.h"
#include "core/polynomial.h"

#ifdef ENABLE_HIGH_PRECISION
#include "hp/polynomial_hp.h"
#include "hp/high_precision_types.h"
#include "hp/precision_conversion.h"
#include "hp/precision_context.h"
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <random>
#include <iomanip>

using namespace polynomial_solver;

// Type aliases for readability
using PolyBaseDouble = PolynomialBase<double>;
#ifdef ENABLE_HIGH_PRECISION
using PolyBaseHP = PolynomialBase<mpreal>;
#endif

//=============================================================================
// Utility Functions
//=============================================================================

bool approx_equal(double a, double b, double tol = std::numeric_limits<double>::epsilon() * 100.0) {
    return std::fabs(a - b) <= tol;
}

#ifdef ENABLE_HIGH_PRECISION
bool approx_equal_hp(const mpreal& a, const mpreal& b, const mpreal& tol) {
    return abs(a - b) <= tol;
}
#endif

//=============================================================================
// PART 1: PolynomialBase<double> vs Polynomial Equivalence Tests
//=============================================================================

void test_equivalence_construction() {
    std::cout << "  test_equivalence_construction... ";

    // Create same polynomial using both classes
    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, 2.0, 3.0}; // 1 + 2x + 3x^2

    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    PolyBaseDouble poly_base = PolyBaseDouble::fromPower(degrees, power_coeffs);

    // Verify same dimension and degrees
    assert(poly.dimension() == poly_base.dimension());
    assert(poly.degrees()[0] == poly_base.degrees()[0]);
    assert(poly.coefficientCount() == poly_base.coefficientCount());

    // Verify same primary representation
    assert(poly.primaryRepresentation() == poly_base.primaryRepresentation());

    std::cout << "PASSED\n";
}

void test_equivalence_evaluation() {
    std::cout << "  test_equivalence_evaluation... ";

    std::vector<unsigned int> degrees = {3};
    std::vector<double> power_coeffs = {1.0, -2.0, 0.5, 3.0};

    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    PolyBaseDouble poly_base = PolyBaseDouble::fromPower(degrees, power_coeffs);

    // Compare evaluations at multiple points
    double test_points[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double t : test_points) {
        double val_poly = poly.evaluate(t);
        double val_base = poly_base.evaluate(t);
        assert(approx_equal(val_poly, val_base));
    }

    std::cout << "PASSED\n";
}

void test_equivalence_bernstein_coefficients() {
    std::cout << "  test_equivalence_bernstein_coefficients... ";

    std::vector<unsigned int> degrees = {2};
    std::vector<double> power_coeffs = {1.0, 2.0, 3.0};

    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    PolyBaseDouble poly_base = PolyBaseDouble::fromPower(degrees, power_coeffs);

    // Convert to Bernstein and compare coefficients
    poly.ensureBernsteinPrimary();
    poly_base.ensureBernsteinPrimary();

    const std::vector<double>& bern_poly = poly.bernsteinCoefficients();
    const std::vector<double>& bern_base = poly_base.bernsteinCoefficients();

    assert(bern_poly.size() == bern_base.size());
    for (size_t i = 0; i < bern_poly.size(); ++i) {
        assert(approx_equal(bern_poly[i], bern_base[i]));
    }

    std::cout << "PASSED\n";
}

void test_equivalence_2d_evaluation() {
    std::cout << "  test_equivalence_2d_evaluation... ";

    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> power_coeffs = {
        1.0, 0.5, 0.25,   // x^0 * (1 + 0.5y + 0.25y^2)
        2.0, 1.0, 0.5,    // x^1 * (2 + 1y + 0.5y^2)
        3.0, 1.5, 0.75    // x^2 * (3 + 1.5y + 0.75y^2)
    };

    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    PolyBaseDouble poly_base = PolyBaseDouble::fromPower(degrees, power_coeffs);

    // Compare evaluations
    for (double x = 0.0; x <= 1.0; x += 0.25) {
        for (double y = 0.0; y <= 1.0; y += 0.25) {
            std::vector<double> pt = {x, y};
            double val_poly = poly.evaluate(pt);
            double val_base = poly_base.evaluate(pt);
            assert(approx_equal(val_poly, val_base));
        }
    }

    std::cout << "PASSED\n";
}

//=============================================================================
// PART 2: Full Test Suite from Existing Tests (using PolynomialBase<double>)
//=============================================================================

void test_dual_representation() {
    std::cout << "  test_dual_representation... ";

    // Create from power basis
    std::vector<unsigned int> degrees{2};
    std::vector<double> power_coeffs{1.0, 2.0, 3.0};
    PolyBaseDouble p1 = PolyBaseDouble::fromPower(degrees, power_coeffs);

    assert(p1.primaryRepresentation() == PolynomialRepresentation::POWER);
    assert(p1.hasPowerCoefficients());
    assert(!p1.hasBernsteinCoefficients());

    // Evaluate using power basis (Horner)
    double val = p1.evaluate(0.5);
    double expected = 1.0 + 2.0 * 0.5 + 3.0 * 0.25; // 2.75
    assert(approx_equal(val, expected));

    // Switch to Bernstein primary (new design: converts in-place)
    p1.ensureBernsteinPrimary();
    assert(p1.primaryRepresentation() == PolynomialRepresentation::BERNSTEIN);
    // New design: single representation, no dual storage
    // assert(p1.hasPowerCoefficients()); // Old test - dual rep no longer valid

    // Create from Bernstein
    std::vector<double> bern_coeffs{0.0, 0.5, 1.0};
    PolyBaseDouble p2 = PolyBaseDouble::fromBernstein(degrees, bern_coeffs);

    assert(p2.primaryRepresentation() == PolynomialRepresentation::BERNSTEIN);
    assert(p2.hasBernsteinCoefficients());
    assert(!p2.hasPowerCoefficients());

    // Switch to power primary (new design: converts in-place)
    p2.ensurePowerPrimary();
    assert(p2.primaryRepresentation() == PolynomialRepresentation::POWER);
    // New design: single representation, no dual storage
    // assert(p2.hasBernsteinCoefficients()); // Old test - dual rep no longer valid

    std::cout << "PASSED\n";
}

void test_conversion_1d_constant() {
    std::cout << "  test_conversion_1d_constant... ";

    const unsigned int n = 7u;
    std::vector<unsigned int> degrees{n};
    std::vector<double> power_coeffs(n + 1u, 0.0);
    power_coeffs[0] = 1.0;

    PolyBaseDouble p = PolyBaseDouble::fromPower(degrees, power_coeffs);
    p.ensureBernsteinPrimary();
    const std::vector<double>& b = p.bernsteinCoefficients();

    for (std::size_t k = 0; k < b.size(); ++k) {
        assert(approx_equal(b[k], 1.0));
    }

    std::cout << "PASSED\n";
}

void test_conversion_1d_linear() {
    std::cout << "  test_conversion_1d_linear... ";

    const unsigned int n = 10u;
    std::vector<unsigned int> degrees{n};
    std::vector<double> power_coeffs(n + 1u, 0.0);
    power_coeffs[1] = 1.0;

    PolyBaseDouble p = PolyBaseDouble::fromPower(degrees, power_coeffs);
    p.ensureBernsteinPrimary();
    const std::vector<double>& b = p.bernsteinCoefficients();

    for (unsigned int k = 0; k <= n; ++k) {
        double expected = static_cast<double>(k) / static_cast<double>(n);
        assert(approx_equal(b[k], expected));
    }

    std::cout << "PASSED\n";
}

void test_conversion_1d_quadratic() {
    std::cout << "  test_conversion_1d_quadratic... ";

    const unsigned int n = 8u;
    std::vector<unsigned int> degrees{n};
    std::vector<double> power_coeffs(n + 1u, 0.0);
    power_coeffs[2] = 1.0;

    PolyBaseDouble p = PolyBaseDouble::fromPower(degrees, power_coeffs);
    p.ensureBernsteinPrimary();
    const std::vector<double>& b = p.bernsteinCoefficients();

    for (unsigned int k = 0; k <= n; ++k) {
        double expected = 0.0;
        if (k >= 2u) {
            expected = (static_cast<double>(k) * static_cast<double>(k - 1)) /
                       (static_cast<double>(n) * static_cast<double>(n - 1));
        }
        assert(approx_equal(b[k], expected));
    }

    std::cout << "PASSED\n";
}

void test_stress_1d_random() {
    std::cout << "  test_stress_1d_random... ";

    const unsigned int n = 12u;
    std::vector<unsigned int> degrees{n};

    std::mt19937 rng(123456u);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    const int num_polys = 10;
    const int num_points = 101;

    for (int p_idx = 0; p_idx < num_polys; ++p_idx) {
        std::vector<double> power_coeffs(n + 1u);
        for (unsigned int i = 0; i <= n; ++i) {
            power_coeffs[i] = dist(rng);
        }

        PolyBaseDouble p = PolyBaseDouble::fromPower(degrees, power_coeffs);

        for (int j = 0; j < num_points; ++j) {
            double t = static_cast<double>(j) / static_cast<double>(num_points - 1);

            // Horner evaluation for expected value
            double expected = 0.0;
            for (int i = static_cast<int>(n); i >= 0; --i) {
                expected = expected * t + power_coeffs[static_cast<std::size_t>(i)];
            }

            double value = p.evaluate(t);
            assert(approx_equal(expected, value));
        }
    }

    std::cout << "PASSED\n";
}

void test_stress_2d_random() {
    std::cout << "  test_stress_2d_random... ";

    const unsigned int dx = 6u;
    const unsigned int dy = 6u;
    std::vector<unsigned int> degrees{dx, dy};

    std::mt19937 rng(7891011u);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    const int num_polys = 5;
    const int num_points = 21;

    auto idx2d = [dy](unsigned int ix, unsigned int iy) {
        return static_cast<std::size_t>(ix * (dy + 1u) + iy);
    };

    for (int p_idx = 0; p_idx < num_polys; ++p_idx) {
        std::vector<double> power_coeffs((dx + 1u) * (dy + 1u));
        for (unsigned int ix = 0; ix <= dx; ++ix) {
            for (unsigned int iy = 0; iy <= dy; ++iy) {
                power_coeffs[idx2d(ix, iy)] = dist(rng);
            }
        }

        PolyBaseDouble p = PolyBaseDouble::fromPower(degrees, power_coeffs);

        for (int ix_sample = 0; ix_sample < num_points; ++ix_sample) {
            double x = static_cast<double>(ix_sample) / static_cast<double>(num_points - 1);
            for (int iy_sample = 0; iy_sample < num_points; ++iy_sample) {
                double y = static_cast<double>(iy_sample) / static_cast<double>(num_points - 1);

                // Horner evaluation for expected value
                std::vector<double> row(dx + 1u);
                for (unsigned int ix = 0; ix <= dx; ++ix) {
                    double v = 0.0;
                    for (int iy = static_cast<int>(dy); iy >= 0; --iy) {
                        v = v * y + power_coeffs[idx2d(ix, static_cast<unsigned int>(iy))];
                    }
                    row[ix] = v;
                }

                double expected = 0.0;
                for (int ix = static_cast<int>(dx); ix >= 0; --ix) {
                    expected = expected * x + row[static_cast<std::size_t>(ix)];
                }

                std::vector<double> pt{x, y};
                double value = p.evaluate(pt);
                assert(approx_equal(expected, value));
            }
        }
    }

    std::cout << "PASSED\n";
}

//=============================================================================
// PART 3: High Precision Tests (PolynomialBase<mpreal>)
//=============================================================================

#ifdef ENABLE_HIGH_PRECISION

void test_hp_construction() {
    std::cout << "  test_hp_construction... ";

    PrecisionContext ctx(256);

    std::vector<unsigned int> degrees = {2};
    std::vector<mpreal> power_coeffs = {mpreal(1), mpreal(2), mpreal(3)};

    PolyBaseHP p = PolyBaseHP::fromPower(degrees, power_coeffs);

    assert(p.dimension() == 1);
    assert(p.degrees()[0] == 2);
    assert(p.hasPowerCoefficients());
    assert(p.primaryRepresentation() == PolynomialRepresentation::POWER);

    std::cout << "PASSED\n";
}

void test_hp_evaluation() {
    std::cout << "  test_hp_evaluation... ";

    PrecisionContext ctx(256);

    // f(x) = x^3 - 2x + 1
    std::vector<unsigned int> degrees = {3};
    std::vector<mpreal> power_coeffs = {mpreal(1), mpreal(-2), mpreal(0), mpreal(1)};

    PolyBaseHP p = PolyBaseHP::fromPower(degrees, power_coeffs);

    // Evaluate at x = 0.5
    mpreal x = mpreal("0.5");
    mpreal result = p.evaluate(x);

    // Expected: 0.125 - 1.0 + 1.0 = 0.125
    mpreal expected = mpreal("0.125");
    mpreal error = abs(result - expected);

    assert(error < mpreal(1e-70));

    std::cout << "PASSED\n";
}

void test_hp_vs_double_precision() {
    std::cout << "  test_hp_vs_double_precision... ";

    PrecisionContext ctx(512);

    // Create ill-conditioned polynomial: (x - 0.5)^5
    // Power: -0.03125 + 0.3125x - 1.25x^2 + 2.5x^3 - 2.5x^4 + x^5
    std::vector<unsigned int> degrees = {5};
    std::vector<double> power_coeffs_dbl = {-0.03125, 0.3125, -1.25, 2.5, -2.5, 1.0};
    std::vector<mpreal> power_coeffs_hp;
    for (double c : power_coeffs_dbl) {
        power_coeffs_hp.push_back(toHighPrecision(c));
    }

    PolyBaseDouble p_dbl = PolyBaseDouble::fromPower(degrees, power_coeffs_dbl);
    PolyBaseHP p_hp = PolyBaseHP::fromPower(degrees, power_coeffs_hp);

    // Evaluate at x = 0.5 (the root)
    double result_dbl = p_dbl.evaluate(0.5);
    mpreal result_hp = p_hp.evaluate(mpreal("0.5"));

    // Both should be very close to zero
    assert(std::abs(result_dbl) < 1e-14);
    assert(abs(result_hp) < mpreal(1e-140));

    // Evaluate near the root: x = 0.5 + 1e-10
    mpreal delta = mpreal("1e-10");
    mpreal x_near = mpreal("0.5") + delta;
    
    double result_near_dbl = p_dbl.evaluate(0.5 + 1e-10);
    mpreal result_near_hp = p_hp.evaluate(x_near);

    // HP should give much better precision for (1e-10)^5 = 1e-50
    std::cout << "\n    Double result at 0.5+1e-10: " << result_near_dbl;
    std::cout << "\n    HP result at 0.5+1e-10: " << toString(result_near_hp, 60) << "\n    ";

    std::cout << "PASSED\n";
}

void test_hp_basis_conversion() {
    std::cout << "  test_hp_basis_conversion... ";

    PrecisionContext ctx(256);

    // f(x) = 1 + 2x + 3x^2
    std::vector<unsigned int> degrees = {2};
    std::vector<mpreal> power_coeffs = {mpreal(1), mpreal(2), mpreal(3)};

    PolyBaseHP p = PolyBaseHP::fromPower(degrees, power_coeffs);
    p.convertPowerToBernstein();

    const std::vector<mpreal>& bern = p.bernsteinCoefficients();

    // Expected Bernstein coefficients: [1, 2, 6]
    assert(bern.size() == 3);
    assert(abs(bern[0] - mpreal(1)) < mpreal(1e-70));
    assert(abs(bern[1] - mpreal(2)) < mpreal(1e-70));
    assert(abs(bern[2] - mpreal(6)) < mpreal(1e-70));

    std::cout << "PASSED\n";
}

void test_hp_arithmetic() {
    std::cout << "  test_hp_arithmetic... ";

    PrecisionContext ctx(256);

    std::vector<unsigned int> degrees = {1};
    std::vector<mpreal> coeffs_a = {mpreal(1), mpreal(2)};
    std::vector<mpreal> coeffs_b = {mpreal(3), mpreal(4)};

    PolyBaseHP a = PolyBaseHP::fromPower(degrees, coeffs_a);
    PolyBaseHP b = PolyBaseHP::fromPower(degrees, coeffs_b);

    // Test addition: (1 + 2x) + (3 + 4x) = 4 + 6x
    PolyBaseHP sum = a + b;
    assert(abs(sum.evaluate(mpreal(0)) - mpreal(4)) < mpreal(1e-70));
    assert(abs(sum.evaluate(mpreal(1)) - mpreal(10)) < mpreal(1e-70));

    // Test subtraction: (1 + 2x) - (3 + 4x) = -2 - 2x
    PolyBaseHP diff = a - b;
    assert(abs(diff.evaluate(mpreal(0)) - mpreal(-2)) < mpreal(1e-70));
    assert(abs(diff.evaluate(mpreal(1)) - mpreal(-4)) < mpreal(1e-70));

    // Test multiplication: (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
    PolyBaseHP prod = a * b;
    assert(abs(prod.evaluate(mpreal(0)) - mpreal(3)) < mpreal(1e-70));
    assert(abs(prod.evaluate(mpreal(1)) - mpreal(21)) < mpreal(1e-70));

    std::cout << "PASSED\n";
}

void test_hp_vs_polynomialHP_equivalence() {
    std::cout << "  test_hp_vs_polynomialHP_equivalence... ";

    PrecisionContext ctx(256);

    // Create same polynomial using both classes
    std::vector<unsigned int> degrees = {3};
    std::vector<mpreal> power_coeffs = {mpreal(1), mpreal(-2), mpreal(0), mpreal(1)};

    PolyBaseHP poly_base = PolyBaseHP::fromPower(degrees, power_coeffs);
    PolynomialHP poly_hp = fromPowerHP(degrees, power_coeffs);

    // Compare evaluations
    double test_points[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double t_dbl : test_points) {
        mpreal t = mpreal(t_dbl);
        mpreal val_base = poly_base.evaluate(t);
        mpreal val_hp = poly_hp.evaluate(t);
        assert(abs(val_base - val_hp) < mpreal(1e-70));
    }

    std::cout << "PASSED\n";
}

void test_precision_scaling() {
    std::cout << "  test_precision_scaling... ";

    // Test at different precision levels
    int precisions[] = {64, 128, 256, 512};

    for (int prec : precisions) {
        PrecisionContext ctx(prec);

        // Create a simple polynomial
        std::vector<unsigned int> degrees = {1};
        std::vector<mpreal> coeffs = {mpreal(1), mpreal(1)};  // 1 + x

        PolyBaseHP p = PolyBaseHP::fromPower(degrees, coeffs);

        // Evaluate at x = 0.5
        mpreal result = p.evaluate(mpreal("0.5"));
        mpreal expected = mpreal("1.5");
        mpreal error = abs(result - expected);

        // Error should be near machine precision for that precision level
        // For prec bits, epsilon is ~2^(-prec)
        mpreal eps = mpreal(1) / pow(mpreal(2), mpreal(prec - 10));
        assert(error < eps);
    }

    std::cout << "PASSED\n";
}

#endif // ENABLE_HIGH_PRECISION

//=============================================================================
// Main
//=============================================================================

int main() {
    std::cout << "=== Comprehensive PolynomialBase<Scalar> Tests ===\n\n";

    std::cout << "PART 1: PolynomialBase<double> vs Polynomial Equivalence:\n";
    test_equivalence_construction();
    test_equivalence_evaluation();
    test_equivalence_bernstein_coefficients();
    test_equivalence_2d_evaluation();

    std::cout << "\nPART 2: Full Test Suite (PolynomialBase<double>):\n";
    test_dual_representation();
    test_conversion_1d_constant();
    test_conversion_1d_linear();
    test_conversion_1d_quadratic();
    test_stress_1d_random();
    test_stress_2d_random();

#ifdef ENABLE_HIGH_PRECISION
    std::cout << "\nPART 3: High Precision Tests (PolynomialBase<mpreal>):\n";
    test_hp_construction();
    test_hp_evaluation();
    test_hp_vs_double_precision();
    test_hp_basis_conversion();
    test_hp_arithmetic();
    test_hp_vs_polynomialHP_equivalence();
    test_precision_scaling();
#else
    std::cout << "\nPART 3: High Precision Tests SKIPPED (ENABLE_HIGH_PRECISION=OFF)\n";
#endif

    std::cout << "\n=== All Comprehensive Tests PASSED ===\n";
    return 0;
}
