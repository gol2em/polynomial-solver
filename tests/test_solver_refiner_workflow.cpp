#include "core/polynomial.h"
#include "solver/solver.h"
#include "refinement/result_refiner.h"
#include "core/differentiation.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace polynomial_solver;

bool approx_equal(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) < tol;
}

int main() {
    std::cout << "=== Test Solver-Refiner Workflow with Dual Representation ===" << std::endl;

    // Test: Create polynomial from power basis
    // p(x) = (x - 0.3)(x - 0.7) = x^2 - x + 0.21
    std::cout << "\n1. Create polynomial from power basis" << std::endl;
    std::vector<unsigned int> degrees{2};
    std::vector<double> power_coeffs{0.21, -1.0, 1.0};
    
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    
    if (poly.primaryRepresentation() != PolynomialRepresentation::POWER) {
        std::cerr << "FAIL: Primary should be POWER after fromPower()" << std::endl;
        return 1;
    }
    std::cout << "  ✓ Polynomial created with POWER as primary" << std::endl;

    // Test: Solver needs Bernstein, so it should call ensureBernsteinPrimary()
    std::cout << "\n2. Solver workflow: switch to Bernstein primary" << std::endl;
    poly.ensureBernsteinPrimary();
    
    if (poly.primaryRepresentation() != PolynomialRepresentation::BERNSTEIN) {
        std::cerr << "FAIL: Primary should be BERNSTEIN after ensureBernsteinPrimary()" << std::endl;
        return 1;
    }
    
    if (!poly.hasPowerCoefficients()) {
        std::cerr << "FAIL: Power coefficients should still be cached" << std::endl;
        return 1;
    }
    std::cout << "  ✓ Primary switched to BERNSTEIN, power still cached" << std::endl;

    // Test: Verify solver can use Bernstein coefficients
    std::cout << "\n3. Solver uses Bernstein coefficients" << std::endl;
    const std::vector<double>& bern = poly.bernsteinCoefficients();
    std::cout << "  Bernstein coeffs: [" << bern[0] << ", " << bern[1] << ", " << bern[2] << "]" << std::endl;
    std::cout << "  ✓ Solver can access Bernstein coefficients" << std::endl;

    // Test: Refiner uses original polynomial (which still has power cached)
    std::cout << "\n4. Refiner workflow: use power basis for Newton" << std::endl;
    
    // Create a fresh polynomial from power (simulating original input)
    Polynomial poly_original = Polynomial::fromPower(degrees, power_coeffs);
    
    if (!poly_original.hasPowerCoefficients()) {
        std::cerr << "FAIL: Original polynomial should have power coefficients" << std::endl;
        return 1;
    }
    
    // Test power-basis differentiation
    Polynomial dpoly = Differentiation::differentiateAxisPower(poly_original, 0);
    
    if (dpoly.primaryRepresentation() != PolynomialRepresentation::POWER) {
        std::cerr << "FAIL: Derivative should have POWER as primary" << std::endl;
        return 1;
    }
    std::cout << "  ✓ Derivative computed in power basis" << std::endl;

    // Test: Verify derivative is correct
    // p(x) = 0.21 - x + x^2, so p'(x) = -1 + 2x
    std::cout << "\n5. Verify power-basis derivative correctness" << std::endl;
    double x_test = 0.5;
    double dp_val = dpoly.evaluate(x_test);
    double expected = -1.0 + 2.0 * x_test; // -1 + 1 = 0
    
    if (!approx_equal(dp_val, expected)) {
        std::cerr << "FAIL: Derivative at x=" << x_test << " should be " << expected 
                  << ", got " << dp_val << std::endl;
        return 1;
    }
    std::cout << "  ✓ Derivative correct: p'(0.5) = " << dp_val << std::endl;

    // Test: Newton iteration using power basis
    std::cout << "\n6. Newton iteration using power basis" << std::endl;
    double x0 = 0.25; // Initial guess near root at 0.3
    double x = x0;
    
    for (int iter = 0; iter < 5; ++iter) {
        double f = poly_original.evaluate(x);
        double df = dpoly.evaluate(x);
        
        if (std::abs(df) < 1e-14) {
            std::cerr << "FAIL: Derivative too small" << std::endl;
            return 1;
        }
        
        double x_new = x - f / df;
        std::cout << "  Iter " << iter << ": x = " << x << ", f(x) = " << f << std::endl;
        
        if (std::abs(f) < 1e-10) {
            break;
        }
        
        x = x_new;
    }
    
    // Check if converged to root at 0.3
    if (!approx_equal(x, 0.3, 1e-8)) {
        std::cerr << "FAIL: Should converge to 0.3, got " << x << std::endl;
        return 1;
    }
    std::cout << "  ✓ Converged to root: x = " << x << std::endl;

    // Test: Verify both representations give same evaluation
    std::cout << "\n7. Verify consistency between representations" << std::endl;
    Polynomial poly_bern = Polynomial::fromPower(degrees, power_coeffs);
    poly_bern.ensureBernsteinPrimary();
    
    Polynomial poly_pow = Polynomial::fromPower(degrees, power_coeffs);
    // poly_pow already has power as primary
    
    for (double t = 0.0; t <= 1.0; t += 0.1) {
        double val_bern = poly_bern.evaluate(t);
        double val_pow = poly_pow.evaluate(t);
        
        if (!approx_equal(val_bern, val_pow, 1e-12)) {
            std::cerr << "FAIL: Evaluations differ at t=" << t 
                      << ": Bernstein=" << val_bern << ", Power=" << val_pow << std::endl;
            return 1;
        }
    }
    std::cout << "  ✓ Both representations give identical evaluations" << std::endl;

    std::cout << "\n=== All workflow tests passed! ===" << std::endl;
    return 0;
}

