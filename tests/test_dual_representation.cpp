#include "core/polynomial.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace polynomial_solver;

bool approx_equal(double a, double b, double tol = 1e-10) {
    return std::fabs(a - b) < tol;
}

int main() {
    std::cout << "=== Test Dual-Representation Polynomial ===" << std::endl;

    // Test 1: Create from power, verify primary is POWER
    std::cout << "\nTest 1: Create from power basis" << std::endl;
    std::vector<unsigned int> degrees{2};
    std::vector<double> power_coeffs{1.0, 2.0, 3.0}; // 1 + 2x + 3x^2
    
    Polynomial p1 = Polynomial::fromPower(degrees, power_coeffs);
    
    if (p1.primaryRepresentation() != PolynomialRepresentation::POWER) {
        std::cerr << "FAIL: Primary should be POWER" << std::endl;
        return 1;
    }
    
    if (!p1.hasPowerCoefficients()) {
        std::cerr << "FAIL: Power coefficients should be valid" << std::endl;
        return 1;
    }
    
    if (p1.hasBernsteinCoefficients()) {
        std::cerr << "FAIL: Bernstein coefficients should not be computed yet" << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ Primary is POWER, power valid, Bernstein not yet computed" << std::endl;

    // Test 2: Evaluate using power basis (should use Horner's method)
    std::cout << "\nTest 2: Evaluate using power basis" << std::endl;
    double val = p1.evaluate(0.5);
    double expected = 1.0 + 2.0 * 0.5 + 3.0 * 0.5 * 0.5; // 1 + 1 + 0.75 = 2.75
    
    if (!approx_equal(val, expected)) {
        std::cerr << "FAIL: Expected " << expected << ", got " << val << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ Evaluation correct: " << val << std::endl;

    // Test 3: Explicitly convert to Bernstein (without changing primary)
    std::cout << "\nTest 3: Explicitly convert to Bernstein (without changing primary)" << std::endl;

    // Create a copy and convert it to Bernstein
    Polynomial p1_copy = p1;
    p1_copy.ensureBernsteinPrimary();

    // Now we can access Bernstein coefficients
    const std::vector<double>& bern = p1_copy.bernsteinCoefficients();

    if (!p1_copy.hasBernsteinCoefficients()) {
        std::cerr << "FAIL: Bernstein coefficients should now be valid" << std::endl;
        return 1;
    }

    if (p1_copy.primaryRepresentation() != PolynomialRepresentation::BERNSTEIN) {
        std::cerr << "FAIL: Primary should now be BERNSTEIN after ensureBernsteinPrimary()" << std::endl;
        return 1;
    }

    std::cout << "  ✓ Bernstein coefficients computed, primary switched to BERNSTEIN" << std::endl;
    std::cout << "  Bernstein coeffs: [" << bern[0] << ", " << bern[1] << ", " << bern[2] << "]" << std::endl;

    // Test 4: Switch primary to Bernstein
    std::cout << "\nTest 4: Switch primary to Bernstein" << std::endl;
    p1.ensureBernsteinPrimary();
    
    if (p1.primaryRepresentation() != PolynomialRepresentation::BERNSTEIN) {
        std::cerr << "FAIL: Primary should now be BERNSTEIN" << std::endl;
        return 1;
    }
    
    if (!p1.hasPowerCoefficients()) {
        std::cerr << "FAIL: Power coefficients should still be cached" << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ Primary switched to BERNSTEIN, power still cached" << std::endl;

    // Test 5: Create from Bernstein, verify primary is BERNSTEIN
    std::cout << "\nTest 5: Create from Bernstein basis" << std::endl;
    std::vector<double> bern_coeffs{0.0, 0.5, 1.0}; // Linear function
    Polynomial p2 = Polynomial::fromBernstein(degrees, bern_coeffs);
    
    if (p2.primaryRepresentation() != PolynomialRepresentation::BERNSTEIN) {
        std::cerr << "FAIL: Primary should be BERNSTEIN" << std::endl;
        return 1;
    }
    
    if (!p2.hasBernsteinCoefficients()) {
        std::cerr << "FAIL: Bernstein coefficients should be valid" << std::endl;
        return 1;
    }
    
    if (p2.hasPowerCoefficients()) {
        std::cerr << "FAIL: Power coefficients should not be computed yet" << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ Primary is BERNSTEIN, Bernstein valid, power not yet computed" << std::endl;

    // Test 6: Switch to power primary
    std::cout << "\nTest 6: Switch primary to power" << std::endl;
    p2.ensurePowerPrimary();
    
    if (p2.primaryRepresentation() != PolynomialRepresentation::POWER) {
        std::cerr << "FAIL: Primary should now be POWER" << std::endl;
        return 1;
    }
    
    if (!p2.hasBernsteinCoefficients()) {
        std::cerr << "FAIL: Bernstein coefficients should still be cached" << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ Primary switched to POWER, Bernstein still cached" << std::endl;

    std::cout << "\n=== All tests passed! ===" << std::endl;
    return 0;
}

