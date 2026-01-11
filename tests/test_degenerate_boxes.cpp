#include "core/polynomial.h"
#include "solver/solver.h"
#include "core/geometry.h"
#include <iostream>
#include <cmath>
#include <limits>

using namespace polynomial_solver;

namespace {

bool approx_equal(double a, double b, double eps = 1e-6) {
    return std::fabs(a - b) <= eps;
}

// Test 1: Single point root (0D degenerate case)
int test_single_point_root() {
    std::cout << "Test 1: Single point root..." << std::endl;

    // 1D system: p(x) = (x - 0.5)
    // Root at x = 0.5
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{-0.5, 1.0};  // -0.5 + x
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    PolynomialSystem system({p});

    // Enable debug output
    getGeometryConfig().debug = false;  // Keep it off for now

    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 0.05;  // Use same tolerance as existing test
    config.max_depth = 20;

    SubdivisionSolverResult result =
        solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);

    // Should find at least one box containing the root
    if (result.boxes.empty()) {
        std::cerr << "  FAIL: No boxes found" << std::endl;
        return 1;
    }

    std::cout << "  Found " << result.boxes.size() << " box(es), "
              << result.num_resolved << " resolved" << std::endl;
    for (size_t i = 0; i < result.boxes.size(); ++i) {
        std::cout << "    Box " << i << ": [" << result.boxes[i].lower[0] << ", "
                  << result.boxes[i].upper[0] << "], width = "
                  << (result.boxes[i].upper[0] - result.boxes[i].lower[0])
                  << ", converged = " << result.boxes[i].converged << std::endl;
    }

    // Check if any box contains x = 0.5
    bool found_root = false;
    for (const auto& box : result.boxes) {
        if (box.lower[0] <= 0.5 && box.upper[0] >= 0.5) {
            found_root = true;

            // Verify center is close to root
            double center = 0.5 * (box.lower[0] + box.upper[0]);
            double value = p.evaluate(center);
            std::cout << "    Center: " << center << ", value: " << value << std::endl;

            if (std::fabs(value) > 0.1) {  // Very loose tolerance
                std::cerr << "  FAIL: Center is not close to root, value = " << value << std::endl;
                return 1;
            }

            break;
        }
    }

    if (!found_root) {
        std::cerr << "  FAIL: Root not found in any box" << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 2: Test isApproximateRoot method
int test_is_approximate_root() {
    std::cout << "Test 2: isApproximateRoot method..." << std::endl;

    // 1D system: p(x) = x - 0.5
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{-0.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    PolynomialSystem system({p});

    // Test point at root
    std::vector<double> root_point{0.5};
    if (!system.isApproximateRoot(root_point, 1e-6)) {
        std::cerr << "  FAIL: Root point not recognized" << std::endl;
        return 1;
    }

    // Test point near root
    std::vector<double> near_root{0.501};
    if (!system.isApproximateRoot(near_root, 0.01)) {
        std::cerr << "  FAIL: Near-root point not recognized with loose tolerance" << std::endl;
        return 1;
    }

    // Test point far from root
    std::vector<double> far_point{0.9};
    if (system.isApproximateRoot(far_point, 0.01)) {
        std::cerr << "  FAIL: Far point incorrectly recognized as root" << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

} // anonymous namespace

// Test 3: Test evaluate method
int test_evaluate() {
    std::cout << "Test 3: evaluate method..." << std::endl;

    // 2D system: p1(x,y) = x - 0.5, p2(x,y) = y - 0.3
    std::vector<unsigned int> degrees{1u, 1u};

    std::vector<double> power1(4, 0.0);
    power1[0] = -0.5;
    power1[2] = 1.0;
    Polynomial p1 = Polynomial::fromPower(degrees, power1);

    std::vector<double> power2(4, 0.0);
    power2[0] = -0.3;
    power2[1] = 1.0;
    Polynomial p2 = Polynomial::fromPower(degrees, power2);

    PolynomialSystem system({p1, p2});

    // Evaluate at root
    std::vector<double> root_point{0.5, 0.3};
    std::vector<double> values;
    system.evaluate(root_point, values);

    if (values.size() != 2) {
        std::cerr << "  FAIL: Expected 2 values, got " << values.size() << std::endl;
        return 1;
    }

    if (std::fabs(values[0]) > 1e-10 || std::fabs(values[1]) > 1e-10) {
        std::cerr << "  FAIL: Values at root not zero: " << values[0] << ", " << values[1] << std::endl;
        return 1;
    }

    // Evaluate at non-root
    std::vector<double> other_point{0.7, 0.8};
    system.evaluate(other_point, values);

    if (std::fabs(values[0] - 0.2) > 1e-10 || std::fabs(values[1] - 0.5) > 1e-10) {
        std::cerr << "  FAIL: Values at (0.7, 0.8) incorrect: " << values[0] << ", " << values[1] << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

int main() {
    std::cout << "Testing degenerate box handling..." << std::endl;
    std::cout << std::endl;

    int failures = 0;
    failures += test_single_point_root();
    failures += test_is_approximate_root();
    failures += test_evaluate();

    std::cout << std::endl;
    if (failures == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << failures << " test(s) failed." << std::endl;
        return 1;
    }
}

