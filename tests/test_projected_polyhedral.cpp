#include "core/polynomial.h"
#include "solver/solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

using namespace polynomial_solver;

// Test 1D linear function with PP method
int test_1d_linear_pp() {
    std::cout << "Test 1: 1D linear function p(x) = x - 0.5 (PP method)" << std::endl;
    std::cout << "Expected root: x = 0.5" << std::endl;
    
    // Create polynomial p(x) = x - 0.5
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{-0.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    
    PolynomialSystem system({p});
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "  Found " << result.boxes.size() << " box(es), " 
              << result.num_resolved << " resolved" << std::endl;
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    
    if (result.num_resolved != 1) {
        std::cerr << "  FAIL: Expected 1 resolved root, got " << result.num_resolved << std::endl;
        return 1;
    }
    
    const SubdivisionBoxResult& box = result.boxes[0];
    std::cout << "  Root box:" << std::endl;
    std::cout << "    Center: " << std::fixed << std::setprecision(16) << box.center[0] << std::endl;
    std::cout << "    Max error: " << std::scientific << std::setprecision(16) << box.max_error[0] << std::endl;
    std::cout << "    Depth: " << box.depth << std::endl;
    
    // Check if error is machine epsilon
    const double machine_eps = std::numeric_limits<double>::epsilon();
    std::cout << "  Machine epsilon: " << std::scientific << std::setprecision(16) << machine_eps << std::endl;
    std::cout << "  Error is machine epsilon: " << (box.max_error[0] == machine_eps ? "YES" : "NO") << std::endl;
    
    // Evaluate at the center
    double eval_value = p.evaluate(box.center[0]);
    std::cout << "  Evaluation: p(" << std::fixed << std::setprecision(16) << box.center[0] << ") = " 
              << std::scientific << std::setprecision(16) << eval_value << std::endl;
    
    if (std::fabs(eval_value) < 1e-10) {
        std::cout << "  PASS: PP method gives exact root!" << std::endl;
        return 0;
    } else {
        std::cerr << "  FAIL: Evaluation error too large" << std::endl;
        return 1;
    }
}

// Test 2D linear system with PP method
int test_2d_linear_pp() {
    std::cout << "Test 2: 2D linear system (PP method)" << std::endl;
    std::cout << "  p1(x,y) = x - 0.5" << std::endl;
    std::cout << "  p2(x,y) = y - 0.3" << std::endl;
    std::cout << "Expected root: (0.5, 0.3)" << std::endl;

    std::vector<unsigned int> degrees{1u, 1u};

    // p1(x,y) = x - 0.5
    std::vector<double> power1(4, 0.0);
    power1[0] = -0.5;  // constant
    power1[2] = 1.0;   // x term
    Polynomial p1 = Polynomial::fromPower(degrees, power1);

    // p2(x,y) = y - 0.3
    std::vector<double> power2(4, 0.0);
    power2[0] = -0.3;  // constant
    power2[1] = 1.0;   // y term
    Polynomial p2 = Polynomial::fromPower(degrees, power2);
    
    PolynomialSystem system({p1, p2});
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "  Found " << result.boxes.size() << " box(es), " 
              << result.num_resolved << " resolved" << std::endl;
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    
    if (result.num_resolved != 1) {
        std::cerr << "  FAIL: Expected 1 resolved root, got " << result.num_resolved << std::endl;
        return 1;
    }
    
    const SubdivisionBoxResult& box = result.boxes[0];
    std::cout << "  Root box:" << std::endl;
    std::cout << "    Center: (" << std::fixed << std::setprecision(16) 
              << box.center[0] << ", " << box.center[1] << ")" << std::endl;
    std::cout << "    Max error: (" << std::scientific << std::setprecision(16) 
              << box.max_error[0] << ", " << box.max_error[1] << ")" << std::endl;
    std::cout << "    Depth: " << box.depth << std::endl;
    
    // Check if error is machine epsilon
    const double machine_eps = std::numeric_limits<double>::epsilon();
    std::cout << "  Machine epsilon: " << std::scientific << std::setprecision(16) << machine_eps << std::endl;
    std::cout << "  Error is machine epsilon: (" 
              << (box.max_error[0] == machine_eps ? "YES" : "NO") << ", "
              << (box.max_error[1] == machine_eps ? "YES" : "NO") << ")" << std::endl;
    
    // Evaluate at the center
    std::vector<double> eval_values;
    system.evaluate(box.center, eval_values);
    std::cout << "  Evaluation:" << std::endl;
    std::cout << "    p1 = " << std::scientific << std::setprecision(16) << eval_values[0] << std::endl;
    std::cout << "    p2 = " << eval_values[1] << std::endl;
    
    if (std::fabs(eval_values[0]) < 1e-10 && std::fabs(eval_values[1]) < 1e-10) {
        std::cout << "  PASS: PP method gives exact root!" << std::endl;
        return 0;
    } else {
        std::cerr << "  FAIL: Evaluation error too large" << std::endl;
        return 1;
    }
}

// Test 2D quadratic system with PP method
int test_2d_quadratic_pp() {
    std::cout << "Test 3: 2D quadratic system (PP method)" << std::endl;
    std::cout << "  p1(x,y) = x^2 - 0.25" << std::endl;
    std::cout << "  p2(x,y) = y - 0.3" << std::endl;
    std::cout << "Expected root: (0.5, 0.3)" << std::endl;

    // f1 = x^2 - 0.25 in power basis with degrees (2, 1)
    std::vector<unsigned int> degrees_f1{2u, 1u};
    std::vector<double> power_coeffs_f1{-0.25, 0.0, 0.0, 0.0, 1.0, 0.0};
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs_f1);

    // f2 = y - 0.3 in power basis with degrees (0, 1)
    std::vector<unsigned int> degrees_f2{0u, 1u};
    std::vector<double> power_coeffs_f2{-0.3, 1.0};
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs_f2);

    PolynomialSystem system({p1, p2});
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 0.05;  // Use same tolerance as GraphHull test
    config.max_depth = 100;
    config.degeneracy_multiplier = 10.0;  // Increase to avoid false degeneracy detection

    SubdivisionSolverResult result = solver.subdivisionSolve(system, config, RootBoundingMethod::ProjectedPolyhedral);

    std::cout << "  Found " << result.boxes.size() << " box(es), "
              << result.num_resolved << " resolved" << std::endl;
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << std::endl;

    if (result.boxes.empty()) {
        std::cerr << "  FAIL: Expected at least one box" << std::endl;
        return 1;
    }

    // Check that at least one box contains the root (0.5, 0.3)
    bool found_root = false;
    for (const SubdivisionBoxResult& box : result.boxes) {
        if (box.lower[0] <= 0.5 && box.upper[0] >= 0.5 &&
            box.lower[1] <= 0.3 && box.upper[1] >= 0.3) {
            found_root = true;
            std::cout << "  Found root box: [" << std::fixed << std::setprecision(4)
                      << box.lower[0] << ", " << box.upper[0] << "] x ["
                      << box.lower[1] << ", " << box.upper[1] << "]" << std::endl;
            std::cout << "    Center: (" << box.center[0] << ", " << box.center[1] << ")" << std::endl;
            std::cout << "    Depth: " << box.depth << std::endl;
            break;
        }
    }

    if (found_root) {
        std::cout << "  PASS: PP method found the root!" << std::endl;
        return 0;
    } else {
        std::cerr << "  FAIL: Root not found in any box" << std::endl;
        return 1;
    }
}

int main() {
    std::cout << "Testing Projected Polyhedral (PP) method..." << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << std::endl;

    int failures = 0;
    failures += test_1d_linear_pp();
    std::cout << std::endl;
    failures += test_2d_linear_pp();
    std::cout << std::endl;
    failures += test_2d_quadratic_pp();

    std::cout << std::endl;
    if (failures == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << failures << " test(s) failed." << std::endl;
        return 1;
    }
}

