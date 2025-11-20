/**
 * Example: Wilkinson Polynomial
 *
 * This example demonstrates solving the famous Wilkinson polynomial,
 * a notoriously ill-conditioned polynomial with 20 roots.
 *
 * Problem:
 *   p(x) = (x-1)(x-2)(x-3)...(x-20)
 *
 * Expected roots: x = 1, 2, 3, ..., 20
 * Domain: [0, 21]
 *
 * This is a challenging test case due to:
 * - Large number of roots (20)
 * - Wide range of root locations
 * - Numerical sensitivity (Wilkinson's polynomial is famously ill-conditioned)
 */

#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using namespace polynomial_solver;

void print_result(const SubdivisionSolverResult& result) {
    std::cout << "Found " << result.num_resolved << " root(s)\n";
    std::cout << "Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";
    std::cout << "Unresolved boxes: " << (result.boxes.size() - result.num_resolved) << "\n\n";

    std::cout << "Resolved roots:\n";
    for (size_t i = 0; i < result.num_resolved; ++i) {
        const SubdivisionBoxResult& box = result.boxes[i];
        std::cout << "  Root " << (i+1) << ":\n";
        std::cout << "    Interval: [" << std::fixed << std::setprecision(10)
                  << box.lower[0] << ", " << box.upper[0] << "]\n";
        std::cout << "    Center: " << box.center[0] << "\n";
        std::cout << "    Width: " << std::scientific << std::setprecision(6)
                  << (box.upper[0] - box.lower[0]) << "\n";
        std::cout << "    Depth: " << box.depth << "\n";
        std::cout << "\n";
    }

    if (result.boxes.size() > result.num_resolved) {
        std::cout << "Unresolved boxes:\n";
        for (size_t i = result.num_resolved; i < result.boxes.size(); ++i) {
            const SubdivisionBoxResult& box = result.boxes[i];
            std::cout << "  Box " << (i - result.num_resolved + 1) << ":\n";
            std::cout << "    Interval: [" << std::fixed << std::setprecision(10)
                      << box.lower[0] << ", " << box.upper[0] << "]\n";
            std::cout << "    Center: " << box.center[0] << "\n";
            std::cout << "    Width: " << std::scientific << std::setprecision(6)
                      << (box.upper[0] - box.lower[0]) << "\n";
            std::cout << "    Depth: " << box.depth << "\n";
            std::cout << "\n";
        }
    }
}

void dump_result(const SubdivisionSolverResult& result, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }

    out << "# Solver Result Dump\n";
    out << "# Dimension: 1\n";
    out << "# Total boxes: " << result.boxes.size() << "\n";
    out << "# Resolved: " << result.num_resolved << "\n";
    out << "# Unresolved: " << (result.boxes.size() - result.num_resolved) << "\n";
    out << "# Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";
    out << "\n";

    // Dump resolved boxes
    out << "## Resolved Boxes\n";
    for (size_t i = 0; i < result.num_resolved; ++i) {
        const SubdivisionBoxResult& box = result.boxes[i];
        out << "Box " << i << "\n";
        out << "  Lower: " << std::setprecision(17) << box.lower[0] << "\n";
        out << "  Upper: " << box.upper[0] << "\n";
        out << "  Center: " << box.center[0] << "\n";
        out << "  Depth: " << box.depth << "\n";
        out << "  Converged: " << (box.converged ? "yes" : "no") << "\n";
        out << "\n";
    }

    // Dump unresolved boxes
    if (result.boxes.size() > result.num_resolved) {
        out << "## Unresolved Boxes\n";
        for (size_t i = result.num_resolved; i < result.boxes.size(); ++i) {
            const SubdivisionBoxResult& box = result.boxes[i];
            out << "Box " << (i - result.num_resolved) << "\n";
            out << "  Lower: " << std::setprecision(17) << box.lower[0] << "\n";
            out << "  Upper: " << box.upper[0] << "\n";
            out << "  Center: " << box.center[0] << "\n";
            out << "  Depth: " << box.depth << "\n";
            out << "  Converged: " << (box.converged ? "yes" : "no") << "\n";
            out << "\n";
        }
    }

    out.close();
    std::cout << "Result dump saved to: " << filename << "\n";
}

// Compute coefficients of (x-1/n)(x-2/n)...(x-(n-1)/n) iteratively
// This scales the Wilkinson polynomial to [0, 1] domain
std::vector<double> wilkinson_coefficients_scaled(int n) {
    // Start with (x - 1/n)
    double scale = 1.0 / static_cast<double>(n);
    std::vector<double> coeffs = {-scale, 1.0};  // -1/n + x

    // Multiply by (x - k/n) for k = 2, 3, ..., n-1
    for (int k = 2; k < n; ++k) {
        std::vector<double> new_coeffs(coeffs.size() + 1, 0.0);
        double root = k * scale;

        // Multiply: coeffs * (x - k/n) = coeffs * x - (k/n) * coeffs
        for (size_t i = 0; i < coeffs.size(); ++i) {
            new_coeffs[i] -= root * coeffs[i];      // -(k/n) * coeffs
            new_coeffs[i + 1] += coeffs[i];         // coeffs * x
        }

        coeffs = new_coeffs;
    }

    return coeffs;
}

int main() {
    std::cout << "Wilkinson Polynomial Root Finding\n";
    std::cout << std::string(60, '=') << "\n\n";

    std::cout << "Problem: p(x) = (x-1/21)(x-2/21)(x-3/21)...(x-20/21)\n";
    std::cout << "Expected roots: x = 1/21, 2/21, 3/21, ..., 20/21\n";
    std::cout << "Domain: [0, 1]\n\n";

    // Compute Wilkinson polynomial coefficients (scaled to [0,1])
    // Roots at 1/21, 2/21, ..., 19/21 (19 roots total)
    std::vector<double> power_coeffs = wilkinson_coefficients_scaled(20);

    std::cout << "Polynomial degree: " << (power_coeffs.size() - 1) << "\n";
    std::cout << "Leading coefficient: " << power_coeffs.back() << "\n";
    std::cout << "Constant term: " << std::scientific << std::setprecision(6)
              << power_coeffs[0] << "\n\n";

    std::vector<unsigned int> degrees{static_cast<unsigned int>(power_coeffs.size() - 1)};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});

    // Verify a few roots
    std::cout << "Verification (sample roots):\n";
    int test_roots_idx[] = {1, 5, 10, 15, 19};
    for (int i = 0; i < 5; ++i) {
        double root = test_roots_idx[i] / 20.0;
        double val = p.evaluate(root);
        std::cout << "  p(" << test_roots_idx[i] << "/20 = " << std::fixed << std::setprecision(4) << root
                  << ") = " << std::scientific << std::setprecision(6) << val << "\n";
    }
    std::cout << "\n";
    
    // Test with ProjectedPolyhedral method (generates geometry dump for visualization)
    Solver solver;

    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;
    config.contraction_threshold = 0.8;
    config.strategy = SubdivisionStrategy::SubdivideFirst;
    config.dump_geometry = true;
    config.dump_prefix = "dumps/wilkinson_1d";

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Solving with SubdivideFirst strategy\n";
    std::cout << std::string(60, '=') << "\n\n";

    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);

    print_result(result);

    // Dump result to file
    dump_result(result, "dumps/wilkinson_1d_result.txt");

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "\nGeometry dump saved to dumps/wilkinson_1d_geometry.txt\n";
    std::cout << "Result dump saved to dumps/wilkinson_1d_result.txt\n";
    std::cout << "Visualize with: python examples/visualize_wilkinson_1d.py\n";

    return 0;
}

