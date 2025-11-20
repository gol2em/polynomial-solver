/**
 * Example: 1D Polynomial with High Multiplicity Root
 *
 * This example demonstrates solving a 1D polynomial with a high multiplicity root.
 *
 * Problem:
 *   p(x) = (x - 0.2)(x - 0.6)^6
 *
 * Expected roots: x = 0.2 (multiplicity 1), x = 0.6 (multiplicity 6)
 *
 * This tests the solver's ability to handle roots with high multiplicity,
 * which can be challenging due to the flat tangency at the root.
 */

#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

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

int main() {
    std::cout << "1D Polynomial with High Multiplicity Root\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    std::cout << "Problem: p(x) = (x - 0.2)(x - 0.6)^6\n";
    std::cout << "Expected roots: x = 0.2 (multiplicity 1), x = 0.6 (multiplicity 6)\n";
    std::cout << "Domain: [0, 1]\n\n";

    // Expand (x - 0.2)(x - 0.6)^6
    // Computed using: (x - 0.6)^6 * (x - 0.2)
    // p(x) = x^7 - 3.8*x^6 + 6.12*x^5 - 5.4*x^4 + 2.808*x^3 - 0.85536*x^2 + 0.139968*x - 0.0093312

    std::vector<unsigned int> degrees{7u};
    std::vector<double> power_coeffs{
        -0.0093312,   // x^0
        0.139968,     // x^1
        -0.85536,     // x^2
        2.808,        // x^3
        -5.4,         // x^4
        6.12,         // x^5
        -3.8,         // x^6
        1.0           // x^7
    };
    
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    // Verify the roots
    std::cout << "Verification:\n";
    double roots[] = {0.2, 0.6};
    for (int i = 0; i < 2; ++i) {
        double val = p.evaluate(roots[i]);
        std::cout << "  p(" << roots[i] << ") = " << std::scientific 
                  << std::setprecision(6) << val << "\n";
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
    config.dump_prefix = "dumps/multiplicity_1d";

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Solving with SubdivideFirst strategy\n";
    std::cout << std::string(60, '=') << "\n\n";

    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);

    print_result(result);

    // Dump result to file
    dump_result(result, "dumps/multiplicity_1d_result.txt");

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "\nGeometry dump saved to dumps/multiplicity_1d_geometry.txt\n";
    std::cout << "Result dump saved to dumps/multiplicity_1d_result.txt\n";
    std::cout << "Visualize with: python examples/visualize_multiplicity_1d.py\n";

    return 0;
}

