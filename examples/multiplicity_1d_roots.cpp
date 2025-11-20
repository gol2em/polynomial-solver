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
#include <cmath>

using namespace polynomial_solver;

void print_result(const SubdivisionSolverResult& result) {
    std::cout << "Found " << result.num_resolved << " root(s)\n";
    std::cout << "Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n\n";

    for (size_t i = 0; i < result.num_resolved; ++i) {
        const SubdivisionBoxResult& box = result.boxes[i];
        std::cout << "Root " << (i+1) << ":\n";
        std::cout << "  Interval: [" << std::fixed << std::setprecision(10)
                  << box.lower[0] << ", " << box.upper[0] << "]\n";
        std::cout << "  Center: " << box.center[0] << "\n";
        std::cout << "  Width: " << std::scientific << std::setprecision(6)
                  << (box.upper[0] - box.lower[0]) << "\n";
        std::cout << "  Depth: " << box.depth << "\n";
        std::cout << "  Converged: " << (box.converged ? "yes" : "no") << "\n";
        std::cout << "\n";
    }
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
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "\nGeometry dump saved to dumps/multiplicity_1d_geometry.txt\n";
    std::cout << "Visualize with: python examples/visualize_multiplicity_1d.py\n";
    
    return 0;
}

