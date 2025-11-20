/**
 * Example: 1D Cubic Polynomial with 3 Roots
 * 
 * This example demonstrates solving a 1D cubic polynomial with 3 known roots.
 * 
 * Problem:
 *   p(x) = (x - 0.2)(x - 0.5)(x - 0.8)
 *        = x^3 - 1.5*x^2 + 0.62*x - 0.08
 * 
 * Expected roots: x = 0.2, 0.5, 0.8
 * 
 * This example tests all three subdivision strategies and generates
 * geometry dumps for 1D visualization.
 */

#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

void print_result(const char* strategy_name, const SubdivisionSolverResult& result) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Strategy: " << strategy_name << "\n";
    std::cout << std::string(60, '=') << "\n";
    
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
    std::cout << "1D Cubic Polynomial Root Finding\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    std::cout << "Problem: p(x) = (x - 0.2)(x - 0.5)(x - 0.8)\n";
    std::cout << "              = x^3 - 1.5*x^2 + 0.66*x - 0.08\n";
    std::cout << "Expected roots: x = 0.2, 0.5, 0.8\n";
    std::cout << "Domain: [0, 1]\n\n";

    // Define p(x) = x^3 - 1.5*x^2 + 0.66*x - 0.08
    // Coefficients in power basis: [c0, c1, c2, c3] for c0 + c1*x + c2*x^2 + c3*x^3
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{-0.08, 0.66, -1.5, 1.0};
    
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    
    // Verify the roots
    std::cout << "Verification:\n";
    double roots[] = {0.2, 0.5, 0.8};
    for (int i = 0; i < 3; ++i) {
        double val = p.evaluate(roots[i]);
        std::cout << "  p(" << roots[i] << ") = " << std::scientific 
                  << std::setprecision(6) << val << "\n";
    }
    std::cout << "\n";
    
    // Test all three strategies
    const char* strategy_names[] = {"ContractFirst", "SubdivideFirst", "Simultaneous"};
    SubdivisionStrategy strategies[] = {
        SubdivisionStrategy::ContractFirst,
        SubdivisionStrategy::SubdivideFirst,
        SubdivisionStrategy::Simultaneous
    };
    
    Solver solver;
    
    for (int s = 0; s < 3; ++s) {
        SubdivisionConfig config;
        config.tolerance = 1e-8;
        config.max_depth = 100;
        config.contraction_threshold = 0.9;
        config.strategy = strategies[s];
        config.dump_geometry = true;
        config.dump_prefix = std::string("dumps/cubic_1d_") + strategy_names[s];
        
        SubdivisionSolverResult result = solver.subdivisionSolve(
            system, config, RootBoundingMethod::ProjectedPolyhedral);
        
        print_result(strategy_names[s], result);
    }
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "\nGeometry dumps saved to dumps/cubic_1d_*.txt\n";
    std::cout << "Visualize with: python examples/visualize_cubic_1d.py\n";
    
    return 0;
}

