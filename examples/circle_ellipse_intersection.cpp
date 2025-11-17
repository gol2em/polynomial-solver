/**
 * Example: Circle-Ellipse Intersection
 * 
 * This example demonstrates solving a system of polynomial equations
 * representing the intersection of a circle and an ellipse.
 * 
 * Problem:
 *   f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)
 *   f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)
 * 
 * Expected root in [0,1]^2: approximately (0.894, 0.447)
 * 
 * This example tests all three subdivision strategies:
 * - ContractFirst: Contract bounds before subdividing
 * - SubdivideFirst: Subdivide before contracting
 * - Simultaneous: Balance both approaches
 */

#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

void print_result(const std::string& strategy_name, 
                  const SubdivisionSolverResult& result,
                  const PolynomialSystem& system) {
    std::cout << "\n" << strategy_name << " Strategy:\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::cout << "  Found " << result.boxes.size() << " box(es), "
              << result.num_resolved << " resolved\n";
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";
    
    if (result.num_resolved > 0) {
        const SubdivisionBoxResult& box = result.boxes[0];
        std::cout << "  Root:\n";
        std::cout << "    Center: (" << std::setprecision(10) << box.center[0]
                  << ", " << box.center[1] << ")\n";
        std::cout << "    Max error: (" << std::scientific << std::setprecision(3)
                  << box.max_error[0] << ", " << box.max_error[1] << ")\n";
        std::cout << "    Depth: " << box.depth << "\n";
        
        // Evaluate at center to verify
        std::vector<double> eval;
        system.evaluate(box.center, eval);
        std::cout << "    Residual: (" << eval[0] << ", " << eval[1] << ")\n";
        
        // Calculate actual distance from expected root
        double expected_x = 2.0 / std::sqrt(5.0);  // ≈ 0.894427191
        double expected_y = 1.0 / std::sqrt(5.0);  // ≈ 0.447213595
        double error_x = std::abs(box.center[0] - expected_x);
        double error_y = std::abs(box.center[1] - expected_y);
        std::cout << "    Error from expected: (" << error_x << ", " << error_y << ")\n";
    }
}

int main() {
    std::cout << "Circle-Ellipse Intersection Example\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    std::cout << "Problem:\n";
    std::cout << "  f1(x,y) = x^2 + y^2 - 1 = 0        (unit circle)\n";
    std::cout << "  f2(x,y) = x^2/4 + 4*y^2 - 1 = 0    (ellipse)\n";
    std::cout << "  Domain: [0, 1] × [0, 1]\n\n";
    
    // Define f1(x,y) = x^2 + y^2 - 1
    std::vector<unsigned int> degrees_f1{2u, 2u};
    std::vector<double> power_coeffs1{
        -1.0,  0.0,  1.0,   // 1, y, y^2
         0.0,  0.0,  0.0,   // x, xy, xy^2
         1.0,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs1);
    
    // Define f2(x,y) = x^2/4 + 4*y^2 - 1
    std::vector<unsigned int> degrees_f2{2u, 2u};
    std::vector<double> power_coeffs2{
        -1.0,   0.0,  4.0,   // 1, y, y^2
         0.0,   0.0,  0.0,   // x, xy, xy^2
         0.25,  0.0,  0.0    // x^2, x^2y, x^2y^2
    };
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs2);
    
    PolynomialSystem system({p1, p2});
    
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
        config.tolerance = 1e-6;
        config.max_depth = 100;
        config.contraction_threshold = 0.9;
        config.strategy = strategies[s];
        config.dump_geometry = true;
        config.dump_prefix = std::string("dumps/example_") + strategy_names[s];
        
        SubdivisionSolverResult result = solver.subdivisionSolve(
            system, config, RootBoundingMethod::ProjectedPolyhedral);
        
        print_result(strategy_names[s], result, system);
    }
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "\nGeometry dumps saved to dumps/example_*.txt\n";
    std::cout << "Visualize with: tools/visualize_solver.py\n";
    
    return 0;
}

