/**
 * Circle Test - Finding the zero set of x² + y² - 1 = 0
 * 
 * This is a simple test case that demonstrates finding the unit circle
 * using the Projected Polyhedral method. The solver finds boxes that
 * contain parts of the circle (the zero set of the polynomial).
 * 
 * Expected result: Boxes distributed around the unit circle
 */

#include <polynomial_solver.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

int main(int argc, char* argv[]) {
    std::cout << "========================================" << std::endl;
    std::cout << "Circle Test: x² + y² - 1 = 0" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Parse command line arguments
    double degeneracy_multiplier = 5.0;
    if (argc > 1) {
        degeneracy_multiplier = std::atof(argv[1]);
    }
    
    std::cout << "Degeneracy multiplier: " << degeneracy_multiplier << "\n" << std::endl;
    
    // Define the polynomial: x² + y² - 1
    // In power basis: f(x,y) = -1 + 0*y + 1*y^2 + 0*x + 0*xy + 0*xy^2 + 1*x^2 + 0*x^2y + 0*x^2y^2
    std::vector<unsigned int> degrees = {2, 2};
    std::vector<double> power_coeffs = {
        -1.0,  0.0,  1.0,   // 1, y, y^2      (x^0)
         0.0,  0.0,  0.0,   // x, xy, xy^2    (x^1)
         1.0,  0.0,  0.0    // x^2, x^2y, x^2y^2  (x^2)
    };

    Polynomial circle_poly = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({circle_poly});
    
    std::cout << "Polynomial: x² + y² - 1" << std::endl;
    std::cout << "  Degree: (" << circle_poly.degrees()[0] << ", " 
              << circle_poly.degrees()[1] << ")" << std::endl;
    std::cout << "  Coefficients: " << circle_poly.coefficientCount() << "\n" << std::endl;
    
    // Solve for zero set
    std::cout << "Finding zero set using PP method..." << std::endl;
    
    SubdivisionConfig config = defaultSolverConfig();
    config.tolerance = 1e-8;
    config.max_depth = 100;
    config.degeneracy_multiplier = degeneracy_multiplier;

#ifdef ENABLE_GEOMETRY_DUMP
    config.dump_geometry = true;
    config.dump_prefix = "dumps/circle";
#endif
    
    Solver solver;
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "\nResults:" << std::endl;
    std::cout << "  Total boxes: " << result.boxes.size() << std::endl;
    std::cout << "  Resolved (converged): " << result.num_resolved << std::endl;
    std::cout << "  Unresolved: " << (result.boxes.size() - result.num_resolved) << std::endl;
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    
    // Write results to file
    std::ofstream result_file("dumps/circle_result_boxes.txt");
    result_file << "# Circle test results: x² + y² - 1 = 0\n";
    result_file << "# Format: x_min x_max y_min y_max converged\n";
    result_file << std::setprecision(17);

    for (std::size_t i = 0; i < result.boxes.size(); ++i) {
        const auto& box = result.boxes[i];
        bool converged = (i < result.num_resolved);
        result_file << box.lower[0] << " " << box.upper[0] << " "
                   << box.lower[1] << " " << box.upper[1] << " "
                   << (converged ? "1" : "0") << "\n";
    }
    result_file.close();
    
    std::cout << "\n✅ Results written to: dumps/circle_result_boxes.txt" << std::endl;
    
    // Compute statistics
    double mean_radius = 0.0;
    double mean_error = 0.0;
    int count = 0;
    
    for (std::size_t i = result.num_resolved; i < result.boxes.size(); ++i) {
        const auto& box = result.boxes[i];
        double x = box.center[0];
        double y = box.center[1];
        double radius = std::sqrt(x * x + y * y);
        mean_radius += radius;
        mean_error += std::abs(radius - 1.0);
        count++;
    }
    
    if (count > 0) {
        mean_radius /= count;
        mean_error /= count;
        
        std::cout << "\nStatistics for unresolved boxes:" << std::endl;
        std::cout << "  Mean radius: " << std::fixed << std::setprecision(6) << mean_radius << std::endl;
        std::cout << "  Mean error from unit circle: " << mean_error << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "To visualize: python3 visualize_circle_boxes_2d.py" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}

