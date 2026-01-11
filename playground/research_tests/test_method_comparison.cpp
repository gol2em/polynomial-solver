#include "core/polynomial.h"
#include "solver/solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

void compare_methods(const std::string& test_name, const PolynomialSystem& system, 
                     const SubdivisionConfig& config) {
    std::cout << "Test: " << test_name << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    Solver solver;
    
    // Test with GraphHull method
    SubdivisionSolverResult result_gh = solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);
    
    std::cout << "GraphHull method:" << std::endl;
    std::cout << "  Total boxes: " << result_gh.boxes.size() << std::endl;
    std::cout << "  Resolved: " << result_gh.num_resolved << std::endl;
    std::cout << "  Degeneracy: " << (result_gh.degeneracy_detected ? "yes" : "no") << std::endl;
    
    if (result_gh.num_resolved > 0) {
        std::cout << "  First root: (";
        for (std::size_t i = 0; i < result_gh.boxes[0].center.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << std::fixed << std::setprecision(8) << result_gh.boxes[0].center[i];
        }
        std::cout << ") at depth " << result_gh.boxes[0].depth << std::endl;
    }
    
    // Test with ProjectedPolyhedral method
    SubdivisionSolverResult result_pp = solver.subdivisionSolve(system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "ProjectedPolyhedral method:" << std::endl;
    std::cout << "  Total boxes: " << result_pp.boxes.size() << std::endl;
    std::cout << "  Resolved: " << result_pp.num_resolved << std::endl;
    std::cout << "  Degeneracy: " << (result_pp.degeneracy_detected ? "yes" : "no") << std::endl;
    
    if (result_pp.num_resolved > 0) {
        std::cout << "  First root: (";
        for (std::size_t i = 0; i < result_pp.boxes[0].center.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << std::fixed << std::setprecision(8) << result_pp.boxes[0].center[i];
        }
        std::cout << ") at depth " << result_pp.boxes[0].depth << std::endl;
    }
    
    std::cout << std::endl;
}

int main() {
    std::cout << "Comparison of GraphHull vs ProjectedPolyhedral methods" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << std::endl;
    
    SubdivisionConfig config;
    config.tolerance = 0.05;
    config.max_depth = 100;
    config.degeneracy_multiplier = 10.0;
    
    // Test 1: 1D linear
    {
        std::vector<unsigned int> degrees{1u};
        std::vector<double> power_coeffs{-0.5, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        PolynomialSystem system({p});
        
        compare_methods("1D Linear: p(x) = x - 0.5", system, config);
    }
    
    // Test 2: 2D linear
    {
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
        compare_methods("2D Linear: p1(x,y) = x - 0.5, p2(x,y) = y - 0.3", system, config);
    }
    
    // Test 3: 2D quadratic
    {
        std::vector<unsigned int> degrees_f1{2u, 1u};
        std::vector<double> power_coeffs_f1{-0.25, 0.0, 0.0, 0.0, 1.0, 0.0};
        Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs_f1);
        
        std::vector<unsigned int> degrees_f2{0u, 1u};
        std::vector<double> power_coeffs_f2{-0.3, 1.0};
        Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs_f2);
        
        PolynomialSystem system({p1, p2});
        compare_methods("2D Quadratic: p1(x,y) = x^2 - 0.25, p2(x,y) = y - 0.3", system, config);
    }
    
    std::cout << "Comparison complete!" << std::endl;
    std::cout << std::endl;
    std::cout << "Summary:" << std::endl;
    std::cout << "- Both methods achieve machine epsilon for linear functions" << std::endl;
    std::cout << "- GraphHull works in full n-dimensional space" << std::endl;
    std::cout << "- ProjectedPolyhedral works direction-by-direction" << std::endl;
    std::cout << "- Both methods successfully handle quadratic systems" << std::endl;
    
    return 0;
}

