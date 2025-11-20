#include "differentiation.h"
#include "polynomial.h"
#include <iostream>
#include <vector>

/**
 * @file test_differentiation_debug.cpp
 * @brief Debug test for 2D polynomial differentiation
 */

using namespace polynomial_solver;

int main()
{
    std::cout << "Debug: Testing f(x,y) = x^2 + y^2" << std::endl;
    
    // For degrees (2,2), we have 3x3 = 9 coefficients
    // Layout: last dimension (y) varies fastest
    // Index order: (0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)
    // Power basis: x^i * y^j
    
    std::vector<unsigned int> degrees{2u, 2u};
    std::vector<double> power_coeffs(9, 0.0);
    
    // x^2 + y^2 means:
    // - coefficient of x^2*y^0 = 1.0  -> index (2,0) = 6
    // - coefficient of x^0*y^2 = 1.0  -> index (0,2) = 2
    
    power_coeffs[6] = 1.0;  // x^2 * y^0
    power_coeffs[2] = 1.0;  // x^0 * y^2
    
    std::cout << "Power coefficients:" << std::endl;
    for (size_t i = 0; i < 9; ++i) {
        size_t ix = i / 3;
        size_t iy = i % 3;
        std::cout << "  [" << i << "] (x^" << ix << "*y^" << iy << ") = " << power_coeffs[i] << std::endl;
    }
    
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    
    std::cout << "\nBernstein coefficients:" << std::endl;
    const std::vector<double>& bern_coeffs = p.bernsteinCoefficients();
    for (size_t i = 0; i < bern_coeffs.size(); ++i) {
        size_t ix = i / 3;
        size_t iy = i % 3;
        std::cout << "  [" << i << "] (i=" << ix << ",j=" << iy << ") = " << bern_coeffs[i] << std::endl;
    }
    
    // Test evaluation
    std::cout << "\nEvaluating f(x,y):" << std::endl;
    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0}, {0.5, 0.5}, {1.0, 1.0}
    };
    
    for (const auto& pt : test_points) {
        double x = pt[0], y = pt[1];
        double expected = x*x + y*y;
        double actual = p.evaluate(pt);
        std::cout << "  f(" << x << "," << y << ") = " << actual 
                  << " (expected " << expected << ")" << std::endl;
    }
    
    // Test gradient
    std::cout << "\nGradient:" << std::endl;
    std::vector<Polynomial> grad = Differentiation::gradient(p);
    
    std::cout << "df/dx Bernstein coefficients:" << std::endl;
    const std::vector<double>& df_dx_coeffs = grad[0].bernsteinCoefficients();
    for (size_t i = 0; i < df_dx_coeffs.size(); ++i) {
        std::cout << "  [" << i << "] = " << df_dx_coeffs[i] << std::endl;
    }
    
    std::cout << "df/dy Bernstein coefficients:" << std::endl;
    const std::vector<double>& df_dy_coeffs = grad[1].bernsteinCoefficients();
    for (size_t i = 0; i < df_dy_coeffs.size(); ++i) {
        std::cout << "  [" << i << "] = " << df_dy_coeffs[i] << std::endl;
    }
    
    std::cout << "\nEvaluating gradient:" << std::endl;
    for (const auto& pt : test_points) {
        double x = pt[0], y = pt[1];
        double expected_df_dx = 2.0 * x;
        double expected_df_dy = 2.0 * y;
        double actual_df_dx = grad[0].evaluate(pt);
        double actual_df_dy = grad[1].evaluate(pt);
        
        std::cout << "  grad f(" << x << "," << y << ") = (" 
                  << actual_df_dx << "," << actual_df_dy << ")"
                  << " (expected (" << expected_df_dx << "," << expected_df_dy << "))" << std::endl;
    }
    
    return 0;
}

