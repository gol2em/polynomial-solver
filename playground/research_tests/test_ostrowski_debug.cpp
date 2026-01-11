#include <iostream>
#include <iomanip>
#include "core/polynomial.h"
#include "hp/polynomial_hp.h"
#include "hp/differentiation_hp.h"
#include "hp/precision_context.h"

using namespace polynomial_solver;

int main() {
    // Set precision
    PrecisionContext ctx(256);
    
    // Test case: (x - 0.5)^3 - triple root
    // Create in power form: x^3 - 1.5x^2 + 0.75x - 0.125
    std::vector<double> power_coeffs = {-0.125, 0.75, -1.5, 1.0};
    std::vector<unsigned int> degrees = {3};
    Polynomial poly_double = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialHP poly(poly_double);
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    
    std::cout << std::setprecision(15) << std::scientific;
    std::cout << "Testing Ostrowski formula for triple root at x=0.5\n";
    std::cout << "=================================================\n\n";
    
    // Perform 3 regular Newton steps from x=0.48
    mpreal x0 = mpreal("0.48");
    mpreal x1, x2, x3;
    
    // Step 1
    {
        mpreal f0 = poly.evaluate(x0);
        mpreal df0 = dpoly.evaluate(x0);
        x1 = x0 - f0 / df0;
        std::cout << "x0 = " << x0 << std::endl;
        std::cout << "f(x0) = " << f0 << ", f'(x0) = " << df0 << std::endl;
        std::cout << "x1 = " << x1 << std::endl;
        std::cout << "step size = " << (x1 - x0) << "\n\n";
    }
    
    // Step 2
    {
        mpreal f1 = poly.evaluate(x1);
        mpreal df1 = dpoly.evaluate(x1);
        x2 = x1 - f1 / df1;
        std::cout << "x1 = " << x1 << std::endl;
        std::cout << "f(x1) = " << f1 << ", f'(x1) = " << df1 << std::endl;
        std::cout << "x2 = " << x2 << std::endl;
        std::cout << "step size = " << (x2 - x1) << "\n\n";
    }
    
    // Step 3
    {
        mpreal f2 = poly.evaluate(x2);
        mpreal df2 = dpoly.evaluate(x2);
        x3 = x2 - f2 / df2;
        std::cout << "x2 = " << x2 << std::endl;
        std::cout << "f(x2) = " << f2 << ", f'(x2) = " << df2 << std::endl;
        std::cout << "x3 = " << x3 << std::endl;
        std::cout << "step size = " << (x3 - x2) << "\n\n";
    }
    
    // Compute error ratios
    mpreal e1 = x1 - x0;  // error at step 1 (approximately)
    mpreal e2 = x2 - x1;  // error at step 2
    mpreal e3 = x3 - x2;  // error at step 3
    
    std::cout << "Error analysis:\n";
    std::cout << "e1 = x1 - x0 = " << e1 << std::endl;
    std::cout << "e2 = x2 - x1 = " << e2 << std::endl;
    std::cout << "e3 = x3 - x2 = " << e3 << std::endl;
    std::cout << "e2/e1 = " << (e2/e1) << " (should be ~(m-1)/m = 2/3 for m=3)" << std::endl;
    std::cout << "e3/e2 = " << (e3/e2) << " (should be ~(m-1)/m = 2/3 for m=3)" << std::endl;
    std::cout << std::endl;
    
    // Apply Ostrowski formula
    mpreal numerator = x1 - x2;
    mpreal denominator = x3 - mpreal("2.0") * x2 + x1;
    mpreal p_est = mpreal("0.5") + numerator / denominator;
    
    std::cout << "Ostrowski formula:\n";
    std::cout << "numerator = x1 - x2 = " << numerator << std::endl;
    std::cout << "denominator = x3 - 2*x2 + x1 = " << denominator << std::endl;
    std::cout << "p = 1/2 + numerator/denominator = " << p_est << std::endl;
    std::cout << "ceil(p) = " << static_cast<int>(ceil(p_est)) << std::endl;
    std::cout << "floor(p) = " << static_cast<int>(floor(p_est)) << std::endl;
    std::cout << "round(p) = " << static_cast<int>(round(p_est)) << std::endl;
    std::cout << std::endl;
    
    // Try alternative formula: m = 1 / (1 - ratio)
    mpreal ratio = e2 / e1;
    mpreal m_from_ratio = mpreal("1.0") / (mpreal("1.0") - ratio);
    std::cout << "Alternative: m = 1/(1 - e2/e1) = " << m_from_ratio << std::endl;
    std::cout << "ceil(m) = " << static_cast<int>(ceil(m_from_ratio)) << std::endl;
    std::cout << "round(m) = " << static_cast<int>(round(m_from_ratio)) << std::endl;
    
    return 0;
}

