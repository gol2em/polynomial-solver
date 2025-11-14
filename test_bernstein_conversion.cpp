#include <iostream>
#include <vector>

// Simple test to verify Bernstein conversion
// For p(x) = x^2 - 0.25, the Bernstein coefficients should be:
// b_0 = p(0) = -0.25
// b_1 = ? (not p(0.5)!)
// b_2 = p(1) = 0.75

// The Bernstein basis for degree 2:
// B_0,2(x) = (1-x)^2
// B_1,2(x) = 2x(1-x)
// B_2,2(x) = x^2

// So p(x) = b_0 * (1-x)^2 + b_1 * 2x(1-x) + b_2 * x^2
//         = b_0 * (1 - 2x + x^2) + b_1 * (2x - 2x^2) + b_2 * x^2
//         = b_0 + (2b_1 - 2b_0)x + (b_0 - 2b_1 + b_2)x^2

// Comparing with p(x) = -0.25 + 0x + 1x^2:
// b_0 = -0.25
// 2b_1 - 2b_0 = 0  =>  b_1 = b_0 = -0.25
// b_0 - 2b_1 + b_2 = 1  =>  -0.25 - 2(-0.25) + b_2 = 1  =>  -0.25 + 0.5 + b_2 = 1  =>  b_2 = 0.75

int main() {
    std::cout << "For p(x) = x^2 - 0.25:" << std::endl;
    std::cout << "Bernstein coefficients should be:" << std::endl;
    std::cout << "b_0 = -0.25" << std::endl;
    std::cout << "b_1 = -0.25" << std::endl;
    std::cout << "b_2 = 0.75" << std::endl;
    
    // Verify by evaluating at x=0.5
    double b0 = -0.25, b1 = -0.25, b2 = 0.75;
    double x = 0.5;
    double B0 = (1-x)*(1-x);
    double B1 = 2*x*(1-x);
    double B2 = x*x;
    double p_bernstein = b0*B0 + b1*B1 + b2*B2;
    double p_power = x*x - 0.25;
    
    std::cout << "\nVerification at x=0.5:" << std::endl;
    std::cout << "p(0.5) using Bernstein = " << p_bernstein << std::endl;
    std::cout << "p(0.5) using power = " << p_power << std::endl;
    std::cout << "Match: " << (std::abs(p_bernstein - p_power) < 1e-10 ? "YES" : "NO") << std::endl;
    
    return 0;
}

