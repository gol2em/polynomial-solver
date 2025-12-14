/**
 * @file test_multiplicity_debug.cpp
 * @brief Debug multiplicity detection for high multiplicities
 */

#ifdef ENABLE_HIGH_PRECISION

#include "result_refiner_hp.h"
#include "polynomial_hp.h"
#include "polynomial.h"
#include "differentiation_hp.h"
#include "high_precision_types.h"
#include "precision_context.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace polynomial_solver;

PolynomialHP createMultipleRootPolynomial(double root, unsigned int multiplicity) {
    std::vector<double> power_coeffs = {-root, 1.0};
    for (unsigned int i = 1; i < multiplicity; ++i) {
        std::vector<double> new_coeffs(power_coeffs.size() + 1, 0.0);
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-root);
            new_coeffs[j+1] += power_coeffs[j];
        }
        power_coeffs = new_coeffs;
    }
    std::vector<unsigned int> degrees = {multiplicity};
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    return PolynomialHP(poly);
}

void debugTaylorMethod(unsigned int true_m, double root, double guess) {
    std::cout << "\n========================================\n";
    std::cout << "Debug Taylor Method for m=" << true_m << "\n";
    std::cout << "========================================\n";
    
    PolynomialHP poly = createMultipleRootPolynomial(root, true_m);
    mpreal x = mpreal(guess);
    
    // Compute derivatives
    std::cout << "\nDerivative values at x=" << x << ":\n";
    std::cout << "Order | Value          | Abs Value      | Ratio to next\n";
    std::cout << "------+----------------+----------------+--------------\n";
    
    std::vector<mpreal> deriv_vals(11);
    for (unsigned int k = 1; k <= 10; ++k) {
        PolynomialHP deriv = DifferentiationHP::derivative(poly, 0, k);
        deriv_vals[k] = deriv.evaluate(x);
    }
    
    for (unsigned int k = 1; k <= 10; ++k) {
        mpreal ratio = mpreal(0);
        if (k < 10 && abs(deriv_vals[k]) > mpreal("1e-100")) {
            ratio = abs(deriv_vals[k+1]) / abs(deriv_vals[k]);
        }
        
        std::cout << std::setw(5) << k << " | ";
        std::cout << std::scientific << std::setprecision(6) << deriv_vals[k] << " | ";
        std::cout << std::scientific << std::setprecision(6) << abs(deriv_vals[k]) << " | ";
        if (k < 10) {
            std::cout << std::scientific << std::setprecision(2) << ratio;
            if (ratio > mpreal(100)) std::cout << " (LARGE)";
        }
        std::cout << "\n";
    }
    
    // Test Taylor method
    mpreal first_nonzero;
    unsigned int m_est = ResultRefinerHP::estimateMultiplicity(x, poly, 10, mpreal("1e-50"), first_nonzero);
    std::cout << "\nTaylor estimate: m=" << m_est << " (true m=" << true_m << ")\n";
    std::cout << "First nonzero deriv: " << first_nonzero << "\n";
}

void debugOstrowskiMethod(unsigned int true_m, double root, double guess) {
    std::cout << "\n========================================\n";
    std::cout << "Debug Ostrowski Method for m=" << true_m << "\n";
    std::cout << "========================================\n";
    
    PolynomialHP poly = createMultipleRootPolynomial(root, true_m);
    PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
    
    mpreal x0 = mpreal(guess);
    mpreal true_root = mpreal(root);
    
    std::cout << "\nRegular Newton iterations:\n";
    std::cout << "Iter | x              | Error          | f(x)           | f'(x)\n";
    std::cout << "-----+----------------+----------------+----------------+----------------\n";
    
    std::vector<mpreal> iterates;
    iterates.push_back(x0);
    
    mpreal x = x0;
    for (int i = 0; i < 5; ++i) {
        mpreal f = poly.evaluate(x);
        mpreal df = dpoly.evaluate(x);
        mpreal error = abs(x - true_root);
        
        std::cout << std::setw(4) << i << " | ";
        std::cout << std::scientific << std::setprecision(6) << x << " | ";
        std::cout << std::scientific << std::setprecision(6) << error << " | ";
        std::cout << std::scientific << std::setprecision(6) << f << " | ";
        std::cout << std::scientific << std::setprecision(6) << df << "\n";
        
        if (abs(df) < mpreal("1e-100")) break;
        x = x - f / df;
        iterates.push_back(x);
    }
    
    // Compute error ratios
    std::cout << "\nError ratios (should be (m-1)/m = " << (true_m-1.0)/true_m << "):\n";
    for (size_t i = 1; i < iterates.size(); ++i) {
        mpreal e_prev = abs(iterates[i-1] - true_root);
        mpreal e_curr = abs(iterates[i] - true_root);
        if (e_prev > mpreal("1e-100")) {
            mpreal ratio = e_curr / e_prev;
            std::cout << "  e[" << i << "]/e[" << i-1 << "] = " << ratio << "\n";
        }
    }
    
    // Test Ostrowski formula
    if (iterates.size() >= 4) {
        mpreal x1 = iterates[1];
        mpreal x2 = iterates[2];
        mpreal x3 = iterates[3];
        
        mpreal num = x1 - x2;
        mpreal denom = x3 - mpreal(2)*x2 + x1;
        mpreal p = mpreal("0.5") + num / denom;
        
        std::cout << "\nOstrowski formula:\n";
        std::cout << "  numerator = " << num << "\n";
        std::cout << "  denominator = " << denom << "\n";
        std::cout << "  p = 0.5 + num/denom = " << p << "\n";
        std::cout << "  round(p) = " << static_cast<int>(round(p)) << "\n";
        std::cout << "  floor(p) = " << static_cast<int>(floor(p)) << "\n";
        std::cout << "  ceil(p) = " << static_cast<int>(ceil(p)) << "\n";
        
        unsigned int m_est = ResultRefinerHP::estimateMultiplicityOstrowski(x1, x2, x3);
        std::cout << "\nOstrowski estimate: m=" << m_est << " (true m=" << true_m << ")\n";
    }
}

int main() {
    PrecisionContext ctx(256);
    
    std::cout << "========================================\n";
    std::cout << "  Multiplicity Detection Debug\n";
    std::cout << "========================================\n";

    // Test multiplicities 2-12 to find limits
    std::cout << "\nSummary Table:\n";
    std::cout << "True m | Taylor | Ostrowski | Status\n";
    std::cout << "-------+--------+-----------+--------\n";

    for (unsigned int m = 2; m <= 12; ++m) {
        PolynomialHP poly = createMultipleRootPolynomial(0.5, m);
        mpreal x = mpreal(0.48);

        // Taylor estimate
        mpreal first_nonzero;
        unsigned int taylor_est = ResultRefinerHP::estimateMultiplicity(x, poly, 15, mpreal("1e-50"), first_nonzero);

        // Ostrowski estimate
        unsigned int ostro_est = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);

        std::cout << std::setw(6) << m << " | ";
        std::cout << std::setw(6) << taylor_est << " | ";
        std::cout << std::setw(9) << ostro_est << " | ";

        if (taylor_est == m && ostro_est == m) {
            std::cout << "✓ BOTH CORRECT";
        } else if (taylor_est == m) {
            std::cout << "Taylor OK, Ostro=" << ostro_est;
        } else if (ostro_est == m) {
            std::cout << "Ostro OK, Taylor=" << taylor_est;
        } else {
            std::cout << "BOTH WRONG";
        }
        std::cout << "\n";
    }

    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

