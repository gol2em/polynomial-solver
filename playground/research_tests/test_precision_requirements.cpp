/**
 * @file test_precision_requirements.cpp
 * @brief Systematically determine precision requirements for multiplicity detection and convergence
 * 
 * This test explores:
 * 1. Minimum precision needed for correct multiplicity detection
 * 2. Minimum precision needed for modified Newton convergence (given correct multiplicity)
 * 3. Optimal workflow design based on these requirements
 */

#ifdef ENABLE_HIGH_PRECISION

#include "hp/result_refiner_hp.h"
#include "hp/polynomial_hp.h"
#include "core/polynomial.h"
#include "hp/differentiation_hp.h"
#include "hp/precision_context.h"
#include "hp/precision_conversion.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace polynomial_solver;

// Create polynomial in HIGH PRECISION
PolynomialHP createMultipleRootPolynomialHP(const mpreal& root, unsigned int multiplicity) {
    std::vector<mpreal> power_coeffs = {-root, mpreal(1)};
    for (unsigned int i = 1; i < multiplicity; ++i) {
        std::vector<mpreal> new_coeffs(power_coeffs.size() + 1, mpreal(0));
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-root);
            new_coeffs[j+1] += power_coeffs[j];
        }
        power_coeffs = new_coeffs;
    }
    std::vector<unsigned int> degrees = {multiplicity};
    return fromPowerHP(degrees, power_coeffs);
}

// Test multiplicity detection at different precisions
void testMultiplicityDetection(unsigned int true_mult) {
    std::cout << "\n=== Multiplicity Detection: m=" << true_mult << " ===\n";
    std::cout << "Precision | Taylor | Ostrowski | Sturm | Status\n";
    std::cout << "----------+--------+-----------+-------+--------\n";
    
    std::vector<unsigned int> precisions = {64, 96, 128, 192, 256, 384, 512, 768, 1024};
    unsigned int min_working_precision = 0;
    
    for (unsigned int prec : precisions) {
        PrecisionContext ctx(prec);
        
        PolynomialHP poly = createMultipleRootPolynomialHP(mpreal("0.5"), true_mult);
        mpreal x = mpreal("0.48");  // Close to root
        
        // Test Taylor method
        mpreal dummy;
        unsigned int taylor = ResultRefinerHP::estimateMultiplicity(
            x, poly, 15, mpreal("1e-50"), dummy, 10.0);
        
        // Test Ostrowski method
        unsigned int ostro = ResultRefinerHP::estimateMultiplicityOstrowskiFromPoint(x, poly);
        
        // Test Sturm method
        unsigned int sturm = ResultRefinerHP::estimateMultiplicitySturm(x, poly, mpreal("1e-10"));
        
        bool taylor_correct = (taylor == true_mult);
        bool ostro_correct = (ostro == true_mult);
        bool sturm_correct = (sturm == true_mult);
        bool all_correct = taylor_correct && ostro_correct && sturm_correct;
        
        std::cout << std::setw(9) << prec << " | ";
        std::cout << std::setw(6) << taylor << " | ";
        std::cout << std::setw(9) << ostro << " | ";
        std::cout << std::setw(5) << sturm << " | ";
        
        if (all_correct) {
            std::cout << "✓ ALL OK";
            if (min_working_precision == 0) {
                min_working_precision = prec;
            }
        } else {
            std::cout << "✗ ";
            if (!taylor_correct) std::cout << "T";
            if (!ostro_correct) std::cout << "O";
            if (!sturm_correct) std::cout << "S";
        }
        std::cout << "\n";
    }
    
    if (min_working_precision > 0) {
        std::cout << "→ Minimum precision for detection: " << min_working_precision << " bits\n";
    } else {
        std::cout << "→ Detection failed at all tested precisions!\n";
    }
}

// Test convergence with known multiplicity at different precisions
void testConvergenceWithKnownMultiplicity(unsigned int true_mult) {
    std::cout << "\n=== Convergence Test: m=" << true_mult << " (multiplicity known) ===\n";
    std::cout << "Precision | Iters | Final Error | Converged | Status\n";
    std::cout << "----------+-------+-------------+-----------+--------\n";
    
    std::vector<unsigned int> precisions = {64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048};
    unsigned int min_working_precision = 0;
    
    for (unsigned int prec : precisions) {
        PrecisionContext ctx(prec);
        
        PolynomialHP poly = createMultipleRootPolynomialHP(mpreal("0.5"), true_mult);
        PolynomialHP dpoly = DifferentiationHP::derivative(poly, 0, 1);
        
        mpreal x = mpreal("0.3");  // Start further away
        mpreal true_root = mpreal("0.5");
        mpreal tolerance = mpreal("1e-50");
        
        unsigned int max_iters = 50;
        unsigned int iter;
        mpreal final_error;
        bool converged = false;
        bool oscillated = false;
        mpreal prev_error = mpreal("1e10");
        
        for (iter = 0; iter < max_iters; ++iter) {
            mpreal error = abs(x - true_root);
            
            if (error < tolerance) {
                converged = true;
                final_error = error;
                break;
            }
            
            // Check for oscillation (error increasing)
            if (iter > 2 && error > prev_error * mpreal("2.0")) {
                oscillated = true;
            }
            
            prev_error = error;
            final_error = error;
            
            // Modified Newton step
            mpreal f = poly.evaluate(x);
            mpreal df = dpoly.evaluate(x);
            
            if (abs(df) < mpreal("1e-200")) {
                break;
            }
            
            mpreal step = mpreal(true_mult) * f / df;
            x = x - step;
        }
        
        std::cout << std::setw(9) << prec << " | ";
        std::cout << std::setw(5) << iter << " | ";
        std::cout << std::scientific << std::setprecision(2) << toDouble(final_error) << " | ";
        std::cout << std::setw(9) << (converged ? "Yes" : "No") << " | ";
        
        if (converged && !oscillated) {
            std::cout << "✓ OK";
            if (min_working_precision == 0) {
                min_working_precision = prec;
            }
        } else if (converged && oscillated) {
            std::cout << "⚠ Oscillated";
        } else if (oscillated) {
            std::cout << "✗ Oscillated";
        } else {
            std::cout << "✗ No conv";
        }
        std::cout << "\n";
    }
    
    if (min_working_precision > 0) {
        std::cout << "→ Minimum precision for convergence: " << min_working_precision << " bits\n";
    } else {
        std::cout << "→ Convergence failed at all tested precisions!\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Precision Requirements Analysis\n";
    std::cout << "========================================\n";
    
    // Test multiplicities 2-10
    for (unsigned int m = 2; m <= 10; ++m) {
        testMultiplicityDetection(m);
        testConvergenceWithKnownMultiplicity(m);
    }
    
    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

