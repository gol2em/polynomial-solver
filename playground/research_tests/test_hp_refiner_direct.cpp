/**
 * @file test_hp_refiner_direct.cpp
 * @brief Test HP refiner directly to verify multiplicity detection
 */

#ifdef ENABLE_HIGH_PRECISION

#include "hp/result_refiner_hp.h"
#include "hp/polynomial_hp.h"
#include "core/polynomial.h"
#include "hp/precision_context.h"
#include "hp/precision_conversion.h"
#include <iostream>
#include <iomanip>

using namespace polynomial_solver;

// Create polynomial in DOUBLE precision (loses precision for high m)
Polynomial createMultipleRootPolynomial(double root, unsigned int multiplicity) {
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
    return Polynomial::fromPower(degrees, power_coeffs);
}

// Create polynomial in HIGH PRECISION (preserves precision for high m)
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

int main() {
    std::cout << "========================================\n";
    std::cout << "  HP Refiner Direct Test\n";
    std::cout << "========================================\n";

    // Test with 512-bit precision for better convergence
    PrecisionContext ctx(512);
    
    std::cout << "\n=== Test 1: Polynomial created in DOUBLE precision ===\n";
    std::cout << "(This simulates current workflow: double→Bernstein→HP)\n";
    std::cout << "LIMITATION: Bernstein coefficients computed in double lose precision!\n";
    std::cout << "True m | Detected m | Converged | Iterations | Error Bound\n";
    std::cout << "-------+------------+-----------+------------+-------------\n";

    for (unsigned int m = 2; m <= 10; ++m) {
        Polynomial poly = createMultipleRootPolynomial(0.5, m);
        PolynomialHP poly_hp(poly);  // Convert double→HP (loses precision!)

        RefinementConfigHP config;
        config.target_tolerance_str = "1e-50";
        config.residual_tolerance_str = "1e-50";
        config.max_newton_iters = 100;
        config.max_multiplicity = 15;
        config.taylor_ratio_threshold = 10.0;
        config.multiplicity_hint = 0;  // No hint - let it auto-detect

        double initial_guess = 0.48;
        RefinedRootHP result = ResultRefinerHP::refineRoot1D(initial_guess, poly_hp, config);

        std::cout << std::setw(6) << m << " | ";
        std::cout << std::setw(10) << result.multiplicity << " | ";
        std::cout << std::setw(9) << (result.converged ? "Yes" : "No") << " | ";
        std::cout << std::setw(10) << result.iterations << " | ";

        if (result.has_guaranteed_bounds) {
            std::cout << std::scientific << std::setprecision(2) << toDouble(result.max_error);
        } else {
            std::cout << "N/A";
        }
        std::cout << "\n";
    }

    std::cout << "\n=== Test 2: Polynomial created in HIGH precision ===\n";
    std::cout << "(IDEAL workflow: power coeffs→HP→Bernstein in HP)\n";
    std::cout << "This preserves full precision of coefficients!\n";
    std::cout << "True m | Detected m | Converged | Iterations | Error Bound\n";
    std::cout << "-------+------------+-----------+------------+-------------\n";

    for (unsigned int m = 2; m <= 10; ++m) {
        PolynomialHP poly_hp = createMultipleRootPolynomialHP(mpreal("0.5"), m);

        RefinementConfigHP config;
        config.target_tolerance_str = "1e-100";  // Stricter tolerance for 512-bit
        config.residual_tolerance_str = "1e-100";
        config.max_newton_iters = 200;  // More iterations
        config.max_multiplicity = 15;
        config.taylor_ratio_threshold = 10.0;
        config.multiplicity_hint = 0;

        double initial_guess = 0.48;
        RefinedRootHP result = ResultRefinerHP::refineRoot1D(initial_guess, poly_hp, config);

        std::cout << std::setw(6) << m << " | ";
        std::cout << std::setw(10) << result.multiplicity << " | ";
        std::cout << std::setw(9) << (result.converged ? "Yes" : "No") << " | ";
        std::cout << std::setw(10) << result.iterations << " | ";

        if (result.has_guaranteed_bounds) {
            std::cout << std::scientific << std::setprecision(2) << toDouble(result.max_error);
        } else {
            std::cout << "N/A";
        }
        std::cout << "\n";
    }
    
    // Test with different ratio thresholds for higher multiplicities
    std::cout << "\n=== Testing with different ratio thresholds (m=8-12) ===\n";
    std::cout << "True m | Thresh=10 | Thresh=50 | Thresh=100\n";
    std::cout << "-------+-----------+-----------+-----------\n";
    
    for (unsigned int m = 8; m <= 12; ++m) {
        Polynomial poly = createMultipleRootPolynomial(0.5, m);
        PolynomialHP poly_hp(poly);
        
        std::cout << std::setw(6) << m << " | ";
        
        for (double thresh : {10.0, 50.0, 100.0}) {
            RefinementConfigHP config;
            config.target_tolerance_str = "1e-50";
            config.max_multiplicity = 20;
            config.taylor_ratio_threshold = thresh;
            config.multiplicity_hint = 0;
            
            RefinedRootHP result = ResultRefinerHP::refineRoot1D(0.48, poly_hp, config);
            std::cout << std::setw(9) << result.multiplicity;
            if (thresh < 100.0) std::cout << " | ";
        }
        std::cout << "\n";
    }
    
    // Test with multiplicity hint
    std::cout << "\n=== Testing with multiplicity hint ===\n";
    std::cout << "True m | Hint | Detected m | Converged\n";
    std::cout << "-------+------+------------+----------\n";
    
    for (unsigned int m = 3; m <= 8; ++m) {
        Polynomial poly = createMultipleRootPolynomial(0.5, m);
        PolynomialHP poly_hp(poly);
        
        // Test with correct hint
        RefinementConfigHP config;
        config.target_tolerance_str = "1e-50";
        config.max_multiplicity = 15;
        config.multiplicity_hint = m;  // Provide correct hint
        
        RefinedRootHP result = ResultRefinerHP::refineRoot1D(0.48, poly_hp, config);
        
        std::cout << std::setw(6) << m << " | ";
        std::cout << std::setw(4) << m << " | ";
        std::cout << std::setw(10) << result.multiplicity << " | ";
        std::cout << std::setw(9) << (result.converged ? "Yes" : "No") << "\n";
    }
    
    return 0;
}

#else
int main() {
    std::cout << "High precision not enabled\n";
    return 1;
}
#endif

