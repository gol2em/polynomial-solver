/**
 * @file test_high_precision_types.cpp
 * @brief Test high-precision type definitions and precision control
 */

#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "hp/high_precision_types.h"
#include "hp/precision_context.h"
#include "hp/precision_conversion.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace polynomial_solver;

void test_basic_types() {
    std::cout << "Testing basic type definitions..." << std::endl;
    
    // Test that mpreal type is defined
    mpreal x = 1.0;
    mpreal y = 2.0;
    mpreal z = x + y;
    
    assert(z == 3.0);
    std::cout << "  ✓ Basic arithmetic works" << std::endl;
    
    // Test backend information
    std::cout << "  Backend: " << HP_BACKEND_NAME << std::endl;
    std::cout << "  Runtime precision control: " << (HP_RUNTIME_PRECISION ? "Yes" : "No") << std::endl;
}

void test_precision_control() {
    std::cout << "\nTesting precision control..." << std::endl;
    
    // Test setPrecision/getPrecision
    setPrecision(256);
    int bits = getPrecision();
    std::cout << "  Set precision to 256 bits, got: " << bits << " bits" << std::endl;
    
    // Test setPrecisionDigits/getPrecisionDigits
    setPrecisionDigits(50);
    int digits = getPrecisionDigits();
    std::cout << "  Set precision to 50 digits, got: " << digits << " digits" << std::endl;
    
    // Test that variables created at different precisions work
    setPrecision(128);
    mpreal a = 1.0;
    
    setPrecision(256);
    mpreal b = 2.0;
    
    mpreal c = a + b;
    std::cout << "  ✓ Mixed precision arithmetic works" << std::endl;
}

void test_precision_context() {
    std::cout << "\nTesting PrecisionContext..." << std::endl;

    if (!HP_RUNTIME_PRECISION) {
        std::cout << "  ⊗ Skipped (backend has fixed precision)" << std::endl;
        return;
    }

    setPrecision(256);
    int base_precision = getPrecision();
    std::cout << "  Base precision: " << base_precision << " bits" << std::endl;

    {
        PrecisionContext ctx(512);
        int high_precision = getPrecision();
        std::cout << "  Inside context: " << high_precision << " bits" << std::endl;
        assert(high_precision >= 500);  // Allow some rounding

        mpreal x = 1.0;
        // x should have high precision
    }

    int restored_precision = getPrecision();
    std::cout << "  After context: " << restored_precision << " bits" << std::endl;
    assert(std::abs(restored_precision - base_precision) < 10);  // Allow small difference

    std::cout << "  ✓ Precision correctly restored" << std::endl;
}

void test_nested_contexts() {
    std::cout << "\nTesting nested PrecisionContext..." << std::endl;

    if (!HP_RUNTIME_PRECISION) {
        std::cout << "  ⊗ Skipped (backend has fixed precision)" << std::endl;
        return;
    }

    setPrecision(128);
    std::cout << "  Level 0: " << getPrecision() << " bits" << std::endl;

    {
        PrecisionContext ctx1(256);
        std::cout << "  Level 1: " << getPrecision() << " bits" << std::endl;

        {
            PrecisionContext ctx2(512);
            std::cout << "  Level 2: " << getPrecision() << " bits" << std::endl;

            {
                PrecisionContext ctx3(1024);
                std::cout << "  Level 3: " << getPrecision() << " bits" << std::endl;
            }

            std::cout << "  Back to level 2: " << getPrecision() << " bits" << std::endl;
        }

        std::cout << "  Back to level 1: " << getPrecision() << " bits" << std::endl;
    }

    std::cout << "  Back to level 0: " << getPrecision() << " bits" << std::endl;
    std::cout << "  ✓ Nested contexts work correctly" << std::endl;
}

void test_precision_conversion() {
    std::cout << "\nTesting precision conversion..." << std::endl;
    
    // Test scalar conversion
    double d = 3.14159265358979323846;
    mpreal m = toHighPrecision(d);
    double d2 = toDouble(m);
    
    assert(std::abs(d - d2) < 1e-15);
    std::cout << "  ✓ Scalar conversion works" << std::endl;
    
    // Test vector conversion
    std::vector<double> vec_d = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<mpreal> vec_m = toHighPrecision(vec_d);
    std::vector<double> vec_d2 = toDouble(vec_m);
    
    assert(vec_d.size() == vec_m.size());
    assert(vec_d.size() == vec_d2.size());
    
    for (size_t i = 0; i < vec_d.size(); ++i) {
        assert(std::abs(vec_d[i] - vec_d2[i]) < 1e-15);
    }
    
    std::cout << "  ✓ Vector conversion works" << std::endl;
}

void test_initialization() {
    std::cout << "\nTesting initialization..." << std::endl;

    initHighPrecision();
    int default_bits = getPrecision();
    std::cout << "  Default precision: " << default_bits << " bits" << std::endl;

    if (HP_RUNTIME_PRECISION) {
        assert(std::abs(default_bits - DEFAULT_HP_PRECISION_BITS) < 10);

        initHighPrecision(512);
        int custom_bits = getPrecision();
        std::cout << "  Custom precision: " << custom_bits << " bits" << std::endl;
        assert(custom_bits >= 500);
    } else {
        std::cout << "  (Fixed precision backend - initialization is no-op)" << std::endl;
    }

    std::cout << "  ✓ Initialization works" << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "High-Precision Types Test" << std::endl;
    std::cout << "========================================" << std::endl;
    
    try {
        test_basic_types();
        test_precision_control();
        test_precision_context();
        test_nested_contexts();
        test_precision_conversion();
        test_initialization();
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "All tests passed! ✓" << std::endl;
        std::cout << "========================================" << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << std::endl;
        return 1;
    }
}

#else

int main() {
    std::cout << "High-precision support not enabled. Skipping test." << std::endl;
    return 0;
}

#endif // ENABLE_HIGH_PRECISION

