#ifndef POLYNOMIAL_SOLVER_HIGH_PRECISION_TYPES_H
#define POLYNOMIAL_SOLVER_HIGH_PRECISION_TYPES_H

/**
 * @file high_precision_types.h
 * @brief Type definitions for high-precision arithmetic
 * 
 * This header defines the multiprecision types used throughout the library.
 * The actual types depend on the backend selected at compile time:
 * 
 * - MPFR Backend: Uses boost::multiprecision::mpfr_float (arbitrary precision, runtime control)
 * - cpp_dec_float Backend: Uses boost::multiprecision::cpp_dec_float (fixed precision)
 * - quadmath Backend: Uses __float128 (fixed 128-bit precision)
 * 
 * @section backends Backend Comparison
 * 
 * | Backend        | Precision      | Runtime Control | Speed    | Dependencies |
 * |----------------|----------------|-----------------|----------|--------------|
 * | MPFR           | Arbitrary      | ✅ Yes          | Fastest  | Boost+MPFR+GMP |
 * | cpp_dec_float  | Fixed (50/100) | ❌ No           | Slower   | Boost only   |
 * | quadmath       | Fixed (128-bit)| ❌ No           | Fast     | libquadmath  |
 * 
 * @section precision_control Runtime Precision Control (MPFR only)
 * 
 * With MPFR backend, you can change precision at runtime:
 * 
 * @code{.cpp}
 * #include "high_precision_types.h"
 * 
 * // Set precision to 256 bits
 * setPrecision(256);
 * 
 * // Create variables with current precision
 * mpreal x = 1.0;
 * mpreal y = 2.0;
 * 
 * // Change precision (affects new variables only)
 * setPrecision(512);
 * mpreal z = 3.0;  // z has 512-bit precision
 * 
 * // Query current precision
 * int bits = getPrecision();
 * @endcode
 */

#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

// Include Boost.Multiprecision headers based on backend
#ifdef USE_MPFR_BACKEND
    #include <boost/multiprecision/mpfr.hpp>
#endif

#ifdef USE_CPP_DEC_FLOAT_BACKEND
    #include <boost/multiprecision/cpp_dec_float.hpp>
#endif

// Use GCC's native quadmath instead of Boost wrapper
#ifdef ENABLE_QUADMATH
    extern "C" {
        #include <quadmath.h>
    }
#endif

namespace polynomial_solver {

// ============================================================================
// Precision Constants
// ============================================================================

/**
 * @brief Default precision for high-precision arithmetic (in bits)
 * 
 * This is used as the default precision for MPFR backend.
 * For cpp_dec_float, this is informational only (precision is fixed at compile time).
 */
constexpr int DEFAULT_HP_PRECISION_BITS = 256;

/**
 * @brief Default precision in decimal digits
 * 
 * Approximately 77 decimal digits for 256 bits.
 * Formula: decimal_digits ≈ bits * log10(2) ≈ bits * 0.30103
 */
constexpr int DEFAULT_HP_DECIMAL_DIGITS = 77;

/**
 * @brief Precision for quadmath (__float128)
 * 
 * Fixed at 113 bits of mantissa (approximately 34 decimal digits).
 */
constexpr int QUADMATH_PRECISION_BITS = 113;
constexpr int QUADMATH_DECIMAL_DIGITS = 34;

// ============================================================================
// Type Definitions
// ============================================================================

#ifdef USE_MPFR_BACKEND
    /**
     * @brief High-precision floating-point type (MPFR backend)
     * 
     * Uses MPFR library for arbitrary precision arithmetic.
     * Precision can be changed at runtime using setPrecision().
     * 
     * @note This is the recommended backend for adaptive precision algorithms.
     */
    using mpreal = boost::multiprecision::mpfr_float;
    
    /**
     * @brief Backend identifier
     */
    constexpr const char* HP_BACKEND_NAME = "MPFR";
    constexpr bool HP_RUNTIME_PRECISION = true;

#elif defined(USE_CPP_DEC_FLOAT_BACKEND)
    /**
     * @brief High-precision floating-point type (cpp_dec_float backend)
     * 
     * Uses Boost's header-only cpp_dec_float for fixed precision arithmetic.
     * Precision is fixed at compile time (100 decimal digits).
     * 
     * @note This backend is slower than MPFR but requires no external libraries.
     */
    using mpreal = boost::multiprecision::cpp_dec_float_100;
    
    /**
     * @brief Backend identifier
     */
    constexpr const char* HP_BACKEND_NAME = "cpp_dec_float";
    constexpr bool HP_RUNTIME_PRECISION = false;

#elif defined(USE_QUADMATH_BACKEND)
    /**
     * @brief High-precision floating-point type (quadmath backend)
     *
     * Uses GCC's native __float128 for fixed 128-bit precision arithmetic.
     * Precision is fixed at 113 bits (approximately 34 decimal digits).
     *
     * @note This backend is fast but limited to 128-bit precision.
     * @note Uses GCC's libquadmath directly, not Boost wrapper.
     */
    using mpreal = __float128;

    /**
     * @brief Backend identifier
     */
    constexpr const char* HP_BACKEND_NAME = "quadmath";
    constexpr bool HP_RUNTIME_PRECISION = false;

#else
    #error "No high-precision backend selected. This should not happen."
#endif

#ifdef ENABLE_QUADMATH
    /**
     * @brief Quadmath type (__float128)
     *
     * Available when quadmath support is enabled.
     * Can be used alongside other backends in template mode.
     * Uses GCC's native __float128 type directly.
     */
    using quad = __float128;
#endif

// ============================================================================
// Runtime Precision Control (MPFR backend only)
// ============================================================================

#ifdef USE_MPFR_BACKEND

/**
 * @brief Set the default precision for new mpreal variables
 *
 * This function sets the default precision (in bits) for all mpreal variables
 * created after this call. Existing variables are not affected.
 *
 * @param bits Precision in bits (must be >= 2)
 *
 * @note This only works with MPFR backend. For other backends, this function
 *       does nothing (precision is fixed at compile time).
 *
 * @warning This is a global setting that affects all threads. For thread-safe
 *          precision control, use PrecisionContext (see precision_context.h).
 *
 * @code{.cpp}
 * setPrecision(256);  // Set to 256 bits
 * mpreal x = 1.0;     // x has 256-bit precision
 *
 * setPrecision(512);  // Change to 512 bits
 * mpreal y = 2.0;     // y has 512-bit precision
 * // x still has 256-bit precision
 * @endcode
 */
inline void setPrecision(int bits) {
    // Convert bits to decimal digits for Boost.Multiprecision API
    // digits10 = floor(bits * log10(2)) ≈ floor(bits * 0.30103)
    int digits10 = static_cast<int>(bits * 0.30103);
    mpreal::default_precision(digits10);
}

/**
 * @brief Get the current default precision
 *
 * @return Current default precision in bits
 *
 * @code{.cpp}
 * setPrecision(256);
 * int bits = getPrecision();  // Returns 256
 * @endcode
 */
inline int getPrecision() {
    // Get precision in decimal digits and convert to bits
    // bits = ceil(digits10 * log2(10)) ≈ ceil(digits10 * 3.32193)
    int digits10 = static_cast<int>(mpreal::default_precision());
    return static_cast<int>(digits10 * 3.32193 + 0.5);
}

/**
 * @brief Set precision in decimal digits
 *
 * Convenience function to set precision in decimal digits instead of bits.
 *
 * @param digits Precision in decimal digits
 *
 * @note Conversion: bits = ceil(digits * log2(10)) ≈ ceil(digits * 3.32193)
 *
 * @code{.cpp}
 * setPrecisionDigits(50);  // Set to 50 decimal digits (≈166 bits)
 * @endcode
 */
inline void setPrecisionDigits(int digits) {
    // Convert decimal digits to bits: bits = ceil(digits * log2(10))
    // log2(10) ≈ 3.32193
    int bits = static_cast<int>(digits * 3.32193 + 0.5);
    setPrecision(bits);
}

/**
 * @brief Get current precision in decimal digits
 *
 * @return Current precision in decimal digits
 *
 * @note Conversion: digits = floor(bits * log10(2)) ≈ floor(bits * 0.30103)
 */
inline int getPrecisionDigits() {
    int bits = getPrecision();
    // Convert bits to decimal digits: digits = floor(bits * log10(2))
    // log10(2) ≈ 0.30103
    return static_cast<int>(bits * 0.30103);
}

#else // !USE_MPFR_BACKEND

/**
 * @brief Set precision (no-op for non-MPFR backends)
 *
 * For cpp_dec_float and quadmath backends, precision is fixed at compile time.
 * This function does nothing but is provided for API compatibility.
 *
 * @param bits Ignored
 */
inline void setPrecision(int bits) {
    (void)bits;  // Suppress unused parameter warning
    // No-op: precision is fixed at compile time
}

/**
 * @brief Get precision (returns fixed precision for non-MPFR backends)
 *
 * @return Fixed precision in bits
 */
inline int getPrecision() {
    #ifdef USE_CPP_DEC_FLOAT_BACKEND
        // cpp_dec_float_100 has 100 decimal digits
        // Convert to bits: 100 * log2(10) ≈ 332 bits
        return 332;
    #elif defined(USE_QUADMATH_BACKEND)
        return QUADMATH_PRECISION_BITS;
    #else
        return 0;
    #endif
}

/**
 * @brief Set precision in digits (no-op for non-MPFR backends)
 */
inline void setPrecisionDigits(int digits) {
    (void)digits;
    // No-op: precision is fixed at compile time
}

/**
 * @brief Get precision in digits (returns fixed precision)
 */
inline int getPrecisionDigits() {
    #ifdef USE_CPP_DEC_FLOAT_BACKEND
        return 100;  // cpp_dec_float_100
    #elif defined(USE_QUADMATH_BACKEND)
        return QUADMATH_DECIMAL_DIGITS;
    #else
        return 0;
    #endif
}

#endif // USE_MPFR_BACKEND

// ============================================================================
// Initialization
// ============================================================================

/**
 * @brief Initialize high-precision arithmetic with default precision
 *
 * This function should be called once at the start of your program.
 * It sets the default precision to DEFAULT_HP_PRECISION_BITS.
 *
 * @note For MPFR backend, this sets the global default precision.
 *       For other backends, this is a no-op.
 *
 * @code{.cpp}
 * int main() {
 *     initHighPrecision();  // Initialize with default precision
 *     // ... your code ...
 * }
 * @endcode
 */
inline void initHighPrecision() {
    #ifdef USE_MPFR_BACKEND
        setPrecision(DEFAULT_HP_PRECISION_BITS);
    #endif
}

/**
 * @brief Initialize high-precision arithmetic with custom precision
 *
 * @param bits Precision in bits
 *
 * @code{.cpp}
 * int main() {
 *     initHighPrecision(512);  // Initialize with 512-bit precision
 *     // ... your code ...
 * }
 * @endcode
 */
inline void initHighPrecision(int bits) {
    setPrecision(bits);
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // POLYNOMIAL_SOLVER_HIGH_PRECISION_TYPES_H

