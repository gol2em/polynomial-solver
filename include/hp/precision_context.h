#ifndef POLYNOMIAL_SOLVER_PRECISION_CONTEXT_H
#define POLYNOMIAL_SOLVER_PRECISION_CONTEXT_H

/**
 * @file precision_context.h
 * @brief RAII-style precision management for temporary precision changes
 * 
 * This header provides PrecisionContext class for safe, temporary precision changes.
 * The precision is automatically restored when the context goes out of scope.
 * 
 * @section usage Usage Example
 * 
 * @code{.cpp}
 * #include "hp/precision_context.h"
 * 
 * void example() {
 *     setPrecision(256);  // Base precision: 256 bits
 *     
 *     mpreal x = 1.0;  // x has 256-bit precision
 *     
 *     {
 *         PrecisionContext ctx(512);  // Temporarily increase to 512 bits
 *         mpreal y = 2.0;             // y has 512-bit precision
 *         
 *         // Do high-precision computation
 *         mpreal result = compute(y);
 *     }  // Precision automatically restored to 256 bits
 *     
 *     mpreal z = 3.0;  // z has 256-bit precision again
 * }
 * @endcode
 * 
 * @section thread_safety Thread Safety
 * 
 * PrecisionContext is NOT thread-safe because it modifies global MPFR state.
 * Do not use PrecisionContext across multiple threads without external synchronization.
 * 
 * @section backends Backend Support
 * 
 * - MPFR backend: Full support, precision is changed and restored
 * - cpp_dec_float backend: No-op, precision is fixed at compile time
 * - quadmath backend: No-op, precision is fixed at 128 bits
 */

#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "hp/high_precision_types.h"

namespace polynomial_solver {

/**
 * @class PrecisionContext
 * @brief RAII wrapper for temporary precision changes
 * 
 * This class provides a safe way to temporarily change the default precision
 * for high-precision arithmetic. The original precision is automatically
 * restored when the PrecisionContext object is destroyed.
 * 
 * @note This only works with MPFR backend. For other backends, this class
 *       does nothing (precision is fixed at compile time).
 * 
 * @warning This class is NOT thread-safe. Do not use across multiple threads.
 * 
 * @section example Example Usage
 * 
 * @code{.cpp}
 * // Base precision: 256 bits
 * setPrecision(256);
 * mpreal x = 1.0;
 * 
 * {
 *     // Temporarily increase precision to 512 bits
 *     PrecisionContext high_prec(512);
 *     mpreal y = 2.0;  // y has 512-bit precision
 *     
 *     // Nested contexts work correctly
 *     {
 *         PrecisionContext ultra_high_prec(1024);
 *         mpreal z = 3.0;  // z has 1024-bit precision
 *     }  // Restored to 512 bits
 *     
 *     mpreal w = 4.0;  // w has 512-bit precision
 * }  // Restored to 256 bits
 * 
 * mpreal v = 5.0;  // v has 256-bit precision
 * @endcode
 */
class PrecisionContext {
public:
    /**
     * @brief Construct a precision context with new precision
     * 
     * Saves the current precision and sets a new default precision.
     * 
     * @param new_precision New precision in bits
     * 
     * @note For non-MPFR backends, this constructor does nothing.
     */
    explicit PrecisionContext(int new_precision)
        : old_precision_(getPrecision())
    {
        setPrecision(new_precision);
    }
    
    /**
     * @brief Destructor - restores the original precision
     * 
     * Automatically restores the precision that was active when this
     * context was created.
     */
    ~PrecisionContext() {
        setPrecision(old_precision_);
    }
    
    /**
     * @brief Get the precision that was active before this context
     * 
     * @return Original precision in bits
     */
    int getOldPrecision() const {
        return old_precision_;
    }
    
    /**
     * @brief Get the current precision
     * 
     * @return Current precision in bits
     */
    int getCurrentPrecision() const {
        return getPrecision();
    }
    
    // Disable copying and moving to prevent precision confusion
    PrecisionContext(const PrecisionContext&) = delete;
    PrecisionContext& operator=(const PrecisionContext&) = delete;
    PrecisionContext(PrecisionContext&&) = delete;
    PrecisionContext& operator=(PrecisionContext&&) = delete;

private:
    int old_precision_;  ///< Precision to restore on destruction
};

/**
 * @class PrecisionContextDigits
 * @brief RAII wrapper for temporary precision changes (decimal digits)
 * 
 * Similar to PrecisionContext, but takes precision in decimal digits
 * instead of bits.
 * 
 * @code{.cpp}
 * setPrecisionDigits(50);  // Base: 50 decimal digits
 * 
 * {
 *     PrecisionContextDigits ctx(100);  // Temporarily: 100 decimal digits
 *     mpreal x = 1.0;  // x has ~332-bit precision
 * }  // Restored to ~166 bits
 * @endcode
 */
class PrecisionContextDigits {
public:
    /**
     * @brief Construct a precision context with new precision (in digits)
     * 
     * @param new_precision_digits New precision in decimal digits
     */
    explicit PrecisionContextDigits(int new_precision_digits)
        : old_precision_(getPrecision())
    {
        setPrecisionDigits(new_precision_digits);
    }
    
    /**
     * @brief Destructor - restores the original precision
     */
    ~PrecisionContextDigits() {
        setPrecision(old_precision_);
    }
    
    // Disable copying and moving
    PrecisionContextDigits(const PrecisionContextDigits&) = delete;
    PrecisionContextDigits& operator=(const PrecisionContextDigits&) = delete;
    PrecisionContextDigits(PrecisionContextDigits&&) = delete;
    PrecisionContextDigits& operator=(PrecisionContextDigits&&) = delete;

private:
    int old_precision_;  ///< Precision (in bits) to restore on destruction
};

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // POLYNOMIAL_SOLVER_PRECISION_CONTEXT_H

