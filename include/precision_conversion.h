#ifndef POLYNOMIAL_SOLVER_PRECISION_CONVERSION_H
#define POLYNOMIAL_SOLVER_PRECISION_CONVERSION_H

/**
 * @file precision_conversion.h
 * @brief Utilities for converting between double and high-precision types
 * 
 * This header provides functions for safe conversion between double precision
 * and high-precision types, with optional warnings for precision loss.
 * 
 * @section usage Usage Example
 * 
 * @code{.cpp}
 * #include "precision_conversion.h"
 * 
 * // Convert double to high precision
 * double x = 3.14159265358979323846;
 * mpreal y = toHighPrecision(x);
 * 
 * // Convert high precision to double (may lose precision)
 * mpreal z = computeHighPrecision();
 * double w = toDouble(z);  // Precision loss warning if enabled
 * 
 * // Convert vectors
 * std::vector<double> vec_d = {1.0, 2.0, 3.0};
 * std::vector<mpreal> vec_hp = toHighPrecision(vec_d);
 * @endcode
 */

#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "high_precision_types.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

// For quadmath backend
#ifdef USE_QUADMATH_BACKEND
    #include <cstdio>
#endif

namespace polynomial_solver {

// ============================================================================
// Configuration
// ============================================================================

/**
 * @brief Enable/disable precision loss warnings
 * 
 * When enabled, converting from high precision to double will print a warning
 * if significant precision loss is detected.
 */
extern bool g_warn_precision_loss;

/**
 * @brief Set whether to warn about precision loss
 * 
 * @param enable True to enable warnings, false to disable
 */
inline void setWarnPrecisionLoss(bool enable) {
    g_warn_precision_loss = enable;
}

// ============================================================================
// Scalar Conversions
// ============================================================================

/**
 * @brief Convert double to high-precision type
 * 
 * @param value Double precision value
 * @return High-precision value
 * 
 * @note No precision is lost in this conversion (double has 53 bits of precision,
 *       high-precision types have much more).
 */
inline mpreal toHighPrecision(double value) {
    return mpreal(value);
}

/**
 * @brief Convert high-precision type to double
 * 
 * @param value High-precision value
 * @return Double precision value
 * 
 * @warning This may lose precision! If the value requires more than 53 bits
 *          of precision, the result will be rounded.
 * 
 * @note If g_warn_precision_loss is true, a warning will be printed if
 *       significant precision loss is detected.
 */
inline double toDouble(const mpreal& value) {
    // TODO: Add precision loss detection and warning
    return static_cast<double>(value);
}

// ============================================================================
// Vector Conversions
// ============================================================================

/**
 * @brief Convert vector of doubles to vector of high-precision values
 * 
 * @param values Vector of double precision values
 * @return Vector of high-precision values
 */
inline std::vector<mpreal> toHighPrecision(const std::vector<double>& values) {
    std::vector<mpreal> result;
    result.reserve(values.size());
    for (double value : values) {
        result.push_back(toHighPrecision(value));
    }
    return result;
}

/**
 * @brief Convert vector of high-precision values to vector of doubles
 * 
 * @param values Vector of high-precision values
 * @return Vector of double precision values
 * 
 * @warning This may lose precision for each element!
 */
inline std::vector<double> toDouble(const std::vector<mpreal>& values) {
    std::vector<double> result;
    result.reserve(values.size());
    for (const mpreal& value : values) {
        result.push_back(toDouble(value));
    }
    return result;
}

// ============================================================================
// String Conversions (for I/O without precision loss)
// ============================================================================

/**
 * @brief Convert high-precision value to string with full precision
 * 
 * @param value High-precision value
 * @param digits Number of decimal digits to output (0 = all available digits)
 * @return String representation
 * 
 * @note This is the recommended way to output high-precision values without
 *       losing precision.
 * 
 * @code{.cpp}
 * mpreal x = computeHighPrecision();
 * std::string str = toString(x, 50);  // 50 decimal digits
 * std::cout << str << std::endl;
 * @endcode
 */
inline std::string toString(const mpreal& value, int digits = 0) {
    if (digits <= 0) {
        digits = getPrecisionDigits();
    }

#ifdef USE_QUADMATH_BACKEND
    // Use quadmath's native formatting
    char buf[256];
    int width = digits + 10;  // Extra space for exponent
    quadmath_snprintf(buf, sizeof(buf), "%.*Qe", digits, value);
    return std::string(buf);
#else
    // Use Boost.Multiprecision's stream operators
    std::ostringstream oss;
    oss.precision(digits);
    oss << std::scientific << value;
    return oss.str();
#endif
}

/**
 * @brief Convert string to high-precision value
 * 
 * @param str String representation of a number
 * @return High-precision value
 * 
 * @note This preserves full precision from the string.
 * 
 * @code{.cpp}
 * std::string str = "3.14159265358979323846264338327950288419716939937510";
 * mpreal x = fromString(str);  // Full precision preserved
 * @endcode
 */
inline mpreal fromString(const std::string& str) {
#ifdef USE_QUADMATH_BACKEND
    // Use quadmath's native parsing
    return strtoflt128(str.c_str(), nullptr);
#else
    // Use Boost.Multiprecision's string constructor
    return mpreal(str);
#endif
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // POLYNOMIAL_SOLVER_PRECISION_CONVERSION_H

