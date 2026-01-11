/**
 * @file precision_conversion.cpp
 * @brief Implementation of precision conversion utilities
 */

#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "hp/precision_conversion.h"

namespace polynomial_solver {

// Global flag for precision loss warnings (default: disabled)
bool g_warn_precision_loss = false;

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

