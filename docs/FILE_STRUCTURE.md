# File Structure After Refactorization

## Overview

This document shows the complete file structure after implementing the three-tier system.

## Directory Structure

```
polynomial-solver/
├── CMakeLists.txt                          # Root build configuration
├── README.md                               # Updated with high-precision info
│
├── cmake/                                  # NEW: CMake modules
│   └── FindHighPrecision.cmake            # NEW: Dependency detection
│
├── include/
│   ├── config.h.in                        # NEW: CMake configuration template
│   ├── config.h                           # GENERATED: Compile-time configuration
│   │
│   ├── high_precision_types.h             # NEW: Type definitions for multiprecision
│   ├── precision_conversion.h             # NEW: Conversion utilities
│   │
│   ├── polynomial.h                       # MODIFIED: Tier 1, optionally delegates to templates
│   ├── de_casteljau.h                     # MODIFIED: Tier 1, optionally delegates to templates
│   ├── differentiation.h                  # MODIFIED: Tier 1, optionally delegates to templates
│   ├── result_refiner.h                   # MODIFIED: Tier 1, optionally delegates to templates
│   │
│   ├── polynomial_hp.h                    # NEW: Tier 2 - Fixed high-precision (no templates)
│   ├── de_casteljau_hp.h                  # NEW: Tier 2 - Fixed high-precision
│   ├── differentiation_hp.h               # NEW: Tier 2 - Fixed high-precision
│   ├── result_refiner_hp.h                # NEW: Tier 2 - Fixed high-precision
│   │
│   ├── polynomial_mp.h                    # NEW: Tier 3 - Template-based mpreal API
│   ├── de_casteljau_mp.h                  # NEW: Tier 3 - Template-based mpreal API
│   ├── differentiation_mp.h               # NEW: Tier 3 - Template-based mpreal API
│   ├── result_refiner_mp.h                # NEW: Tier 3 - Template-based mpreal API
│   │
│   ├── polynomial_quad.h                  # NEW: Tier 3 - Template-based quadmath API (optional)
│   ├── de_casteljau_quad.h                # NEW: Tier 3 - Template-based quadmath API (optional)
│   ├── differentiation_quad.h             # NEW: Tier 3 - Template-based quadmath API (optional)
│   ├── result_refiner_quad.h              # NEW: Tier 3 - Template-based quadmath API (optional)
│   │
│   ├── polynomial_solver.h                # NEW: Unified header (includes appropriate headers)
│   │
│   └── detail/                            # NEW: Private template implementations
│       ├── polynomial_impl.h              # NEW: Template implementation for Polynomial
│       ├── de_casteljau_impl.h            # NEW: Template implementation for DeCasteljau
│       ├── differentiation_impl.h         # NEW: Template implementation for Differentiation
│       └── result_refiner_impl.h          # NEW: Template implementation for ResultRefiner
│
├── src/
│   ├── polynomial.cpp                     # MODIFIED: Optionally delegates to templates
│   ├── de_casteljau.cpp                   # MODIFIED: Optionally delegates to templates
│   ├── differentiation.cpp                # MODIFIED: Optionally delegates to templates
│   ├── result_refiner.cpp                 # MODIFIED: Optionally delegates to templates
│   │
│   ├── precision_conversion.cpp           # NEW: Conversion utilities implementation
│   │
│   ├── polynomial_hp.cpp                  # NEW: Tier 2 implementation
│   ├── de_casteljau_hp.cpp                # NEW: Tier 2 implementation
│   ├── differentiation_hp.cpp             # NEW: Tier 2 implementation
│   ├── result_refiner_hp.cpp              # NEW: Tier 2 implementation
│   │
│   ├── polynomial_mp.cpp                  # NEW: Tier 3 - Delegates to templates
│   ├── de_casteljau_mp.cpp                # NEW: Tier 3 - Delegates to templates
│   ├── differentiation_mp.cpp             # NEW: Tier 3 - Delegates to templates
│   ├── result_refiner_mp.cpp              # NEW: Tier 3 - Delegates to templates
│   │
│   ├── polynomial_quad.cpp                # NEW: Tier 3 - Quadmath (optional)
│   ├── de_casteljau_quad.cpp              # NEW: Tier 3 - Quadmath (optional)
│   ├── differentiation_quad.cpp           # NEW: Tier 3 - Quadmath (optional)
│   └── result_refiner_quad.cpp            # NEW: Tier 3 - Quadmath (optional)
│
├── tests/
│   ├── test_polynomial.cpp                # EXISTING: Tests for Tier 1
│   ├── test_de_casteljau.cpp              # EXISTING: Tests for Tier 1
│   ├── test_differentiation.cpp           # EXISTING: Tests for Tier 1
│   ├── test_result_refiner.cpp            # EXISTING: Tests for Tier 1
│   │
│   ├── test_high_precision_fallback.cpp   # NEW: Tests for Tier 2
│   ├── test_high_precision_templates.cpp  # NEW: Tests for Tier 3
│   ├── test_precision_conversion.cpp      # NEW: Tests for conversion utilities
│   └── test_all_configurations.cpp        # NEW: Comprehensive tests for all tiers
│
├── examples/
│   ├── example_double_precision.cpp       # NEW: Example for Tier 1
│   ├── example_high_precision_fallback.cpp # NEW: Example for Tier 2
│   ├── example_high_precision_templates.cpp # NEW: Example for Tier 3
│   └── example_wilkinson_polynomial.cpp   # NEW: Challenging example using all tiers
│
├── docs/
│   ├── REFACTORIZATION_PLAN.md            # This plan
│   ├── FILE_STRUCTURE.md                  # This document
│   ├── BUILD.md                           # NEW: Build instructions
│   ├── API.md                             # NEW: API documentation
│   ├── HIGH_PRECISION_MINIMAL_DEPENDENCIES.md
│   ├── HIGH_PRECISION_SETUP_GUIDE.md
│   ├── HIGH_PRECISION_MULTIPLE_TYPES.md
│   └── ...                                # Other existing docs
│
└── external/                              # OPTIONAL: For bundled dependencies
    └── boost/                             # OPTIONAL: Boost headers (if not system-installed)
```

## File Count Summary

| Category | Tier 1 | Tier 2 | Tier 3 | Total |
|----------|--------|--------|--------|-------|
| **Headers** | 4 (existing) | +4 new | +8 new (+4 detail) | 20 |
| **Sources** | 4 (existing) | +5 new | +8 new | 17 |
| **Tests** | 4 (existing) | +1 new | +3 new | 8 |
| **Examples** | 0 | +1 new | +3 new | 4 |
| **CMake** | 1 (existing) | +1 new | 0 | 2 |
| **Docs** | ~10 (existing) | +2 new | +1 new | ~13 |

**Total new files**: ~40 files

## Compilation Matrix

| Configuration | Files Compiled | Tier |
|---------------|----------------|------|
| Default (no flags) | Tier 1 only (4 cpp files) | Tier 1 |
| `-DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON` | Tier 1 + Tier 2 (9 cpp files) | Tier 2 |
| `-DENABLE_HIGH_PRECISION=ON` | Tier 1 + Tier 3 (12 cpp files) | Tier 3 |
| `-DENABLE_HIGH_PRECISION=ON -DENABLE_QUADMATH=ON` | Tier 1 + Tier 3 + Quad (16 cpp files) | Tier 3 |

## Key Files Explained

### CMake Configuration

**`cmake/FindHighPrecision.cmake`**
- Detects Boost headers
- Detects MPFR/GMP libraries
- Detects quadmath library
- Sets configuration variables
- Provides helpful error messages

**`include/config.h.in`** → **`include/config.h`**
- CMake template that generates compile-time configuration
- Defines macros: `ENABLE_HIGH_PRECISION`, `USE_TEMPLATES`, etc.

### Type Definitions

**`include/high_precision_types.h`**
- Defines `mpreal` type (MPFR or cpp_dec_float backend)
- Defines `quad` type (if quadmath available)
- Defines precision constants

**`include/precision_conversion.h`**
- Utilities to convert between double and high-precision types
- Handles precision loss warnings

### Tier 1: Minimum Version (Always Available)

**`include/polynomial.h`**, **`src/polynomial.cpp`**
- Original double-precision implementation
- Optionally delegates to templates if `USE_TEMPLATES` is defined
- API remains unchanged

### Tier 2: Fallback Version (No Templates)

**`include/polynomial_hp.h`**, **`src/polynomial_hp.cpp`**
- Separate class for high-precision
- Copy of algorithms with `double` → `mpreal`
- Only compiled if `ENABLE_HIGH_PRECISION=ON` and `DISABLE_TEMPLATES=ON`

### Tier 3: Full Version (Templates)

**`include/detail/polynomial_impl.h`**
- Private template implementation
- Header-only (in `detail` namespace)
- Used by both double and high-precision APIs

**`include/polynomial_mp.h`**, **`src/polynomial_mp.cpp`**
- Public API for mpreal
- Delegates to `detail::PolynomialImpl<mpreal>`
- Only compiled if `ENABLE_HIGH_PRECISION=ON` and `USE_TEMPLATES=ON`

### Unified Header

**`include/polynomial_solver.h`**
- Single include for users
- Includes appropriate headers based on configuration
- Provides type aliases for convenience

## Usage Examples

### Tier 1: Double Precision Only

```cpp
#include "polynomial.h"

Polynomial poly(degrees, coeffs);
double result = poly.evaluate(params);
```

### Tier 2: Fallback High Precision

```cpp
#include "polynomial.h"
#include "polynomial_hp.h"
#include "precision_conversion.h"

// Start with double
Polynomial poly_d(degrees, coeffs);

// Convert to high precision
PolynomialHP poly_hp = PolynomialHP::fromPolynomial(poly_d);

// Evaluate with high precision
mpreal result_hp = poly_hp.evaluate(toHighPrecision(params));
```

### Tier 3: Template High Precision

```cpp
#include "polynomial.h"
#include "polynomial_mp.h"

// Double precision
Polynomial poly_d(degrees, coeffs_d);
double result_d = poly_d.evaluate(params_d);

// High precision (mpreal)
PolynomialMP poly_mp(degrees, coeffs_mp);
mpreal result_mp = poly_mp.evaluate(params_mp);

// Both use the same template implementation!
```

## Next Steps

This file structure will be created incrementally following the phases in `REFACTORIZATION_PLAN.md`.

