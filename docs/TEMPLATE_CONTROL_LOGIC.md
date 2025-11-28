# Template Control Logic - Complete Explanation

## Your Question

> "How is template disabled for fallback and default version? The implementation is still written with templates. Explain the logic of the template control."

## The Answer: Conditional Compilation + Separate Files

**Key insight**: Templates are NOT disabled in the template code itself. Instead, we use **conditional compilation** to choose between:
1. **Template implementation** (in `detail/` headers)
2. **Non-template implementation** (in regular `.cpp` files)

---

## Three-Tier Architecture

### Tier 1: Default (Double Precision Only)

**Build**: `cmake ..` (no flags)

**What gets compiled**:
```
src/polynomial.cpp          ✅ Compiled (double precision)
src/de_casteljau.cpp        ✅ Compiled (double precision)
include/detail/*.h          ❌ NOT included (templates not used)
src/polynomial_hp.cpp       ❌ NOT compiled (high precision disabled)
src/polynomial_mp.cpp       ❌ NOT compiled (high precision disabled)
```

**Result**: Only double-precision code, no templates, no high-precision.

---

### Tier 2: Fallback (No Templates)

**Build**: `cmake .. -DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON`

**What gets compiled**:
```
src/polynomial.cpp          ✅ Compiled (double precision)
src/polynomial_hp.cpp       ✅ Compiled (mpreal, NO templates)
include/detail/*.h          ❌ NOT included (templates disabled)
src/polynomial_mp.cpp       ❌ NOT compiled (templates disabled)
```

**Result**: Two separate implementations (double and mpreal), no templates, code duplication.

---

### Tier 3: Full (Templates)

**Build**: `cmake .. -DENABLE_HIGH_PRECISION=ON`

**What gets compiled**:
```
src/polynomial.cpp          ✅ Compiled (delegates to templates)
src/polynomial_mp.cpp       ✅ Compiled (delegates to templates)
include/detail/*.h          ✅ Included (template implementation)
src/polynomial_hp.cpp       ❌ NOT compiled (fallback not needed)
```

**Result**: Template implementation shared by both double and mpreal, no code duplication.

---

## How It Works: File-Level Control

### CMakeLists.txt Controls Which Files to Compile

```cmake
# Always compile core files
add_library(polynomial_solver
    src/polynomial.cpp
    src/de_casteljau.cpp
    src/differentiation.cpp
    # ...
)

# Conditionally compile high-precision files
if(ENABLE_HIGH_PRECISION)
    if(DISABLE_TEMPLATES)
        # Tier 2: Fallback version (no templates)
        target_sources(polynomial_solver PRIVATE
            src/polynomial_hp.cpp      # Non-template high-precision
            src/de_casteljau_hp.cpp
        )
    else()
        # Tier 3: Template version
        target_sources(polynomial_solver PRIVATE
            src/polynomial_mp.cpp      # Template-based high-precision
            src/de_casteljau_mp.cpp
        )
    endif()
endif()
```

**Key point**: Different `.cpp` files are compiled based on flags!

---

## How It Works: Code-Level Control

### Tier 1: Default Implementation (src/polynomial.cpp)

```cpp
// src/polynomial.cpp
#include "polynomial.h"
#include "de_casteljau.h"

// NO templates, just plain double precision
double Polynomial::evaluate(const std::vector<double>& parameters) const {
    return DeCasteljau::evaluateTensorProduct(degrees_, bernstein_coeffs_, parameters);
}

// This file NEVER includes template headers
// This file is ALWAYS compiled
```

**No templates here!** Just regular C++ code with `double`.

---

### Tier 2: Fallback Implementation (src/polynomial_hp.cpp)

```cpp
// src/polynomial_hp.cpp
// Only compiled if: ENABLE_HIGH_PRECISION=ON and DISABLE_TEMPLATES=ON

#include "polynomial_hp.h"
#include "high_precision_types.h"

// NO templates, just plain mpreal
mpreal PolynomialHP::evaluate(const std::vector<mpreal>& parameters) const {
    // Same algorithm as Polynomial::evaluate, but with mpreal instead of double
    return evaluateTensorProductHP(degrees_, bernstein_coeffs_, parameters);
}

// This is a COPY of the algorithm with double → mpreal
// Code duplication, but no templates!
```

**No templates here either!** Just regular C++ code with `mpreal`.

**This file is only compiled when templates are disabled.**

---

### Tier 3: Template Implementation

#### Step 1: Template Code (include/detail/polynomial_impl.h)

```cpp
// include/detail/polynomial_impl.h
// Header-only template implementation
// Only included when USE_TEMPLATES is defined

#ifndef POLYNOMIAL_IMPL_H
#define POLYNOMIAL_IMPL_H

#include "config.h"

#ifdef USE_TEMPLATES  // Only compile if templates enabled

namespace detail {

// This is the template implementation
template<typename T>
class PolynomialImpl {
    std::vector<unsigned int> degrees_;
    std::vector<T> bernstein_coeffs_;
    
public:
    PolynomialImpl(const std::vector<unsigned int>& degrees,
                   const std::vector<T>& coeffs)
        : degrees_(degrees), bernstein_coeffs_(coeffs) {}
    
    T evaluate(const std::vector<T>& parameters) const {
        // Template algorithm works for any T
        return evaluateTensorProductImpl(degrees_, bernstein_coeffs_, parameters);
    }
    
private:
    template<typename U>
    static U evaluateTensorProductImpl(const std::vector<unsigned int>& degrees,
                                       const std::vector<U>& coeffs,
                                       const std::vector<U>& params) {
        // Generic algorithm using type U
        // Works for double, mpreal, __float128, etc.
    }
};

} // namespace detail

#endif // USE_TEMPLATES

#endif // POLYNOMIAL_IMPL_H
```

**This file contains templates**, but it's **only included when `USE_TEMPLATES` is defined**.

#### Step 2: Public API Delegates to Templates (src/polynomial_mp.cpp)

```cpp
// src/polynomial_mp.cpp
// Only compiled if: ENABLE_HIGH_PRECISION=ON and USE_TEMPLATES=ON

#include "polynomial_mp.h"
#include "high_precision_types.h"
#include "detail/polynomial_impl.h"  // Include template implementation

// Public API for mpreal
class PolynomialMP {
    detail::PolynomialImpl<mpreal> impl_;  // Use template with mpreal
    
public:
    PolynomialMP(const std::vector<unsigned int>& degrees,
                 const std::vector<mpreal>& coeffs)
        : impl_(degrees, coeffs) {}
    
    mpreal evaluate(const std::vector<mpreal>& parameters) const {
        return impl_.evaluate(parameters);  // Delegate to template
    }
};
```

**This file uses templates**, but it's **only compiled when templates are enabled**.

#### Step 3: Default API Can Also Use Templates (src/polynomial.cpp)

```cpp
// src/polynomial.cpp
#include "polynomial.h"
#include "config.h"

#ifdef USE_TEMPLATES
    #include "detail/polynomial_impl.h"  // Include template if enabled
#endif

double Polynomial::evaluate(const std::vector<double>& parameters) const {
    #ifdef USE_TEMPLATES
        // Use template implementation
        detail::PolynomialImpl<double> impl(degrees_, bernstein_coeffs_);
        return impl.evaluate(parameters);
    #else
        // Use original implementation
        return DeCasteljau::evaluateTensorProduct(degrees_, bernstein_coeffs_, parameters);
    #endif
}
```

**This file can use templates OR not**, depending on `USE_TEMPLATES` flag.

---

## Preprocessor Flags Control Everything

### config.h (Generated by CMake)

```cpp
// config.h (generated from config.h.in)

// Tier 1: Neither flag defined
// (nothing)

// Tier 2: Only ENABLE_HIGH_PRECISION defined
#define ENABLE_HIGH_PRECISION
// USE_TEMPLATES is NOT defined

// Tier 3: Both flags defined
#define ENABLE_HIGH_PRECISION
#define USE_TEMPLATES
```

### How Flags Control Compilation

```cpp
// Example: polynomial.cpp

#include "config.h"

#ifdef USE_TEMPLATES
    // This code only compiled in Tier 3
    #include "detail/polynomial_impl.h"
#endif

double Polynomial::evaluate(...) const {
    #ifdef USE_TEMPLATES
        // Tier 3: Use template
        return detail::PolynomialImpl<double>(...).evaluate(...);
    #else
        // Tier 1 & 2: Use original code
        return originalImplementation(...);
    #endif
}
```

---

## Summary: How Templates Are "Disabled"

### Method 1: Don't Compile Template Files

```cmake
# CMakeLists.txt
if(DISABLE_TEMPLATES)
    # Don't add template-based files
    # target_sources(... src/polynomial_mp.cpp)  # NOT added
else()
    # Add template-based files
    target_sources(... src/polynomial_mp.cpp)    # Added
endif()
```

### Method 2: Don't Include Template Headers

```cpp
// src/polynomial.cpp
#ifdef USE_TEMPLATES
    #include "detail/polynomial_impl.h"  // Only if templates enabled
#endif
```

### Method 3: Use Preprocessor Guards in Template Headers

```cpp
// include/detail/polynomial_impl.h
#ifdef USE_TEMPLATES
    // Template code here
    template<typename T> class PolynomialImpl { ... };
#endif
```

---

## Key Insight

**Templates are not "disabled" in the template code itself.**

Instead:
1. **Template files are not compiled** when templates are disabled
2. **Template headers are not included** when templates are disabled
3. **Fallback files are compiled instead** when templates are disabled

**Result**: The compiler never sees the template code when templates are disabled!

---

## Visual Diagram: File Compilation Flow

### Tier 1: Default (No High Precision)

```
CMake Configuration:
  ENABLE_HIGH_PRECISION = OFF
  USE_TEMPLATES = (not set)

Files Compiled:
  ✅ src/polynomial.cpp (double only, no templates)
  ✅ src/de_casteljau.cpp (double only)
  ❌ src/polynomial_hp.cpp (not compiled)
  ❌ src/polynomial_mp.cpp (not compiled)
  ❌ include/detail/*.h (not included)

Result: Double precision only, no templates
```

### Tier 2: Fallback (No Templates)

```
CMake Configuration:
  ENABLE_HIGH_PRECISION = ON
  DISABLE_TEMPLATES = ON
  → USE_TEMPLATES = (not defined)

Files Compiled:
  ✅ src/polynomial.cpp (double, no templates)
  ✅ src/polynomial_hp.cpp (mpreal, NO templates, code duplication)
  ✅ src/de_casteljau_hp.cpp (mpreal, NO templates)
  ❌ src/polynomial_mp.cpp (not compiled)
  ❌ include/detail/*.h (not included)

Result: Two separate implementations, no templates
```

### Tier 3: Full (Templates)

```
CMake Configuration:
  ENABLE_HIGH_PRECISION = ON
  DISABLE_TEMPLATES = OFF (default)
  → USE_TEMPLATES = ON

Files Compiled:
  ✅ src/polynomial.cpp (delegates to template)
  ✅ src/polynomial_mp.cpp (delegates to template)
  ✅ include/detail/polynomial_impl.h (template, included)
  ❌ src/polynomial_hp.cpp (not compiled, not needed)

Result: Template implementation shared, no duplication
```

---

## Code Flow Comparison

### Tier 1: Direct Implementation

```
User calls: Polynomial::evaluate(params)
    ↓
src/polynomial.cpp: double Polynomial::evaluate(...)
    ↓
Direct implementation with double
    ↓
Return double result
```

### Tier 2: Separate Implementations

```
User calls: Polynomial::evaluate(params)          PolynomialHP::evaluate(params)
    ↓                                                  ↓
src/polynomial.cpp                                 src/polynomial_hp.cpp
    ↓                                                  ↓
Direct implementation with double                  Direct implementation with mpreal
    ↓                                                  ↓
Return double result                               Return mpreal result

Note: Same algorithm, different types, CODE DUPLICATION
```

### Tier 3: Template Implementation

```
User calls: Polynomial::evaluate(params)          PolynomialMP::evaluate(params)
    ↓                                                  ↓
src/polynomial.cpp                                 src/polynomial_mp.cpp
    ↓                                                  ↓
detail::PolynomialImpl<double>                     detail::PolynomialImpl<mpreal>
    ↓                                                  ↓
    └──────────────────┬──────────────────────────────┘
                       ↓
        include/detail/polynomial_impl.h (SHARED template)
                       ↓
        Template algorithm works for any type T
                       ↓
        Return T result (double or mpreal)

Note: Same algorithm, different instantiations, NO DUPLICATION
```

---

## Why This Design?

### Problem: Some Users Don't Want Templates

**Reasons**:
1. Longer compile times (template instantiation)
2. Larger binary size (multiple instantiations)
3. Harder to debug (template error messages)
4. Compatibility issues (older compilers)

### Solution: Provide Both Options

**Tier 2 (Fallback)**:
- ✅ No templates
- ✅ Faster compilation
- ✅ Smaller binary
- ❌ Code duplication
- ❌ Less flexible

**Tier 3 (Templates)**:
- ✅ No code duplication
- ✅ Flexible (supports any type)
- ✅ Better architecture
- ⚠️ Slightly longer compile time
- ⚠️ Slightly larger binary

### User Choice

```bash
# Don't want templates? Use Tier 2
cmake .. -DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON

# Want templates? Use Tier 3 (default)
cmake .. -DENABLE_HIGH_PRECISION=ON
```

---

## Implementation Status

### Phase 1: ✅ COMPLETE
- CMake configuration
- Preprocessor flags
- Conditional compilation logic

### Phase 2-7: 🚧 TO BE IMPLEMENTED
- Actual template implementations
- Fallback implementations
- Public APIs

**Current status**: Infrastructure is ready, implementation follows the plan in `REFACTORIZATION_PLAN.md`.

---

## Conclusion

**Templates are "disabled" by**:
1. ✅ Not compiling template-based `.cpp` files
2. ✅ Not including template headers
3. ✅ Compiling fallback `.cpp` files instead
4. ✅ Using preprocessor guards (`#ifdef USE_TEMPLATES`)

**The template code itself is not modified** - it's simply not compiled when templates are disabled.

**This gives users the choice**: templates (flexible, no duplication) vs no templates (simpler, faster compilation).

