# Supporting Multiple High-Precision Types

## Your Question

Can private templates support:
1. **double** (always available)
2. **__float128** (quadmath, 128-bit)
3. **mpreal** (Boost.Multiprecision, arbitrary precision: 128, 256, 512, 1024 bits)

**Answer: YES! Templates work with ANY numeric type!**

## The Magic: Templates Are Type-Agnostic

The template implementation works with **any type** that supports arithmetic operations (`+`, `-`, `*`, `/`).

### Single Template Implementation

```cpp
// src/de_casteljau.cpp

namespace detail {

// This template works with ANY numeric type T!
template<typename T>
T evaluate1D_impl(const std::vector<T>& bernstein_coeffs, T t) {
    if (bernstein_coeffs.empty()) {
        return T(0);
    }

    std::vector<T> work = bernstein_coeffs;
    const std::size_t n = work.size();

    for (std::size_t k = 1; k < n; ++k) {
        for (std::size_t i = 0; i < n - k; ++i) {
            work[i] = (T(1) - t) * work[i] + t * work[i + 1];
        }
    }

    return work[0];
}

} // namespace detail
```

**This ONE template works with:**
- `double`
- `__float128`
- `mpreal` (at any precision: 128, 256, 512, 1024 bits)
- `long double`
- `float`
- Any custom numeric type!

## Architecture: Multiple APIs, One Implementation

### Header: `include/de_casteljau.h`

```cpp
#ifndef DE_CASTELJAU_H
#define DE_CASTELJAU_H

#include <vector>

namespace polynomial_solver {

// API 1: Double precision (always available)
class DeCasteljau {
public:
    static double evaluate1D(const std::vector<double>& coeffs, double t);
};

#ifdef ENABLE_QUADMATH
    // API 2: Quad precision (128-bit float)
    #include <quadmath.h>
    
    class DeCasteljauQuad {
    public:
        static __float128 evaluate1D(const std::vector<__float128>& coeffs, __float128 t);
    };
#endif

#ifdef ENABLE_HIGH_PRECISION
    // API 3: Arbitrary precision (Boost.Multiprecision)
    #include <boost/multiprecision/mpfr.hpp>
    using mpreal = boost::multiprecision::mpfr_float;
    
    class DeCasteljauMP {
    public:
        static mpreal evaluate1D(const std::vector<mpreal>& coeffs, const mpreal& t);
    };
#endif

} // namespace polynomial_solver

#endif // DE_CASTELJAU_H
```

### Implementation: `src/de_casteljau.cpp`

```cpp
#include "de_casteljau.h"

namespace polynomial_solver {

namespace detail {

// Template implementation (WRITTEN ONCE!)
template<typename T>
T evaluate1D_impl(const std::vector<T>& bernstein_coeffs, T t) {
    if (bernstein_coeffs.empty()) {
        return T(0);
    }

    std::vector<T> work = bernstein_coeffs;
    const std::size_t n = work.size();

    for (std::size_t k = 1; k < n; ++k) {
        for (std::size_t i = 0; i < n - k; ++i) {
            work[i] = (T(1) - t) * work[i] + t * work[i + 1];
        }
    }

    return work[0];
}

} // namespace detail

// API 1: Double precision
double DeCasteljau::evaluate1D(const std::vector<double>& coeffs, double t) {
    return detail::evaluate1D_impl<double>(coeffs, t);
}

#ifdef ENABLE_QUADMATH
    // API 2: Quad precision
    __float128 DeCasteljauQuad::evaluate1D(const std::vector<__float128>& coeffs, __float128 t) {
        return detail::evaluate1D_impl<__float128>(coeffs, t);
    }
#endif

#ifdef ENABLE_HIGH_PRECISION
    // API 3: Arbitrary precision
    mpreal DeCasteljauMP::evaluate1D(const std::vector<mpreal>& coeffs, const mpreal& t) {
        return detail::evaluate1D_impl<mpreal>(coeffs, t);
    }
#endif

} // namespace polynomial_solver
```

**Result**: ONE template, THREE APIs (or more)!

## Usage Examples

### Example 1: Double Precision (Always Available)

```cpp
#include "de_casteljau.h"

std::vector<double> coeffs = {1.0, 2.0, 3.0};
double result = DeCasteljau::evaluate1D(coeffs, 0.5);
// result ≈ 2.0 (15-17 decimal digits)
```

### Example 2: Quad Precision (128-bit)

```cpp
#include "de_casteljau.h"

#ifdef ENABLE_QUADMATH
    std::vector<__float128> coeffs_quad = {1.0Q, 2.0Q, 3.0Q};
    __float128 result_quad = DeCasteljauQuad::evaluate1D(coeffs_quad, 0.5Q);
    // result_quad ≈ 2.0 (33-36 decimal digits)
    
    // Print with quadmath
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.36Qe", result_quad);
    std::cout << "Quad result: " << buf << std::endl;
#endif
```

### Example 3: Arbitrary Precision (Boost.Multiprecision)

```cpp
#include "de_casteljau.h"

#ifdef ENABLE_HIGH_PRECISION
    // 128-bit precision
    mpreal::set_default_prec(128);
    std::vector<mpreal> coeffs_128 = {mpreal(1), mpreal(2), mpreal(3)};
    mpreal result_128 = DeCasteljauMP::evaluate1D(coeffs_128, mpreal("0.5"));
    // result_128 ≈ 2.0 (38 decimal digits)
    
    // 256-bit precision
    mpreal::set_default_prec(256);
    std::vector<mpreal> coeffs_256 = {mpreal(1), mpreal(2), mpreal(3)};
    mpreal result_256 = DeCasteljauMP::evaluate1D(coeffs_256, mpreal("0.5"));
    // result_256 ≈ 2.0 (77 decimal digits)
    
    // 1024-bit precision
    mpreal::set_default_prec(1024);
    std::vector<mpreal> coeffs_1024 = {mpreal(1), mpreal(2), mpreal(3)};
    mpreal result_1024 = DeCasteljauMP::evaluate1D(coeffs_1024, mpreal("0.5"));
    // result_1024 ≈ 2.0 (308 decimal digits)
#endif
```

### Example 4: Using All Types Together

```cpp
#include "de_casteljau.h"
#include <iostream>
#include <iomanip>

void compare_precisions() {
    // Test polynomial: (x-1)^3 = x^3 - 3x^2 + 3x - 1
    // Bernstein form: [-1, 0, 0, 0]
    
    // Double precision
    std::vector<double> coeffs_d = {-1.0, 0.0, 0.0, 0.0};
    double t_d = 1.0 + 1e-10;  // Slightly perturbed
    double result_d = DeCasteljau::evaluate1D(coeffs_d, t_d);
    std::cout << "Double:  " << std::setprecision(17) << result_d << std::endl;
    
#ifdef ENABLE_QUADMATH
    // Quad precision
    std::vector<__float128> coeffs_q = {-1.0Q, 0.0Q, 0.0Q, 0.0Q};
    __float128 t_q = 1.0Q + 1e-10Q;
    __float128 result_q = DeCasteljauQuad::evaluate1D(coeffs_q, t_q);
    
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.36Qe", result_q);
    std::cout << "Quad:    " << buf << std::endl;
#endif
    
#ifdef ENABLE_HIGH_PRECISION
    // Arbitrary precision (256-bit)
    mpreal::set_default_prec(256);
    std::vector<mpreal> coeffs_mp = {mpreal(-1), mpreal(0), mpreal(0), mpreal(0)};
    mpreal t_mp = mpreal(1) + mpreal("1e-10");
    mpreal result_mp = DeCasteljauMP::evaluate1D(coeffs_mp, t_mp);
    std::cout << "Mpreal:  " << std::setprecision(50) << result_mp << std::endl;
#endif
}
```

## Build Configuration

### CMakeLists.txt

```cmake
# Option 1: Double precision only (default)
# cmake ..

# Option 2: Double + Quadmath
# cmake -DENABLE_QUADMATH=ON ..

# Option 3: Double + Boost.Multiprecision
# cmake -DENABLE_HIGH_PRECISION=ON ..

# Option 4: All three!
# cmake -DENABLE_QUADMATH=ON -DENABLE_HIGH_PRECISION=ON ..

option(ENABLE_QUADMATH "Enable quadmath (__float128) support" OFF)
option(ENABLE_HIGH_PRECISION "Enable Boost.Multiprecision support" OFF)

if(ENABLE_QUADMATH)
    find_library(QUADMATH_LIBRARY quadmath)
    if(QUADMATH_LIBRARY)
        add_definitions(-DENABLE_QUADMATH)
        target_link_libraries(polynomial_solver_lib ${QUADMATH_LIBRARY})
        message(STATUS "Quadmath support enabled")
    else()
        message(WARNING "Quadmath library not found")
    endif()
endif()

if(ENABLE_HIGH_PRECISION)
    find_package(Boost REQUIRED)
    find_library(MPFR_LIBRARY mpfr)
    find_library(GMP_LIBRARY gmp)
    
    if(MPFR_LIBRARY AND GMP_LIBRARY)
        add_definitions(-DENABLE_HIGH_PRECISION)
        target_link_libraries(polynomial_solver_lib ${MPFR_LIBRARY} ${GMP_LIBRARY})
        message(STATUS "Boost.Multiprecision support enabled")
    else()
        message(WARNING "MPFR/GMP libraries not found")
    endif()
endif()
```

## Precision Comparison Table

| Type | Bits | Decimal Digits | Speed | Use Case |
|------|------|----------------|-------|----------|
| `double` | 64 | 15-17 | Fastest | Default, 99% of cases |
| `__float128` | 128 | 33-36 | Fast | Fixed high precision |
| `mpreal(128)` | 128 | 38 | Medium | Flexible, same as quad |
| `mpreal(256)` | 256 | 77 | Slower | Very high precision |
| `mpreal(512)` | 512 | 154 | Slow | Extreme precision |
| `mpreal(1024)` | 1024 | 308 | Very slow | Research/validation |

## Key Advantages

### 1. Single Implementation

```cpp
// Write ONCE:
template<typename T>
T evaluate1D_impl(const std::vector<T>& coeffs, T t) {
    // Algorithm
}

// Works with ALL types:
// - double
// - __float128
// - mpreal (at any precision)
// - long double
// - float
// - Any custom type!
```

### 2. Flexible Precision with mpreal

```cpp
// Same API, different precisions:
mpreal::set_default_prec(128);   // 38 digits
mpreal result_128 = DeCasteljauMP::evaluate1D(coeffs, t);

mpreal::set_default_prec(256);   // 77 digits
mpreal result_256 = DeCasteljauMP::evaluate1D(coeffs, t);

mpreal::set_default_prec(1024);  // 308 digits
mpreal result_1024 = DeCasteljauMP::evaluate1D(coeffs, t);
```

### 3. Mix and Match

```cpp
// Start with double
double root_approx = solve_with_double();

#ifdef ENABLE_QUADMATH
    // Refine with quadmath
    __float128 root_quad = refine_with_quad(root_approx);
#endif

#ifdef ENABLE_HIGH_PRECISION
    // Further refine with mpreal
    mpreal::set_default_prec(512);
    mpreal root_hp = refine_with_mpreal(root_quad);
#endif
```

## Summary

### Question: Can private templates support multiple high-precision types?

**Answer: YES! Templates work with ANY numeric type!**

✅ **One template implementation** works with:
- `double` (always)
- `__float128` (if `ENABLE_QUADMATH=ON`)
- `mpreal` at 128, 256, 512, 1024 bits (if `ENABLE_HIGH_PRECISION=ON`)
- Any other numeric type you want to add later!

✅ **No code duplication** - Algorithm written once

✅ **Flexible precision** - Choose precision at runtime with `mpreal::set_default_prec()`

✅ **Compile-time control** - Enable only what you need

✅ **Zero overhead** - Disabled types don't affect binary size or performance

**This is the ultimate flexibility!**

