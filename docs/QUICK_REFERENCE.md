# Quick Reference: Three-Tier System

## Build Configurations

### Tier 1: Minimum Version (Default)
```bash
cmake ..
make
```
- ✅ Double precision only
- ✅ Zero dependencies
- ✅ Always available

### Tier 2: Fallback Version (No Templates)
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON
make
```
- ✅ Fixed high precision (mpreal)
- ✅ No templates
- ⚠️ Code duplication

### Tier 3: Full Version (Templates, Recommended)
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON
make
```
- ✅ Template-based
- ✅ Flexible precision
- ✅ No code duplication
- ✅ Best architecture

### Tier 3 with Quadmath
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DENABLE_QUADMATH=ON
make
```
- ✅ All features
- ✅ Supports double, mpreal, __float128

## API Quick Reference

### Tier 1: Double Precision

```cpp
#include "polynomial.h"

// Create polynomial
std::vector<unsigned int> degrees = {3, 2};
std::vector<double> coeffs = {1.0, 2.0, 3.0, ...};
Polynomial poly(degrees, coeffs);

// Evaluate
std::vector<double> params = {0.5, 0.3};
double result = poly.evaluate(params);

// Differentiate
Polynomial dpoly = Differentiation::derivative(poly, 0);

// Refine root
std::vector<double> initial_guess = {1.0, 0.5};
RefinedRoot result = ResultRefiner::refineRoot(poly, initial_guess);
```

### Tier 2: Fallback High Precision

```cpp
#include "polynomial.h"
#include "polynomial_hp.h"
#include "precision_conversion.h"

// Convert from double to high precision
Polynomial poly_d(degrees, coeffs_d);
PolynomialHP poly_hp = PolynomialHP::fromPolynomial(poly_d);

// Or create directly with mpreal
std::vector<mpreal> coeffs_hp = toHighPrecision(coeffs_d);
PolynomialHP poly_hp(degrees, coeffs_hp);

// Evaluate
std::vector<mpreal> params_hp = toHighPrecision(params_d);
mpreal result_hp = poly_hp.evaluate(params_hp);

// Differentiate
PolynomialHP dpoly_hp = DifferentiationHP::derivative(poly_hp, 0);

// Refine root
std::vector<mpreal> initial_guess_hp = toHighPrecision(initial_guess_d);
RefinedRootHP result_hp = ResultRefinerHP::refineRoot(poly_hp, initial_guess_hp);
```

### Tier 3: Template High Precision

```cpp
#include "polynomial.h"
#include "polynomial_mp.h"

// Set precision (runtime configurable!)
mpreal::set_default_prec(256);  // 77 decimal digits

// Create high-precision polynomial
std::vector<mpreal> coeffs_mp = {mpreal(1), mpreal(2), mpreal(3), ...};
PolynomialMP poly_mp(degrees, coeffs_mp);

// Evaluate
std::vector<mpreal> params_mp = {mpreal("0.5"), mpreal("0.3")};
mpreal result_mp = poly_mp.evaluate(params_mp);

// Differentiate
PolynomialMP dpoly_mp = DifferentiationMP::derivative(poly_mp, 0);

// Refine root
std::vector<mpreal> initial_guess_mp = {mpreal(1), mpreal("0.5")};
RefinedRootMP result_mp = ResultRefinerMP::refineRoot(poly_mp, initial_guess_mp);

// Change precision dynamically
mpreal::set_default_prec(512);  // 154 decimal digits
// ... continue with higher precision
```

### Tier 3: Quadmath (Optional)

```cpp
#include "polynomial_quad.h"

// Create quadmath polynomial
std::vector<__float128> coeffs_quad = {1.0Q, 2.0Q, 3.0Q, ...};
PolynomialQuad poly_quad(degrees, coeffs_quad);

// Evaluate
std::vector<__float128> params_quad = {0.5Q, 0.3Q};
__float128 result_quad = poly_quad.evaluate(params_quad);

// Print (requires quadmath_snprintf)
char buf[128];
quadmath_snprintf(buf, sizeof(buf), "%.36Qe", result_quad);
std::cout << "Result: " << buf << std::endl;
```

## Precision Comparison

| Type | Bits | Decimal Digits | Speed | Tier |
|------|------|----------------|-------|------|
| `double` | 64 | 15-17 | ⚡⚡⚡⚡⚡ | 1, 2, 3 |
| `__float128` | 128 | 33-36 | ⚡⚡⚡⚡ | 3 (optional) |
| `mpreal(128)` | 128 | 38 | ⚡⚡⚡ | 2, 3 |
| `mpreal(256)` | 256 | 77 | ⚡⚡ | 2, 3 |
| `mpreal(512)` | 512 | 154 | ⚡ | 2, 3 |
| `mpreal(1024)` | 1024 | 308 | 🐌 | 2, 3 |

## Dependencies

### Tier 1
- **None** (always available)

### Tier 2 & 3
- **Boost headers** (~50 MB, header-only)
- **MPFR library** (~1 MB, optional, can build from source)
- **GMP library** (~1 MB, optional, can build from source)
- **Fallback**: cpp_dec_float (header-only, no MPFR/GMP needed)

### Quadmath (Optional)
- **libquadmath** (usually included with GCC)

## When to Use Each Tier

### Use Tier 1 (Default) When:
- ✅ Standard precision (15-17 digits) is sufficient
- ✅ You want zero dependencies
- ✅ You want maximum performance
- ✅ You're solving well-conditioned problems

### Use Tier 2 (Fallback) When:
- ✅ You need high precision
- ✅ You don't want templates in your codebase
- ✅ You only need one fixed precision level
- ✅ You want simple integration

### Use Tier 3 (Templates) When:
- ✅ You need flexible precision (change at runtime)
- ✅ You want to support multiple numeric types
- ✅ You want the best long-term architecture
- ✅ You want to avoid code duplication
- ✅ You're solving ill-conditioned problems

## Common Workflows

### Workflow 1: Cascading Precision

```cpp
// Start with double (fast)
Polynomial poly_d(degrees, coeffs_d);
RefinedRoot result_d = ResultRefiner::refineRoot(poly_d, guess_d);

if (result_d.residual < 1e-14) {
    // Success with double precision
    return result_d;
}

#ifdef ENABLE_HIGH_PRECISION
    // Refine with high precision
    mpreal::set_default_prec(256);
    PolynomialMP poly_mp = convertToHighPrecision(poly_d);
    RefinedRootMP result_mp = ResultRefinerMP::refineRoot(
        poly_mp, convertToHighPrecision(result_d.root)
    );
    return result_mp;
#endif
```

### Workflow 2: Direct High Precision

```cpp
#ifdef ENABLE_HIGH_PRECISION
    // Start with high precision for ill-conditioned problems
    mpreal::set_default_prec(512);
    PolynomialMP poly_mp(degrees, coeffs_mp);
    RefinedRootMP result_mp = ResultRefinerMP::refineRoot(poly_mp, guess_mp);
#endif
```

### Workflow 3: Adaptive Precision

```cpp
#ifdef ENABLE_HIGH_PRECISION
    // Try increasing precisions until convergence
    for (int prec : {128, 256, 512, 1024}) {
        mpreal::set_default_prec(prec);
        RefinedRootMP result = ResultRefinerMP::refineRoot(poly_mp, guess_mp);
        
        if (result.converged) {
            std::cout << "Converged at " << prec << " bits" << std::endl;
            return result;
        }
    }
#endif
```

## Troubleshooting

### Error: "Boost headers not found"
```bash
# Option 1: Install system-wide
sudo apt install libboost-dev

# Option 2: Copy to project
wget https://boostorg.jfrog.io/.../boost_1_83_0.tar.gz
tar xzf boost_1_83_0.tar.gz
mkdir -p external
mv boost_1_83_0/boost external/
```

### Error: "MPFR library not found"
```bash
# Option 1: Install system-wide
sudo apt install libmpfr-dev libgmp-dev

# Option 2: Build locally (no admin rights)
# See docs/HIGH_PRECISION_SETUP_GUIDE.md

# Option 3: Use cpp_dec_float (automatic fallback)
# CMake will automatically use cpp_dec_float if MPFR not found
```

### Warning: "Templates disabled, using fallback"
```bash
# This is expected if you set DISABLE_TEMPLATES=ON
# To use templates (recommended):
cmake .. -DENABLE_HIGH_PRECISION=ON
# (without -DDISABLE_TEMPLATES=ON)
```

## Performance Tips

1. **Use Tier 1 first**: Always try double precision first, it's 10-100× faster
2. **Cascade precision**: Start low, increase only if needed
3. **Static linking**: Use static MPFR/GMP libraries to avoid runtime dependencies
4. **Minimize conversions**: Convert to high precision once, not repeatedly
5. **Reuse polynomials**: Don't recreate polynomial objects unnecessarily

## Further Reading

- **`docs/REFACTORIZATION_PLAN.md`** - Detailed implementation plan
- **`docs/FILE_STRUCTURE.md`** - Complete file structure
- **`docs/HIGH_PRECISION_SETUP_GUIDE.md`** - Setup instructions
- **`docs/HIGH_PRECISION_MINIMAL_DEPENDENCIES.md`** - Dependency options
- **`docs/BUILD.md`** - Build instructions (to be created)
- **`docs/API.md`** - API documentation (to be created)

