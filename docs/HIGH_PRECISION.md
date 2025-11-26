# High-Precision Arithmetic Support

## Overview

The polynomial solver supports three tiers of precision:

1. **Tier 1: Minimum Version** (default) - Double precision only, zero dependencies
2. **Tier 2: Fallback Version** - Fixed high precision, no templates
3. **Tier 3: Full Version** - Template-based, flexible precision, supports multiple types

## Quick Start

### Default Build (Tier 1 - Double Precision)
```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### High-Precision Build (Tier 3 - Recommended)
```bash
mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

### Fallback Build (Tier 2 - No Templates)
```bash
mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON
make -j$(nproc)
```

## Dependencies

### Tier 1 (Default)
- **None** - Always available

### Tier 2 & 3 (High Precision)
- **Boost headers** (~50 MB, header-only)
- **MPFR library** (~1 MB, optional, can build from source)
- **GMP library** (~1 MB, optional, can build from source)
- **Fallback**: cpp_dec_float (header-only, no MPFR/GMP needed)

### Optional
- **libquadmath** (for __float128 support)

## Precision Comparison

| Type | Bits | Decimal Digits | Speed | Tier |
|------|------|----------------|-------|------|
| `double` | 64 | 15-17 | ⚡⚡⚡⚡⚡ | 1, 2, 3 |
| `__float128` | 128 | 33-36 | ⚡⚡⚡⚡ | 3 (optional) |
| `mpreal(128)` | 128 | 38 | ⚡⚡⚡ | 2, 3 |
| `mpreal(256)` | 256 | 77 | ⚡⚡ | 2, 3 |
| `mpreal(512)` | 512 | 154 | ⚡ | 2, 3 |
| `mpreal(1024)` | 1024 | 308 | 🐌 | 2, 3 |

## API Examples

### Tier 1: Double Precision
```cpp
#include "polynomial.h"
#include "result_refiner.h"

Polynomial poly(degrees, coeffs);
RefinedRoot result = ResultRefiner::refineRoot(poly, initial_guess);
// Precision: ~15-17 decimal digits
```

### Tier 2: Fallback High Precision
```cpp
#include "polynomial_hp.h"
#include "result_refiner_hp.h"

PolynomialHP poly_hp(degrees, coeffs_hp);
RefinedRootHP result_hp = ResultRefinerHP::refineRoot(poly_hp, initial_guess_hp);
// Precision: Fixed at compile time (e.g., 100 decimal digits)
```

### Tier 3: Template High Precision
```cpp
#include "polynomial_mp.h"
#include "result_refiner_mp.h"

// Set precision at runtime
mpreal::set_default_prec(256);  // 77 decimal digits

PolynomialMP poly_mp(degrees, coeffs_mp);
RefinedRootMP result_mp = ResultRefinerMP::refineRoot(poly_mp, initial_guess_mp);

// Change precision dynamically
mpreal::set_default_prec(512);  // 154 decimal digits
// ... continue with higher precision
```

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
- ✅ You're solving ill-conditioned problems

## Build Configuration

See [SETUP.md](../SETUP.md) for detailed build instructions including:
- Installing dependencies
- Building without admin rights
- Specifying custom library paths
- Using static vs dynamic linking

## Implementation Details

See [REFACTORIZATION_PLAN.md](REFACTORIZATION_PLAN.md) for the complete implementation plan.

See [FILE_STRUCTURE.md](FILE_STRUCTURE.md) for the file organization.

## Minimal Dependencies Setup

See [HIGH_PRECISION_MINIMAL_DEPENDENCIES.md](HIGH_PRECISION_MINIMAL_DEPENDENCIES.md) for:
- Building MPFR/GMP from source without admin rights
- Using header-only cpp_dec_float backend
- Minimal Boost header installation

## Multiple Precision Types

See [HIGH_PRECISION_MULTIPLE_TYPES.md](HIGH_PRECISION_MULTIPLE_TYPES.md) for:
- Supporting double, __float128, and mpreal simultaneously
- Using templates to avoid code duplication
- Flexible precision with mpreal

## Setup Guide

See [HIGH_PRECISION_SETUP_GUIDE.md](HIGH_PRECISION_SETUP_GUIDE.md) for:
- Step-by-step setup on dedicated servers
- Building without admin rights
- Troubleshooting common issues

## Quick Reference

See [QUICK_REFERENCE.md](QUICK_REFERENCE.md) for:
- Build commands
- API quick reference
- Common workflows
- Performance tips

