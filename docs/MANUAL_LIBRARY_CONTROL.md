# Manual Library Control

## Overview

The build system now provides fine-grained control over which multiprecision libraries to use. This allows users to:

- **Disable specific libraries** they don't want to use
- **Force specific backends** for testing or compatibility
- **Minimize dependencies** for deployment
- **Avoid conflicts** with existing library versions

## Available Flags

All flags default to `ON` (enabled):

| Flag | Description | Default | Dependencies |
|------|-------------|---------|--------------|
| `USE_BOOST` | Enable Boost multiprecision | ON | None |
| `USE_MPFR` | Enable MPFR backend | ON | Boost |
| `USE_GMP` | Enable GMP library | ON | Boost, MPFR |
| `USE_QUADMATH` | Enable quadmath support | ON | None |

## Backend Selection Logic

The build system automatically selects the best available backend based on enabled libraries:

```
1. If USE_BOOST=ON and USE_MPFR=ON and USE_GMP=ON
   → MPFR backend (fastest, arbitrary precision)

2. If USE_BOOST=ON but MPFR/GMP disabled or not found
   → cpp_dec_float backend (header-only, fixed precision)

3. If USE_BOOST=OFF and USE_QUADMATH=ON
   → quadmath backend (standalone, fixed precision)

4. If all disabled or not found
   → Graceful degradation to double precision
```

## Usage Examples

### Example 1: Use Only Boost (cpp_dec_float)

Disable MPFR/GMP to use header-only cpp_dec_float backend:

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DUSE_MPFR=OFF \
  -DUSE_GMP=OFF
```

**Result**:
- Backend: cpp_dec_float
- Dependencies: Boost headers only (no libraries)
- Precision: Fixed (50 or 100 decimal digits)
- Performance: Good (2-5× slower than MPFR)

### Example 2: Use Only Quadmath

Disable Boost to use standalone quadmath:

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DUSE_BOOST=OFF
```

**Result**:
- Backend: quadmath
- Dependencies: libquadmath (usually included with GCC)
- Precision: Fixed (33-36 decimal digits)
- Performance: Very good

### Example 3: Disable Quadmath

Use Boost but disable quadmath support:

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DUSE_QUADMATH=OFF
```

**Result**:
- Backend: MPFR or cpp_dec_float (depending on MPFR/GMP availability)
- No __float128 support in templates
- Useful if quadmath causes compatibility issues

### Example 4: Force MPFR Backend

Ensure MPFR backend is used (fail if not available):

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DUSE_BOOST=ON \
  -DUSE_MPFR=ON \
  -DUSE_GMP=ON \
  -DUSE_QUADMATH=OFF
```

If MPFR/GMP not found, build will gracefully degrade with warning.

### Example 5: Minimal Dependencies

Use only what's available, prefer quadmath:

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DUSE_BOOST=OFF \
  -DUSE_QUADMATH=ON
```

**Result**:
- Smallest dependency footprint
- Only requires libquadmath (usually pre-installed)

## CMake Output

The build system shows your preferences and the detection results:

```
========================================
Detecting High-Precision Dependencies
========================================

User preferences:
  USE_BOOST:     ON
  USE_MPFR:      OFF
  USE_GMP:       OFF
  USE_QUADMATH:  ON

Step 1: Checking for Boost headers...
  ✓ Boost found: /usr/include

Step 2: Checking for MPFR and GMP libraries...
  ⊗ MPFR disabled by user (USE_MPFR=OFF)
  ⊗ GMP disabled by user (USE_GMP=OFF)

Step 3: Checking for quadmath library...
  ✓ Quadmath library found: /usr/lib/x86_64-linux-gnu/libquadmath.so.0

========================================
High-Precision Configuration Summary
========================================
  Status:           Available
  Backend:          CPP_DEC_FLOAT
  Boost headers:    Found
  MPFR/GMP:         Disabled by user (using cpp_dec_float)
  Performance:      Good (2-5× slower than MPFR)
  Precision:        Fixed (50, 100 decimal digits)
  Quadmath:         Available (can be used with templates)
========================================
```

## Use Cases

### 1. Deployment with Minimal Dependencies

**Scenario**: Deploy to servers where you can't install Boost libraries

**Solution**: Use quadmath only
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF
```

### 2. Header-Only Build

**Scenario**: Want high precision but can't link against libraries

**Solution**: Use cpp_dec_float backend
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF
```

### 3. Testing Different Backends

**Scenario**: Compare performance of different backends

**Solution**: Build multiple times with different flags
```bash
# Test MPFR
cmake .. -DENABLE_HIGH_PRECISION=ON

# Test cpp_dec_float
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF

# Test quadmath
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF
```

