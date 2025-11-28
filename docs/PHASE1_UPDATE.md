# Phase 1 Update: Enhanced Dependency Detection

## Summary

Updated the high-precision build system based on user feedback to improve flexibility and usability.

## Changes Made

### 1. Auto-Detection of High-Precision Libraries

**Before**: High-precision was always OFF by default, user had to explicitly enable it.

**After**: 
- CMake always checks for high-precision libraries
- If libraries are detected, shows a helpful message suggesting how to enable
- User can still explicitly enable/disable with `-DENABLE_HIGH_PRECISION=ON/OFF`

**Example Output**:
```
High-precision libraries detected!
  To enable: cmake .. -DENABLE_HIGH_PRECISION=ON
  (Currently disabled by default)
```

### 2. Quadmath as Standalone Backend

**Before**: Quadmath was only an optional add-on to Boost-based backends.

**After**:
- Quadmath can now be used as a standalone high-precision backend
- If Boost is not available but quadmath is, quadmath becomes the backend
- Fallback version (Tier 2) now supports quadmath
- Three possible backends:
  1. **MPFR** (Boost + MPFR + GMP) - Fastest, arbitrary precision
  2. **CPP_DEC_FLOAT** (Boost only) - Good, fixed precision (50/100 digits)
  3. **QUADMATH** (quadmath only) - Very good, fixed precision (33-36 digits)

### 3. Graceful Degradation

**Before**: If high-precision was enabled but libraries not found, CMake would error out.

**After**:
- If `-DENABLE_HIGH_PRECISION=ON` but no libraries found:
  - Shows a **warning** (not error)
  - Automatically disables high-precision
  - Continues with double precision
  - Provides clear installation instructions

**Example Output**:
```
CMake Warning:
  High-precision support requested but no suitable libraries found.
  Boost or quadmath is required.
  Disabling high-precision support and continuing with double precision.
  
  To enable high-precision, install dependencies:
    Ubuntu/Debian: sudo apt-get install libboost-dev libmpfr-dev libgmp-dev
    ...
  Or specify custom paths:
    -DBOOST_ROOT=/path/to/boost
    -DMPFR_ROOT=/path/to/mpfr -DGMP_ROOT=/path/to/gmp
```

### 4. Improved Quadmath Detection

**Before**: Quadmath was not reliably detected on Linux systems.

**After**:
- Searches standard GCC library locations
- Checks for `.so.0` versioned libraries
- Manually checks common paths if find_library fails
- Supports both x86_64 and aarch64 architectures
- Supports GCC versions 11, 12, 13

**Search Locations**:
- `/usr/lib/x86_64-linux-gnu/libquadmath.so.0`
- `/usr/lib/gcc/x86_64-linux-gnu/{11,12,13}/libquadmath.so`
- And more...

### 5. Enhanced Configuration Summary

**Before**: Basic summary showing only backend type.

**After**: Comprehensive summary showing:
- Status (Available / Not available)
- Backend type (MPFR / CPP_DEC_FLOAT / QUADMATH)
- Which libraries were found
- Performance characteristics
- Precision capabilities
- Installation instructions if not available

**Example Output**:
```
========================================
High-Precision Configuration Summary
========================================
  Status:           Available
  Backend:          MPFR
  Boost headers:    Found
  MPFR library:     Found
  GMP library:      Found
  Performance:      Fastest
  Precision:        Arbitrary (runtime configurable)
  Quadmath:         Available (can be used with templates)
========================================
```

### 6. Comprehensive Documentation for Custom Paths

Added detailed section to `SETUP.md` covering:

**Method 1**: Individual library paths (recommended)
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DBOOST_ROOT=$HOME/local \
  -DMPFR_ROOT=$HOME/local \
  -DGMP_ROOT=$HOME/local
```

**Method 2**: CMake prefix path
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DCMAKE_PREFIX_PATH="$HOME/local;/opt/custom"
```

**Method 3**: Environment variables
```bash
export CMAKE_PREFIX_PATH=$HOME/local
cmake .. -DENABLE_HIGH_PRECISION=ON
```

**Method 4**: Project-local libraries
```bash
# Copy to external/ directory
cp -r /path/to/boost external/
cmake .. -DENABLE_HIGH_PRECISION=ON
```

**Troubleshooting**: 4 solutions for common path issues

### 7. Updated Configuration Header

Added `USE_QUADMATH_BACKEND` definition to `config.h.in`:
```cpp
// Multiprecision backend (mutually exclusive)
#cmakedefine USE_MPFR_BACKEND
#cmakedefine USE_CPP_DEC_FLOAT_BACKEND
#cmakedefine USE_QUADMATH_BACKEND

// Quadmath support (can be used alongside other backends in template mode)
#cmakedefine ENABLE_QUADMATH
```

## Testing

### Test 1: Default Build
```bash
cmake ..
```
✅ Shows detection summary and suggestion to enable if libraries found
✅ Continues with double precision

### Test 2: Enable with Libraries Available
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON
```
✅ Detects MPFR backend
✅ Detects quadmath
✅ Shows comprehensive summary
✅ Enables template support

### Test 3: Enable without Libraries
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DBOOST_ROOT=/nonexistent
```
✅ Shows warning (not error)
✅ Provides installation instructions
✅ Continues with double precision

### Test 4: Quadmath Detection
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON
```
✅ Detects quadmath in `/usr/lib/x86_64-linux-gnu/libquadmath.so.0`
✅ Shows quadmath as available

## Files Modified

- `CMakeLists.txt` - Auto-detection, graceful degradation
- `cmake/FindHighPrecision.cmake` - Quadmath backend, improved detection
- `include/config.h.in` - Added USE_QUADMATH_BACKEND
- `SETUP.md` - Comprehensive custom path documentation

## Benefits

1. **Better User Experience**: Auto-detection with helpful messages
2. **More Flexible**: Quadmath as standalone option
3. **More Robust**: Graceful degradation instead of errors
4. **Better Documentation**: Comprehensive path configuration guide
5. **Wider Compatibility**: Works on more systems out of the box

## Next Steps

Continue with Phase 2: Type Definition Headers

