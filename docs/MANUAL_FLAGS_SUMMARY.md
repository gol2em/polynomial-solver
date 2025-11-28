# Manual Library Control - Implementation Summary

## Overview

Implemented fine-grained control over multiprecision library usage through CMake flags.

## What Was Implemented

### 1. New CMake Flags

Added four new option flags in `CMakeLists.txt`:

```cmake
option(USE_BOOST "Use Boost multiprecision library" ON)
option(USE_MPFR "Use MPFR library (requires Boost)" ON)
option(USE_GMP "Use GMP library (requires Boost and MPFR)" ON)
option(USE_QUADMATH "Use quadmath library for __float128 support" ON)
```

**All flags default to ON** - auto-detect and use if available.

### 2. Updated Detection Logic

Modified `cmake/FindHighPrecision.cmake` to:

- **Show user preferences** at the start of detection
- **Skip detection** if user disabled a library
- **Show clear messages** indicating when libraries are disabled by user vs not found
- **Respect dependencies** (e.g., MPFR requires Boost)

### 3. Enhanced CMake Output

The build system now shows:

```
User preferences:
  USE_BOOST:     ON/OFF
  USE_MPFR:      ON/OFF
  USE_GMP:       ON/OFF
  USE_QUADMATH:  ON/OFF
```

And detection results:
- `✓` - Found and enabled
- `✗` - Not found
- `⊗` - Disabled by user

### 4. Backend Selection Logic

```
Priority 1: MPFR backend
  - Requires: USE_BOOST=ON, USE_MPFR=ON, USE_GMP=ON
  - All libraries must be found
  - Result: Fastest, arbitrary precision

Priority 2: cpp_dec_float backend
  - Requires: USE_BOOST=ON
  - MPFR/GMP disabled or not found
  - Result: Header-only, fixed precision

Priority 3: quadmath backend
  - Requires: USE_QUADMATH=ON
  - Boost disabled or not found
  - Result: Standalone, fixed precision

Priority 4: Graceful degradation
  - All libraries disabled or not found
  - Result: Double precision with warning
```

## Testing Results

### Test 1: All Enabled (Default)

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON
```

**Result**: ✅ MPFR backend selected
- All libraries detected and used
- Fastest performance
- Arbitrary precision

### Test 2: Disable MPFR/GMP

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF
```

**Result**: ✅ cpp_dec_float backend selected
- Boost headers used
- MPFR/GMP explicitly disabled
- Header-only, no library linking
- Fixed precision (50/100 digits)

### Test 3: Disable Boost

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF
```

**Result**: ✅ quadmath backend selected
- Boost explicitly disabled
- Quadmath detected and used
- Standalone backend
- Fixed precision (33-36 digits)

### Test 4: Disable All

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF -DUSE_QUADMATH=OFF
```

**Result**: ✅ Graceful degradation
- All libraries disabled
- Warning message shown
- Falls back to double precision
- Build continues successfully

## Use Cases

### 1. Header-Only Build

**Need**: High precision without linking libraries

**Solution**:
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF
```

### 2. Minimal Dependencies

**Need**: Smallest dependency footprint

**Solution**:
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF
```

### 3. Force Specific Backend

**Need**: Test specific backend performance

**Solution**:
```bash
# Force cpp_dec_float
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF

# Force quadmath
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF

# Force MPFR (or fail gracefully)
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=ON -DUSE_MPFR=ON -DUSE_GMP=ON
```

### 4. Avoid Library Conflicts

**Need**: System has conflicting Boost version

**Solution**:
```bash
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF
```

## Documentation

### Files Updated

1. **CMakeLists.txt** - Added USE_* option flags
2. **cmake/FindHighPrecision.cmake** - Respect flags, enhanced output
3. **SETUP.md** - Added "Manual Library Control" section
4. **docs/MANUAL_LIBRARY_CONTROL.md** - Comprehensive guide with examples

### Quick Reference

| Want | Command |
|------|---------|
| MPFR backend | `cmake .. -DENABLE_HIGH_PRECISION=ON` |
| cpp_dec_float | `cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF` |
| quadmath only | `cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF` |
| No quadmath | `cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_QUADMATH=OFF` |

## Benefits

1. **Flexibility**: Users control exactly which libraries to use
2. **Testing**: Easy to test different backends
3. **Deployment**: Minimize dependencies for production
4. **Compatibility**: Avoid conflicts with existing libraries
5. **Transparency**: Clear output shows what's enabled/disabled
6. **Graceful**: Never fails, always provides working build

## Commits

- **9f5ec18**: Updated check_prerequisites.sh
- **ec83b04**: Added manual library control flags (current)

Both pushed to gitee and github.

