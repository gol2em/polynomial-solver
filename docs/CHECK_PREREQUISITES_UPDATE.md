# check_prerequisites.sh Update

## Summary

Updated `check_prerequisites.sh` to align with the new CMake-based high-precision build system.

## Question

**User asked**: "Check if check_prerequisite is still needed."

## Answer

**Yes, the script is still useful and has been updated!**

## What Changed

### 1. Added High-Precision Library Detection

**New function**: `check_high_precision()`

Checks for:
- ✅ **Boost headers** (`/usr/include/boost`, `/usr/local/include/boost`)
- ✅ **MPFR library** (via `ldconfig -p` or direct file check)
- ✅ **GMP library** (via `ldconfig -p` or direct file check)
- ✅ **Quadmath library** (via `ldconfig -p` or direct file check)

**Intelligent backend detection**:
- If Boost + MPFR + GMP found → "MPFR backend - fastest"
- If Boost only found → "cpp_dec_float backend"
- If quadmath only found → "quadmath backend"
- If none found → Shows installation instructions

### 2. Updated Build Instructions

**Before**: Referenced deleted `build.sh` script

**After**: Shows modern CMake-based workflow:
```bash
# Standard build (double precision)
mkdir -p build && cd build
cmake ..
make -j$(nproc)

# High-precision build (if libraries available)
mkdir -p build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)

# Run tests
cd build && ctest --output-on-failure
```

### 3. Enhanced Error Messages

Added high-precision installation instructions to the failure message:
```bash
For high-precision support (optional):
  Ubuntu/Debian: sudo apt-get install libboost-dev libmpfr-dev libgmp-dev
  Fedora/RHEL:   sudo dnf install boost-devel mpfr-devel gmp-devel
  macOS:         brew install boost mpfr gmp
```

## Why Keep the Script?

### Benefits

1. **Quick Pre-Check**: Verify environment before running CMake
2. **System Overview**: Shows OS, architecture, kernel, distribution
3. **Comprehensive Check**: Checks all tools in one place:
   - C++ compiler (g++)
   - Build system (cmake, make)
   - Version control (git)
   - Python tools (python3, numpy, matplotlib, uv)
   - High-precision libraries (Boost, MPFR, GMP, quadmath)
4. **User-Friendly**: Color-coded output, clear summary
5. **Helpful Guidance**: Shows next steps and installation commands

### Complementary to CMake

- **Script**: Quick overview, checks all tools
- **CMake**: Detailed detection, configuration, error messages

They work together:
1. Run `./check_prerequisites.sh` - Get overview
2. Run `cmake ..` - Get detailed configuration
3. If issues, both provide helpful error messages

## Example Output

```
========================================
Polynomial Solver - Prerequisite Check
========================================

System Information...
✓ OS: Linux
✓ Architecture: x86_64
✓ Kernel: 6.6.87.2-microsoft-standard-WSL2
✓ Distribution: Ubuntu 22.04.5 LTS

Checking C++ Compiler...
✓ GCC C++ Compiler found: g++ (GCC) 13.2.0
  → GCC version 13.2.0 supports C++11

Checking Build System...
✓ CMake found: cmake version 4.1.3
  → CMake version 4.1.3 meets minimum requirement (3.15)
✓ GNU Make found: GNU Make 4.3

Checking Python (optional for visualization)...
✓ Python 3 found: Python 3.10.12
  → NumPy is installed
  → Matplotlib not found (optional for visualization)
✓ UV (Python package manager) found: uv 0.9.9

Checking Version Control...
✓ Git found: git version 2.34.1

Checking High-Precision Libraries (optional)...
✓ Boost headers found
✓ MPFR library found
✓ GMP library found
✓ Quadmath library found
  → High-precision support available (MPFR backend - fastest)
  → Enable with: cmake .. -DENABLE_HIGH_PRECISION=ON

========================================
Summary
========================================
Passed:   10
Warnings: 0
Failed:   0

✓ All required prerequisites are met!
  You can proceed with building the project.

Next steps:
  # Standard build (double precision)
  mkdir -p build && cd build
  cmake ..
  make -j$(nproc)

  # High-precision build (if libraries available)
  mkdir -p build && cd build
  cmake .. -DENABLE_HIGH_PRECISION=ON
  make -j$(nproc)

  # Run tests
  cd build && ctest --output-on-failure

See SETUP.md for more build options and configurations.
```

## Conclusion

✅ **Keep the script** - It's still useful and now better than ever!

The script provides a quick, user-friendly overview that complements CMake's detailed configuration.

