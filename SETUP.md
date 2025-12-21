# Complete Setup Guide

This guide provides comprehensive instructions for setting up the Polynomial Solver project in any environment.

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Installation Steps](#installation-steps)
3. [Verification](#verification)
4. [Troubleshooting](#troubleshooting)

## System Requirements

### Required Software

| Software | Minimum Version | Recommended | Purpose |
|----------|----------------|-------------|---------|
| GCC/Clang | 4.8+ / 3.4+ | 13.2+ / 15+ | C++11 compiler |
| CMake | 3.15 | 3.20+ | Build system |
| Make | 4.0+ | 4.3+ | Build automation |

### Optional Software

| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.11+ | Visualization |
| NumPy | Latest | Numerical computation |
| Matplotlib | Latest | Plotting |
| UV | Latest | Python package manager |

### Supported Platforms

- **Linux**: Ubuntu 20.04+, Fedora 30+, Arch Linux, Debian 10+
- **macOS**: 10.14+ (with Xcode Command Line Tools)
- **Windows**: WSL2 (Ubuntu 20.04+) or MinGW-w64

## Installation Steps

### Step 1: Install Prerequisites

```bash
# Ubuntu/Debian
sudo apt-get update && sudo apt-get install -y build-essential cmake

# Fedora/RHEL/CentOS
sudo dnf install -y gcc-c++ cmake make

# Arch Linux
sudo pacman -S base-devel cmake

# macOS
xcode-select --install && brew install cmake
```

### Step 2: Extract and Build

```bash
tar -xzf polynomial-solver.tar.gz
cd polynomial-solver
./check_prerequisites.sh        # Verify prerequisites
mkdir build && cd build
cmake .. && make -j$(nproc)
ctest                           # Run tests
```

For high-precision support, see [High-Precision Build Options](#high-precision-build-options).

## Verification

### Test the Installation

```bash
# Run all tests
cd build
ctest

# Expected output:
# 100% tests passed, 0 tests failed out of 21
```

### Run Standard Verification Workflow

The circle-ellipse intersection serves as the standard verification:

```bash
# 1. Run the test with all 3 strategies
./build/bin/test_strategies

# 2. Set up Python environment (one-time)
uv venv .venv
source .venv/bin/activate
uv pip install numpy matplotlib

# 3. Visualize results
python examples/visualize_circle_ellipse.py
```

**Expected Results:**
- All 21 tests pass
- 3 geometry dumps generated in `dumps/`
- Visualizations generated in `visualizations/`
- All strategies converge to root: (0.894427, 0.447214)

### Run Other Example Programs

```bash
# Method comparison
./build/bin/test_method_comparison

# Linear function test (machine epsilon precision)
./build/bin/test_linear_graph_hull

# ProjectedPolyhedral method test
./build/bin/test_projected_polyhedral
```

### Verify Build Artifacts

```bash
# List built executables
ls -lh build/bin/

# List built libraries
ls -lh build/lib/
```

Expected executables:
- `example_*` - Example programs (4 examples)
- `test_*` - Test executables (21 tests)

Expected libraries:
- `libpolynomial_solver.a` - Unified library (link against this)
- `libpolynomial.a` - Polynomial module
- `libsolver.a` - Solver module
- `libresult_refiner.a` - Result refiner module
- `libdifferentiation.a` - Differentiation module
- `libgeometry.a` - Geometry module
- `libde_casteljau.a` - De Casteljau module

## High-Precision Build Options

### Install Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install libboost-dev libmpfr-dev libgmp-dev

# Fedora/RHEL
sudo dnf install boost-devel mpfr-devel gmp-devel

# macOS
brew install boost mpfr gmp
```

### Build

```bash
mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

### Verify

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON 2>&1 | grep -i "backend"
# Expected: "-- Using MPFR backend" or "-- Using cpp_dec_float backend"
```

For advanced options (custom paths, building from source, backend selection), see [docs/HIGH_PRECISION.md](docs/HIGH_PRECISION.md).

## Troubleshooting

| Issue | Solution |
|-------|----------|
| CMake too old | Download from https://cmake.org/download/ |
| C++11 not supported | Update compiler: `sudo apt-get install gcc-9 g++-9` |
| Linking errors | Clean rebuild: `rm -rf build && mkdir build && cd build && cmake .. && make` |
| Tests fail | Run verbose: `ctest --output-on-failure --verbose` |
| Python visualization | Install: `pip3 install numpy matplotlib` |

## Summary

After following this guide, you should have:

1. ✅ All prerequisites installed
2. ✅ Project extracted and built
3. ✅ All tests passing
4. ✅ Example programs running

For quick reference, see [README.md](README.md).

For algorithm details, see documentation in [docs/](docs/).

## Support

- **Email**: wenyd@lsec.cc.ac.cn


