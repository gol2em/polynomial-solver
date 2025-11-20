# Complete Setup Guide

This guide provides comprehensive instructions for setting up the Polynomial Solver project in any environment.

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Installation Steps](#installation-steps)
3. [Verification](#verification)
4. [Troubleshooting](#troubleshooting)
5. [Development Setup](#development-setup)

## System Requirements

### Required Software

| Software | Minimum Version | Recommended | Purpose |
|----------|----------------|-------------|---------|
| GCC/Clang | 4.8+ / 3.4+ | 13.2+ / 15+ | C++11 compiler |
| CMake | 3.15 | 3.20+ | Build system |
| Make | 4.0+ | 4.3+ | Build automation |
| Git | 2.0+ | 2.34+ | Version control |

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

#### Ubuntu/Debian

```bash
# Update package list
sudo apt-get update

# Install required packages
sudo apt-get install -y build-essential cmake git

# Optional: Install Python tools (for visualization)
curl -LsSf https://astral.sh/uv/install.sh | sh
# Then use: uv venv .venv && source .venv/bin/activate && uv pip install numpy matplotlib
```

#### Fedora/RHEL/CentOS

```bash
# Install required packages
sudo dnf install -y gcc-c++ cmake git make

# Optional: Install Python tools (for visualization)
curl -LsSf https://astral.sh/uv/install.sh | sh
```

#### Arch Linux

```bash
# Install required packages
sudo pacman -S base-devel cmake git

# Optional: Install Python tools (for visualization)
curl -LsSf https://astral.sh/uv/install.sh | sh
```

#### macOS

```bash
# Install Xcode Command Line Tools
xcode-select --install

# Install Homebrew (if not already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install required packages
brew install cmake git

# Optional: Install Python tools (for visualization)
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Step 2: Clone Repository

```bash
# Clone the repository
git clone https://github.com/gol2em/polynomial-solver.git

# Navigate to project directory
cd polynomial-solver
```

### Step 3: Verify Prerequisites

```bash
# Run the prerequisite checker
./check_prerequisites.sh
```

Expected output:
```
✓ All required prerequisites are met!
```

If any prerequisites are missing, the script will provide installation hints.

### Step 4: Build Project

```bash
# Build with tests
./build.sh --test
```

Or manually:

```bash
# Create build directory
mkdir -p build
cd build

# Configure
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build
make -j$(nproc)

# Run tests
ctest --output-on-failure
```

## Verification

### Test the Installation

```bash
# Run all tests
cd build
ctest

# Expected output:
# 100% tests passed, 0 tests failed out of 12
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
- All 12 tests pass
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
- `polynomial_solver_app` - Main application
- `test_*` - Test executables (12 tests)
- `test_strategies` - Strategy comparison with geometry dumps
- `example_circle_ellipse` - Circle-ellipse intersection example

Expected libraries:
- `libpolynomial.a` - Polynomial module
- `libsolver.a` - Solver module
- `libgeometry.a` - Geometry module
- `libde_casteljau.a` - De Casteljau module

## Troubleshooting

### Issue: CMake version too old

**Error**: `CMake 3.15 or higher is required`

**Solution**:
```bash
# Ubuntu/Debian: Install from Kitware repository
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc | sudo apt-key add -
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ focal main'
sudo apt-get update
sudo apt-get install cmake

# Or download from https://cmake.org/download/
```

### Issue: Compiler doesn't support C++11

**Error**: `error: unrecognized command line option '-std=c++11'`

**Solution**:
```bash
# Update GCC
sudo apt-get install gcc-9 g++-9
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 90
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 90
```

### Issue: Build fails with linking errors

**Solution**:
```bash
# Clean build and rebuild
./build.sh --clean
```

### Issue: Tests fail

**Solution**:
```bash
# Run tests with verbose output
cd build
ctest --output-on-failure --verbose

# Run specific test
./bin/test_polynomial_conversion
```

### Issue: Python visualization doesn't work

**Solution**:
```bash
# Install Python dependencies
pip3 install numpy matplotlib

# Or use UV
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync
source .venv/bin/activate
```

## Development Setup

### Setting Up Development Environment

```bash
# Clone repository
git clone https://github.com/gol2em/polynomial-solver.git
cd polynomial-solver

# Configure Git
git config user.name "gol2em"
git config user.email "wenyd@lsec.cc.ac.cn"

# Build in debug mode
./build.sh --debug --test
```

### Project Layout

```
polynomial-solver/
├── include/              # Public headers
├── src/                  # Implementation
├── tests/                # Test suite
├── python/               # Visualization tools
├── docs/                 # Documentation
├── build/                # Build artifacts (generated)
└── .venv/                # Python virtual env (optional)
```

### Build Options

```bash
# Release build (optimized)
./build.sh

# Debug build (with symbols)
./build.sh --debug

# Clean build
./build.sh --clean

# Build and test
./build.sh --test

# Verbose build
./build.sh --verbose

# Custom parallel jobs
./build.sh -j 8
```

### Running Individual Tests

```bash
cd build/bin

# Polynomial tests
./test_polynomial_conversion
./test_polynomial_graph
./test_polynomial_system

# Solver tests
./test_solver_subdivision
./test_linear_graph_hull
./test_projected_polyhedral

# Geometry tests
./test_geometry_convex
./test_2d_convex_hull_robust
./test_2d_hyperplane_intersection
./test_2d_hull_vertex_order

# Degeneracy tests
./test_degenerate_boxes

# Comparison tool
./test_method_comparison
```

## Moving to Another Environment

### Export Project

```bash
# Create a tarball
cd /path/to/parent/directory
tar -czf polynomial-solver.tar.gz polynomial-solver/

# Or use Git
cd polynomial-solver
git bundle create polynomial-solver.bundle --all
```

### Import to New Environment

```bash
# From tarball
tar -xzf polynomial-solver.tar.gz
cd polynomial-solver

# From Git bundle
git clone polynomial-solver.bundle polynomial-solver
cd polynomial-solver

# Verify and build
./check_prerequisites.sh
./build.sh --test
```

### Minimal Files Needed

If you want to minimize the transfer size, you only need:

**Required**:
- `include/` - Header files
- `src/` - Source files
- `tests/` - Test suite
- `CMakeLists.txt` - Build configuration
- `tests/CMakeLists.txt` - Test configuration
- `build.sh` - Build script
- `check_prerequisites.sh` - Prerequisite checker

**Optional**:
- `python/` - Visualization tools
- `docs/` - Documentation
- `README.md`, `QUICKSTART.md`, `SETUP.md` - Documentation
- `.git/` - Git history

**Not needed** (will be regenerated):
- `build/` - Build artifacts
- `.venv/` - Python virtual environment
- `*.png`, `*.csv` - Generated files

## Summary

After following this guide, you should have:

1. ✅ All prerequisites installed
2. ✅ Project cloned and built
3. ✅ All 11 tests passing
4. ✅ Example programs running
5. ✅ Ready for development or deployment

For quick reference, see [QUICKSTART.md](QUICKSTART.md).

For algorithm details, see documentation in [docs/](docs/).

## Support

- **Email**: wenyd@lsec.cc.ac.cn
- **GitHub**: https://github.com/gol2em/polynomial-solver
- **Issues**: https://github.com/gol2em/polynomial-solver/issues


