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

#### Basic Build (Double Precision Only)

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

#### Build with High-Precision Support

See [High-Precision Build Options](#high-precision-build-options) section below for detailed instructions.

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

The polynomial solver supports three tiers of precision. See [docs/HIGH_PRECISION.md](docs/HIGH_PRECISION.md) for complete details.

### Tier 1: Default (Double Precision Only)

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**Features**:
- ✅ Double precision (15-17 decimal digits)
- ✅ Zero dependencies
- ✅ Maximum performance

### Tier 2: High-Precision Version (Recommended)

#### Prerequisites

Install dependencies:

```bash
# Ubuntu/Debian
sudo apt-get install libboost-dev libmpfr-dev libgmp-dev

# Fedora/RHEL
sudo dnf install boost-devel mpfr-devel gmp-devel

# macOS
brew install boost mpfr gmp
```

#### Build

```bash
mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

**Features**:
- ✅ Fixed high-precision types (PolynomialHP, ResultRefinerHP)
- ✅ Runtime-configurable precision (MPFR backend)
- ✅ Supports MPFR, cpp_dec_float, __float128 backends
- ✅ Production-ready (Tier 2 complete)

### Manual Library Control

You can explicitly control which multiprecision libraries to use:

```bash
# Use only Boost (cpp_dec_float backend, no MPFR/GMP)
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF

# Use only quadmath (no Boost)
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF

# Disable quadmath even if available
cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_QUADMATH=OFF

# Custom combination
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DUSE_BOOST=ON \
  -DUSE_MPFR=OFF \
  -DUSE_GMP=OFF \
  -DUSE_QUADMATH=ON
```

**Available flags**:
- `USE_BOOST=ON/OFF` - Enable/disable Boost multiprecision (default: ON)
- `USE_MPFR=ON/OFF` - Enable/disable MPFR backend (default: ON, requires Boost)
- `USE_GMP=ON/OFF` - Enable/disable GMP (default: ON, requires Boost and MPFR)
- `USE_QUADMATH=ON/OFF` - Enable/disable quadmath support (default: ON)

**Backend selection**:
- If `USE_BOOST=ON` and `USE_MPFR=ON` and `USE_GMP=ON` → MPFR backend (fastest)
- If `USE_BOOST=ON` but MPFR/GMP disabled → cpp_dec_float backend
- If `USE_BOOST=OFF` and `USE_QUADMATH=ON` → quadmath backend (standalone)

### Future: Tier 3 (Template-Based, Not Yet Implemented)

Tier 3 will provide template-based high-precision support for maximum flexibility.
Currently, Tier 2 (fixed high-precision types) is the recommended approach.

To enable templates in the future (when implemented):

```bash
mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON -DENABLE_TEMPLATES=ON
make -j$(nproc)
```

**Planned Features**:
- ⏳ Template-based (supports any numeric type)
- ⏳ Flexible precision (change at compile time)
- ⏳ No code duplication

### Custom Library Paths

If libraries are not installed globally (e.g., on dedicated servers without admin rights), you can specify custom paths.

#### Method 1: Individual Library Paths (Recommended)

Specify the root directory for each library:

```bash
mkdir build && cd build
cmake .. \
  -DENABLE_HIGH_PRECISION=ON \
  -DBOOST_ROOT=$HOME/local \
  -DGMP_ROOT=$HOME/local \
  -DMPFR_ROOT=$HOME/local
make -j$(nproc)
```

**Expected directory structure**:
```
$HOME/local/
├── include/
│   ├── boost/           # Boost headers
│   ├── gmp.h            # GMP header
│   └── mpfr.h           # MPFR header
└── lib/
    ├── libgmp.a         # GMP library
    └── libmpfr.a        # MPFR library
```

#### Method 2: CMake Prefix Path

Use a common prefix for all libraries:

```bash
cmake .. \
  -DENABLE_HIGH_PRECISION=ON \
  -DCMAKE_PREFIX_PATH="$HOME/local;/opt/custom"
make -j$(nproc)
```

**Note**: CMake will search `<prefix>/include` and `<prefix>/lib` for each prefix.

#### Method 3: Environment Variables

Set environment variables before running CMake:

```bash
export CMAKE_PREFIX_PATH=$HOME/local
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

**Tip**: Add these to your `~/.bashrc` or `~/.bash_profile` for persistence.

#### Method 4: Project-Local Libraries

Copy libraries directly to your project:

```bash
# Create external directory
mkdir -p external

# Copy Boost headers
cp -r /path/to/boost external/

# CMake will automatically detect external/boost/
cmake .. -DENABLE_HIGH_PRECISION=ON
```

**Expected structure**:
```
polynomial-solver/
├── external/
│   └── boost/
│       └── version.hpp
├── include/
├── src/
└── CMakeLists.txt
```

#### Troubleshooting Custom Paths

**Problem**: CMake can't find libraries even with custom paths

**Solution 1**: Verify directory structure
```bash
# Check if headers exist
ls $HOME/local/include/boost/version.hpp
ls $HOME/local/include/mpfr.h
ls $HOME/local/include/gmp.h

# Check if libraries exist
ls $HOME/local/lib/libmpfr.*
ls $HOME/local/lib/libgmp.*
```

**Solution 2**: Use absolute paths
```bash
cmake .. \
  -DENABLE_HIGH_PRECISION=ON \
  -DBOOST_ROOT=/home/username/local \
  -DMPFR_ROOT=/home/username/local \
  -DGMP_ROOT=/home/username/local
```

**Solution 3**: Check CMake output
```bash
# Run CMake with verbose output
cmake .. -DENABLE_HIGH_PRECISION=ON --debug-find 2>&1 | grep -i "boost\|mpfr\|gmp"
```

**Solution 4**: Manually specify library files
```bash
cmake .. \
  -DENABLE_HIGH_PRECISION=ON \
  -DBOOST_ROOT=$HOME/local \
  -DMPFR_LIBRARY=$HOME/local/lib/libmpfr.a \
  -DGMP_LIBRARY=$HOME/local/lib/libgmp.a \
  -DMPFR_INCLUDE_DIR=$HOME/local/include \
  -DGMP_INCLUDE_DIR=$HOME/local/include
```

### Building Dependencies from Source (No Admin Rights)

If you don't have admin rights, build dependencies locally:

#### Build GMP

```bash
cd ~
mkdir -p local/src && cd local/src

wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz && cd gmp-6.3.0

./configure --prefix=$HOME/local --enable-static --disable-shared
make -j$(nproc)
make install
```

#### Build MPFR

```bash
cd ~/local/src

wget https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz
tar xf mpfr-4.2.1.tar.xz && cd mpfr-4.2.1

./configure --prefix=$HOME/local --with-gmp=$HOME/local --enable-static --disable-shared
make -j$(nproc)
make install
```

#### Get Boost Headers

```bash
cd ~/local/src

wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz
tar xzf boost_1_83_0.tar.gz

# Copy headers to local installation
cp -r boost_1_83_0/boost ~/local/include/
```

#### Build Project with Local Libraries

```bash
cd ~/polynomial-solver
mkdir build && cd build

cmake .. \
  -DENABLE_HIGH_PRECISION=ON \
  -DCMAKE_PREFIX_PATH=$HOME/local \
  -DBOOST_ROOT=$HOME/local

make -j$(nproc)
```

**Result**: Static linking, no runtime dependencies!

### Header-Only Backend (Zero Dependencies)

If you can't install MPFR/GMP, use the header-only cpp_dec_float backend:

```bash
# Just copy Boost headers to your project
mkdir -p external
cd external
wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz
tar xzf boost_1_83_0.tar.gz
mv boost_1_83_0/boost .

# Build (CMake will automatically use cpp_dec_float if MPFR not found)
cd ..
mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

**Trade-off**: 2-5× slower than MPFR, but works everywhere with zero external dependencies.

### Verify High-Precision Build

```bash
# Check which backend is being used
cd build
cmake .. -DENABLE_HIGH_PRECISION=ON 2>&1 | grep -i "precision\|mpfr\|backend"

# Expected output (MPFR backend):
# -- Using MPFR backend
# -- MPFR: /usr/lib/x86_64-linux-gnu/libmpfr.so
# -- GMP:  /usr/lib/x86_64-linux-gnu/libgmp.so

# Or (cpp_dec_float backend):
# -- MPFR not found, using cpp_dec_float backend
```

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
cd build
make clean
make -j$(nproc)

# Or remove build directory entirely
cd ..
rm -rf build
mkdir build && cd build
cmake ..
make -j$(nproc)
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

# Build in debug mode
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)

# Run tests
ctest --output-on-failure
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

#### Build Types

```bash
# Release build (optimized, default)
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Debug build (with symbols and debug info)
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)

# RelWithDebInfo (optimized with debug symbols)
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make -j$(nproc)
```

#### Parallel Build

```bash
# Use all available cores
make -j$(nproc)

# Use specific number of cores
make -j4

# Verbose build output
make VERBOSE=1 -j$(nproc)
```

#### Clean Build

```bash
# Clean and rebuild
cd build
make clean
make -j$(nproc)

# Or remove build directory entirely
cd ..
rm -rf build
mkdir build && cd build
cmake ..
make -j$(nproc)
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
- `check_prerequisites.sh` - Prerequisite checker

**Optional**:
- `examples/` - Example programs
- `tools/` - Visualization and analysis tools
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
3. ✅ All 21 tests passing
4. ✅ Example programs running
5. ✅ Ready for development or deployment

For quick reference, see [README.md](README.md).

For algorithm details, see documentation in [docs/](docs/).

## Support

- **Email**: wenyd@lsec.cc.ac.cn
- **GitHub**: https://github.com/gol2em/polynomial-solver
- **Issues**: https://github.com/gol2em/polynomial-solver/issues


