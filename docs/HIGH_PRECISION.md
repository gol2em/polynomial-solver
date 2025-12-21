# High-Precision Arithmetic Support

## Quick Start

```bash
# Install dependencies
sudo apt-get install libboost-dev libmpfr-dev libgmp-dev  # Ubuntu/Debian

# Build with high precision
mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

## Precision Tiers

| Tier | Description | Dependencies | Use Case |
|------|-------------|--------------|----------|
| 1 | Double precision (15-17 digits) | None | Default, maximum performance |
| 2 | Fixed high precision | Boost + MPFR/GMP | Production, simple integration |
| 3 | Template-based (planned) | Boost + MPFR/GMP | Flexible precision at runtime |

## Backend Options

| Backend | Requires | Precision | Speed | Notes |
|---------|----------|-----------|-------|-------|
| MPFR | Boost + MPFR + GMP | Configurable | ⚡⚡⚡ | Recommended |
| cpp_dec_float | Boost only | 100 digits | ⚡⚡ | Header-only, no library deps |
| quadmath | libquadmath | 34 digits | ⚡⚡⚡⚡ | GCC native |

## API Examples

```cpp
// Tier 1: Double precision
#include "polynomial.h"
Polynomial poly(degrees, coeffs);

// Tier 2: High precision
#include "polynomial_hp.h"
PolynomialHP poly_hp(degrees, coeffs_hp);
```

## CMake Flags

| Flag | Default | Description |
|------|---------|-------------|
| `ENABLE_HIGH_PRECISION` | OFF | Enable high-precision support |
| `USE_BOOST` | ON | Enable Boost multiprecision |
| `USE_MPFR` | ON | Enable MPFR backend (requires Boost) |
| `USE_GMP` | ON | Enable GMP (requires Boost and MPFR) |
| `USE_QUADMATH` | ON | Enable quadmath support |

**Backend selection logic**:
- `USE_BOOST=ON` + `USE_MPFR=ON` + `USE_GMP=ON` → MPFR backend
- `USE_BOOST=ON` + MPFR/GMP disabled → cpp_dec_float backend
- `USE_BOOST=OFF` + `USE_QUADMATH=ON` → quadmath backend

## Custom Library Paths

```bash
# Method 1: Individual paths
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DBOOST_ROOT=$HOME/local \
  -DGMP_ROOT=$HOME/local \
  -DMPFR_ROOT=$HOME/local

# Method 2: CMake prefix
cmake .. -DENABLE_HIGH_PRECISION=ON \
  -DCMAKE_PREFIX_PATH="$HOME/local"

# Method 3: Environment variables
export CMAKE_PREFIX_PATH=$HOME/local
cmake .. -DENABLE_HIGH_PRECISION=ON
```

## Building Dependencies from Source

For systems without admin rights:

```bash
# Build GMP
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz && cd gmp-6.3.0
./configure --prefix=$HOME/local --enable-static --disable-shared
make -j$(nproc) && make install

# Build MPFR
wget https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz
tar xf mpfr-4.2.1.tar.xz && cd mpfr-4.2.1
./configure --prefix=$HOME/local --with-gmp=$HOME/local --enable-static --disable-shared
make -j$(nproc) && make install

# Get Boost headers
wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz
tar xzf boost_1_83_0.tar.gz
cp -r boost_1_83_0/boost $HOME/local/include/
```

## Header-Only Backend (Zero Dependencies)

```bash
# Copy Boost headers to project
mkdir -p external && cd external
wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz
tar xzf boost_1_83_0.tar.gz && mv boost_1_83_0/boost .

# Build (auto-detects cpp_dec_float if MPFR not found)
cd .. && mkdir build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

**Trade-off**: 2-5× slower than MPFR, but zero external dependencies.

## Related Documentation

- [HIGH_PRECISION_MINIMAL_DEPENDENCIES.md](HIGH_PRECISION_MINIMAL_DEPENDENCIES.md) - Minimal setup guide
- [HIGH_PRECISION_MULTIPLE_TYPES.md](HIGH_PRECISION_MULTIPLE_TYPES.md) - Multiple precision types
- [HIGH_PRECISION_SETUP_GUIDE.md](HIGH_PRECISION_SETUP_GUIDE.md) - Server setup guide
- [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - Quick reference
