# Step-by-Step Setup Guide for Dedicated Server

This guide shows you how to set up high-precision arithmetic on a dedicated server **without admin rights**.

## Quick Decision Tree

```
Can you install packages? (sudo apt install)
├─ YES → Use MPFR backend (fastest)
└─ NO → Can you build from source?
    ├─ YES → Build MPFR locally (fastest)
    └─ NO → Use cpp_dec_float (header-only, slower but works)
```

## Option 1: Build MPFR Locally (Recommended)

**Time**: 30 minutes  
**Result**: Fastest high-precision arithmetic, no admin rights needed

### Step 1: Build GMP and MPFR

```bash
# On your dedicated server (no admin rights needed)
cd ~
mkdir -p local/src
cd local/src

# Download and build GMP
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --prefix=$HOME/local --enable-static --disable-shared
make -j4
make install

# Download and build MPFR
cd ~/local/src
wget https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz
tar xf mpfr-4.2.1.tar.xz
cd mpfr-4.2.1
./configure --prefix=$HOME/local --with-gmp=$HOME/local --enable-static --disable-shared
make -j4
make install

# Verify installation
ls ~/local/lib/libgmp.a    # Should exist
ls ~/local/lib/libmpfr.a   # Should exist
ls ~/local/include/gmp.h   # Should exist
ls ~/local/include/mpfr.h  # Should exist
```

### Step 2: Get Boost Headers

```bash
cd ~/local/src

# Download Boost (we only need headers)
wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz
tar xzf boost_1_83_0.tar.gz

# Copy headers to your project (or use in-place)
cp -r boost_1_83_0/boost ~/local/include/

# Or create a symlink
ln -s ~/local/src/boost_1_83_0/boost ~/local/include/boost
```

### Step 3: Configure Your Project

```bash
cd ~/myproject

# Create build directory
mkdir build
cd build

# Configure with local libraries
cmake .. \
  -DENABLE_HIGH_PRECISION=ON \
  -DCMAKE_PREFIX_PATH=$HOME/local \
  -DBOOST_ROOT=$HOME/local

# Build
make -j4

# Test
./bin/test_high_precision
```

### Step 4: Update CMakeLists.txt

```cmake
option(ENABLE_HIGH_PRECISION "Enable Boost.Multiprecision support" OFF)

if(ENABLE_HIGH_PRECISION)
    # Find Boost headers
    set(BOOST_ROOT "$ENV{HOME}/local" CACHE PATH "Boost root directory")
    find_package(Boost)
    
    if(NOT Boost_FOUND)
        # Try local installation
        include_directories($ENV{HOME}/local/include)
    endif()
    
    # Find MPFR and GMP (prefer static libraries)
    find_library(MPFR_LIBRARY 
        NAMES libmpfr.a mpfr
        PATHS $ENV{HOME}/local/lib
        NO_DEFAULT_PATH
    )
    
    find_library(GMP_LIBRARY 
        NAMES libgmp.a gmp
        PATHS $ENV{HOME}/local/lib
        NO_DEFAULT_PATH
    )
    
    if(MPFR_LIBRARY AND GMP_LIBRARY)
        message(STATUS "Using MPFR backend")
        message(STATUS "  MPFR: ${MPFR_LIBRARY}")
        message(STATUS "  GMP:  ${GMP_LIBRARY}")
        
        add_definitions(-DENABLE_HIGH_PRECISION -DUSE_MPFR_BACKEND)
        include_directories($ENV{HOME}/local/include)
        target_link_libraries(polynomial_solver_lib ${MPFR_LIBRARY} ${GMP_LIBRARY})
    else()
        message(WARNING "MPFR/GMP not found, falling back to cpp_dec_float")
        add_definitions(-DENABLE_HIGH_PRECISION -DUSE_CPP_DEC_FLOAT_BACKEND)
    endif()
endif()
```

**Result**: Static linking, no runtime dependencies, fastest performance!

## Option 2: Use cpp_dec_float (Zero Dependencies)

**Time**: 5 minutes  
**Result**: Good high-precision arithmetic, works everywhere, 2-5× slower than MPFR

### Step 1: Get Boost Headers

```bash
cd ~/myproject

# Download Boost
wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz
tar xzf boost_1_83_0.tar.gz

# Create external directory
mkdir -p external
mv boost_1_83_0/boost external/

# Clean up
rm -rf boost_1_83_0 boost_1_83_0.tar.gz
```

### Step 2: Update CMakeLists.txt

```cmake
option(ENABLE_HIGH_PRECISION "Enable Boost.Multiprecision support" OFF)

if(ENABLE_HIGH_PRECISION)
    # Use local Boost headers
    include_directories(${CMAKE_SOURCE_DIR}/external/boost)
    
    # Use cpp_dec_float (header-only, no external libraries)
    add_definitions(-DENABLE_HIGH_PRECISION -DUSE_CPP_DEC_FLOAT_BACKEND)
    
    message(STATUS "Using cpp_dec_float backend (header-only)")
endif()
```

### Step 3: Update Code

```cpp
// include/high_precision_types.h

#ifdef ENABLE_HIGH_PRECISION
    #ifdef USE_MPFR_BACKEND
        #include <boost/multiprecision/mpfr.hpp>
        using mpreal = boost::multiprecision::mpfr_float;
    #else
        // cpp_dec_float backend (header-only)
        #include <boost/multiprecision/cpp_dec_float.hpp>
        
        // Use 100 decimal digits (equivalent to ~332 bits)
        using mpreal = boost::multiprecision::cpp_dec_float_100;
        
        // Or define custom precision
        // using mpreal = boost::multiprecision::number<
        //     boost::multiprecision::cpp_dec_float<150>
        // >;
    #endif
#endif
```

### Step 4: Build

```bash
cd ~/myproject/build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j4
```

**Result**: Zero external dependencies, works everywhere!

## Comparison

| Aspect | MPFR Backend | cpp_dec_float Backend |
|--------|--------------|----------------------|
| **Setup time** | 30 minutes | 5 minutes |
| **Dependencies** | GMP + MPFR (~2 MB) | None |
| **Performance** | Fastest | 2-5× slower |
| **Precision** | Arbitrary (runtime) | Fixed (compile-time) |
| **Portability** | Need to build | Works everywhere |
| **Recommendation** | If you can build | If you can't install anything |

## Testing Your Setup

Create a test file to verify everything works:

```cpp
// test_precision.cpp

#include <iostream>
#include <iomanip>

#ifdef ENABLE_HIGH_PRECISION
    #ifdef USE_MPFR_BACKEND
        #include <boost/multiprecision/mpfr.hpp>
        using mpreal = boost::multiprecision::mpfr_float;
    #else
        #include <boost/multiprecision/cpp_dec_float.hpp>
        using mpreal = boost::multiprecision::cpp_dec_float_100;
    #endif
#endif

int main() {
#ifdef ENABLE_HIGH_PRECISION
    #ifdef USE_MPFR_BACKEND
        std::cout << "Using MPFR backend" << std::endl;
        mpreal::set_default_prec(256);  // 77 decimal digits
    #else
        std::cout << "Using cpp_dec_float backend" << std::endl;
    #endif
    
    // Compute pi with high precision
    mpreal pi = boost::multiprecision::atan(mpreal(1)) * 4;
    
    std::cout << "Pi with high precision:" << std::endl;
    std::cout << std::setprecision(50) << pi << std::endl;
    
    // Expected: 3.1415926535897932384626433832795028841971693993751...
#else
    std::cout << "High precision not enabled" << std::endl;
#endif
    
    return 0;
}
```

Build and run:

```bash
cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make test_precision
./test_precision
```

## Troubleshooting

### Problem: "boost/multiprecision/mpfr.hpp: No such file"

**Solution**: Boost headers not found

```bash
# Check if Boost headers exist
ls ~/local/include/boost/multiprecision/mpfr.hpp

# If not, copy them
cp -r ~/local/src/boost_1_83_0/boost ~/local/include/
```

### Problem: "undefined reference to __gmpz_init"

**Solution**: GMP library not linked

```cmake
# Make sure you're linking the libraries
target_link_libraries(polynomial_solver_lib 
    ${MPFR_LIBRARY} 
    ${GMP_LIBRARY}
)
```

### Problem: Build takes too long

**Solution**: Use parallel build

```bash
make -j4  # Use 4 cores
# or
make -j$(nproc)  # Use all available cores
```

## Summary

### Minimal Dependencies

**For MPFR backend** (fastest):
- Boost headers: ~50 MB (header-only)
- GMP library: ~1 MB (can build from source)
- MPFR library: ~1 MB (can build from source)
- **Total**: ~52 MB, 30 minutes to set up

**For cpp_dec_float backend** (easiest):
- Boost headers: ~50 MB (header-only)
- **Total**: ~50 MB, 5 minutes to set up

### Recommendation

1. **Try MPFR first**: Build from source in `~/local`, best performance
2. **Fallback to cpp_dec_float**: If build fails, use header-only backend
3. **Both work**: Your code is compatible with both backends!

**No admin rights needed for either option!**

