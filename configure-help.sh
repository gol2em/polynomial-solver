#!/bin/bash
# Configuration help for Polynomial Solver

cat << 'EOF'
==========================================
Polynomial Solver Configuration Options
==========================================

High-Precision Options:
  -DENABLE_HIGH_PRECISION=ON/OFF    Enable high-precision arithmetic (default: OFF)
  -DDISABLE_TEMPLATES=ON/OFF        Disable template-based implementation (default: OFF)

Library Control (manual override):
  -DUSE_BOOST=ON/OFF                Use Boost multiprecision (default: ON)
  -DUSE_MPFR=ON/OFF                 Use MPFR library (default: ON)
  -DUSE_GMP=ON/OFF                  Use GMP library (default: ON)
  -DUSE_QUADMATH=ON/OFF             Use quadmath for __float128 (default: ON)

Library Path Hints:
  -DBOOST_ROOT=/path/to/boost       Boost installation directory
  -DMPFR_ROOT=/path/to/mpfr         MPFR installation directory
  -DGMP_ROOT=/path/to/gmp           GMP installation directory
  -DCMAKE_PREFIX_PATH=/path         Search path for all libraries

Build Options:
  -DBUILD_EXAMPLES=ON/OFF           Build example programs (default: ON)
  -DCMAKE_BUILD_TYPE=Debug/Release  Build type (default: Debug)

==========================================
Common Configurations
==========================================

1. Default (double precision only):
   mkdir build && cd build
   cmake ..
   make -j$(nproc)

2. High-precision with MPFR (best performance):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON
   make -j$(nproc)

3. High-precision with cpp_dec_float (header-only):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF
   make -j$(nproc)

4. High-precision without templates (fallback):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON
   make -j$(nproc)

5. Release build with high-precision:
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DCMAKE_BUILD_TYPE=Release
   make -j$(nproc)

6. Custom library paths:
   mkdir build && cd build
   cmake .. \
     -DENABLE_HIGH_PRECISION=ON \
     -DBOOST_ROOT=$HOME/local \
     -DMPFR_ROOT=$HOME/local \
     -DGMP_ROOT=$HOME/local
   make -j$(nproc)

==========================================
Checking Configuration
==========================================

To see current configuration:
  cd build
  cmake -L ..                    # Show project-specific options
  cmake -LA ..                   # Show all options including advanced

To see detailed configuration summary:
  cd build
  cmake ..                       # Configuration summary is printed

To reconfigure:
  cd build
  rm -rf *                       # Clean build directory
  cmake .. [options]             # Reconfigure with new options

==========================================
Three-Tier System
==========================================

Tier 1: Minimum Version (Default)
  - No multiprecision dependencies
  - Only double precision
  - Zero overhead
  - Always available

Tier 2: Fallback Version (Partial Multiprecision)
  - Fixed high-precision support
  - No templates (separate implementation)
  - Enable with: -DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON
  - Supports one or two fixed precision types

Tier 3: Full Version (Template-Based)
  - Flexible multiprecision support
  - Template-based (supports any numeric type)
  - Enable with: -DENABLE_HIGH_PRECISION=ON
  - Supports double, __float128, mpreal at any precision

==========================================
For More Information
==========================================

Documentation:
  docs/HIGH_PRECISION.md              - High-precision overview
  docs/MANUAL_LIBRARY_CONTROL.md      - Manual flag documentation
  docs/CONFIGURATION_TEST_RESULTS.md  - Tested configurations
  docs/REFACTORIZATION_PROGRESS.md    - Implementation progress
  SETUP.md                            - Build instructions

Test all configurations:
  ./test_all_configurations.sh

==========================================
EOF

