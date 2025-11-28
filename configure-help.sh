#!/bin/bash
# Configuration help for Polynomial Solver

cat << 'EOF'
==========================================
Polynomial Solver Configuration Options
==========================================

High-Precision Options:
  -DENABLE_HIGH_PRECISION=ON/OFF    Enable high-precision arithmetic (default: OFF)
  -DENABLE_TEMPLATES=ON/OFF         Enable template-based implementation (default: ON)
                                      ON  = Tier 3 (flexible, any numeric type)
                                      OFF = Tier 2 (fixed backend, no templates)

Backend Selection (Tier 2 only - when ENABLE_TEMPLATES=OFF):
  -DHP_BACKEND=AUTO|MPFR|CPP_DEC_FLOAT|QUADMATH
                                    Choose backend for Tier 2 (default: AUTO)
                                      AUTO          = Use detected backend
                                      MPFR          = Use MPFR (best performance)
                                      CPP_DEC_FLOAT = Use cpp_dec_float (header-only)
                                      QUADMATH      = Use __float128 (128-bit only)

Library Control (manual override for detection):
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

1. Tier 1 - Double precision only (default):
   mkdir build && cd build
   cmake ..
   make -j$(nproc)

2. Tier 3 - Template-based with MPFR (best performance, flexible):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON
   make -j$(nproc)

3. Tier 3 - Template-based with cpp_dec_float (header-only):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF
   make -j$(nproc)

4. Tier 2 - Fixed high-precision with MPFR (no templates):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DENABLE_TEMPLATES=OFF -DHP_BACKEND=MPFR
   make -j$(nproc)

5. Tier 2 - Fixed high-precision with cpp_dec_float (no templates):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DENABLE_TEMPLATES=OFF -DHP_BACKEND=CPP_DEC_FLOAT
   make -j$(nproc)

6. Tier 2 - Fixed high-precision with quadmath (no templates):
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DENABLE_TEMPLATES=OFF -DHP_BACKEND=QUADMATH
   make -j$(nproc)

7. Release build with high-precision:
   mkdir build && cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON -DCMAKE_BUILD_TYPE=Release
   make -j$(nproc)

8. Custom library paths:
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

Tier 1: Double Precision Only (Default)
  - No multiprecision dependencies
  - Only double precision
  - Zero overhead
  - Always available
  - Enable with: (default, no flags needed)

Tier 2: Fixed High-Precision (No Templates)
  - Single high-precision backend (MPFR, cpp_dec_float, or quadmath)
  - No templates, simpler implementation
  - Choose backend at configure time with -DHP_BACKEND=<backend>
  - Enable with: -DENABLE_HIGH_PRECISION=ON -DENABLE_TEMPLATES=OFF
  - Use when: Need high-precision with single backend

Tier 3: Template-Based (Flexible)
  - Flexible multiprecision support
  - Template-based (supports any numeric type)
  - Can switch between backends at runtime (MPFR)
  - Enable with: -DENABLE_HIGH_PRECISION=ON -DENABLE_TEMPLATES=ON (default)
  - Use when: Need to switch between different numeric types

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

