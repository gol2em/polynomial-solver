# FindHighPrecision.cmake
#
# Detects available multiprecision libraries and configures high-precision support
#
# This module sets the following variables:
#   HIGH_PRECISION_FOUND          - TRUE if all required components are found
#   HIGH_PRECISION_BACKEND         - Backend type: "MPFR" or "CPP_DEC_FLOAT"
#   HIGH_PRECISION_INCLUDE_DIRS    - Include directories for high-precision libraries
#   HIGH_PRECISION_LIBRARIES       - Libraries to link for high-precision support
#   QUADMATH_FOUND                 - TRUE if quadmath is available
#   QUADMATH_LIBRARY               - Path to quadmath library
#
# The following cache variables are also set:
#   BOOST_ROOT                     - Root directory of Boost installation
#   GMP_ROOT                       - Root directory of GMP installation
#   MPFR_ROOT                      - Root directory of MPFR installation

# Set policy for FindBoost module (CMake 3.30+)
# The old FindBoost module is deprecated, but we only need Boost headers
# so we can safely use the old behavior or manual detection
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

message(STATUS "")
message(STATUS "========================================")
message(STATUS "Detecting High-Precision Dependencies")
message(STATUS "========================================")
message(STATUS "")
message(STATUS "User preferences:")
message(STATUS "  USE_BOOST:     ${USE_BOOST}")
message(STATUS "  USE_MPFR:      ${USE_MPFR}")
message(STATUS "  USE_GMP:       ${USE_GMP}")
message(STATUS "  USE_QUADMATH:  ${USE_QUADMATH}")

# Initialize variables
set(HIGH_PRECISION_FOUND FALSE)
set(HIGH_PRECISION_BACKEND "NONE")
set(HIGH_PRECISION_INCLUDE_DIRS "")
set(HIGH_PRECISION_LIBRARIES "")
set(QUADMATH_FOUND FALSE)
set(QUADMATH_LIBRARY "")
set(BOOST_FOUND_LOCAL FALSE)

# ============================================
# Step 1: Find Boost Headers (Optional)
# ============================================

message(STATUS "")
message(STATUS "Step 1: Checking for Boost headers...")

if(NOT USE_BOOST)
    message(STATUS "  ⊗ Boost disabled by user (USE_BOOST=OFF)")
    set(BOOST_FOUND_LOCAL FALSE)
else()
    # Try to find Boost
    find_package(Boost QUIET)

    if(Boost_FOUND)
        message(STATUS "  ✓ Boost found: ${Boost_INCLUDE_DIRS}")
        list(APPEND HIGH_PRECISION_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
        set(BOOST_FOUND_LOCAL TRUE)
    else()
        # Try common locations
        set(BOOST_SEARCH_PATHS
            ${BOOST_ROOT}
            $ENV{BOOST_ROOT}
            ${CMAKE_PREFIX_PATH}
            $ENV{HOME}/local
            /usr/local
            /usr
            /opt/local
            /opt
        )

        foreach(SEARCH_PATH ${BOOST_SEARCH_PATHS})
            if(EXISTS "${SEARCH_PATH}/include/boost/version.hpp")
                message(STATUS "  ✓ Boost found: ${SEARCH_PATH}/include")
                list(APPEND HIGH_PRECISION_INCLUDE_DIRS "${SEARCH_PATH}/include")
                set(BOOST_FOUND_LOCAL TRUE)
                break()
            endif()
        endforeach()

        # Check in project's external directory
        if(NOT BOOST_FOUND_LOCAL AND EXISTS "${CMAKE_SOURCE_DIR}/external/boost/version.hpp")
            message(STATUS "  ✓ Boost found: ${CMAKE_SOURCE_DIR}/external/boost")
            list(APPEND HIGH_PRECISION_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/external")
            set(BOOST_FOUND_LOCAL TRUE)
        endif()
    endif()

    if(NOT BOOST_FOUND_LOCAL)
        message(STATUS "  ✗ Boost headers not found")
        message(STATUS "    (Will check for quadmath as alternative)")
    endif()
endif()

# ============================================
# Step 2: Find MPFR and GMP (Optional)
# ============================================

message(STATUS "")
message(STATUS "Step 2: Checking for MPFR and GMP libraries...")

if(NOT USE_MPFR OR NOT USE_GMP OR NOT BOOST_FOUND_LOCAL)
    if(NOT USE_MPFR)
        message(STATUS "  ⊗ MPFR disabled by user (USE_MPFR=OFF)")
    endif()
    if(NOT USE_GMP)
        message(STATUS "  ⊗ GMP disabled by user (USE_GMP=OFF)")
    endif()
    if(NOT BOOST_FOUND_LOCAL)
        message(STATUS "  ⊗ MPFR/GMP require Boost (Boost not available)")
    endif()
    set(MPFR_LIBRARY "")
    set(GMP_LIBRARY "")
    set(MPFR_INCLUDE_DIR "")
    set(GMP_INCLUDE_DIR "")
else()
    # Search paths for MPFR and GMP
    set(LIB_SEARCH_PATHS
        ${MPFR_ROOT}
        ${GMP_ROOT}
        $ENV{MPFR_ROOT}
        $ENV{GMP_ROOT}
        ${CMAKE_PREFIX_PATH}
        $ENV{HOME}/local
        /usr/local
        /usr
        /opt/local
        /opt
    )

    # Find MPFR library
    find_library(MPFR_LIBRARY
        NAMES mpfr libmpfr libmpfr.a
        PATHS ${LIB_SEARCH_PATHS}
        PATH_SUFFIXES lib lib64
        DOC "Path to MPFR library"
    )

    # Find GMP library
    find_library(GMP_LIBRARY
        NAMES gmp libgmp libgmp.a
        PATHS ${LIB_SEARCH_PATHS}
        PATH_SUFFIXES lib lib64
        DOC "Path to GMP library"
    )

    # Find MPFR header
    find_path(MPFR_INCLUDE_DIR
        NAMES mpfr.h
        PATHS ${LIB_SEARCH_PATHS}
        PATH_SUFFIXES include
        DOC "Path to MPFR include directory"
    )

    # Find GMP header
    find_path(GMP_INCLUDE_DIR
        NAMES gmp.h
        PATHS ${LIB_SEARCH_PATHS}
        PATH_SUFFIXES include
        DOC "Path to GMP include directory"
    )
endif()

if(MPFR_LIBRARY AND GMP_LIBRARY AND MPFR_INCLUDE_DIR AND GMP_INCLUDE_DIR AND BOOST_FOUND_LOCAL)
    message(STATUS "  ✓ MPFR library found: ${MPFR_LIBRARY}")
    message(STATUS "  ✓ GMP library found:  ${GMP_LIBRARY}")
    message(STATUS "  ✓ MPFR headers found: ${MPFR_INCLUDE_DIR}")
    message(STATUS "  ✓ GMP headers found:  ${GMP_INCLUDE_DIR}")

    set(HIGH_PRECISION_BACKEND "MPFR")
    list(APPEND HIGH_PRECISION_LIBRARIES ${MPFR_LIBRARY} ${GMP_LIBRARY})
    list(APPEND HIGH_PRECISION_INCLUDE_DIRS ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    set(HIGH_PRECISION_FOUND TRUE)
elseif(BOOST_FOUND_LOCAL)
    message(STATUS "  ✗ MPFR or GMP not found")
    message(STATUS "  → Will use cpp_dec_float backend (header-only, slower)")
    message(STATUS "")
    message(STATUS "  To use MPFR backend (recommended for performance):")
    message(STATUS "    Ubuntu/Debian: sudo apt-get install libmpfr-dev libgmp-dev")
    message(STATUS "    Fedora/RHEL:   sudo dnf install mpfr-devel gmp-devel")
    message(STATUS "    macOS:         brew install mpfr gmp")
    message(STATUS "    Or build from source and specify:")
    message(STATUS "      -DMPFR_ROOT=/path/to/mpfr -DGMP_ROOT=/path/to/gmp")

    set(HIGH_PRECISION_BACKEND "CPP_DEC_FLOAT")
    set(HIGH_PRECISION_FOUND TRUE)
else()
    message(STATUS "  ✗ MPFR or GMP not found")
    message(STATUS "  → Boost not available, checking for quadmath...")
endif()

# ============================================
# Step 3: Find Quadmath (Always check)
# ============================================

message(STATUS "")
message(STATUS "Step 3: Checking for quadmath library...")

if(NOT USE_QUADMATH)
    message(STATUS "  ⊗ Quadmath disabled by user (USE_QUADMATH=OFF)")
    set(QUADMATH_LIBRARY "")
else()
    # Try to find quadmath library
    # First try standard find_library
    find_library(QUADMATH_LIBRARY
        NAMES quadmath
        PATHS
            ${LIB_SEARCH_PATHS}
            /usr/lib
            /usr/local/lib
        PATH_SUFFIXES
            lib lib64
            x86_64-linux-gnu
            aarch64-linux-gnu
            gcc/x86_64-linux-gnu/11 gcc/x86_64-linux-gnu/12 gcc/x86_64-linux-gnu/13
            gcc/aarch64-linux-gnu/11 gcc/aarch64-linux-gnu/12 gcc/aarch64-linux-gnu/13
        DOC "Path to quadmath library"
    )

    # If not found, check common locations manually
    if(NOT QUADMATH_LIBRARY)
        set(QUADMATH_COMMON_PATHS
            "/usr/lib/x86_64-linux-gnu/libquadmath.so.0"
            "/usr/lib/x86_64-linux-gnu/libquadmath.so"
            "/usr/lib/aarch64-linux-gnu/libquadmath.so.0"
            "/usr/lib/aarch64-linux-gnu/libquadmath.so"
            "/usr/lib64/libquadmath.so.0"
            "/usr/lib64/libquadmath.so"
            "/usr/lib/gcc/x86_64-linux-gnu/12/libquadmath.so"
            "/usr/lib/gcc/x86_64-linux-gnu/11/libquadmath.so"
            "/usr/lib/gcc/x86_64-linux-gnu/13/libquadmath.so"
        )

        foreach(QUAD_PATH ${QUADMATH_COMMON_PATHS})
            if(EXISTS "${QUAD_PATH}")
                set(QUADMATH_LIBRARY "${QUAD_PATH}")
                break()
            endif()
        endforeach()
    endif()
endif()

if(QUADMATH_LIBRARY)
    message(STATUS "  ✓ Quadmath library found: ${QUADMATH_LIBRARY}")
    set(QUADMATH_FOUND TRUE)

    # If no other backend found, use quadmath as the backend
    if(NOT HIGH_PRECISION_FOUND)
        message(STATUS "  → Using quadmath as high-precision backend")
        set(HIGH_PRECISION_BACKEND "QUADMATH")
        list(APPEND HIGH_PRECISION_LIBRARIES ${QUADMATH_LIBRARY})
        set(HIGH_PRECISION_FOUND TRUE)
    endif()
else()
    message(STATUS "  ✗ Quadmath library not found")
    if(NOT HIGH_PRECISION_FOUND)
        message(STATUS "  → No high-precision libraries available")
    endif()
    message(STATUS "")
    message(STATUS "  To enable quadmath support:")
    message(STATUS "    Ubuntu/Debian: sudo apt-get install libquadmath0")
    message(STATUS "    (Usually included with GCC)")
endif()

# ============================================
# Summary
# ============================================

message(STATUS "")
message(STATUS "========================================")
message(STATUS "High-Precision Configuration Summary")
message(STATUS "========================================")

if(HIGH_PRECISION_FOUND)
    message(STATUS "  Status:           Available")
    message(STATUS "  Backend:          ${HIGH_PRECISION_BACKEND}")

    if(BOOST_FOUND_LOCAL)
        message(STATUS "  Boost headers:    Found")
    elseif(NOT USE_BOOST)
        message(STATUS "  Boost headers:    Disabled by user")
    else()
        message(STATUS "  Boost headers:    Not found")
    endif()

    if(HIGH_PRECISION_BACKEND STREQUAL "MPFR")
        message(STATUS "  MPFR library:     Found")
        message(STATUS "  GMP library:      Found")
        message(STATUS "  Performance:      Fastest")
        message(STATUS "  Precision:        Arbitrary (runtime configurable)")
    elseif(HIGH_PRECISION_BACKEND STREQUAL "CPP_DEC_FLOAT")
        if(NOT USE_MPFR OR NOT USE_GMP)
            message(STATUS "  MPFR/GMP:         Disabled by user (using cpp_dec_float)")
        else()
            message(STATUS "  MPFR library:     Not found (using cpp_dec_float)")
        endif()
        message(STATUS "  Performance:      Good (2-5× slower than MPFR)")
        message(STATUS "  Precision:        Fixed (50, 100 decimal digits)")
    elseif(HIGH_PRECISION_BACKEND STREQUAL "QUADMATH")
        message(STATUS "  Quadmath:         Found (standalone backend)")
        message(STATUS "  Performance:      Very good")
        message(STATUS "  Precision:        Fixed (33-36 decimal digits)")
    endif()

    if(QUADMATH_FOUND AND NOT HIGH_PRECISION_BACKEND STREQUAL "QUADMATH")
        message(STATUS "  Quadmath:         Available (can be used with templates)")
    elseif(NOT QUADMATH_FOUND AND USE_QUADMATH)
        message(STATUS "  Quadmath:         Not found")
    elseif(NOT USE_QUADMATH)
        message(STATUS "  Quadmath:         Disabled by user")
    endif()
else()
    message(STATUS "  Status:           Not available")

    if(NOT USE_BOOST)
        message(STATUS "  Boost headers:    Disabled by user")
    else()
        message(STATUS "  Boost headers:    Not found")
    endif()

    if(NOT USE_MPFR OR NOT USE_GMP)
        message(STATUS "  MPFR/GMP:         Disabled by user")
    else()
        message(STATUS "  MPFR library:     Not found")
    endif()

    if(NOT USE_QUADMATH)
        message(STATUS "  Quadmath:         Disabled by user")
    else()
        message(STATUS "  Quadmath:         Not found")
    endif()

    message(STATUS "")
    message(STATUS "  To enable libraries, use:")
    message(STATUS "    -DUSE_BOOST=ON -DUSE_MPFR=ON -DUSE_GMP=ON -DUSE_QUADMATH=ON")
    message(STATUS "")
    message(STATUS "  Install at least one of:")
    message(STATUS "    • Boost (for mpreal with MPFR or cpp_dec_float)")
    message(STATUS "    • Quadmath (for __float128)")
    message(STATUS "")
    message(STATUS "  Installation:")
    message(STATUS "    Ubuntu/Debian: sudo apt-get install libboost-dev libmpfr-dev libgmp-dev")
    message(STATUS "    Fedora/RHEL:   sudo dnf install boost-devel mpfr-devel gmp-devel")
    message(STATUS "    macOS:         brew install boost mpfr gmp")
    message(STATUS "")
    message(STATUS "  Or specify custom paths:")
    message(STATUS "    -DBOOST_ROOT=/path/to/boost")
    message(STATUS "    -DMPFR_ROOT=/path/to/mpfr -DGMP_ROOT=/path/to/gmp")
endif()

message(STATUS "========================================")
message(STATUS "")

