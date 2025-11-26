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

message(STATUS "")
message(STATUS "========================================")
message(STATUS "Detecting High-Precision Dependencies")
message(STATUS "========================================")

# Initialize variables
set(HIGH_PRECISION_FOUND FALSE)
set(HIGH_PRECISION_BACKEND "NONE")
set(HIGH_PRECISION_INCLUDE_DIRS "")
set(HIGH_PRECISION_LIBRARIES "")
set(QUADMATH_FOUND FALSE)
set(QUADMATH_LIBRARY "")

# ============================================
# Step 1: Find Boost Headers (Required)
# ============================================

message(STATUS "")
message(STATUS "Step 1: Checking for Boost headers...")

# Try to find Boost
find_package(Boost QUIET)

if(Boost_FOUND)
    message(STATUS "  ✓ Boost found: ${Boost_INCLUDE_DIRS}")
    list(APPEND HIGH_PRECISION_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
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
            set(Boost_FOUND TRUE)
            break()
        endif()
    endforeach()
    
    # Check in project's external directory
    if(NOT Boost_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/external/boost/version.hpp")
        message(STATUS "  ✓ Boost found: ${CMAKE_SOURCE_DIR}/external/boost")
        list(APPEND HIGH_PRECISION_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/external")
        set(Boost_FOUND TRUE)
    endif()
endif()

if(NOT Boost_FOUND)
    message(STATUS "  ✗ Boost headers not found")
    message(STATUS "")
    message(FATAL_ERROR 
        "Boost headers are required for high-precision support.\n"
        "Please install Boost or copy headers to your project:\n"
        "  Ubuntu/Debian: sudo apt-get install libboost-dev\n"
        "  Fedora/RHEL:   sudo dnf install boost-devel\n"
        "  macOS:         brew install boost\n"
        "  Or copy to:    ${CMAKE_SOURCE_DIR}/external/boost/\n"
        "  Or specify:    -DBOOST_ROOT=/path/to/boost\n"
    )
endif()

# ============================================
# Step 2: Find MPFR and GMP (Optional)
# ============================================

message(STATUS "")
message(STATUS "Step 2: Checking for MPFR and GMP libraries...")

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

if(MPFR_LIBRARY AND GMP_LIBRARY AND MPFR_INCLUDE_DIR AND GMP_INCLUDE_DIR)
    message(STATUS "  ✓ MPFR library found: ${MPFR_LIBRARY}")
    message(STATUS "  ✓ GMP library found:  ${GMP_LIBRARY}")
    message(STATUS "  ✓ MPFR headers found: ${MPFR_INCLUDE_DIR}")
    message(STATUS "  ✓ GMP headers found:  ${GMP_INCLUDE_DIR}")
    
    set(HIGH_PRECISION_BACKEND "MPFR")
    list(APPEND HIGH_PRECISION_LIBRARIES ${MPFR_LIBRARY} ${GMP_LIBRARY})
    list(APPEND HIGH_PRECISION_INCLUDE_DIRS ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    set(HIGH_PRECISION_FOUND TRUE)
else()
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
endif()

# ============================================
# Step 3: Find Quadmath (Optional)
# ============================================

if(ENABLE_QUADMATH)
    message(STATUS "")
    message(STATUS "Step 3: Checking for quadmath library...")
    
    find_library(QUADMATH_LIBRARY
        NAMES quadmath libquadmath
        PATHS ${LIB_SEARCH_PATHS}
        PATH_SUFFIXES lib lib64
        DOC "Path to quadmath library"
    )
    
    if(QUADMATH_LIBRARY)
        message(STATUS "  ✓ Quadmath library found: ${QUADMATH_LIBRARY}")
        set(QUADMATH_FOUND TRUE)
    else()
        message(STATUS "  ✗ Quadmath library not found")
        message(STATUS "  → __float128 support will be disabled")
        message(STATUS "")
        message(STATUS "  To enable quadmath support:")
        message(STATUS "    Ubuntu/Debian: sudo apt-get install libquadmath0")
        message(STATUS "    (Usually included with GCC)")
    endif()
endif()

# ============================================
# Summary
# ============================================

message(STATUS "")
message(STATUS "========================================")
message(STATUS "High-Precision Configuration Summary")
message(STATUS "========================================")
message(STATUS "  Backend:          ${HIGH_PRECISION_BACKEND}")
message(STATUS "  Boost headers:    Found")
if(HIGH_PRECISION_BACKEND STREQUAL "MPFR")
    message(STATUS "  MPFR library:     Found")
    message(STATUS "  GMP library:      Found")
    message(STATUS "  Performance:      Fastest")
else()
    message(STATUS "  MPFR library:     Not found (using cpp_dec_float)")
    message(STATUS "  Performance:      Good (2-5× slower than MPFR)")
endif()
if(ENABLE_QUADMATH)
    if(QUADMATH_FOUND)
        message(STATUS "  Quadmath:         Enabled")
    else()
        message(STATUS "  Quadmath:         Not found")
    endif()
endif()
message(STATUS "========================================")
message(STATUS "")

