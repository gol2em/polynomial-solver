cmake_minimum_required(VERSION 3.15)
project(PolynomialSolver VERSION 1.0.0 LANGUAGES CXX)

# Always export compile_commands.json for IDE support (VS Code, CLion, etc.)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add cmake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#==============================================================================
# Configuration Options
#==============================================================================

# Geometry dump macro control
# In Debug mode: macro is defined, feature is enabled and controlled by runtime flag
# In Release mode: macro is not defined, all dump code is compiled out
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "")
    add_compile_definitions(ENABLE_GEOMETRY_DUMP)
    message(STATUS "Geometry dump: ENABLED (Debug mode - controlled by runtime flag)")
else()
    message(STATUS "Geometry dump: DISABLED (Release mode - code compiled out)")
endif()

# High-precision arithmetic options
option(ENABLE_HIGH_PRECISION "Enable high-precision arithmetic support" OFF)
option(ENABLE_TEMPLATES "Enable template-based implementation (Tier 3)" ON)

# Manual control over which libraries to use
option(USE_BOOST "Use Boost multiprecision library" ON)
option(USE_MPFR "Use MPFR library (requires Boost)" ON)
option(USE_GMP "Use GMP library (requires Boost and MPFR)" ON)
option(USE_QUADMATH "Use quadmath library for __float128 support" ON)

# Backend selection for Tier 2 (no templates)
# Only used when ENABLE_HIGH_PRECISION=ON and ENABLE_TEMPLATES=OFF
set(HP_BACKEND "AUTO" CACHE STRING "High-precision backend for Tier 2: AUTO, MPFR, CPP_DEC_FLOAT, QUADMATH")
set_property(CACHE HP_BACKEND PROPERTY STRINGS AUTO MPFR CPP_DEC_FLOAT QUADMATH)

# Build options
option(BUILD_EXAMPLES "Build example programs" ON)

# Note: For configuration help, run: ./configure-help.sh
# Or use: cmake -L .. to list all options

#==============================================================================
# High-Precision Configuration
#==============================================================================

# Always check for high-precision dependencies to enable auto-detection
include(FindHighPrecision)

# Auto-enable high-precision if libraries are detected and not explicitly disabled
if(NOT ENABLE_HIGH_PRECISION AND HIGH_PRECISION_FOUND)
    message(STATUS "")
    message(STATUS "High-precision libraries detected!")
    message(STATUS "  To enable: cmake .. -DENABLE_HIGH_PRECISION=ON")
    message(STATUS "  (Currently disabled by default)")
    message(STATUS "")
endif()

if(ENABLE_HIGH_PRECISION)
    if(NOT HIGH_PRECISION_FOUND)
        message(WARNING
            "High-precision support requested but no suitable libraries found.\n"
            "  Boost or quadmath is required.\n"
            "  Disabling high-precision support and continuing with double precision.\n"
            "\n"
            "  To enable high-precision, install dependencies:\n"
            "    Ubuntu/Debian: sudo apt-get install libboost-dev libmpfr-dev libgmp-dev\n"
            "    Fedora/RHEL:   sudo dnf install boost-devel mpfr-devel gmp-devel\n"
            "    macOS:         brew install boost mpfr gmp\n"
            "  Or specify custom paths:\n"
            "    -DBOOST_ROOT=/path/to/boost\n"
            "    -DMPFR_ROOT=/path/to/mpfr -DGMP_ROOT=/path/to/gmp\n"
        )
        set(ENABLE_HIGH_PRECISION OFF CACHE BOOL "Enable high-precision arithmetic support" FORCE)
    else()
        # Define ENABLE_HIGH_PRECISION for both compilation and IDE support
        # Also defined in config.h with #ifndef guard to prevent redefinition
        add_compile_definitions(ENABLE_HIGH_PRECISION)

        # Determine if templates should be used
        if(ENABLE_TEMPLATES)
            # Set CMake variable for config.h generation
            set(USE_TEMPLATES ON)
            # Define USE_TEMPLATES for both compilation and IDE support
            # Also defined in config.h with #ifndef guard to prevent redefinition
            add_compile_definitions(USE_TEMPLATES)
            message(STATUS "Implementation: Tier 3 (template-based, flexible)")
        else()
            message(STATUS "Implementation: Tier 2 (fixed high-precision, no templates)")

            # For Tier 2, enforce single backend selection
            if(HP_BACKEND STREQUAL "AUTO")
                # Auto-select based on available libraries
                message(STATUS "Backend selection: AUTO (using detected backend: ${HIGH_PRECISION_BACKEND})")
            else()
                # User specified backend - validate it's available
                if(HP_BACKEND STREQUAL "MPFR" AND NOT HIGH_PRECISION_BACKEND STREQUAL "MPFR")
                    message(FATAL_ERROR "HP_BACKEND=MPFR requested but MPFR not available. Install MPFR or choose different backend.")
                elseif(HP_BACKEND STREQUAL "CPP_DEC_FLOAT" AND NOT BOOST_FOUND_LOCAL)
                    message(FATAL_ERROR "HP_BACKEND=CPP_DEC_FLOAT requested but Boost not available. Install Boost or choose different backend.")
                elseif(HP_BACKEND STREQUAL "QUADMATH" AND NOT QUADMATH_FOUND)
                    message(FATAL_ERROR "HP_BACKEND=QUADMATH requested but quadmath not available. Use GCC compiler or choose different backend.")
                endif()

                # Override detected backend with user choice
                set(HIGH_PRECISION_BACKEND ${HP_BACKEND})
                message(STATUS "Backend selection: ${HP_BACKEND} (user specified)")
            endif()
        endif()

        # Set backend type
        # Set CMake variables for config.h generation AND compile definitions for IDE support
        if(HIGH_PRECISION_BACKEND STREQUAL "MPFR")
            set(USE_MPFR_BACKEND ON)
            add_compile_definitions(USE_MPFR_BACKEND)
            message(STATUS "Backend: MPFR (best performance, runtime precision)")
        elseif(HIGH_PRECISION_BACKEND STREQUAL "CPP_DEC_FLOAT")
            set(USE_CPP_DEC_FLOAT_BACKEND ON)
            add_compile_definitions(USE_CPP_DEC_FLOAT_BACKEND)
            message(STATUS "Backend: cpp_dec_float (header-only, fixed precision)")
        elseif(HIGH_PRECISION_BACKEND STREQUAL "QUADMATH")
            set(USE_QUADMATH_BACKEND ON)
            add_compile_definitions(USE_QUADMATH_BACKEND)
            message(STATUS "Backend: quadmath (native __float128, 128-bit)")
        endif()

        # Add quadmath support if available
        # Set CMake variable for config.h generation AND compile definition for IDE support
        if(QUADMATH_FOUND)
            set(ENABLE_QUADMATH ON)
            add_compile_definitions(ENABLE_QUADMATH)
            message(STATUS "Quadmath support: ENABLED")
        elseif(ENABLE_QUADMATH AND NOT QUADMATH_FOUND)
            message(WARNING "Quadmath support requested but library not found")
        endif()

        # Add include directories
        include_directories(${HIGH_PRECISION_INCLUDE_DIRS})
    endif()
endif()

if(NOT ENABLE_HIGH_PRECISION)
    message(STATUS "High-precision support: DISABLED (using double precision only)")
endif()

# Generate config.h from config.h.in
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/config.h"
    @ONLY
)

# Add generated include directory
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

#==============================================================================
# Library: polynomial_solver (unified library for users)
#==============================================================================

# Create a unified library that users can link against
add_library(polynomial_solver STATIC
    src/polynomial.cpp
    src/geometry.cpp
    src/de_casteljau.cpp
    src/differentiation.cpp
    src/result_refiner.cpp
    src/solver.cpp
)

# Add high-precision source files if enabled
if(ENABLE_HIGH_PRECISION)
    target_sources(polynomial_solver PRIVATE
        src/precision_conversion.cpp
        src/polynomial_hp.cpp
        src/differentiation_hp.cpp
        src/result_refiner_hp.cpp
    )
endif()

# Set include directories for the library
target_include_directories(polynomial_solver PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link high-precision libraries if enabled
if(ENABLE_HIGH_PRECISION)
    # Add include directories for high-precision libraries
    if(HIGH_PRECISION_INCLUDE_DIRS)
        target_include_directories(polynomial_solver PUBLIC ${HIGH_PRECISION_INCLUDE_DIRS})
    endif()

    # Link libraries
    target_link_libraries(polynomial_solver PUBLIC ${HIGH_PRECISION_LIBRARIES})

    if(ENABLE_QUADMATH AND QUADMATH_FOUND)
        target_link_libraries(polynomial_solver PUBLIC ${QUADMATH_LIBRARY})
    endif()
endif()

# Install library and headers
install(TARGETS polynomial_solver
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(FILES
    include/polynomial_solver.h
    include/polynomial.h
    include/solver.h
    include/result_refiner.h
    include/differentiation.h
    include/geometry.h
    include/de_casteljau.h
    DESTINATION include
)

#==============================================================================
# Legacy individual libraries (for backward compatibility)
#==============================================================================

# Build polynomial library
add_library(polynomial STATIC
    src/polynomial.cpp
)

# Link HP libraries if enabled (polynomial.cpp uses HP for conversion)
if(ENABLE_HIGH_PRECISION)
    if(HIGH_PRECISION_INCLUDE_DIRS)
        target_include_directories(polynomial PUBLIC ${HIGH_PRECISION_INCLUDE_DIRS})
    endif()
    target_link_libraries(polynomial PUBLIC ${HIGH_PRECISION_LIBRARIES})

    # Also need polynomial_hp for fromPowerHP function
    target_sources(polynomial PRIVATE
        src/polynomial_hp.cpp
        src/precision_conversion.cpp
    )
endif()

# Build geometry library
add_library(geometry STATIC
    src/geometry.cpp
)

# Build de_casteljau library
add_library(de_casteljau STATIC
    src/de_casteljau.cpp
)
target_link_libraries(de_casteljau polynomial geometry)

# Build differentiation library
add_library(differentiation STATIC
    src/differentiation.cpp
)
target_link_libraries(differentiation polynomial)

# Build result_refiner library
add_library(result_refiner STATIC
    src/result_refiner.cpp
)
target_link_libraries(result_refiner polynomial differentiation)

# Build solver library
add_library(solver STATIC
    src/solver.cpp
)
target_link_libraries(solver polynomial geometry de_casteljau)

#==============================================================================
# Executables and Tools
#==============================================================================

# Build main executable
add_executable(polynomial_solver_app
    src/main.cpp
)
target_link_libraries(polynomial_solver_app polynomial_solver)

# Build refinement from dumps tool (main tool for users)
add_executable(refine_from_dumps
    tools/refine_from_dumps.cpp
)
target_link_libraries(refine_from_dumps polynomial_solver)

# Legacy tools (for backward compatibility)
add_executable(test_refiner_wilkinson
    tools/test_refiner_wilkinson.cpp
)
target_link_libraries(test_refiner_wilkinson polynomial_solver)

add_executable(test_multiplicity_detection
    tools/test_multiplicity_detection.cpp
)
target_link_libraries(test_multiplicity_detection polynomial_solver)

add_executable(test_1d_refinement
    tools/test_1d_refinement.cpp
)
target_link_libraries(test_1d_refinement polynomial_solver)

# Optional: Enable testing
enable_testing()
add_subdirectory(tests)

# Optional: Build examples
option(BUILD_EXAMPLES "Build example programs" ON)
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
    message(STATUS "Examples: ENABLED")
else()
    message(STATUS "Examples: DISABLED")
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "========== Build Configuration ==========")
message(STATUS "Project: ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "=========================================")
message(STATUS "")

#==============================================================================
# IDE Support - Trigger VS Code IntelliSense Reload
#==============================================================================

# Touch .vscode/c_cpp_properties.json to trigger VS Code file watcher
# This forces IntelliSense to reload compile_commands.json without manual reload
if(EXISTS "${CMAKE_SOURCE_DIR}/.vscode/c_cpp_properties.json")
    file(TOUCH "${CMAKE_SOURCE_DIR}/.vscode/c_cpp_properties.json")
endif()

#==============================================================================
# Generate Playground Makefiles
#==============================================================================
# Generate simple Makefiles in playground/ and playground/research_tests/
# so developers can quickly compile test files with: make test.cpp
# These Makefiles automatically stay in sync with library dependencies.

# Collect all link flags needed for playground executables
# Note: We link statically against polynomial_solver.a which already contains
# all necessary object files, so we only need external libraries (MPFR, GMP, quadmath)
set(PLAYGROUND_LINK_FLAGS "-L${CMAKE_BINARY_DIR}/lib -lpolynomial_solver -lm")

# Add high-precision library flags if enabled
if(ENABLE_HIGH_PRECISION)
    # Note: Boost is header-only for multiprecision, no need to link boost_system
    # Add all libraries that were found during high-precision detection
    foreach(LIB ${HIGH_PRECISION_LIBRARIES})
        get_filename_component(LIB_NAME ${LIB} NAME_WE)
        string(REGEX REPLACE "^lib" "" LIB_NAME ${LIB_NAME})
        set(PLAYGROUND_LINK_FLAGS "${PLAYGROUND_LINK_FLAGS} -l${LIB_NAME}")
    endforeach()

    # Add quadmath if available
    if(QUADMATH_FOUND AND QUADMATH_LIBRARY)
        set(PLAYGROUND_LINK_FLAGS "${PLAYGROUND_LINK_FLAGS} -lquadmath")
    endif()

    set(PLAYGROUND_HP_FLAGS "CXXFLAGS += -DENABLE_HIGH_PRECISION")
else()
    set(PLAYGROUND_HP_FLAGS "# High precision disabled")
endif()

# Generate playground/Makefile
configure_file(
    "${CMAKE_SOURCE_DIR}/playground/Makefile.in"
    "${CMAKE_SOURCE_DIR}/playground/Makefile"
    @ONLY
)

# Generate playground/research_tests/Makefile
configure_file(
    "${CMAKE_SOURCE_DIR}/playground/research_tests/Makefile.in"
    "${CMAKE_SOURCE_DIR}/playground/research_tests/Makefile"
    @ONLY
)

message(STATUS "Generated playground Makefiles:")
message(STATUS "  - playground/Makefile")
message(STATUS "  - playground/research_tests/Makefile")
message(STATUS "Developers can now use: cd playground && make <name>.cpp")
