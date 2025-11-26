cmake_minimum_required(VERSION 3.15)
project(PolynomialSolver VERSION 1.0.0 LANGUAGES CXX)

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
option(DISABLE_TEMPLATES "Disable template-based implementation (use fallback)" OFF)
option(ENABLE_QUADMATH "Enable quadmath (__float128) support" OFF)

#==============================================================================
# High-Precision Configuration
#==============================================================================

if(ENABLE_HIGH_PRECISION)
    # Find high-precision dependencies
    include(FindHighPrecision)

    if(NOT HIGH_PRECISION_FOUND)
        message(FATAL_ERROR "High-precision support requested but dependencies not found")
    endif()

    # Set preprocessor definitions
    add_compile_definitions(ENABLE_HIGH_PRECISION)

    # Determine if templates should be used
    if(DISABLE_TEMPLATES)
        message(STATUS "Template support: DISABLED (using fallback implementation)")
    else()
        add_compile_definitions(USE_TEMPLATES)
        message(STATUS "Template support: ENABLED")
    endif()

    # Set backend type
    if(HIGH_PRECISION_BACKEND STREQUAL "MPFR")
        add_compile_definitions(USE_MPFR_BACKEND)
    elseif(HIGH_PRECISION_BACKEND STREQUAL "CPP_DEC_FLOAT")
        add_compile_definitions(USE_CPP_DEC_FLOAT_BACKEND)
    endif()

    # Add quadmath support if requested and available
    if(ENABLE_QUADMATH AND QUADMATH_FOUND)
        add_compile_definitions(ENABLE_QUADMATH)
        message(STATUS "Quadmath support: ENABLED")
    elseif(ENABLE_QUADMATH AND NOT QUADMATH_FOUND)
        message(WARNING "Quadmath support requested but library not found")
    endif()

    # Add include directories
    include_directories(${HIGH_PRECISION_INCLUDE_DIRS})
else()
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

# Set include directories for the library
target_include_directories(polynomial_solver PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link high-precision libraries if enabled
if(ENABLE_HIGH_PRECISION)
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

