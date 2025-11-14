cmake_minimum_required(VERSION 3.15)
project(PolynomialSolver VERSION 1.0.0 LANGUAGES CXX)

# Set C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

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

# Build solver library
add_library(solver STATIC
    src/solver.cpp
)
target_link_libraries(solver polynomial geometry de_casteljau)

# Build main executable
add_executable(polynomial_solver_app
    src/main.cpp
)
target_link_libraries(polynomial_solver_app solver polynomial geometry de_casteljau)

# Optional: Enable testing
enable_testing()
add_subdirectory(tests)

# Print configuration summary
message(STATUS "")
message(STATUS "========== Build Configuration ==========")
message(STATUS "Project: ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "=========================================")
message(STATUS "")

