# Examples for polynomial solver
# All examples use the unified polynomial_solver library

# Simple Cubic Example (2-line workflow demonstration)
add_executable(example_simple_cubic simple_cubic.cpp)
target_link_libraries(example_simple_cubic polynomial_solver)

# Circle-Ellipse Intersection Example (2D)
add_executable(example_circle_ellipse circle_ellipse_intersection.cpp)
target_link_libraries(example_circle_ellipse polynomial_solver)

# 1D High Multiplicity Root Example
add_executable(example_multiplicity_1d multiplicity_1d_roots.cpp)
target_link_libraries(example_multiplicity_1d polynomial_solver)

# 1D Wilkinson Polynomial Example (ill-conditioned)
add_executable(example_wilkinson_1d wilkinson_1d_roots.cpp)
target_link_libraries(example_wilkinson_1d polynomial_solver)

# Set output directory
set_target_properties(
    example_simple_cubic
    example_circle_ellipse
    example_multiplicity_1d
    example_wilkinson_1d
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add to examples target (optional, for building all examples at once)
add_custom_target(examples)
add_dependencies(examples
    example_simple_cubic
    example_circle_ellipse
    example_multiplicity_1d
    example_wilkinson_1d
)

