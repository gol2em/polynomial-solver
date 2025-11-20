# Examples for polynomial solver

# Circle-Ellipse Intersection Example
add_executable(example_circle_ellipse circle_ellipse_intersection.cpp)
target_link_libraries(example_circle_ellipse polynomial solver geometry de_casteljau)
target_include_directories(example_circle_ellipse PRIVATE ${CMAKE_SOURCE_DIR}/include)

# 1D Cubic Polynomial Example
add_executable(example_cubic_1d cubic_1d_roots.cpp)
target_link_libraries(example_cubic_1d polynomial solver geometry de_casteljau)
target_include_directories(example_cubic_1d PRIVATE ${CMAKE_SOURCE_DIR}/include)

# 1D High Multiplicity Root Example
add_executable(example_multiplicity_1d multiplicity_1d_roots.cpp)
target_link_libraries(example_multiplicity_1d polynomial solver geometry de_casteljau)
target_include_directories(example_multiplicity_1d PRIVATE ${CMAKE_SOURCE_DIR}/include)

# 1D Wilkinson Polynomial Example
add_executable(example_wilkinson_1d wilkinson_1d_roots.cpp)
target_link_libraries(example_wilkinson_1d polynomial solver geometry de_casteljau)
target_include_directories(example_wilkinson_1d PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Set output directory
set_target_properties(example_circle_ellipse example_cubic_1d example_multiplicity_1d example_wilkinson_1d PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add to examples target (optional, for building all examples at once)
add_custom_target(examples)
add_dependencies(examples example_circle_ellipse example_cubic_1d example_multiplicity_1d example_wilkinson_1d)

