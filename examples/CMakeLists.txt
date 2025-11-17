# Examples for polynomial solver

# Circle-Ellipse Intersection Example
add_executable(example_circle_ellipse circle_ellipse_intersection.cpp)
target_link_libraries(example_circle_ellipse polynomial solver geometry de_casteljau)
target_include_directories(example_circle_ellipse PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Set output directory
set_target_properties(example_circle_ellipse PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add to examples target (optional, for building all examples at once)
add_custom_target(examples)
add_dependencies(examples example_circle_ellipse)

