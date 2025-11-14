# Solver module (main interface)
set(SOLVER_SOURCES
    solver.cpp
)

set(SOLVER_HEADERS
    ${PROJECT_SOURCE_DIR}/include/solver/solver.h
)

add_library(solver ${SOLVER_SOURCES} ${SOLVER_HEADERS})

target_include_directories(solver PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

# Link with all other modules
target_link_libraries(solver PUBLIC
    polynomial
    geometry
    de_casteljau
)

# Set properties
set_target_properties(solver PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)

# Create main executable
add_executable(polynomial_solver_app main.cpp)
target_link_libraries(polynomial_solver_app PRIVATE solver)

