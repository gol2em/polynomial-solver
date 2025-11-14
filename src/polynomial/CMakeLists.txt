# Polynomial module
set(POLYNOMIAL_SOURCES
    polynomial.cpp
)

set(POLYNOMIAL_HEADERS
    ${PROJECT_SOURCE_DIR}/include/polynomial/polynomial.h
)

add_library(polynomial ${POLYNOMIAL_SOURCES} ${POLYNOMIAL_HEADERS})

target_include_directories(polynomial PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

# Set properties
set_target_properties(polynomial PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)

