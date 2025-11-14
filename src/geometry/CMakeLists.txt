# Geometry module
set(GEOMETRY_SOURCES
    geometry.cpp
)

set(GEOMETRY_HEADERS
    ${PROJECT_SOURCE_DIR}/include/geometry/geometry.h
)

add_library(geometry ${GEOMETRY_SOURCES} ${GEOMETRY_HEADERS})

target_include_directories(geometry PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

# Set properties
set_target_properties(geometry PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)

