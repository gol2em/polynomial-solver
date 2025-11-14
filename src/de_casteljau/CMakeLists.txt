# De Casteljau subdivision module
set(DE_CASTELJAU_SOURCES
    de_casteljau.cpp
)

set(DE_CASTELJAU_HEADERS
    ${PROJECT_SOURCE_DIR}/include/de_casteljau/de_casteljau.h
)

add_library(de_casteljau ${DE_CASTELJAU_SOURCES} ${DE_CASTELJAU_HEADERS})

target_include_directories(de_casteljau PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

# Link with polynomial and geometry modules
target_link_libraries(de_casteljau PUBLIC
    polynomial
    geometry
)

# Set properties
set_target_properties(de_casteljau PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)

