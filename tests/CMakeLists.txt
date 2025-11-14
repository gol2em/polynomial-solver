# Tests directory

add_executable(test_polynomial_conversion test_polynomial_conversion.cpp)
target_link_libraries(test_polynomial_conversion PRIVATE polynomial de_casteljau)
add_test(NAME PolynomialConversionTest COMMAND test_polynomial_conversion)

add_executable(test_polynomial_graph test_polynomial_graph.cpp)
target_link_libraries(test_polynomial_graph PRIVATE solver polynomial de_casteljau geometry)
add_test(NAME PolynomialGraphTest COMMAND test_polynomial_graph)

message(STATUS "Tests directory configured (polynomial conversion and graph tests enabled)")

