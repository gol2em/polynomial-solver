# Tests directory

# High-precision types test (only if high precision is enabled)
if(ENABLE_HIGH_PRECISION)
    add_executable(test_high_precision_types test_high_precision_types.cpp)
    target_link_libraries(test_high_precision_types PRIVATE polynomial_solver)
    add_test(NAME HighPrecisionTypesTest COMMAND test_high_precision_types)

    add_executable(test_polynomial_hp test_polynomial_hp.cpp)
    target_link_libraries(test_polynomial_hp PRIVATE polynomial_solver)
    add_test(NAME PolynomialHPTest COMMAND test_polynomial_hp)

    add_executable(test_differentiation_hp test_differentiation_hp.cpp)
    target_link_libraries(test_differentiation_hp PRIVATE polynomial_solver)
    add_test(NAME DifferentiationHPTest COMMAND test_differentiation_hp)

    add_executable(test_result_refiner_hp test_result_refiner_hp.cpp)
    target_link_libraries(test_result_refiner_hp PRIVATE polynomial_solver)
    add_test(NAME ResultRefinerHPTest COMMAND test_result_refiner_hp)
endif()

add_executable(test_polynomial_conversion test_polynomial_conversion.cpp)
target_link_libraries(test_polynomial_conversion PRIVATE polynomial de_casteljau)
add_test(NAME PolynomialConversionTest COMMAND test_polynomial_conversion)

add_executable(test_polynomial_graph test_polynomial_graph.cpp)
target_link_libraries(test_polynomial_graph PRIVATE solver polynomial de_casteljau geometry)
add_test(NAME PolynomialGraphTest COMMAND test_polynomial_graph)

add_executable(test_polynomial_system test_polynomial_system.cpp)
target_link_libraries(test_polynomial_system PRIVATE solver polynomial de_casteljau geometry)
add_test(NAME PolynomialSystemExampleTest COMMAND test_polynomial_system)

add_executable(test_solver_subdivision test_solver_subdivision.cpp)
target_link_libraries(test_solver_subdivision PRIVATE solver polynomial de_casteljau geometry)
add_test(NAME SolverSubdivisionTest COMMAND test_solver_subdivision)
add_executable(test_geometry_convex test_geometry_convex.cpp)
target_link_libraries(test_geometry_convex PRIVATE geometry)
add_test(NAME GeometryConvexTest COMMAND test_geometry_convex)

add_executable(test_2d_hyperplane_intersection test_2d_hyperplane_intersection.cpp)
target_link_libraries(test_2d_hyperplane_intersection PRIVATE geometry)
add_test(NAME 2DHyperplaneIntersectionTest COMMAND test_2d_hyperplane_intersection)

add_executable(test_2d_convex_hull_robust test_2d_convex_hull_robust.cpp)
target_link_libraries(test_2d_convex_hull_robust PRIVATE geometry)
add_test(NAME 2DConvexHullRobustTest COMMAND test_2d_convex_hull_robust)

add_executable(test_2d_hull_vertex_order test_2d_hull_vertex_order.cpp)
target_link_libraries(test_2d_hull_vertex_order PRIVATE geometry)
add_test(NAME 2DHullVertexOrderTest COMMAND test_2d_hull_vertex_order)

add_executable(test_degenerate_boxes test_degenerate_boxes.cpp)
target_link_libraries(test_degenerate_boxes PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME DegenerateBoxesTest COMMAND test_degenerate_boxes)

add_executable(test_linear_graph_hull test_linear_graph_hull.cpp)
target_link_libraries(test_linear_graph_hull PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME LinearGraphHullTest COMMAND test_linear_graph_hull)

add_executable(test_projected_polyhedral test_projected_polyhedral.cpp)
target_link_libraries(test_projected_polyhedral PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME ProjectedPolyhedralTest COMMAND test_projected_polyhedral)

add_executable(test_method_comparison test_method_comparison.cpp)
target_link_libraries(test_method_comparison PRIVATE solver polynomial geometry de_casteljau)
# Note: This is a comparison tool, not added to ctest

add_executable(test_ellipse_dump test_ellipse_dump.cpp)
target_link_libraries(test_ellipse_dump PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME EllipseDumpTest COMMAND test_ellipse_dump)

add_executable(test_strategies test_strategies.cpp)
target_link_libraries(test_strategies PRIVATE solver polynomial geometry de_casteljau)
# Note: This is a comparison tool for different strategies, not added to ctest

add_executable(test_differentiation test_differentiation.cpp)
target_link_libraries(test_differentiation PRIVATE differentiation polynomial de_casteljau geometry)
add_test(NAME DifferentiationTest COMMAND test_differentiation)

add_executable(test_differentiation_debug test_differentiation_debug.cpp)
target_link_libraries(test_differentiation_debug PRIVATE differentiation polynomial de_casteljau geometry)

add_executable(test_degree_zero_derivative test_degree_zero_derivative.cpp)
target_link_libraries(test_degree_zero_derivative PRIVATE differentiation polynomial de_casteljau geometry)

add_executable(test_result_refiner test_result_refiner.cpp)
target_link_libraries(test_result_refiner PRIVATE result_refiner solver differentiation polynomial de_casteljau geometry)
add_test(NAME ResultRefinerTest COMMAND test_result_refiner)

message(STATUS "Tests directory configured (conversion, graph, system, differentiation, and result refiner tests enabled)")

