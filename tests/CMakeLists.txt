#==============================================================================
# Tests Directory
#==============================================================================
# This file defines all unit tests for the polynomial solver library.
# All tests are added to CTest and run with `ctest` or `make test`.
#
# Test organization:
# 1. Polynomial Module Tests
# 2. Differentiation Module Tests
# 3. Geometry Module Tests
# 4. Solver Module Tests
# 5. Result Refiner Module Tests
# 6. High-Precision Module Tests (if ENABLE_HIGH_PRECISION)
# 7. Integration Tests
#==============================================================================

#==============================================================================
# 1. Polynomial Module Tests
#==============================================================================

# Test: Power ↔ Bernstein conversion
add_executable(test_polynomial_conversion test_polynomial_conversion.cpp)
target_link_libraries(test_polynomial_conversion PRIVATE polynomial de_casteljau)
add_test(NAME PolynomialConversionTest COMMAND test_polynomial_conversion)

# Test: Dual representation (power + Bernstein coefficients)
add_executable(test_dual_representation test_dual_representation.cpp)
target_link_libraries(test_dual_representation PRIVATE polynomial_solver)
add_test(NAME DualRepresentationTest COMMAND test_dual_representation)

# Test: Graph control points generation
add_executable(test_polynomial_graph test_polynomial_graph.cpp)
target_link_libraries(test_polynomial_graph PRIVATE solver polynomial de_casteljau geometry)
add_test(NAME PolynomialGraphTest COMMAND test_polynomial_graph)

# Test: PolynomialSystem class and De Casteljau subdivision
add_executable(test_polynomial_system test_polynomial_system.cpp)
target_link_libraries(test_polynomial_system PRIVATE solver polynomial de_casteljau geometry)
add_test(NAME PolynomialSystemExampleTest COMMAND test_polynomial_system)

#==============================================================================
# 2. Differentiation Module Tests
#==============================================================================

# Test: Bernstein differentiation (1D, 2D, gradient, Hessian)
add_executable(test_differentiation test_differentiation.cpp)
target_link_libraries(test_differentiation PRIVATE differentiation polynomial de_casteljau geometry)
add_test(NAME DifferentiationTest COMMAND test_differentiation)

# Test: Edge case - degree 0 after differentiation
add_executable(test_degree_zero_derivative test_degree_zero_derivative.cpp)
target_link_libraries(test_degree_zero_derivative PRIVATE differentiation polynomial de_casteljau geometry)
add_test(NAME DegreeZeroDerivativeTest COMMAND test_degree_zero_derivative)

#==============================================================================
# 3. Geometry Module Tests
#==============================================================================

# Test: Convex hull and bounding box
add_executable(test_geometry_convex test_geometry_convex.cpp)
target_link_libraries(test_geometry_convex PRIVATE geometry)
add_test(NAME GeometryConvexTest COMMAND test_geometry_convex)

# Test: Polygon-hyperplane intersection
add_executable(test_2d_hyperplane_intersection test_2d_hyperplane_intersection.cpp)
target_link_libraries(test_2d_hyperplane_intersection PRIVATE geometry)
add_test(NAME 2DHyperplaneIntersectionTest COMMAND test_2d_hyperplane_intersection)

# Test: Robustness (collinear, duplicate points)
add_executable(test_2d_convex_hull_robust test_2d_convex_hull_robust.cpp)
target_link_libraries(test_2d_convex_hull_robust PRIVATE geometry)
add_test(NAME 2DConvexHullRobustTest COMMAND test_2d_convex_hull_robust)

# Test: CCW vertex ordering
add_executable(test_2d_hull_vertex_order test_2d_hull_vertex_order.cpp)
target_link_libraries(test_2d_hull_vertex_order PRIVATE geometry)
add_test(NAME 2DHullVertexOrderTest COMMAND test_2d_hull_vertex_order)


#==============================================================================
# 4. Solver Module Tests
#==============================================================================

# Test: Subdivision solver with GraphHull method
add_executable(test_solver_subdivision test_solver_subdivision.cpp)
target_link_libraries(test_solver_subdivision PRIVATE solver polynomial de_casteljau geometry)
add_test(NAME SolverSubdivisionTest COMMAND test_solver_subdivision)

# Test: ProjectedPolyhedral root bounding method
add_executable(test_projected_polyhedral test_projected_polyhedral.cpp)
target_link_libraries(test_projected_polyhedral PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME ProjectedPolyhedralTest COMMAND test_projected_polyhedral)

# Test: GraphHull method on linear functions (1D, 2D)
add_executable(test_linear_graph_hull test_linear_graph_hull.cpp)
target_link_libraries(test_linear_graph_hull PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME LinearGraphHullTest COMMAND test_linear_graph_hull)

# Test: Degeneracy detection (0D, 1D boxes)
add_executable(test_degenerate_boxes test_degenerate_boxes.cpp)
target_link_libraries(test_degenerate_boxes PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME DegenerateBoxesTest COMMAND test_degenerate_boxes)

# Test: Geometry dump feature
add_executable(test_ellipse_dump test_ellipse_dump.cpp)
target_link_libraries(test_ellipse_dump PRIVATE solver polynomial geometry de_casteljau)
add_test(NAME EllipseDumpTest COMMAND test_ellipse_dump)

#==============================================================================
# 5. Result Refiner Module Tests
#==============================================================================

# Test: Newton refinement (double precision)
add_executable(test_result_refiner test_result_refiner.cpp)
target_link_libraries(test_result_refiner PRIVATE polynomial_solver)
add_test(NAME ResultRefinerTest COMMAND test_result_refiner)

#==============================================================================
# 6. High-Precision Module Tests (only if ENABLE_HIGH_PRECISION)
#==============================================================================

if(ENABLE_HIGH_PRECISION)
    # Test: High-precision type definitions and precision control
    add_executable(test_high_precision_types test_high_precision_types.cpp)
    target_link_libraries(test_high_precision_types PRIVATE polynomial_solver)
    add_test(NAME HighPrecisionTypesTest COMMAND test_high_precision_types)

    # Test: High-precision polynomial operations
    add_executable(test_polynomial_hp test_polynomial_hp.cpp)
    target_link_libraries(test_polynomial_hp PRIVATE polynomial_solver)
    add_test(NAME PolynomialHPTest COMMAND test_polynomial_hp)

    # Test: High-precision differentiation
    add_executable(test_differentiation_hp test_differentiation_hp.cpp)
    target_link_libraries(test_differentiation_hp PRIVATE polynomial_solver)
    add_test(NAME DifferentiationHPTest COMMAND test_differentiation_hp)

    # Test: High-precision Newton refinement
    add_executable(test_result_refiner_hp test_result_refiner_hp.cpp)
    target_link_libraries(test_result_refiner_hp PRIVATE polynomial_solver)
    add_test(NAME ResultRefinerHPTest COMMAND test_result_refiner_hp)

    # Test: Interval arithmetic for rigorous bounds
    add_executable(test_interval_arithmetic test_interval_arithmetic.cpp)
    target_link_libraries(test_interval_arithmetic PRIVATE polynomial_solver)
    add_test(NAME IntervalArithmeticTest COMMAND test_interval_arithmetic)
endif()

#==============================================================================
# 7. Template-Based Implementation Tests (Tier 3)
#==============================================================================
# These tests verify the templated implementations (PolynomialBase<T>, etc.)
# They are always built when high-precision is enabled (USE_TEMPLATES is ON by default)

if(ENABLE_HIGH_PRECISION)
    # Test: PolynomialBase<double> basic functionality
    add_executable(test_polynomial_base test_polynomial_base.cpp)
    target_link_libraries(test_polynomial_base PRIVATE polynomial_solver)
    add_test(NAME PolynomialBaseTest COMMAND test_polynomial_base)

    # Test: PolynomialBase comprehensive (vs Polynomial, HP types)
    add_executable(test_polynomial_base_comprehensive test_polynomial_base_comprehensive.cpp)
    target_link_libraries(test_polynomial_base_comprehensive PRIVATE polynomial_solver)
    add_test(NAME PolynomialBaseComprehensiveTest COMMAND test_polynomial_base_comprehensive)

    # Test: PolynomialBase solver compatibility (shadow tests)
    add_executable(test_polynomial_base_solver_compat test_polynomial_base_solver_compat.cpp)
    target_link_libraries(test_polynomial_base_solver_compat PRIVATE polynomial_solver)
    add_test(NAME PolynomialBaseSolverCompatTest COMMAND test_polynomial_base_solver_compat)

    # Test: Templated differentiation
    add_executable(test_differentiation_base test_differentiation_base.cpp)
    target_link_libraries(test_differentiation_base PRIVATE polynomial_solver)
    add_test(NAME DifferentiationBaseTest COMMAND test_differentiation_base)

    # Test: Templated result refiner
    add_executable(test_result_refiner_base test_result_refiner_base.cpp)
    target_link_libraries(test_result_refiner_base PRIVATE polynomial_solver)
    add_test(NAME ResultRefinerBaseTest COMMAND test_result_refiner_base)

    # Test: Region tracking for subdivision
    add_executable(test_region_tracking test_region_tracking.cpp)
    target_link_libraries(test_region_tracking PRIVATE polynomial_solver)
    add_test(NAME RegionTrackingTest COMMAND test_region_tracking)

    # Test: Templated solver
    add_executable(test_solver_base test_solver_base.cpp)
    target_link_libraries(test_solver_base PRIVATE polynomial_solver)
    add_test(NAME SolverBaseTest COMMAND test_solver_base)

    # Test: Templated geometry
    add_executable(test_geometry_base test_geometry_base.cpp)
    target_link_libraries(test_geometry_base PRIVATE polynomial_solver)
    add_test(NAME GeometryBaseTest COMMAND test_geometry_base)
endif()

#==============================================================================
# 8. Integration Tests
#==============================================================================

# Test: Full workflow (solve → refine)
add_executable(test_solver_refiner_workflow test_solver_refiner_workflow.cpp)
target_link_libraries(test_solver_refiner_workflow PRIVATE polynomial_solver)
add_test(NAME SolverRefinerWorkflowTest COMMAND test_solver_refiner_workflow)

# Test: Multidimensional polynomials (3D+)
add_executable(test_multidim test_multidim.cpp)
target_link_libraries(test_multidim PRIVATE polynomial_solver)
add_test(NAME MultidimTest COMMAND test_multidim)

# Test: Solver vs SolverBase comparison (verify identical results)
add_executable(test_solver_comparison test_solver_comparison.cpp)
target_link_libraries(test_solver_comparison PRIVATE polynomial_solver)
add_test(NAME SolverComparisonTest COMMAND test_solver_comparison)

# Test: Multidimensional Newton refinement with convergence analysis
if(ENABLE_HIGH_PRECISION)
    add_executable(test_newton_multidim test_newton_multidim.cpp)
    target_link_libraries(test_newton_multidim PRIVATE polynomial_solver)
    add_test(NAME NewtonMultidimTest COMMAND test_newton_multidim)
endif()

#==============================================================================
# Summary
#==============================================================================

if(ENABLE_HIGH_PRECISION)
    message(STATUS "Tests configured: 30 tests (all added to CTest)")
else()
    message(STATUS "Tests configured: 17 tests (all added to CTest)")
endif()
message(STATUS "  - Polynomial module: 4 tests")
message(STATUS "  - Differentiation module: 2 tests")
message(STATUS "  - Geometry module: 4 tests")
message(STATUS "  - Solver module: 5 tests")
message(STATUS "  - Result refiner module: 1 test")
if(ENABLE_HIGH_PRECISION)
    message(STATUS "  - High-precision module: 5 tests")
    message(STATUS "  - Tier 3 template module: 8 tests")
else()
    message(STATUS "  - High-precision module: 0 tests (ENABLE_HIGH_PRECISION=OFF)")
    message(STATUS "  - Tier 3 template module: 0 tests (ENABLE_HIGH_PRECISION=OFF)")
endif()
message(STATUS "  - Integration tests: 1 test")
message(STATUS "Research/debug tests moved to playground/research_tests/")
