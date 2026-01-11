/**
 * @file test_geometry_base.cpp
 * @brief Tests for templated geometry module
 *
 * Tests that:
 * 1. Double precision gives identical results to original geometry module
 * 2. High-precision types work correctly
 * 3. Multidimensional operations work
 */

#include "core/geometry_base.h"
#include "core/geometry.h"  // Original for comparison

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef ENABLE_HIGH_PRECISION
#include "hp/high_precision_types.h"
#endif

using namespace polynomial_solver;

// Tolerance for comparisons
constexpr double DOUBLE_TOL = 1e-12;

bool approx_equal(double a, double b, double tol = DOUBLE_TOL) {
    return std::fabs(a - b) < tol;
}

//=============================================================================
// Test 1: Bounding box - double vs original
//=============================================================================

void test_bounding_box_double() {
    std::cout << "  test_bounding_box_double... ";

    // Create points
    std::vector<std::vector<double>> points = {
        {0.0, 0.0}, {1.0, 0.0}, {0.5, 1.0}, {0.2, 0.3}
    };

    // Original
    ConvexPolyhedron orig_poly;
    orig_poly.vertices = points;
    orig_poly.intrinsic_dim = 2;
    ConvexPolyhedronBox orig_box = bounding_box(orig_poly);

    // Templated
    ConvexPolyhedronBase<double> tmpl_poly;
    tmpl_poly.vertices = points;
    tmpl_poly.intrinsic_dim = 2;
    ConvexPolyhedronBoxBase<double> tmpl_box = bounding_box_impl(tmpl_poly);

    // Compare
    assert(orig_box.dimension() == tmpl_box.dimension());
    for (std::size_t i = 0; i < orig_box.dimension(); ++i) {
        assert(approx_equal(orig_box.min_coords[i], tmpl_box.min_coords[i]));
        assert(approx_equal(orig_box.max_coords[i], tmpl_box.max_coords[i]));
    }

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 2: 2D Convex hull - double vs original
//=============================================================================

void test_convex_hull_2d_double() {
    std::cout << "  test_convex_hull_2d_double... ";

    std::vector<std::vector<double>> points = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}
    };

    // Original
    ConvexPolyhedron orig_hull = convex_hull(points);

    // Templated
    ConvexPolyhedronBase<double> tmpl_hull = convex_hull_impl(points);

    // Both should have 4 vertices (the corners)
    assert(orig_hull.vertices.size() == 4);
    assert(tmpl_hull.vertices.size() == 4);
    assert(orig_hull.intrinsic_dim == tmpl_hull.intrinsic_dim);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 3: Hyperplane intersection - double vs original
//=============================================================================

void test_hyperplane_intersection_double() {
    std::cout << "  test_hyperplane_intersection_double... ";

    // Square straddling y=0
    std::vector<std::vector<double>> points = {
        {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}
    };

    ConvexPolyhedron orig_poly = convex_hull(points);
    ConvexPolyhedronBase<double> tmpl_poly = convex_hull_impl(points);

    ConvexPolyhedron orig_inter;
    ConvexPolyhedronBase<double> tmpl_inter;

    bool orig_has = intersect_convex_polyhedron_with_last_coordinate_zero(orig_poly, orig_inter);
    bool tmpl_has = intersect_convex_polyhedron_with_last_coordinate_zero_impl(tmpl_poly, tmpl_inter);

    assert(orig_has == tmpl_has);
    assert(orig_has);

    // Both should produce a line segment from (-1, 0) to (1, 0)
    assert(orig_inter.vertices.size() >= 2);
    assert(tmpl_inter.vertices.size() >= 2);

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 4: Tolerance traits
//=============================================================================

void test_tolerance_traits() {
    std::cout << "  test_tolerance_traits... ";

    double eps_double = GeometryToleranceTraits<double>::epsilon();
    assert(eps_double > 0);
    assert(eps_double < 1e-10);

    float eps_float = GeometryToleranceTraits<float>::epsilon();
    assert(eps_float > 0);
    assert(eps_float > eps_double);  // float has less precision

    std::cout << "PASSED (double eps=" << eps_double << ", float eps=" << eps_float << ")\n";
}

//=============================================================================
// Test 5: Multidimensional bounding box
//=============================================================================

void test_multidim_bounding_box() {
    std::cout << "  test_multidim_bounding_box... ";

    // 4D points
    std::vector<std::vector<double>> points = {
        {0.0, 0.0, 0.0, 0.0},
        {1.0, 2.0, 3.0, 4.0},
        {0.5, 1.5, 2.5, 3.5}
    };

    ConvexPolyhedronBase<double> poly;
    poly.vertices = points;
    auto box = bounding_box_impl(poly);

    assert(box.dimension() == 4);
    assert(approx_equal(box.min_coords[0], 0.0));
    assert(approx_equal(box.max_coords[0], 1.0));
    assert(approx_equal(box.min_coords[3], 0.0));
    assert(approx_equal(box.max_coords[3], 4.0));

    std::cout << "PASSED\n";
}

//=============================================================================
// Test 6: High-precision geometry (if available)
//=============================================================================

#ifdef USE_MPFR_BACKEND
void test_high_precision_geometry() {
    std::cout << "  test_high_precision_geometry... ";

    using mpreal = boost::multiprecision::mpfr_float;

    // Set high precision
    mpreal::default_precision(256);

    // Check tolerance is appropriately small
    mpreal eps = GeometryToleranceTraits<mpreal>::epsilon();
    assert(eps > mpreal(0));
    assert(eps < mpreal("1e-50"));  // Should be very small with 256 bits

    // 2D points in high precision
    std::vector<std::vector<mpreal>> points = {
        {mpreal(0), mpreal(0)},
        {mpreal(1), mpreal(0)},
        {mpreal(1), mpreal(1)},
        {mpreal(0), mpreal(1)}
    };

    // Convex hull
    ConvexPolyhedronBase<mpreal> hull = convex_hull_impl(points);
    assert(hull.vertices.size() == 4);
    assert(hull.intrinsic_dim == 2);

    // Bounding box
    auto box = bounding_box_impl(hull);
    assert(box.dimension() == 2);
    assert(box.min_coords[0] == mpreal(0));
    assert(box.max_coords[0] == mpreal(1));

    std::cout << "PASSED (mpfr eps=" << eps << ")\n";
}
#endif

#ifdef ENABLE_QUADMATH
void test_quadmath_geometry() {
    std::cout << "  test_quadmath_geometry... ";

    // Check tolerance
    __float128 eps = GeometryToleranceTraits<__float128>::epsilon();
    // __float128 has ~34 decimal digits precision
    __float128 zero = strtoflt128("0", nullptr);
    __float128 eps_limit = strtoflt128("1e-25", nullptr);
    assert(eps > zero);
    assert(eps < eps_limit);

    // Simple 2D test
    __float128 v0 = strtoflt128("0", nullptr);
    __float128 v1 = strtoflt128("1", nullptr);
    __float128 v05 = strtoflt128("0.5", nullptr);

    std::vector<std::vector<__float128>> points = {
        {v0, v0},
        {v1, v0},
        {v05, v1}
    };

    ConvexPolyhedronBase<__float128> hull = convex_hull_impl(points);
    assert(hull.vertices.size() == 3);

    std::cout << "PASSED\n";
}
#endif

//=============================================================================
// Main
//=============================================================================

int main() {
    std::cout << "Running geometry_base tests...\n\n";

    std::cout << "Double precision (vs original):\n";
    test_bounding_box_double();
    test_convex_hull_2d_double();
    test_hyperplane_intersection_double();

    std::cout << "\nTolerance traits:\n";
    test_tolerance_traits();

    std::cout << "\nMultidimensional:\n";
    test_multidim_bounding_box();

#ifdef USE_MPFR_BACKEND
    std::cout << "\nHigh-precision (MPFR):\n";
    test_high_precision_geometry();
#endif

#ifdef ENABLE_QUADMATH
    std::cout << "\nQuadmath:\n";
    test_quadmath_geometry();
#endif

    std::cout << "\n=== All geometry_base tests PASSED ===\n";
    return 0;
}

