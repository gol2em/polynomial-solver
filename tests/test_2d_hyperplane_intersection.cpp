#include "core/geometry.h"
#include <iostream>
#include <cmath>

using namespace polynomial_solver;

static bool approx_equal(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

// Test 1: Polygon intersecting y=0 to get a line segment
int test_polygon_line_segment_intersection() {
    std::cout << "Test 1: Polygon intersecting y=0 to get line segment..." << std::endl;
    
    // Square from (-1, -1) to (1, 1) in 2D
    // Should intersect y=0 at segment from (-1, 0) to (1, 0)
    std::vector<std::vector<double>> pts;
    pts.push_back({-1.0, -1.0});
    pts.push_back({ 1.0, -1.0});
    pts.push_back({ 1.0,  1.0});
    pts.push_back({-1.0,  1.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;
    
    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }
    
    // Should get 2 points: (-1, 0) and (1, 0)
    if (inter.vertices.size() != 2) {
        std::cerr << "  FAIL: Expected 2 vertices, got " << inter.vertices.size() << std::endl;
        return 1;
    }
    
    // Check intrinsic dimension is 1 (line segment)
    if (inter.intrinsic_dim != 1) {
        std::cerr << "  FAIL: Expected intrinsic_dim=1, got " << inter.intrinsic_dim << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 2: Polygon intersecting y=0 at a single point
int test_polygon_point_intersection() {
    std::cout << "Test 2: Polygon intersecting y=0 at single point..." << std::endl;
    
    // Triangle with one vertex at (0, 0) and others above y=0
    std::vector<std::vector<double>> pts;
    pts.push_back({-1.0,  1.0});
    pts.push_back({ 1.0,  1.0});
    pts.push_back({ 0.0,  0.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;
    
    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }
    
    // Should get 1 point: (0, 0)
    if (inter.vertices.size() != 1) {
        std::cerr << "  FAIL: Expected 1 vertex, got " << inter.vertices.size() << std::endl;
        return 1;
    }
    
    // Check intrinsic dimension is 0 (point)
    if (inter.intrinsic_dim != 0) {
        std::cerr << "  FAIL: Expected intrinsic_dim=0, got " << inter.intrinsic_dim << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 3: Polygon entirely above y=0 (no intersection)
int test_polygon_no_intersection() {
    std::cout << "Test 3: Polygon entirely above y=0 (no intersection)..." << std::endl;
    
    // Square from (0, 1) to (2, 3)
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 1.0});
    pts.push_back({2.0, 1.0});
    pts.push_back({2.0, 3.0});
    pts.push_back({0.0, 3.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;
    
    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (ok) {
        std::cerr << "  FAIL: Expected no intersection" << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 4: Line segment crossing y=0
int test_segment_point_intersection() {
    std::cout << "Test 4: Line segment crossing y=0..." << std::endl;
    
    // Segment from (0, -1) to (0, 1)
    // Should intersect at (0, 0)
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, -1.0});
    pts.push_back({0.0,  1.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;
    
    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }
    
    // Should get 1 point: (0, 0)
    if (inter.vertices.size() != 1) {
        std::cerr << "  FAIL: Expected 1 vertex, got " << inter.vertices.size() << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 5: Line segment on y=0
int test_segment_on_line() {
    std::cout << "Test 5: Line segment on y=0..." << std::endl;
    
    // Segment from (-1, 0) to (1, 0)
    std::vector<std::vector<double>> pts;
    pts.push_back({-1.0, 0.0});
    pts.push_back({ 1.0, 0.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;
    
    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }

    // Should get 2 points: (-1, 0) and (1, 0)
    if (inter.vertices.size() != 2) {
        std::cerr << "  FAIL: Expected 2 vertices, got " << inter.vertices.size() << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 6: Point on y=0
int test_point_on_line() {
    std::cout << "Test 6: Point on y=0..." << std::endl;

    // Single point at (0, 0)
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});

    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;

    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }

    // Should get 1 point: (0, 0)
    if (inter.vertices.size() != 1) {
        std::cerr << "  FAIL: Expected 1 vertex, got " << inter.vertices.size() << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 7: Point above y=0 (no intersection)
int test_point_no_intersection() {
    std::cout << "Test 7: Point above y=0 (no intersection)..." << std::endl;

    // Single point at (0, 1)
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 1.0});

    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;

    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (ok) {
        std::cerr << "  FAIL: Expected no intersection" << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 8: Polygon with edge on y=0
int test_polygon_edge_on_line() {
    std::cout << "Test 8: Polygon with edge on y=0..." << std::endl;

    // Triangle with one edge on y=0
    std::vector<std::vector<double>> pts;
    pts.push_back({-1.0, 0.0});
    pts.push_back({ 1.0, 0.0});
    pts.push_back({ 0.0, 1.0});

    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;

    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }

    // Should get 2 points: (-1, 0) and (1, 0)
    if (inter.vertices.size() != 2) {
        std::cerr << "  FAIL: Expected 2 vertices, got " << inter.vertices.size() << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 9: Complex polygon (hexagon) intersecting y=0
int test_complex_polygon() {
    std::cout << "Test 9: Hexagon intersecting y=0..." << std::endl;

    // Regular hexagon centered at origin with radius 2
    // Vertices at angles 0, 60, 120, 180, 240, 300 degrees
    std::vector<std::vector<double>> pts;
    const double pi = 3.14159265358979323846;
    for (int i = 0; i < 6; ++i) {
        double angle = i * pi / 3.0;
        pts.push_back({2.0 * std::cos(angle), 2.0 * std::sin(angle)});
    }

    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;

    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }

    // Should get 2 points where hexagon crosses y=0
    // At angles 0 and 180 degrees: (2, 0) and (-2, 0)
    if (inter.vertices.size() != 2) {
        std::cerr << "  FAIL: Expected 2 vertices, got " << inter.vertices.size() << std::endl;
        for (const auto& v : inter.vertices) {
            std::cerr << "    (" << v[0] << ", " << v[1] << ")" << std::endl;
        }
        return 1;
    }

    // Check the points are approximately (-2, 0) and (2, 0)
    bool found_left = false, found_right = false;
    for (const auto& v : inter.vertices) {
        if (approx_equal(v[0], -2.0) && approx_equal(v[1], 0.0)) {
            found_left = true;
        }
        if (approx_equal(v[0], 2.0) && approx_equal(v[1], 0.0)) {
            found_right = true;
        }
    }

    if (!found_left || !found_right) {
        std::cerr << "  FAIL: Expected points at (-2, 0) and (2, 0)" << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 10: Thin polygon (nearly degenerate)
int test_thin_polygon() {
    std::cout << "Test 10: Thin polygon crossing y=0..." << std::endl;

    // Very thin rectangle crossing y=0
    std::vector<std::vector<double>> pts;
    pts.push_back({-1.0, -0.5});
    pts.push_back({ 1.0, -0.5});
    pts.push_back({ 1.0,  0.5});
    pts.push_back({-1.0,  0.5});

    ConvexPolyhedron poly = convex_hull(pts);
    ConvexPolyhedron inter;

    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter);
    if (!ok) {
        std::cerr << "  FAIL: Expected intersection" << std::endl;
        return 1;
    }

    // Should get 2 points: (-1, 0) and (1, 0)
    if (inter.vertices.size() != 2) {
        std::cerr << "  FAIL: Expected 2 vertices, got " << inter.vertices.size() << std::endl;
        return 1;
    }

    std::cout << "  PASS" << std::endl;
    return 0;
}

int main(int argc, char** argv) {
    std::cout << "Running 2D hyperplane intersection tests..." << std::endl;
    std::cout << std::endl;

    int failures = 0;
    failures += test_polygon_line_segment_intersection();
    failures += test_polygon_point_intersection();
    failures += test_polygon_no_intersection();
    failures += test_segment_point_intersection();
    failures += test_segment_on_line();
    failures += test_point_on_line();
    failures += test_point_no_intersection();
    failures += test_polygon_edge_on_line();
    failures += test_complex_polygon();
    failures += test_thin_polygon();

    std::cout << std::endl;
    if (failures == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << failures << " test(s) failed." << std::endl;
        return 1;
    }
}
