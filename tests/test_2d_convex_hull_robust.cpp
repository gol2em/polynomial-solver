#include "core/geometry.h"
#include <iostream>
#include <cmath>

using namespace polynomial_solver;

static bool approx_equal(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

// Test 1: Collinear points (should produce line segment)
int test_collinear_points() {
    std::cout << "Test 1: Collinear points..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({1.0, 1.0});
    pts.push_back({2.0, 2.0});
    pts.push_back({0.5, 0.5});
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    // Should get 2 points: (0, 0) and (2, 2)
    if (poly.vertices.size() != 2) {
        std::cerr << "  FAIL: Expected 2 vertices, got " << poly.vertices.size() << std::endl;
        return 1;
    }
    
    // Check intrinsic dimension is 1
    if (poly.intrinsic_dim != 1) {
        std::cerr << "  FAIL: Expected intrinsic_dim=1, got " << poly.intrinsic_dim << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 2: Duplicate points
int test_duplicate_points() {
    std::cout << "Test 2: Duplicate points..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({1.0, 0.0});
    pts.push_back({0.0, 1.0});
    pts.push_back({0.0, 0.0});  // duplicate
    pts.push_back({1.0, 0.0});  // duplicate
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    // Should get 3 unique points
    if (poly.vertices.size() != 3) {
        std::cerr << "  FAIL: Expected 3 vertices, got " << poly.vertices.size() << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 3: Single point
int test_single_point() {
    std::cout << "Test 3: Single point..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({1.0, 2.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    if (poly.vertices.size() != 1) {
        std::cerr << "  FAIL: Expected 1 vertex, got " << poly.vertices.size() << std::endl;
        return 1;
    }
    
    // Check intrinsic dimension is 0
    if (poly.intrinsic_dim != 0) {
        std::cerr << "  FAIL: Expected intrinsic_dim=0, got " << poly.intrinsic_dim << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 4: Two points
int test_two_points() {
    std::cout << "Test 4: Two points..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({1.0, 1.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    if (poly.vertices.size() != 2) {
        std::cerr << "  FAIL: Expected 2 vertices, got " << poly.vertices.size() << std::endl;
        return 1;
    }
    
    // Check intrinsic dimension is 1
    if (poly.intrinsic_dim != 1) {
        std::cerr << "  FAIL: Expected intrinsic_dim=1, got " << poly.intrinsic_dim << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 5: Nearly collinear points (numerical stability)
int test_nearly_collinear() {
    std::cout << "Test 5: Nearly collinear points..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({1.0, 0.0});
    pts.push_back({2.0, 0.0});
    pts.push_back({1.0, 1e-13});  // Very small deviation
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    // Should be treated as collinear (2 points) or as a very thin triangle (4 points)
    // Either is acceptable depending on epsilon tolerance
    if (poly.vertices.size() < 2) {
        std::cerr << "  FAIL: Expected at least 2 vertices, got " << poly.vertices.size() << std::endl;
        return 1;
    }
    
    std::cout << "  PASS (got " << poly.vertices.size() << " vertices)" << std::endl;
    return 0;
}

// Test 6: Square (proper 2D polygon)
int test_square() {
    std::cout << "Test 6: Square..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({1.0, 0.0});
    pts.push_back({1.0, 1.0});
    pts.push_back({0.0, 1.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    if (poly.vertices.size() != 4) {
        std::cerr << "  FAIL: Expected 4 vertices, got " << poly.vertices.size() << std::endl;
        return 1;
    }
    
    // Check intrinsic dimension is 2
    if (poly.intrinsic_dim != 2) {
        std::cerr << "  FAIL: Expected intrinsic_dim=2, got " << poly.intrinsic_dim << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

int main() {
    std::cout << "Running 2D convex hull robustness tests..." << std::endl;
    std::cout << std::endl;
    
    int failures = 0;
    failures += test_collinear_points();
    failures += test_duplicate_points();
    failures += test_single_point();
    failures += test_two_points();
    failures += test_nearly_collinear();
    failures += test_square();
    
    std::cout << std::endl;
    if (failures == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << failures << " test(s) failed." << std::endl;
        return 1;
    }
}

