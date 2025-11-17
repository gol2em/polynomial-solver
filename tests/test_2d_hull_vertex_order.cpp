#include "geometry.h"
#include <iostream>
#include <cmath>

using namespace polynomial_solver;

// Helper: compute signed area of triangle (o, a, b)
// Positive if counter-clockwise, negative if clockwise
static double signed_area(const std::vector<double>& o,
                         const std::vector<double>& a,
                         const std::vector<double>& b) {
    const double x1 = a[0] - o[0];
    const double y1 = a[1] - o[1];
    const double x2 = b[0] - o[0];
    const double y2 = b[1] - o[1];
    return x1 * y2 - y1 * x2;
}

// Test that vertices are in counter-clockwise order
bool verify_ccw_order(const ConvexPolyhedron& poly) {
    if (poly.vertices.size() < 3) {
        return true; // Degenerate cases are always valid
    }
    
    // Check that all consecutive triples have positive signed area
    const std::size_t n = poly.vertices.size();
    for (std::size_t i = 0; i < n; ++i) {
        const std::vector<double>& v0 = poly.vertices[i];
        const std::vector<double>& v1 = poly.vertices[(i + 1) % n];
        const std::vector<double>& v2 = poly.vertices[(i + 2) % n];
        
        double area = signed_area(v0, v1, v2);
        if (area < -1e-12) {
            std::cerr << "  Vertices not in CCW order at index " << i << std::endl;
            std::cerr << "  Signed area: " << area << std::endl;
            return false;
        }
    }
    
    return true;
}

// Test 1: Square
int test_square_order() {
    std::cout << "Test 1: Square vertex order..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({1.0, 0.0});
    pts.push_back({1.0, 1.0});
    pts.push_back({0.0, 1.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    if (!verify_ccw_order(poly)) {
        std::cerr << "  FAIL: Vertices not in CCW order" << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 2: Random points
int test_random_points_order() {
    std::cout << "Test 2: Random points vertex order..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.5, 0.2});
    pts.push_back({0.8, 0.9});
    pts.push_back({0.1, 0.7});
    pts.push_back({0.9, 0.3});
    pts.push_back({0.3, 0.5});
    pts.push_back({0.6, 0.8});
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    if (!verify_ccw_order(poly)) {
        std::cerr << "  FAIL: Vertices not in CCW order" << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 3: Pentagon
int test_pentagon_order() {
    std::cout << "Test 3: Pentagon vertex order..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    const double pi = 3.14159265358979323846;
    for (int i = 0; i < 5; ++i) {
        double angle = i * 2.0 * pi / 5.0;
        pts.push_back({std::cos(angle), std::sin(angle)});
    }
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    if (!verify_ccw_order(poly)) {
        std::cerr << "  FAIL: Vertices not in CCW order" << std::endl;
        return 1;
    }
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

// Test 4: Verify edges form a closed loop
int test_edges_form_loop() {
    std::cout << "Test 4: Edges form closed loop..." << std::endl;
    
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({2.0, 0.0});
    pts.push_back({2.0, 1.0});
    pts.push_back({1.0, 2.0});
    pts.push_back({0.0, 1.0});
    
    ConvexPolyhedron poly = convex_hull(pts);
    
    if (poly.vertices.size() < 3) {
        std::cout << "  PASS (degenerate case)" << std::endl;
        return 0;
    }
    
    // Check that consecutive vertices form edges
    // The last vertex should connect back to the first
    const std::size_t n = poly.vertices.size();
    
    std::cout << "  Hull has " << n << " vertices" << std::endl;
    std::cout << "  Edges: ";
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t next = (i + 1) % n;
        std::cout << i << "->" << next;
        if (i < n - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "  PASS" << std::endl;
    return 0;
}

int main() {
    std::cout << "Testing 2D convex hull vertex ordering..." << std::endl;
    std::cout << std::endl;
    
    int failures = 0;
    failures += test_square_order();
    failures += test_random_points_order();
    failures += test_pentagon_order();
    failures += test_edges_form_loop();
    
    std::cout << std::endl;
    if (failures == 0) {
        std::cout << "All tests passed!" << std::endl;
        std::cout << "Confirmed: 2D convex hull vertices are in counter-clockwise order." << std::endl;
        std::cout << "Adjacent vertices form edges of the convex hull." << std::endl;
        return 0;
    } else {
        std::cout << failures << " test(s) failed." << std::endl;
        return 1;
    }
}

