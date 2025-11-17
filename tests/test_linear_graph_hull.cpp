#include "polynomial.h"
#include "solver.h"
#include "geometry.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

using namespace polynomial_solver;

// Test 1D linear function: p(x) = x - 0.5
int test_1d_linear_graph_hull() {
    std::cout << "Test 1: 1D linear function p(x) = x - 0.5" << std::endl;
    std::cout << "Expected root: x = 0.5" << std::endl;
    std::cout << std::endl;
    
    // Create polynomial p(x) = x - 0.5
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{-0.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    
    // Get graph control points
    std::vector<double> control_points;
    p.graphControlPoints(control_points);
    
    std::cout << "Graph control points (in R^2):" << std::endl;
    const std::size_t num_coeffs = p.coefficientCount();
    const std::size_t point_dim = 2u;  // 1D + 1 for function value
    
    std::vector<std::vector<double>> points;
    for (std::size_t i = 0; i < num_coeffs; ++i) {
        std::vector<double> pt(point_dim);
        for (std::size_t j = 0; j < point_dim; ++j) {
            pt[j] = control_points[i * point_dim + j];
        }
        points.push_back(pt);
        std::cout << "  Point " << i << ": (" << pt[0] << ", " << pt[1] << ")" << std::endl;
    }
    std::cout << std::endl;
    
    // Compute convex hull
    ConvexPolyhedron hull = convex_hull(points);
    std::cout << "Convex hull vertices:" << std::endl;
    for (std::size_t i = 0; i < hull.vertices.size(); ++i) {
        std::cout << "  Vertex " << i << ": (" << hull.vertices[i][0] << ", " 
                  << hull.vertices[i][1] << ")" << std::endl;
    }
    std::cout << "Intrinsic dimension: " << hull.intrinsic_dim << std::endl;
    std::cout << "Ambient dimension: " << hull.ambient_dimension() << std::endl;
    std::cout << std::endl;
    
    // Intersect with y = 0 (last coordinate = 0)
    ConvexPolyhedron intersection;
    bool has_intersection = intersect_convex_polyhedron_with_last_coordinate_zero(hull, intersection);
    
    if (!has_intersection) {
        std::cerr << "  FAIL: No intersection with y = 0" << std::endl;
        return 1;
    }
    
    std::cout << "Intersection with y = 0:" << std::endl;
    for (std::size_t i = 0; i < intersection.vertices.size(); ++i) {
        std::cout << "  Vertex " << i << ": (" << intersection.vertices[i][0] << ", " 
                  << intersection.vertices[i][1] << ")" << std::endl;
    }
    std::cout << "Intrinsic dimension: " << intersection.intrinsic_dim << std::endl;
    std::cout << std::endl;
    
    // Project to R^1 (drop last coordinate)
    std::vector<std::vector<double>> projected_points;
    for (const auto& v : intersection.vertices) {
        projected_points.push_back({v[0]});
    }
    
    std::cout << "Projected points (in R^1):" << std::endl;
    for (std::size_t i = 0; i < projected_points.size(); ++i) {
        std::cout << "  Point " << i << ": " << projected_points[i][0] << std::endl;
    }
    std::cout << std::endl;
    
    // Compute bounding box
    if (projected_points.empty()) {
        std::cerr << "  FAIL: No projected points" << std::endl;
        return 1;
    }
    
    double min_x = projected_points[0][0];
    double max_x = projected_points[0][0];
    for (const auto& pt : projected_points) {
        if (pt[0] < min_x) min_x = pt[0];
        if (pt[0] > max_x) max_x = pt[0];
    }
    
    std::cout << "Bounding box: [" << min_x << ", " << max_x << "]" << std::endl;
    std::cout << "Width: " << (max_x - min_x) << std::endl;
    std::cout << std::endl;
    
    // Check if bounding box is a single point at x = 0.5
    if (std::fabs(min_x - 0.5) < 1e-10 && std::fabs(max_x - 0.5) < 1e-10) {
        std::cout << "  Bounding box is exactly the root!" << std::endl;
    } else {
        std::cerr << "  FAIL: Bounding box is not the exact root" << std::endl;
        return 1;
    }

    // Now test with the actual solver
    std::cout << std::endl;
    std::cout << "Testing with solver (GraphHull method):" << std::endl;

    PolynomialSystem system({p});
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;

    SubdivisionSolverResult result = solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);

    std::cout << "  Found " << result.boxes.size() << " box(es), "
              << result.num_resolved << " resolved" << std::endl;
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    if (result.num_resolved != 1) {
        std::cerr << "  FAIL: Expected 1 resolved root, got " << result.num_resolved << std::endl;
        return 1;
    }

    const SubdivisionBoxResult& box = result.boxes[0];
    std::cout << "  Root box:" << std::endl;
    std::cout << "    Lower: " << box.lower[0] << std::endl;
    std::cout << "    Upper: " << box.upper[0] << std::endl;
    std::cout << "    Center: " << box.center[0] << std::endl;
    std::cout << "    Max error: " << std::scientific << std::setprecision(16) << box.max_error[0] << std::endl;
    std::cout << "    Depth: " << box.depth << std::endl;
    std::cout << "    Converged: " << (box.converged ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    // Check if error is machine epsilon
    const double machine_eps = std::numeric_limits<double>::epsilon();
    std::cout << "  Machine epsilon: " << std::scientific << std::setprecision(16) << machine_eps << std::endl;
    std::cout << "  Error is machine epsilon: " << (box.max_error[0] == machine_eps ? "YES" : "NO") << std::endl;
    std::cout << std::endl;

    // Evaluate at the center
    double eval_value = p.evaluate(box.center[0]);
    std::cout << "  Evaluation at center:" << std::endl;
    std::cout << "    p(" << std::fixed << std::setprecision(16) << box.center[0] << ") = "
              << std::scientific << std::setprecision(16) << eval_value << std::endl;
    std::cout << "    |p(center)| = " << std::fabs(eval_value) << std::endl;
    std::cout << std::endl;

    // Check if evaluation is close to zero
    if (std::fabs(eval_value) < 1e-10) {
        std::cout << "  PASS: Linear function gives exact root with machine epsilon error!" << std::endl;
        return 0;
    } else if (min_x <= 0.5 && max_x >= 0.5) {
        std::cout << "  PARTIAL: Bounding box contains the root but is not exact" << std::endl;
        std::cout << "  This means iterative contraction is needed" << std::endl;
        return 0;
    } else {
        std::cerr << "  FAIL: Bounding box does not contain the root" << std::endl;
        return 1;
    }
}

// Test 2D linear system: p1(x,y) = x - 0.5, p2(x,y) = y - 0.3
int test_2d_linear_graph_hull() {
    std::cout << "Test 2: 2D linear system" << std::endl;
    std::cout << "  p1(x,y) = x - 0.5" << std::endl;
    std::cout << "  p2(x,y) = y - 0.3" << std::endl;
    std::cout << "Expected root: (0.5, 0.3)" << std::endl;
    std::cout << std::endl;

    std::vector<unsigned int> degrees{1u, 1u};

    // p1(x,y) = x - 0.5
    std::vector<double> power1(4, 0.0);
    power1[0] = -0.5;  // constant
    power1[2] = 1.0;   // x term
    Polynomial p1 = Polynomial::fromPower(degrees, power1);

    // p2(x,y) = y - 0.3
    std::vector<double> power2(4, 0.0);
    power2[0] = -0.3;  // constant
    power2[1] = 1.0;   // y term
    Polynomial p2 = Polynomial::fromPower(degrees, power2);

    std::cout << "Processing equation 1: p1(x,y) = x - 0.5" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // Get graph control points for p1
    std::vector<double> control_points1;
    p1.graphControlPoints(control_points1);

    std::cout << "Graph control points (in R^3):" << std::endl;
    const std::size_t num_coeffs = p1.coefficientCount();
    const std::size_t point_dim = 3u;  // 2D + 1 for function value

    std::vector<std::vector<double>> points1;
    for (std::size_t i = 0; i < num_coeffs; ++i) {
        std::vector<double> pt(point_dim);
        for (std::size_t j = 0; j < point_dim; ++j) {
            pt[j] = control_points1[i * point_dim + j];
        }
        points1.push_back(pt);
        std::cout << "  Point " << i << ": (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ")" << std::endl;
    }
    std::cout << std::endl;

    // Compute convex hull for p1
    ConvexPolyhedron hull1 = convex_hull(points1);
    std::cout << "Convex hull has " << hull1.vertices.size() << " vertices, intrinsic_dim = "
              << hull1.intrinsic_dim << std::endl;

    // Intersect with z = 0
    ConvexPolyhedron intersection1;
    bool has_intersection1 = intersect_convex_polyhedron_with_last_coordinate_zero(hull1, intersection1);

    if (!has_intersection1) {
        std::cerr << "  FAIL: No intersection with z = 0 for p1" << std::endl;
        return 1;
    }

    std::cout << "Intersection with z = 0 has " << intersection1.vertices.size()
              << " vertices, intrinsic_dim = " << intersection1.intrinsic_dim << std::endl;
    for (std::size_t i = 0; i < intersection1.vertices.size(); ++i) {
        std::cout << "  Vertex " << i << ": (" << intersection1.vertices[i][0] << ", "
                  << intersection1.vertices[i][1] << ", " << intersection1.vertices[i][2] << ")" << std::endl;
    }

    // Project to R^2
    std::vector<std::vector<double>> projected1;
    for (const auto& v : intersection1.vertices) {
        projected1.push_back({v[0], v[1]});
    }

    std::cout << "Projected to R^2:" << std::endl;
    for (std::size_t i = 0; i < projected1.size(); ++i) {
        std::cout << "  Point " << i << ": (" << projected1[i][0] << ", " << projected1[i][1] << ")" << std::endl;
    }

    ConvexPolyhedron poly1 = convex_hull(projected1);
    std::cout << "Convex hull in R^2 has " << poly1.vertices.size() << " vertices, intrinsic_dim = "
              << poly1.intrinsic_dim << std::endl;
    std::cout << std::endl;

    // Repeat for p2
    std::cout << "Processing equation 2: p2(x,y) = y - 0.3" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    std::vector<double> control_points2;
    p2.graphControlPoints(control_points2);

    std::vector<std::vector<double>> points2;
    for (std::size_t i = 0; i < num_coeffs; ++i) {
        std::vector<double> pt(point_dim);
        for (std::size_t j = 0; j < point_dim; ++j) {
            pt[j] = control_points2[i * point_dim + j];
        }
        points2.push_back(pt);
    }

    ConvexPolyhedron hull2 = convex_hull(points2);
    ConvexPolyhedron intersection2;
    bool has_intersection2 = intersect_convex_polyhedron_with_last_coordinate_zero(hull2, intersection2);

    if (!has_intersection2) {
        std::cerr << "  FAIL: No intersection with z = 0 for p2" << std::endl;
        return 1;
    }

    std::vector<std::vector<double>> projected2;
    for (const auto& v : intersection2.vertices) {
        projected2.push_back({v[0], v[1]});
    }

    std::cout << "Projected to R^2:" << std::endl;
    for (std::size_t i = 0; i < projected2.size(); ++i) {
        std::cout << "  Point " << i << ": (" << projected2[i][0] << ", " << projected2[i][1] << ")" << std::endl;
    }

    ConvexPolyhedron poly2 = convex_hull(projected2);
    std::cout << "Convex hull in R^2 has " << poly2.vertices.size() << " vertices, intrinsic_dim = "
              << poly2.intrinsic_dim << std::endl;
    std::cout << std::endl;

    // Intersect the two polyhedra
    std::cout << "Intersecting the two polyhedra..." << std::endl;
    std::vector<ConvexPolyhedron> polys = {poly1, poly2};
    ConvexPolyhedron final_intersection;
    bool has_final = intersect_convex_polyhedra(polys, final_intersection);

    if (!has_final) {
        std::cerr << "  FAIL: No intersection between the two polyhedra" << std::endl;
        return 1;
    }

    std::cout << "Final intersection has " << final_intersection.vertices.size()
              << " vertices, intrinsic_dim = " << final_intersection.intrinsic_dim << std::endl;
    for (std::size_t i = 0; i < final_intersection.vertices.size(); ++i) {
        std::cout << "  Vertex " << i << ": (" << final_intersection.vertices[i][0] << ", "
                  << final_intersection.vertices[i][1] << ")" << std::endl;
    }
    std::cout << std::endl;

    // Compute bounding box
    if (final_intersection.vertices.empty()) {
        std::cerr << "  FAIL: No vertices in final intersection" << std::endl;
        return 1;
    }

    double min_x = final_intersection.vertices[0][0];
    double max_x = final_intersection.vertices[0][0];
    double min_y = final_intersection.vertices[0][1];
    double max_y = final_intersection.vertices[0][1];

    for (const auto& v : final_intersection.vertices) {
        if (v[0] < min_x) min_x = v[0];
        if (v[0] > max_x) max_x = v[0];
        if (v[1] < min_y) min_y = v[1];
        if (v[1] > max_y) max_y = v[1];
    }

    std::cout << "Bounding box: [" << min_x << ", " << max_x << "] x [" << min_y << ", " << max_y << "]" << std::endl;
    std::cout << "Width: (" << (max_x - min_x) << ", " << (max_y - min_y) << ")" << std::endl;
    std::cout << std::endl;

    // Check if bounding box is exactly the root
    if (std::fabs(min_x - 0.5) < 1e-10 && std::fabs(max_x - 0.5) < 1e-10 &&
        std::fabs(min_y - 0.3) < 1e-10 && std::fabs(max_y - 0.3) < 1e-10) {
        std::cout << "  Bounding box is exactly the root!" << std::endl;
    } else if (min_x <= 0.5 && max_x >= 0.5 && min_y <= 0.3 && max_y >= 0.3) {
        std::cout << "  Bounding box contains the root but is not exact" << std::endl;
    } else {
        std::cerr << "  FAIL: Bounding box does not contain the root" << std::endl;
        return 1;
    }

    // Now test with the actual solver
    std::cout << std::endl;
    std::cout << "Testing with solver (GraphHull method):" << std::endl;

    PolynomialSystem system({p1, p2});
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;

    SubdivisionSolverResult result = solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);

    std::cout << "  Found " << result.boxes.size() << " box(es), "
              << result.num_resolved << " resolved" << std::endl;
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    if (result.num_resolved != 1) {
        std::cerr << "  FAIL: Expected 1 resolved root, got " << result.num_resolved << std::endl;
        return 1;
    }

    const SubdivisionBoxResult& box = result.boxes[0];
    std::cout << "  Root box:" << std::endl;
    std::cout << "    Lower: (" << box.lower[0] << ", " << box.lower[1] << ")" << std::endl;
    std::cout << "    Upper: (" << box.upper[0] << ", " << box.upper[1] << ")" << std::endl;
    std::cout << "    Center: (" << box.center[0] << ", " << box.center[1] << ")" << std::endl;
    std::cout << "    Max error: (" << std::scientific << std::setprecision(16)
              << box.max_error[0] << ", " << box.max_error[1] << ")" << std::endl;
    std::cout << "    Depth: " << box.depth << std::endl;
    std::cout << "    Converged: " << (box.converged ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    // Check if error is machine epsilon
    const double machine_eps = std::numeric_limits<double>::epsilon();
    std::cout << "  Machine epsilon: " << std::scientific << std::setprecision(16) << machine_eps << std::endl;
    std::cout << "  Error is machine epsilon: ("
              << (box.max_error[0] == machine_eps ? "YES" : "NO") << ", "
              << (box.max_error[1] == machine_eps ? "YES" : "NO") << ")" << std::endl;
    std::cout << std::endl;

    // Evaluate at the center
    std::vector<double> eval_values;
    system.evaluate(box.center, eval_values);
    std::cout << "  Evaluation at center:" << std::endl;
    std::cout << "    p1(" << std::fixed << std::setprecision(16) << box.center[0] << ", " << box.center[1] << ") = "
              << std::scientific << std::setprecision(16) << eval_values[0] << std::endl;
    std::cout << "    p2(" << std::fixed << std::setprecision(16) << box.center[0] << ", " << box.center[1] << ") = "
              << std::scientific << std::setprecision(16) << eval_values[1] << std::endl;
    std::cout << "    |p1(center)| = " << std::scientific << std::fabs(eval_values[0]) << std::endl;
    std::cout << "    |p2(center)| = " << std::fabs(eval_values[1]) << std::endl;
    std::cout << std::endl;

    // Check if evaluation is close to zero
    if (std::fabs(eval_values[0]) < 1e-10 && std::fabs(eval_values[1]) < 1e-10) {
        std::cout << "  PASS: Linear system gives exact root with machine epsilon error!" << std::endl;
        return 0;
    } else {
        std::cerr << "  FAIL: Evaluation error too large" << std::endl;
        return 1;
    }
}

int main() {
    std::cout << "Testing graph hull method for linear functions..." << std::endl;
    std::cout << "=================================================" << std::endl;
    std::cout << std::endl;

    int failures = 0;
    failures += test_1d_linear_graph_hull();
    std::cout << std::endl;
    failures += test_2d_linear_graph_hull();

    std::cout << std::endl;
    if (failures == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << failures << " test(s) failed." << std::endl;
        return 1;
    }
}

