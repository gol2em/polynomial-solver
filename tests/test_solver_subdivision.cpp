#include "polynomial.h"
#include "solver.h"
#include "geometry.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

using namespace polynomial_solver;

namespace {

bool approx_equal(double a, double b,
                  double eps = std::numeric_limits<double>::epsilon() * 1000.0)
{
    return std::fabs(a - b) <= eps;
}

int test_subdivision_solver_uniform_grid()
{
    // 1D polynomial p(t) = t on [0,1]. The current subdivision solver does not
    // yet use function values in its decisions, but we attach a simple
    // normalized system so that the workflow is exercised end-to-end.
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{0.0, 1.0}; // p(t) = t
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    PolynomialSystem system({p});
    if (system.dimension() != 1u || system.equationCount() != 1u) {
        std::cerr << "Subdivision solver: unexpected system shape" << '\n';
        return 1;
    }

    Solver solver;

    SubdivisionConfig config;
    config.crit_shrink_factor = 0.8;   // Bounding step is a no-op for now.
    config.min_box_width = 0.25;       // Stop once boxes reach width 1/4.
    config.max_depth = 10u;            // Large enough to not be active here.

    const std::vector<SubdivisionBoxResult> boxes =
        solver.subdivisionSolve(system, config, RootBoundingMethod::None);

    if (boxes.size() != 4u) {
        std::cerr << "Subdivision solver: expected 4 final boxes, got "
                  << boxes.size() << '\n';
        return 1;
    }

    const double expected_width = 0.25;
    const double eps = std::numeric_limits<double>::epsilon() * 1000.0;

    for (std::size_t i = 0; i < boxes.size(); ++i) {
        const SubdivisionBoxResult& box = boxes[i];

        if (box.lower.size() != 1u || box.upper.size() != 1u) {
            std::cerr << "Subdivision solver: box dimension mismatch" << '\n';
            return 1;
        }

        const double width = box.upper[0] - box.lower[0];
        if (!approx_equal(width, expected_width, eps)) {
            std::cerr << "Subdivision solver: unexpected box width " << width
                      << ", expected " << expected_width << '\n';
            return 1;
        }

        if (!box.converged) {
            std::cerr << "Subdivision solver: box not marked as converged" << '\n';
            return 1;
        }

        if (box.depth != 2u) {
            std::cerr << "Subdivision solver: expected depth 2, got "
                      << box.depth << '\n';
            return 1;
        }

        if (box.lower[0] < -eps || box.upper[0] > 1.0 + eps) {
            std::cerr << "Subdivision solver: box outside [0,1]" << '\n';
            return 1;
        }
    }

    return 0;
}

int test_subdivision_solver_graph_hull_1d()
{
    // 1D polynomial p(t) = t - 0.5 on [0,1], which has a root at t = 0.5.
    // Using GraphHull method should contract the box around the root.
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{-0.5, 1.0}; // p(t) = t - 0.5
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    PolynomialSystem system({p});

    Solver solver;

    SubdivisionConfig config;
    config.crit_shrink_factor = 0.8;
    config.min_box_width = 0.05;
    config.max_depth = 10u;

    // First solve with None method to get baseline.
    const std::vector<SubdivisionBoxResult> boxes_none =
        solver.subdivisionSolve(system, config, RootBoundingMethod::None);

    // Then solve with GraphHull method.
    const std::vector<SubdivisionBoxResult> boxes_hull =
        solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);

    // With GraphHull, we should get fewer or equal number of boxes due to contraction.
    if (boxes_hull.empty()) {
        std::cerr << "GraphHull 1D: expected at least one box" << '\n';
        return 1;
    }

    if (boxes_hull.size() > boxes_none.size()) {
        std::cerr << "GraphHull 1D: expected fewer boxes than None method, got "
                  << boxes_hull.size() << " vs " << boxes_none.size() << '\n';
        return 1;
    }

    // Check that at least one box contains the root t = 0.5.
    bool found_root = false;
    for (const SubdivisionBoxResult& box : boxes_hull) {
        if (box.lower[0] <= 0.5 && box.upper[0] >= 0.5) {
            found_root = true;
            break;
        }
    }

    if (!found_root) {
        std::cerr << "GraphHull 1D: no box contains the root at t = 0.5" << '\n';
        return 1;
    }

    std::cout << "GraphHull 1D: None method produced " << boxes_none.size()
              << " boxes, GraphHull produced " << boxes_hull.size() << " boxes" << std::endl;

    return 0;
}

int test_subdivision_solver_graph_hull_1d_debug()
{
    // 1D polynomial p(t) = t - 0.5 on [0,1], which has a root at t = 0.5.
    // Dump all intermediate data for visualization.
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{-0.5, 1.0}; // p(t) = t - 0.5
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    // Get graph control points.
    std::vector<double> control_points;
    p.graphControlPoints(control_points);

    // Open CSV file for dumping data.
    std::ofstream csv("graph_hull_1d_debug.csv");
    csv << "step,type,dim,index,coord0,coord1,coord2\n";

    // Dump control points.
    const std::size_t num_coeffs = p.coefficientCount();
    const std::size_t point_dim = 2u; // 1D + 1 for function value
    for (std::size_t i = 0; i < num_coeffs; ++i) {
        csv << "0,control_point,2," << i << ","
            << control_points[i * point_dim + 0] << ","
            << control_points[i * point_dim + 1] << ",0\n";
    }

    // Build convex hull.
    std::vector<std::vector<double>> points;
    for (std::size_t i = 0; i < num_coeffs; ++i) {
        std::vector<double> pt(point_dim, 0.0);
        for (std::size_t j = 0; j < point_dim; ++j) {
            pt[j] = control_points[i * point_dim + j];
        }
        points.push_back(pt);
    }

    ConvexPolyhedron hull = convex_hull(points);

    // Dump hull vertices.
    for (std::size_t i = 0; i < hull.vertices.size(); ++i) {
        csv << "1,hull_vertex,2," << i << ","
            << hull.vertices[i][0] << ","
            << hull.vertices[i][1] << ",0\n";
    }

    // Intersect with hyperplane x_1 = 0.
    ConvexPolyhedron hyperplane_intersection;
    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(hull, hyperplane_intersection);

    if (ok) {
        for (std::size_t i = 0; i < hyperplane_intersection.vertices.size(); ++i) {
            csv << "2,hyperplane_intersection,2," << i << ","
                << hyperplane_intersection.vertices[i][0] << ","
                << hyperplane_intersection.vertices[i][1] << ",0\n";
        }

        // Compute bounding box.
        ConvexPolyhedronBox bbox = bounding_box(hyperplane_intersection);
        csv << "3,bbox_min,2,0,"
            << bbox.min_coords[0] << ","
            << bbox.min_coords[1] << ",0\n";
        csv << "3,bbox_max,2,0,"
            << bbox.max_coords[0] << ","
            << bbox.max_coords[1] << ",0\n";
    }

    csv.close();

    std::cout << "1D debug data written to graph_hull_1d_debug.csv" << std::endl;

    return 0;
}

int test_subdivision_solver_graph_hull_2d()
{
    // 2D system:
    //   f1(x, y) = x - 0.5
    //   f2(x, y) = y - 0.3
    // Root at (0.5, 0.3).
    std::vector<unsigned int> degrees_x{1u, 0u};
    std::vector<double> power_coeffs_x{-0.5, 1.0}; // f1 = x - 0.5
    Polynomial p1 = Polynomial::fromPower(degrees_x, power_coeffs_x);

    std::vector<unsigned int> degrees_y{0u, 1u};
    std::vector<double> power_coeffs_y{-0.3, 1.0}; // f2 = y - 0.3
    Polynomial p2 = Polynomial::fromPower(degrees_y, power_coeffs_y);

    PolynomialSystem system({p1, p2});

    Solver solver;

    SubdivisionConfig config;
    config.crit_shrink_factor = 0.8;
    config.min_box_width = 0.05;
    config.max_depth = 10u;

    // First solve with None method to get baseline.
    const std::vector<SubdivisionBoxResult> boxes_none =
        solver.subdivisionSolve(system, config, RootBoundingMethod::None);

    // Then solve with GraphHull method.
    const std::vector<SubdivisionBoxResult> boxes_hull =
        solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);

    if (boxes_hull.empty()) {
        std::cerr << "GraphHull 2D: expected at least one box" << '\n';
        return 1;
    }

    // With GraphHull, we should get fewer or equal number of boxes.
    if (boxes_hull.size() > boxes_none.size()) {
        std::cerr << "GraphHull 2D: expected fewer boxes than None method, got "
                  << boxes_hull.size() << " vs " << boxes_none.size() << '\n';
        return 1;
    }

    // Check that at least one box contains the root (0.5, 0.3).
    bool found_root = false;
    for (const SubdivisionBoxResult& box : boxes_hull) {
        if (box.lower[0] <= 0.5 && box.upper[0] >= 0.5 &&
            box.lower[1] <= 0.3 && box.upper[1] >= 0.3) {
            found_root = true;
            break;
        }
    }

    if (!found_root) {
        std::cerr << "GraphHull 2D: no box contains the root at (0.5, 0.3)" << '\n';
        return 1;
    }

    std::cout << "GraphHull 2D: None method produced " << boxes_none.size()
              << " boxes, GraphHull produced " << boxes_hull.size() << " boxes" << std::endl;

    return 0;
}

int test_subdivision_solver_graph_hull_2d_debug()
{
    // 2D system:
    //   f1(x, y) = x - 0.5
    //   f2(x, y) = y - 0.3
    // Root at (0.5, 0.3).
    std::vector<unsigned int> degrees_x{1u, 0u};
    std::vector<double> power_coeffs_x{-0.5, 1.0}; // f1 = x - 0.5
    Polynomial p1 = Polynomial::fromPower(degrees_x, power_coeffs_x);

    std::vector<unsigned int> degrees_y{0u, 1u};
    std::vector<double> power_coeffs_y{-0.3, 1.0}; // f2 = y - 0.3
    Polynomial p2 = Polynomial::fromPower(degrees_y, power_coeffs_y);

    // Open CSV file for dumping data.
    std::ofstream csv("graph_hull_2d_debug.csv");
    csv << "equation,step,type,dim,index,coord0,coord1,coord2\n";

    // Process equation 1.
    std::vector<double> control_points_1;
    p1.graphControlPoints(control_points_1);

    const std::size_t num_coeffs_1 = p1.coefficientCount();
    const std::size_t point_dim = 3u; // 2D + 1 for function value

    for (std::size_t i = 0; i < num_coeffs_1; ++i) {
        csv << "1,0,control_point,3," << i << ","
            << control_points_1[i * point_dim + 0] << ","
            << control_points_1[i * point_dim + 1] << ","
            << control_points_1[i * point_dim + 2] << "\n";
    }

    std::vector<std::vector<double>> points_1;
    for (std::size_t i = 0; i < num_coeffs_1; ++i) {
        std::vector<double> pt(point_dim, 0.0);
        for (std::size_t j = 0; j < point_dim; ++j) {
            pt[j] = control_points_1[i * point_dim + j];
        }
        points_1.push_back(pt);
    }

    ConvexPolyhedron hull_1 = convex_hull(points_1);

    for (std::size_t i = 0; i < hull_1.vertices.size(); ++i) {
        csv << "1,1,hull_vertex,3," << i << ","
            << hull_1.vertices[i][0] << ","
            << hull_1.vertices[i][1] << ","
            << hull_1.vertices[i][2] << "\n";
    }

    // Process equation 2.
    std::vector<double> control_points_2;
    p2.graphControlPoints(control_points_2);

    const std::size_t num_coeffs_2 = p2.coefficientCount();

    for (std::size_t i = 0; i < num_coeffs_2; ++i) {
        csv << "2,0,control_point,3," << i << ","
            << control_points_2[i * point_dim + 0] << ","
            << control_points_2[i * point_dim + 1] << ","
            << control_points_2[i * point_dim + 2] << "\n";
    }

    std::vector<std::vector<double>> points_2;
    for (std::size_t i = 0; i < num_coeffs_2; ++i) {
        std::vector<double> pt(point_dim, 0.0);
        for (std::size_t j = 0; j < point_dim; ++j) {
            pt[j] = control_points_2[i * point_dim + j];
        }
        points_2.push_back(pt);
    }

    ConvexPolyhedron hull_2 = convex_hull(points_2);

    for (std::size_t i = 0; i < hull_2.vertices.size(); ++i) {
        csv << "2,1,hull_vertex,3," << i << ","
            << hull_2.vertices[i][0] << ","
            << hull_2.vertices[i][1] << ","
            << hull_2.vertices[i][2] << "\n";
    }

    // NEW WORKFLOW: Intersect each hull with hyperplane first, then intersect the results.

    // Intersect hull_1 with hyperplane x_2 = 0.
    ConvexPolyhedron hyperplane_1;
    bool ok1 = intersect_convex_polyhedron_with_last_coordinate_zero(hull_1, hyperplane_1);

    if (ok1) {
        for (std::size_t i = 0; i < hyperplane_1.vertices.size(); ++i) {
            csv << "1,2,hyperplane_intersection,3," << i << ","
                << hyperplane_1.vertices[i][0] << ","
                << hyperplane_1.vertices[i][1] << ","
                << hyperplane_1.vertices[i][2] << "\n";
        }
    }

    // Intersect hull_2 with hyperplane x_2 = 0.
    ConvexPolyhedron hyperplane_2;
    bool ok2 = intersect_convex_polyhedron_with_last_coordinate_zero(hull_2, hyperplane_2);

    if (ok2) {
        for (std::size_t i = 0; i < hyperplane_2.vertices.size(); ++i) {
            csv << "2,2,hyperplane_intersection,3," << i << ","
                << hyperplane_2.vertices[i][0] << ","
                << hyperplane_2.vertices[i][1] << ","
                << hyperplane_2.vertices[i][2] << "\n";
        }
    }

    // Project to 2D (drop last coordinate).
    ConvexPolyhedron projected_1, projected_2;
    if (ok1) {
        for (const auto& v : hyperplane_1.vertices) {
            projected_1.vertices.push_back({v[0], v[1]});
            csv << "1,3,projected,2," << (projected_1.vertices.size()-1) << ","
                << v[0] << "," << v[1] << ",0\n";
        }
    }
    if (ok2) {
        for (const auto& v : hyperplane_2.vertices) {
            projected_2.vertices.push_back({v[0], v[1]});
            csv << "2,3,projected,2," << (projected_2.vertices.size()-1) << ","
                << v[0] << "," << v[1] << ",0\n";
        }
    }

    // Now intersect the two projected polyhedra in 2D.
    ConvexPolyhedron intersection;
    bool ok = false;
    if (ok1 && ok2) {
        ok = intersect_convex_polyhedra({projected_1, projected_2}, intersection);
    }

    if (ok) {
        for (std::size_t i = 0; i < intersection.vertices.size(); ++i) {
            csv << "0,4,final_intersection,2," << i << ","
                << intersection.vertices[i][0] << ","
                << intersection.vertices[i][1] << ",0\n";
        }

        // Compute bounding box.
        ConvexPolyhedronBox bbox = bounding_box(intersection);
        csv << "0,5,bbox_min,2,0,"
            << bbox.min_coords[0] << ","
            << bbox.min_coords[1] << ",0\n";
        csv << "0,5,bbox_max,2,0,"
            << bbox.max_coords[0] << ","
            << bbox.max_coords[1] << ",0\n";
    }

    csv.close();

    std::cout << "2D debug data written to graph_hull_2d_debug.csv" << std::endl;

    return 0;
}

} // namespace

int test_subdivision_solver_graph_hull_2d_quadratic()
{
    // 2D system with one quadratic equation:
    //   f1(x, y) = x^2 - 0.25  (root at x = 0.5)
    //   f2(x, y) = y - 0.3     (root at y = 0.3)
    // Root at (0.5, 0.3).

    // f1 = x^2 - 0.25 in power basis with degrees (2, 1)
    // Layout (last dimension fastest): [x^0y^0, x^0y^1, x^1y^0, x^1y^1, x^2y^0, x^2y^1]
    // f1(x,y) = -0.25 + 0*y + 0*x + 0*xy + 1*x^2 + 0*x^2y
    std::vector<unsigned int> degrees_f1{2u, 1u};
    std::vector<double> power_coeffs_f1{-0.25, 0.0, 0.0, 0.0, 1.0, 0.0};
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs_f1);

    // f2 = y - 0.3 in power basis: [-0.3, 1] with degrees (0, 1)
    std::vector<unsigned int> degrees_f2{0u, 1u};
    std::vector<double> power_coeffs_f2{-0.3, 1.0};
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs_f2);

    PolynomialSystem system({p1, p2});

    Solver solver;

    SubdivisionConfig config;
    config.crit_shrink_factor = 0.8;
    config.min_box_width = 0.05;
    config.max_depth = 10u;

    // First solve with None method to get baseline.
    const std::vector<SubdivisionBoxResult> boxes_none =
        solver.subdivisionSolve(system, config, RootBoundingMethod::None);

    // Then solve with GraphHull method.
    const std::vector<SubdivisionBoxResult> boxes_hull =
        solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);

    if (boxes_hull.empty()) {
        std::cerr << "GraphHull 2D quadratic: expected at least one box" << '\n';
        return 1;
    }

    // With GraphHull, we should get fewer or equal number of boxes.
    if (boxes_hull.size() > boxes_none.size()) {
        std::cerr << "GraphHull 2D quadratic: expected fewer boxes than None method, got "
                  << boxes_hull.size() << " vs " << boxes_none.size() << '\n';
        return 1;
    }

    // Check that at least one box contains the root (0.5, 0.3).
    bool found_root = false;
    for (const SubdivisionBoxResult& box : boxes_hull) {
        if (box.lower[0] <= 0.5 && box.upper[0] >= 0.5 &&
            box.lower[1] <= 0.3 && box.upper[1] >= 0.3) {
            found_root = true;
            break;
        }
    }

    if (!found_root) {
        std::cerr << "GraphHull 2D quadratic: root (0.5, 0.3) not found in any box" << '\n';
        return 1;
    }

    std::cout << "GraphHull 2D quadratic: None method produced " << boxes_none.size()
              << " boxes, GraphHull produced " << boxes_hull.size() << " boxes" << std::endl;

    return 0;
}

int test_subdivision_solver_graph_hull_2d_quadratic_debug()
{
    // Same as quadratic test but with debug output
    std::vector<unsigned int> degrees_f1{2u, 1u};
    std::vector<double> power_coeffs_f1{-0.25, 0.0, 0.0, 0.0, 1.0, 0.0};
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs_f1);

    std::vector<unsigned int> degrees_f2{0u, 1u};
    std::vector<double> power_coeffs_f2{-0.3, 1.0};
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs_f2);

    // Open CSV file for dumping data.
    std::ofstream csv("graph_hull_2d_quadratic_debug.csv");
    csv << "equation,step,type,dim,index,coord0,coord1,coord2\n";

    // For each equation, dump control points and compute convex hull
    std::vector<ConvexPolyhedron> hyperplane_intersections;
    for (std::size_t eq_idx = 0; eq_idx < 2; ++eq_idx) {
        const Polynomial& poly = (eq_idx == 0) ? p1 : p2;

        // Step 0: Control points
        std::vector<double> control_points;
        poly.graphControlPoints(control_points);
        const std::size_t num_coeffs = poly.coefficientCount();

        std::vector<std::vector<double>> points;
        for (std::size_t i = 0; i < num_coeffs; ++i) {
            csv << (eq_idx + 1) << ",0,control_point,3," << i << ","
                << control_points[i * 3 + 0] << ","
                << control_points[i * 3 + 1] << ","
                << control_points[i * 3 + 2] << "\n";

            std::vector<double> pt(3);
            pt[0] = control_points[i * 3 + 0];
            pt[1] = control_points[i * 3 + 1];
            pt[2] = control_points[i * 3 + 2];
            points.push_back(pt);
        }

        // Step 1: Convex hull
        ConvexPolyhedron hull = convex_hull(points);
        for (std::size_t i = 0; i < hull.vertices.size(); ++i) {
            csv << (eq_idx + 1) << ",1,hull_vertex,3," << i << ","
                << hull.vertices[i][0] << ","
                << hull.vertices[i][1] << ","
                << hull.vertices[i][2] << "\n";
        }

        // Step 2: Intersect with hyperplane z=0
        ConvexPolyhedron hyperplane_intersection;
        bool has_intersection = intersect_convex_polyhedron_with_last_coordinate_zero(
            hull, hyperplane_intersection);

        if (has_intersection) {
            for (std::size_t i = 0; i < hyperplane_intersection.vertices.size(); ++i) {
                csv << (eq_idx + 1) << ",2,hyperplane_intersection,3," << i << ","
                    << hyperplane_intersection.vertices[i][0] << ","
                    << hyperplane_intersection.vertices[i][1] << ","
                    << hyperplane_intersection.vertices[i][2] << "\n";
            }

            // Step 3: Project to 2D
            std::vector<std::vector<double>> projected_points;
            projected_points.reserve(hyperplane_intersection.vertices.size());
            for (const auto& v : hyperplane_intersection.vertices) {
                std::vector<double> proj_v(2);
                proj_v[0] = v[0];
                proj_v[1] = v[1];
                projected_points.push_back(proj_v);
            }

            // Recompute convex hull in 2D to ensure proper vertex ordering
            ConvexPolyhedron projected = convex_hull(projected_points);

            for (std::size_t i = 0; i < projected.vertices.size(); ++i) {
                csv << (eq_idx + 1) << ",3,projected,2," << i << ","
                    << projected.vertices[i][0] << "," << projected.vertices[i][1] << ",0\n";
            }

            hyperplane_intersections.push_back(projected);
        }
    }

    // Step 4: Intersect all projected polyhedra
    if (!hyperplane_intersections.empty()) {
        ConvexPolyhedron final_intersection;
        bool success = intersect_convex_polyhedra(hyperplane_intersections, final_intersection);

        if (success && !final_intersection.vertices.empty()) {
            for (std::size_t i = 0; i < final_intersection.vertices.size(); ++i) {
                csv << "0,4,final_intersection,2," << i << ","
                    << final_intersection.vertices[i][0] << ","
                    << final_intersection.vertices[i][1] << ",0\n";
            }

            // Step 5: Bounding box
            ConvexPolyhedronBox bbox = bounding_box(final_intersection);
            csv << "0,5,bbox_min,2,0,"
                << bbox.min_coords[0] << ","
                << bbox.min_coords[1] << ",0\n";
            csv << "0,5,bbox_max,2,0,"
                << bbox.max_coords[0] << ","
                << bbox.max_coords[1] << ",0\n";
        }
    }

    csv.close();
    std::cout << "Quadratic debug data written to graph_hull_2d_quadratic_debug.csv" << std::endl;

    return 0;
}

int main(int argc, char* argv[])
{
    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--debug") == 0) {
            getGeometryConfig().debug = true;
            std::cout << "Debug mode enabled" << std::endl;
        }
    }

    int status = 0;
    status |= test_subdivision_solver_uniform_grid();
    status |= test_subdivision_solver_graph_hull_1d();
    status |= test_subdivision_solver_graph_hull_2d();
    status |= test_subdivision_solver_graph_hull_2d_quadratic();
    status |= test_subdivision_solver_graph_hull_1d_debug();
    status |= test_subdivision_solver_graph_hull_2d_quadratic_debug();

    if (status == 0) {
        std::cout << "Subdivision solver tests passed." << std::endl;
    }

    return status;
}

