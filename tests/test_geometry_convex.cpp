#include "geometry.h"

#include <cmath>
#include <algorithm>

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

int test_convex_hull_box_basic()
{
    // Three points in R^3.
    std::vector<std::vector<double>> points;
    points.push_back({0.0, 1.0, -1.0});
    points.push_back({2.0, -1.0, 3.0});
    points.push_back({-1.0, 0.5, 0.0});

    ConvexPolyhedron poly = convex_hull(points);
    ConvexPolyhedronBox box = bounding_box(poly);

    if (box.dimension() != 3u) {
        std::cerr << "convex_hull/bounding_box: unexpected dimension "
                  << box.dimension() << '\n';
        return 1;
    }

    if (!approx_equal(box.min_coords[0], -1.0) ||
        !approx_equal(box.max_coords[0],  2.0) ||
        !approx_equal(box.min_coords[1], -1.0) ||
        !approx_equal(box.max_coords[1],  1.0) ||
        !approx_equal(box.min_coords[2], -1.0) ||
        !approx_equal(box.max_coords[2],  3.0)) {
        std::cerr << "convex_hull/bounding_box: incorrect bounds" << '\n';
        return 1;
    }

    return 0;
}

int test_convex_hull_2d_graham()
{
    // Points forming a square with interior and boundary points.
    std::vector<std::vector<double>> points;
    points.push_back({0.0, 0.0});
    points.push_back({2.0, 0.0});
    points.push_back({2.0, 2.0});
    points.push_back({0.0, 2.0});
    points.push_back({1.0, 1.0}); // interior
    points.push_back({1.0, 0.0}); // boundary
    points.push_back({2.0, 1.0}); // boundary

    ConvexPolyhedron poly = convex_hull(points);

    if (poly.dimension() != 2u) {
        std::cerr << "convex_hull (2D): unexpected dimension "
                  << poly.dimension() << '\n';
        return 1;
    }

    if (poly.vertices.size() != 4u) {
        std::cerr << "convex_hull (2D): expected 4 hull vertices, got "
                  << poly.vertices.size() << '\n';
        return 1;
    }

    // Sort hull vertices lexicographically to compare with expected corners.
    std::vector<std::vector<double>> verts = poly.vertices;
    std::sort(verts.begin(), verts.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
                  if (a[0] < b[0]) return true;
                  if (a[0] > b[0]) return false;
                  return a[1] < b[1];
              });

    const std::vector<std::vector<double>> expected = {
        {0.0, 0.0},
        {0.0, 2.0},
        {2.0, 0.0},
        {2.0, 2.0}
    };

    if (verts.size() != expected.size()) {
        std::cerr << "convex_hull (2D): incorrect hull size after sort" << '\n';
        return 1;
    }

    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (!approx_equal(verts[i][0], expected[i][0]) ||
            !approx_equal(verts[i][1], expected[i][1])) {
            std::cerr << "convex_hull (2D): incorrect hull vertex at index "
                      << i << '\n';
            return 1;
        }
    }

    return 0;
}



int test_convex_hull_3d_box()
{
    // Points forming a unit cube in R^3 with additional interior / face points.
    std::vector<std::vector<double>> points;
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                points.push_back({static_cast<double>(ix),
                                  static_cast<double>(iy),
                                  static_cast<double>(iz)});
            }
        }
    }

    // Interior point and points on faces that must not become hull vertices.
    points.push_back({0.5, 0.5, 0.5});
    points.push_back({0.0, 0.5, 0.5});
    points.push_back({1.0, 0.5, 0.5});

    ConvexPolyhedron poly = convex_hull(points);

    if (poly.dimension() != 3u) {
        std::cerr << "convex_hull (3D): unexpected dimension "
                  << poly.dimension() << '\n';
        return 1;
    }

    if (poly.vertices.size() != 8u) {
        std::cerr << "convex_hull (3D): expected 8 hull vertices, got "
                  << poly.vertices.size() << '\n';
        return 1;
    }

    std::vector<std::vector<double>> verts = poly.vertices;
    std::sort(verts.begin(), verts.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
                  if (a[0] < b[0]) return true;
                  if (a[0] > b[0]) return false;
                  if (a[1] < b[1]) return true;
                  if (a[1] > b[1]) return false;
                  return a[2] < b[2];
              });

    std::vector<std::vector<double>> expected;
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                expected.push_back({static_cast<double>(ix),
                                    static_cast<double>(iy),
                                    static_cast<double>(iz)});
            }
        }
    }

    if (verts.size() != expected.size()) {
        std::cerr << "convex_hull (3D): incorrect hull size after sort" << '\n';
        return 1;
    }

    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (!approx_equal(verts[i][0], expected[i][0]) ||
            !approx_equal(verts[i][1], expected[i][1]) ||
            !approx_equal(verts[i][2], expected[i][2])) {
            std::cerr << "convex_hull (3D): incorrect hull vertex at index "
                      << i << '\n';
            return 1;
        }
    }

    return 0;
}


int test_intersect_convex_polyhedra()
{
    // Define two axis-aligned rectangles in R^2 via their corners and build
    // convex polyhedra from them.
    std::vector<std::vector<double>> a_pts;
    a_pts.push_back({0.0, 0.0});
    a_pts.push_back({2.0, 0.0});
    a_pts.push_back({0.0, 2.0});
    a_pts.push_back({2.0, 2.0});
    ConvexPolyhedron a_poly = convex_hull(a_pts);

    std::vector<std::vector<double>> b_pts;
    b_pts.push_back({1.0, -1.0});
    b_pts.push_back({3.0, -1.0});
    b_pts.push_back({1.0,  1.0});
    b_pts.push_back({3.0,  1.0});
    ConvexPolyhedron b_poly = convex_hull(b_pts);

    ConvexPolyhedron inter_poly;
    bool ok = intersect_convex_polyhedra({a_poly, b_poly}, inter_poly);
    if (!ok) {
        std::cerr << "intersect_convex_polyhedra: empty intersection" << '\n';
        return 1;
    }

    ConvexPolyhedronBox inter_box = bounding_box(inter_poly);

    if (inter_box.dimension() != 2u) {
        std::cerr << "intersect_convex_polyhedra: unexpected dimension" << '\n';
        return 1;
    }

    if (!approx_equal(inter_box.min_coords[0], 1.0) ||
        !approx_equal(inter_box.max_coords[0], 2.0) ||
        !approx_equal(inter_box.min_coords[1], 0.0) ||
        !approx_equal(inter_box.max_coords[1], 1.0)) {
        std::cerr << "intersect_convex_polyhedra: incorrect bounds" << '\n';
        return 1;
    }

    // A disjoint case.
    std::vector<std::vector<double>> c_pts;
    c_pts.push_back({3.0, 0.0});
    c_pts.push_back({4.0, 0.0});
    c_pts.push_back({3.0, 1.0});
    c_pts.push_back({4.0, 1.0});
    ConvexPolyhedron c_poly = convex_hull(c_pts);

    ok = intersect_convex_polyhedra({a_poly, c_poly}, inter_poly);
    if (ok) {
        std::cerr << "intersect_convex_polyhedra: expected empty intersection" << '\n';
        return 1;
    }

    return 0;
}

int test_intersect_polyhedron_with_last_coordinate_zero()
{
    // Start from an axis-aligned box in R^3, represented via its corners.
    ConvexPolyhedronBox box;
    box.min_coords = {-1.0, -2.0, -0.5};
    box.max_coords = { 1.0,  3.0,  2.0};

    std::vector<std::vector<double>> pts;
    for (std::size_t mask = 0; mask < 8u; ++mask) {
        std::vector<double> p(3u, 0.0);
        p[0] = (mask & 1u) ? box.max_coords[0] : box.min_coords[0];
        p[1] = (mask & 2u) ? box.max_coords[1] : box.min_coords[1];
        p[2] = (mask & 4u) ? box.max_coords[2] : box.min_coords[2];
        pts.push_back(p);
    }
    ConvexPolyhedron poly = convex_hull(pts);

    ConvexPolyhedron inter_poly;
    bool ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly, inter_poly);
    if (!ok) {
        std::cerr << "intersect_convex_polyhedron_with_last_coordinate_zero: unexpected empty" << '\n';
        return 1;
    }

    ConvexPolyhedronBox inter_box = bounding_box(inter_poly);

    if (inter_box.dimension() != 3u) {
        std::cerr << "intersect_convex_polyhedron_with_last_coordinate_zero: unexpected dim" << '\n';
        return 1;
    }

    if (!approx_equal(inter_box.min_coords[2], 0.0) ||
        !approx_equal(inter_box.max_coords[2], 0.0)) {
        std::cerr << "intersect_convex_polyhedron_with_last_coordinate_zero: last coord not zero" << '\n';
        return 1;
    }

    // Case where the last coordinate interval does not contain 0.
    ConvexPolyhedronBox box2;
    box2.min_coords = {-1.0, -2.0, 0.1};
    box2.max_coords = { 1.0,  3.0, 2.0};

    pts.clear();
    for (std::size_t mask = 0; mask < 8u; ++mask) {
        std::vector<double> p(3u, 0.0);
        p[0] = (mask & 1u) ? box2.max_coords[0] : box2.min_coords[0];
        p[1] = (mask & 2u) ? box2.max_coords[1] : box2.min_coords[1];
        p[2] = (mask & 4u) ? box2.max_coords[2] : box2.min_coords[2];
        pts.push_back(p);
    }
    ConvexPolyhedron poly2 = convex_hull(pts);

    ok = intersect_convex_polyhedron_with_last_coordinate_zero(poly2, inter_poly);
    if (ok) {
        std::cerr << "intersect_convex_polyhedron_with_last_coordinate_zero: expected empty" << '\n';
        return 1;
    }

    return 0;
}

int test_intersect_2d_polygons()
{
    // Test 1: Two overlapping squares.
    std::vector<std::vector<double>> square1_pts;
    square1_pts.push_back({0.0, 0.0});
    square1_pts.push_back({2.0, 0.0});
    square1_pts.push_back({2.0, 2.0});
    square1_pts.push_back({0.0, 2.0});
    ConvexPolyhedron square1 = convex_hull(square1_pts);

    std::vector<std::vector<double>> square2_pts;
    square2_pts.push_back({1.0, 1.0});
    square2_pts.push_back({3.0, 1.0});
    square2_pts.push_back({3.0, 3.0});
    square2_pts.push_back({1.0, 3.0});
    ConvexPolyhedron square2 = convex_hull(square2_pts);

    ConvexPolyhedron inter;
    bool ok = intersect_convex_polyhedra({square1, square2}, inter);
    if (!ok) {
        std::cerr << "intersect_2d_polygons: expected non-empty intersection" << '\n';
        return 1;
    }

    // The intersection should be the square [1,2] x [1,2].
    ConvexPolyhedronBox inter_box = bounding_box(inter);
    if (!approx_equal(inter_box.min_coords[0], 1.0) ||
        !approx_equal(inter_box.max_coords[0], 2.0) ||
        !approx_equal(inter_box.min_coords[1], 1.0) ||
        !approx_equal(inter_box.max_coords[1], 2.0)) {
        std::cerr << "intersect_2d_polygons: incorrect intersection bounds" << '\n';
        return 1;
    }

    // Test 2: Disjoint squares.
    std::vector<std::vector<double>> square3_pts;
    square3_pts.push_back({5.0, 5.0});
    square3_pts.push_back({6.0, 5.0});
    square3_pts.push_back({6.0, 6.0});
    square3_pts.push_back({5.0, 6.0});
    ConvexPolyhedron square3 = convex_hull(square3_pts);

    ok = intersect_convex_polyhedra({square1, square3}, inter);
    if (ok) {
        std::cerr << "intersect_2d_polygons: expected empty intersection for disjoint" << '\n';
        return 1;
    }

    // Test 3: Triangle intersecting square.
    std::vector<std::vector<double>> tri_pts;
    tri_pts.push_back({0.5, 0.5});
    tri_pts.push_back({2.5, 0.5});
    tri_pts.push_back({1.5, 2.5});
    ConvexPolyhedron tri = convex_hull(tri_pts);

    ok = intersect_convex_polyhedra({square1, tri}, inter);
    if (!ok) {
        std::cerr << "intersect_2d_polygons: expected non-empty for triangle-square" << '\n';
        return 1;
    }

    // The intersection should be a quadrilateral inside [0,2] x [0,2].
    inter_box = bounding_box(inter);
    if (inter_box.min_coords[0] < -0.01 || inter_box.max_coords[0] > 2.01 ||
        inter_box.min_coords[1] < -0.01 || inter_box.max_coords[1] > 2.01) {
        std::cerr << "intersect_2d_polygons: triangle-square intersection out of bounds" << '\n';
        return 1;
    }

    return 0;
}

int test_intersect_2d_line_segments()
{
    // Test intersecting two line segments in 2D (embedded in R^2).
    // Segment 1: from (0, 0) to (2, 2).
    // Segment 2: from (0, 2) to (2, 0).
    // They should intersect at (1, 1).

    std::vector<std::vector<double>> seg1_pts;
    seg1_pts.push_back({0.0, 0.0});
    seg1_pts.push_back({2.0, 2.0});
    ConvexPolyhedron seg1 = convex_hull(seg1_pts);

    std::vector<std::vector<double>> seg2_pts;
    seg2_pts.push_back({0.0, 2.0});
    seg2_pts.push_back({2.0, 0.0});
    ConvexPolyhedron seg2 = convex_hull(seg2_pts);

    ConvexPolyhedron inter;
    bool ok = intersect_convex_polyhedra({seg1, seg2}, inter);
    if (!ok) {
        std::cerr << "intersect_2d_line_segments: expected intersection" << '\n';
        return 1;
    }

    // The intersection should be a single point (1, 1).
    if (inter.vertices.size() != 1u) {
        std::cerr << "intersect_2d_line_segments: expected 1 vertex, got "
                  << inter.vertices.size() << '\n';
        return 1;
    }

    if (!approx_equal(inter.vertices[0][0], 1.0) ||
        !approx_equal(inter.vertices[0][1], 1.0)) {
        std::cerr << "intersect_2d_line_segments: incorrect intersection point" << '\n';
        return 1;
    }

    // Test parallel non-overlapping segments.
    std::vector<std::vector<double>> seg3_pts;
    seg3_pts.push_back({0.0, 0.0});
    seg3_pts.push_back({1.0, 0.0});
    ConvexPolyhedron seg3 = convex_hull(seg3_pts);

    std::vector<std::vector<double>> seg4_pts;
    seg4_pts.push_back({0.0, 1.0});
    seg4_pts.push_back({1.0, 1.0});
    ConvexPolyhedron seg4 = convex_hull(seg4_pts);

    ok = intersect_convex_polyhedra({seg3, seg4}, inter);
    if (ok) {
        std::cerr << "intersect_2d_line_segments: expected no intersection for parallel" << '\n';
        return 1;
    }

    // Test overlapping collinear segments.
    std::vector<std::vector<double>> seg5_pts;
    seg5_pts.push_back({0.0, 0.0});
    seg5_pts.push_back({2.0, 0.0});
    ConvexPolyhedron seg5 = convex_hull(seg5_pts);

    std::vector<std::vector<double>> seg6_pts;
    seg6_pts.push_back({1.0, 0.0});
    seg6_pts.push_back({3.0, 0.0});
    ConvexPolyhedron seg6 = convex_hull(seg6_pts);

    ok = intersect_convex_polyhedra({seg5, seg6}, inter);
    if (!ok) {
        std::cerr << "intersect_2d_line_segments: expected intersection for collinear" << '\n';
        return 1;
    }

    // The intersection should be the segment [1, 2] on the x-axis.
    ConvexPolyhedronBox inter_box = bounding_box(inter);
    if (!approx_equal(inter_box.min_coords[0], 1.0) ||
        !approx_equal(inter_box.max_coords[0], 2.0) ||
        !approx_equal(inter_box.min_coords[1], 0.0) ||
        !approx_equal(inter_box.max_coords[1], 0.0)) {
        std::cerr << "intersect_2d_line_segments: incorrect collinear intersection" << '\n';
        return 1;
    }

    return 0;
}

} // namespace

int main()
{
    int status = 0;
    status |= test_convex_hull_box_basic();
    status |= test_convex_hull_2d_graham();
    status |= test_convex_hull_3d_box();
    status |= test_intersect_convex_polyhedra();
    status |= test_intersect_polyhedron_with_last_coordinate_zero();
    status |= test_intersect_2d_polygons();
    status |= test_intersect_2d_line_segments();

    if (status == 0) {
        std::cout << "Geometry convex tests passed." << std::endl;
    }

    return status;
}

