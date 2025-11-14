#include "geometry.h"

#include <algorithm>
#include <limits>
#include <map>
#include <cmath>
#include <iostream>


/**
 * @file geometry.cpp
 * @brief Implementation of geometry tools
 */

namespace polynomial_solver {

// Global geometry configuration
static GeometryConfig g_geometry_config;

GeometryConfig& getGeometryConfig() {
    return g_geometry_config;
}

// Point2D implementation
Point2D::Point2D() : x_(0.0), y_(0.0) {
}

Point2D::Point2D(double x, double y) : x_(x), y_(y) {
}

Point2D::~Point2D() {
}

// Interval implementation
Interval::Interval() : min_(0.0), max_(0.0) {
}

Interval::Interval(double min, double max) : min_(min), max_(max) {
}

Interval::~Interval() {
}

namespace {

/**
 * @brief Compute squared distance between two points.
 */
double squared_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    const std::size_t dim = a.size();
    for (std::size_t i = 0; i < dim; ++i) {
        const double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum;
}

/**
 * @brief Check if points are collinear (lie on a line).
 */
bool are_collinear(const std::vector<std::vector<double>>& points, double eps = GEOMETRY_EPSILON) {
    if (points.size() <= 2) {
        return true;
    }

    const std::size_t dim = points[0].size();
    const std::size_t n = points.size();

    // Find two points that are not identical
    std::size_t i0 = 0, i1 = 1;
    while (i1 < n && squared_distance(points[i0], points[i1]) < eps * eps) {
        ++i1;
    }
    if (i1 >= n) {
        return true; // All points are identical
    }

    // Direction vector from points[i0] to points[i1]
    std::vector<double> dir(dim);
    double dir_len_sq = 0.0;
    for (std::size_t j = 0; j < dim; ++j) {
        dir[j] = points[i1][j] - points[i0][j];
        dir_len_sq += dir[j] * dir[j];
    }

    // Check if all other points lie on the line through points[i0] and points[i1]
    for (std::size_t i = 0; i < n; ++i) {
        if (i == i0 || i == i1) continue;

        // Vector from points[i0] to points[i]
        std::vector<double> v(dim);
        for (std::size_t j = 0; j < dim; ++j) {
            v[j] = points[i][j] - points[i0][j];
        }

        // Compute cross product magnitude (for 2D/3D) or perpendicular distance
        double cross_sq = 0.0;
        if (dim == 2) {
            // 2D cross product: |v × dir|
            double cross = v[0] * dir[1] - v[1] * dir[0];
            cross_sq = cross * cross;
        } else if (dim == 3) {
            // 3D cross product: |v × dir|²
            double cx = v[1] * dir[2] - v[2] * dir[1];
            double cy = v[2] * dir[0] - v[0] * dir[2];
            double cz = v[0] * dir[1] - v[1] * dir[0];
            cross_sq = cx * cx + cy * cy + cz * cz;
        } else {
            // General case: compute perpendicular distance
            // Project v onto dir: proj = (v·dir / |dir|²) * dir
            double dot = 0.0;
            for (std::size_t j = 0; j < dim; ++j) {
                dot += v[j] * dir[j];
            }
            double proj_len = dot / dir_len_sq;

            // Perpendicular component: v - proj
            double perp_sq = 0.0;
            for (std::size_t j = 0; j < dim; ++j) {
                double perp_j = v[j] - proj_len * dir[j];
                perp_sq += perp_j * perp_j;
            }
            cross_sq = perp_sq;
        }

        // If perpendicular distance is too large, not collinear
        if (cross_sq > eps * eps * dir_len_sq) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Check if points are coplanar (lie on a plane) in 3D.
 */
bool are_coplanar_3d(const std::vector<std::vector<double>>& points, double eps = GEOMETRY_EPSILON) {
    if (points.size() <= 3) {
        return true;
    }

    const std::size_t n = points.size();

    // Find three non-collinear points
    std::size_t i0 = 0, i1 = 1, i2 = 2;

    // Find i1 != i0
    while (i1 < n && squared_distance(points[i0], points[i1]) < eps * eps) {
        ++i1;
    }
    if (i1 >= n) return true;

    // Find i2 not collinear with i0, i1
    while (i2 < n) {
        if (i2 == i0 || i2 == i1) {
            ++i2;
            continue;
        }

        // Check if i2 is collinear with i0, i1
        std::vector<std::vector<double>> three_pts = {points[i0], points[i1], points[i2]};
        if (!are_collinear(three_pts, eps)) {
            break;
        }
        ++i2;
    }
    if (i2 >= n) return true; // All points are collinear

    // Compute plane normal: (p1 - p0) × (p2 - p0)
    std::vector<double> v1(3), v2(3);
    for (std::size_t j = 0; j < 3; ++j) {
        v1[j] = points[i1][j] - points[i0][j];
        v2[j] = points[i2][j] - points[i0][j];
    }

    std::vector<double> normal(3);
    normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
    normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
    normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

    double normal_len_sq = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
    if (normal_len_sq < eps * eps) {
        return true; // Degenerate normal
    }

    // Check if all other points lie on the plane
    for (std::size_t i = 0; i < n; ++i) {
        if (i == i0 || i == i1 || i == i2) continue;

        // Vector from points[i0] to points[i]
        std::vector<double> v(3);
        for (std::size_t j = 0; j < 3; ++j) {
            v[j] = points[i][j] - points[i0][j];
        }

        // Distance to plane: |v · normal| / |normal|
        double dot = v[0] * normal[0] + v[1] * normal[1] + v[2] * normal[2];
        double dist_sq = (dot * dot) / normal_len_sq;

        if (dist_sq > eps * eps) {
            return false;
        }
    }

    return true;
}

// Helper: compute an axis-aligned bounding box of a set of points.
ConvexPolyhedronBox
convex_hull_box_from_points(const std::vector<std::vector<double>>& points)
{
    ConvexPolyhedronBox box;
    if (points.empty()) {
        return box;
    }

    const std::size_t dim = points.front().size();
    box.min_coords.assign(dim, std::numeric_limits<double>::infinity());
    box.max_coords.assign(dim, -std::numeric_limits<double>::infinity());

    for (const std::vector<double>& p : points) {
        if (p.size() != dim) {
            // Ignore points of mismatched dimension for robustness.
            continue;
        }
        for (std::size_t j = 0; j < dim; ++j) {
            box.min_coords[j] = std::min(box.min_coords[j], p[j]);
            box.max_coords[j] = std::max(box.max_coords[j], p[j]);
        }
    }

    return box;
}

// Helper: intersect a collection of axis-aligned boxes.
bool intersect_boxes(
    const std::vector<ConvexPolyhedronBox>& boxes,
    ConvexPolyhedronBox& intersection)
{
    if (boxes.empty()) {
        return false;
    }

    const std::size_t dim = boxes.front().dimension();
    intersection = boxes.front();

    for (std::size_t k = 1; k < boxes.size(); ++k) {
        const ConvexPolyhedronBox& box = boxes[k];
        if (box.dimension() != dim) {
            // Dimension mismatch: treat as empty intersection.
            return false;
        }
        for (std::size_t j = 0; j < dim; ++j) {
            intersection.min_coords[j] =
                std::max(intersection.min_coords[j], box.min_coords[j]);
            intersection.max_coords[j] =
                std::min(intersection.max_coords[j], box.max_coords[j]);
        }
    }

    for (std::size_t j = 0; j < dim; ++j) {
        if (intersection.min_coords[j] > intersection.max_coords[j]) {
            return false;
        }
    }

    return true;
}

// Helper: intersect an axis-aligned box with the hyperplane x_{n+1} = 0.
bool intersect_box_with_last_coordinate_zero_impl(
    const ConvexPolyhedronBox& box,
    ConvexPolyhedronBox& intersection)
{
    const std::size_t dim = box.dimension();
    if (dim == 0u) {
        return false;
    }

    intersection = box;

    const std::size_t last = dim - 1u;
    const double min_v = intersection.min_coords[last];
    const double max_v = intersection.max_coords[last];

    if (min_v > 0.0 || max_v < 0.0) {
        // The interval in the last coordinate does not contain 0.
        return false;
    }

    intersection.min_coords[last] = 0.0;
    intersection.max_coords[last] = 0.0;

    return true;
}

// Helper: enumerate the corners of an axis-aligned box as vertices.
void box_corners(const ConvexPolyhedronBox& box,
                 std::vector<std::vector<double>>& corners)
{
    const std::size_t dim = box.dimension();
    corners.clear();
    if (dim == 0u) {
        return;
    }

    const std::size_t num_corners = static_cast<std::size_t>(1u) << dim;
    corners.reserve(num_corners);

    for (std::size_t mask = 0; mask < num_corners; ++mask) {
        std::vector<double> vertex(dim, 0.0);
        for (std::size_t j = 0; j < dim; ++j) {
            const bool use_max = ((mask >> j) & 1u) != 0u;
            vertex[j] = use_max ? box.max_coords[j] : box.min_coords[j];
        }
        corners.push_back(vertex);
    }
}

struct Face3D {
    int a;
    int b;
    int c;
    std::vector<double> normal;
    bool valid;
};

static std::vector<double> subtract3(const std::vector<double>& u,
                                     const std::vector<double>& v)
{
    std::vector<double> r(3u, 0.0);
    r[0] = u[0] - v[0];
    r[1] = u[1] - v[1];
    r[2] = u[2] - v[2];
    return r;
}

static std::vector<double> cross3(const std::vector<double>& u,
                                  const std::vector<double>& v)
{
    std::vector<double> r(3u, 0.0);
    r[0] = u[1] * v[2] - u[2] * v[1];
    r[1] = u[2] * v[0] - u[0] * v[2];
    r[2] = u[0] * v[1] - u[1] * v[0];
    return r;
}

static double dot3(const std::vector<double>& u,
                   const std::vector<double>& v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

static std::vector<double> triangle_normal3D(
    const std::vector<std::vector<double>>& pts,
    int ia, int ib, int ic)
{
    const std::vector<double>& a = pts[ia];
    const std::vector<double>& b = pts[ib];
    const std::vector<double>& c = pts[ic];
    const std::vector<double> ab = subtract3(b, a);
    const std::vector<double> ac = subtract3(c, a);
    return cross3(ab, ac);
}

static double signed_distance_to_face(
    const Face3D& face,
    const std::vector<std::vector<double>>& pts,
    const std::vector<double>& p)
{
    const std::vector<double>& a = pts[face.a];
    const std::vector<double> ap = subtract3(p, a);
    return dot3(face.normal, ap);
}



static std::vector<std::vector<double>>
convex_hull_3d(const std::vector<std::vector<double>>& pts_in)
{
    const double eps = GEOMETRY_EPSILON;
    std::vector<std::vector<double>> hull;

    if (pts_in.empty()) {
        return hull;
    }

    // Work on a local copy so we can reorder if needed.
    std::vector<std::vector<double>> pts = pts_in;
    const std::size_t n = pts.size();

    if (n <= 3u) {
        // Remove exact duplicates.
        std::sort(pts.begin(), pts.end(),
                  [](const std::vector<double>& a,
                     const std::vector<double>& b) {
                      if (a[0] < b[0]) return true;
                      if (a[0] > b[0]) return false;
                      if (a[1] < b[1]) return true;
                      if (a[1] > b[1]) return false;
                      return a[2] < b[2];
                  });
        pts.erase(std::unique(pts.begin(), pts.end(),
                              [](const std::vector<double>& a,
                                 const std::vector<double>& b) {
                                  return a[0] == b[0] && a[1] == b[1] &&
                                         a[2] == b[2];
                              }),
                  pts.end());
        return pts;
    }

    // Step 1: pick two farthest points to define a line.
    int i0 = 0;
    int i1 = -1;
    double max_dist2 = 0.0;
    for (std::size_t i = 1; i < n; ++i) {
        const std::vector<double> diff = subtract3(pts[i], pts[i0]);
        const double d2 = dot3(diff, diff);
        if (d2 > max_dist2) {
            max_dist2 = d2;
            i1 = static_cast<int>(i);
        }
    }
    if (i1 < 0 || max_dist2 <= eps * eps) {
        // All points coincide.
        hull.push_back(pts[i0]);
        return hull;
    }

    // Step 2: find a third point maximizing area (to avoid collinearity).
    int i2 = -1;
    double max_area2 = 0.0;
    const std::vector<double> dir = subtract3(pts[i1], pts[i0]);
    for (std::size_t i = 0; i < n; ++i) {
        if (static_cast<int>(i) == i0 || static_cast<int>(i) == i1) {
            continue;
        }
        const std::vector<double> vi = subtract3(pts[i], pts[i0]);
        const std::vector<double> cr = cross3(dir, vi);
        const double area2 = dot3(cr, cr);
        if (area2 > max_area2) {
            max_area2 = area2;
            i2 = static_cast<int>(i);
        }
    }
    if (i2 < 0 || max_area2 <= eps * eps) {
        // Points are (almost) collinear: take extreme points along dir.
        int min_idx = i0;
        int max_idx = i0;
        double min_proj = 0.0;
        double max_proj = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::vector<double> vi = subtract3(pts[i], pts[i0]);
            const double t = dot3(vi, dir);
            if (t < min_proj) {
                min_proj = t;
                min_idx = static_cast<int>(i);
            }
            if (t > max_proj) {
                max_proj = t;
                max_idx = static_cast<int>(i);
            }
        }
        hull.push_back(pts[min_idx]);
        if (max_idx != min_idx) {
            hull.push_back(pts[max_idx]);
        }
        return hull;
    }

    // Step 3: find a fourth point giving non-zero volume (to avoid coplanarity).
    int i3 = -1;
    double max_volume = 0.0;
    const std::vector<double> base_u = subtract3(pts[i1], pts[i0]);
    const std::vector<double> base_v = subtract3(pts[i2], pts[i0]);
    const std::vector<double> base_normal = cross3(base_u, base_v);
    for (std::size_t i = 0; i < n; ++i) {
        if (static_cast<int>(i) == i0 || static_cast<int>(i) == i1 ||
            static_cast<int>(i) == i2) {
            continue;
        }
        const std::vector<double> w = subtract3(pts[i], pts[i0]);
        const double vol = dot3(base_normal, w);
        const double abs_vol = (vol >= 0.0) ? vol : -vol;
        if (abs_vol > max_volume) {
            max_volume = abs_vol;
            i3 = static_cast<int>(i);
        }
    }

    if (i3 < 0 || max_volume <= eps) {
        // All points are (almost) coplanar: project to 2D and apply a 2D hull.
        struct ProjectedPoint {
            double x;
            double y;
            int index;
        };

        std::vector<ProjectedPoint> proj;
        proj.reserve(n);

        const std::vector<double> origin = pts[i0];
        const std::vector<double> e1 = base_u;
        const std::vector<double> nvec = base_normal;
        const std::vector<double> e2 = cross3(nvec, e1);

        for (std::size_t i = 0; i < n; ++i) {
            const std::vector<double> w = subtract3(pts[i], origin);
            const double x = dot3(w, e1);
            const double y = dot3(w, e2);
            ProjectedPoint p;
            p.x = x;
            p.y = y;
            p.index = static_cast<int>(i);
            proj.push_back(p);
        }

        std::sort(proj.begin(), proj.end(),
                  [](const ProjectedPoint& a, const ProjectedPoint& b) {
                      if (a.x < b.x) return true;
                      if (a.x > b.x) return false;
                      return a.y < b.y;
                  });
        proj.erase(std::unique(proj.begin(), proj.end(),
                               [](const ProjectedPoint& a,
                                  const ProjectedPoint& b) {
                                   return a.x == b.x && a.y == b.y;
                               }),
                   proj.end());

        if (proj.empty()) {
            return hull;
        }
        if (proj.size() == 1u) {
            hull.push_back(pts[proj.front().index]);
            return hull;
        }

        auto cross2d = [](const ProjectedPoint& o,
                          const ProjectedPoint& a,
                          const ProjectedPoint& b) {
            const double x1 = a.x - o.x;
            const double y1 = a.y - o.y;
            const double x2 = b.x - o.x;
            const double y2 = b.y - o.y;
            return x1 * y2 - y1 * x2;
        };

        std::vector<ProjectedPoint> lower;
        std::vector<ProjectedPoint> upper;

        for (const ProjectedPoint& p : proj) {
            while (lower.size() >= 2u &&
                   cross2d(lower[lower.size() - 2u],
                           lower[lower.size() - 1u], p) <= 0.0) {
                lower.pop_back();
            }
            lower.push_back(p);
        }

        for (std::size_t i = proj.size(); i-- > 0u;) {
            const ProjectedPoint& p = proj[i];
            while (upper.size() >= 2u &&
                   cross2d(upper[upper.size() - 2u],
                           upper[upper.size() - 1u], p) <= 0.0) {
                upper.pop_back();
            }
            upper.push_back(p);
        }

        std::vector<int> hull_indices;
        for (std::size_t i = 0; i < lower.size(); ++i) {
            hull_indices.push_back(lower[i].index);
        }
        for (std::size_t i = 1; i + 1 < upper.size(); ++i) {
            hull_indices.push_back(upper[i].index);
        }

        std::sort(hull_indices.begin(), hull_indices.end());
        hull_indices.erase(
            std::unique(hull_indices.begin(), hull_indices.end()),
            hull_indices.end());

        for (int idx : hull_indices) {
            hull.push_back(pts[static_cast<std::size_t>(idx)]);
        }
        return hull;
    }

    // Proper 3D case: incremental convex hull starting from a tetrahedron.
    std::vector<double> interior(3u, 0.0);
    for (int k : {i0, i1, i2, i3}) {
        interior[0] += pts[static_cast<std::size_t>(k)][0];
        interior[1] += pts[static_cast<std::size_t>(k)][1];
        interior[2] += pts[static_cast<std::size_t>(k)][2];
    }
    interior[0] *= 0.25;
    interior[1] *= 0.25;
    interior[2] *= 0.25;

    std::vector<Face3D> faces;
    faces.reserve(16u);

    auto add_face = [&](int a, int b, int c) {
        Face3D f;
        f.a = a;
        f.b = b;
        f.c = c;
        f.valid = true;
        f.normal = triangle_normal3D(pts, a, b, c);
        const double dotc =
            dot3(f.normal, subtract3(interior, pts[static_cast<std::size_t>(a)]));
        if (dotc > 0.0) {
            std::swap(f.b, f.c);
            f.normal = triangle_normal3D(pts, a, f.b, f.c);
        }
        faces.push_back(f);
    };

    add_face(i0, i1, i2);
    add_face(i0, i3, i1);
    add_face(i0, i2, i3);
    add_face(i1, i3, i2);

    // Process all remaining points.
    for (std::size_t i = 0; i < n; ++i) {
        const int idx = static_cast<int>(i);
        if (idx == i0 || idx == i1 || idx == i2 || idx == i3) {
            continue;
        }
        const std::vector<double>& p = pts[i];

        // Determine which faces are visible from p.
        std::vector<bool> visible(faces.size(), false);
        bool any_visible = false;
        for (std::size_t fi = 0; fi < faces.size(); ++fi) {
            if (!faces[fi].valid) {
                continue;
            }
            const double dist = signed_distance_to_face(faces[fi], pts, p);
            if (dist > eps) {
                visible[fi] = true;
                any_visible = true;
            }
        }
        if (!any_visible) {
            continue;
        }

        // Collect horizon edges (edges of visible faces that are not shared by
        // two visible faces).
        std::map<std::pair<int, int>, int> edge_count;
        auto add_edge = [&](int u, int v) {
            if (u > v) {
                std::swap(u, v);
            }
            ++edge_count[std::make_pair(u, v)];
        };

        for (std::size_t fi = 0; fi < faces.size(); ++fi) {
            if (!visible[fi]) {
                continue;
            }
            const Face3D& f = faces[fi];
            add_edge(f.a, f.b);
            add_edge(f.b, f.c);
            add_edge(f.c, f.a);
        }

        // Remove visible faces from the hull.
        for (std::size_t fi = 0; fi < faces.size(); ++fi) {
            if (visible[fi]) {
                faces[fi].valid = false;
            }
        }

        // Create new faces by connecting the new point to each horizon edge.
        for (const auto& kv : edge_count) {
            if (kv.second != 1) {
                continue;
            }
            int u = kv.first.first;
            int v = kv.first.second;
            Face3D f;
            f.a = u;
            f.b = v;
            f.c = idx;
            f.valid = true;
            f.normal = triangle_normal3D(pts, f.a, f.b, f.c);
            const double dotc = dot3(
                f.normal,
                subtract3(interior, pts[static_cast<std::size_t>(f.a)]));
            if (dotc > 0.0) {
                std::swap(f.b, f.c);
                f.normal = triangle_normal3D(pts, f.a, f.b, f.c);
            }
            faces.push_back(f);
        }
    }

    // Collect unique vertex indices from all remaining faces.
    std::vector<bool> is_vertex(n, false);
    for (const Face3D& f : faces) {
        if (!f.valid) {
            continue;
        }
        is_vertex[static_cast<std::size_t>(f.a)] = true;
        is_vertex[static_cast<std::size_t>(f.b)] = true;
        is_vertex[static_cast<std::size_t>(f.c)] = true;
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (is_vertex[i]) {
            hull.push_back(pts[i]);
        }
    }

    return hull;
}

} // namespace

std::size_t compute_intrinsic_dimension(
    const std::vector<std::vector<double>>& points)
{
    if (points.empty()) {
        return 0;
    }

    const std::size_t ambient_dim = points.front().size();
    const std::size_t n = points.size();

    // Filter out duplicate points
    std::vector<std::vector<double>> unique_pts;
    unique_pts.reserve(n);
    const double eps = GEOMETRY_EPSILON;

    for (const std::vector<double>& p : points) {
        if (p.size() != ambient_dim) continue;

        bool is_duplicate = false;
        for (const std::vector<double>& q : unique_pts) {
            if (squared_distance(p, q) < eps * eps) {
                is_duplicate = true;
                break;
            }
        }
        if (!is_duplicate) {
            unique_pts.push_back(p);
        }
    }

    if (unique_pts.empty()) {
        return 0;
    }
    if (unique_pts.size() == 1) {
        return 0; // Single point
    }

    // Check if collinear
    if (are_collinear(unique_pts, eps)) {
        return 1; // Line segment or line
    }

    // For 3D, check if coplanar
    if (ambient_dim == 3 && are_coplanar_3d(unique_pts, eps)) {
        return 2; // Planar
    }

    // For 2D, if not collinear, must be 2D
    if (ambient_dim == 2) {
        return 2;
    }

    // For 3D, if not coplanar, must be 3D
    if (ambient_dim == 3) {
        return 3;
    }

    // For higher dimensions, assume full-dimensional
    return ambient_dim;
}

ConvexPolyhedron
convex_hull(const std::vector<std::vector<double>>& points)
{
    ConvexPolyhedron poly;
    if (points.empty()) {
        return poly;
    }

    const std::size_t dim = points.front().size();

    // Filter out points with mismatched dimension.
    std::vector<std::vector<double>> pts;
    pts.reserve(points.size());
    for (const std::vector<double>& p : points) {
        if (p.size() == dim) {
            pts.push_back(p);
        }
    }

    if (pts.empty()) {
        poly.intrinsic_dim = 0;
        return poly;
    }

    // Compute intrinsic dimension
    poly.intrinsic_dim = compute_intrinsic_dimension(pts);

    if (dim == 3u) {
        poly.vertices = convex_hull_3d(pts);
        return poly;
    }

    // For dimensions other than 2, we currently just store the input points as
    // the vertex set. This keeps the API generic while the exact convex hull is
    // implemented explicitly for 2D using a Graham-scan/monotone-chain variant.
    if (dim != 2u) {
        poly.vertices = pts;
        return poly;
    }

    // 2D convex hull via monotone chain (a variant of Graham's scan).
    // Sort points lexicographically by (x, y) and remove duplicates.
    std::sort(pts.begin(), pts.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
                  if (a[0] < b[0]) return true;
                  if (a[0] > b[0]) return false;
                  return a[1] < b[1];
              });

    pts.erase(std::unique(pts.begin(), pts.end(),
                          [](const std::vector<double>& a,
                             const std::vector<double>& b) {
                              return a[0] == b[0] && a[1] == b[1];
                          }),
              pts.end());

    if (pts.size() <= 1u) {
        poly.vertices = pts;
        return poly;
    }

    // Special case for 2 points (line segment)
    if (pts.size() == 2u) {
        poly.vertices = pts;
        return poly;
    }

    auto cross = [](const std::vector<double>& o,
                    const std::vector<double>& a,
                    const std::vector<double>& b) {
        const double x1 = a[0] - o[0];
        const double y1 = a[1] - o[1];
        const double x2 = b[0] - o[0];
        const double y2 = b[1] - o[1];
        return x1 * y2 - y1 * x2;
    };

    std::vector<std::vector<double>> lower;
    std::vector<std::vector<double>> upper;

    // Build lower hull.
    for (const std::vector<double>& p : pts) {
        while (lower.size() >= 2u &&
               cross(lower[lower.size() - 2u],
                     lower[lower.size() - 1u], p) <= 0.0) {
            lower.pop_back();
        }
        lower.push_back(p);
    }

    // Build upper hull.
    for (std::size_t i = pts.size(); i-- > 0u; ) {
        const std::vector<double>& p = pts[i];
        while (upper.size() >= 2u &&
               cross(upper[upper.size() - 2u],
                     upper[upper.size() - 1u], p) <= 0.0) {
            upper.pop_back();
        }
        upper.push_back(p);
    }

    // Concatenate lower and upper to get full hull. The last point of each list
    // is the starting point of the other list and is omitted to avoid
    // duplication.
    poly.vertices.clear();
    poly.vertices.reserve(lower.size() + upper.size() - 2u);

    for (std::size_t i = 0; i < lower.size(); ++i) {
        poly.vertices.push_back(lower[i]);
    }
    for (std::size_t i = 1; i + 1 < upper.size(); ++i) {
        poly.vertices.push_back(upper[i]);
    }

    return poly;
}

ConvexPolyhedronBox
bounding_box(const ConvexPolyhedron& polyhedron)
{
    return convex_hull_box_from_points(polyhedron.vertices);
}

// Helper: intersect two line segments in 2D (both embedded in R^2 with last coord = 0).
// Each segment is given by two endpoints.
// Returns true if intersection is non-empty, storing result in intersection.
static bool intersect_line_segments_2d(
    const std::vector<double>& a1, const std::vector<double>& a2,
    const std::vector<double>& b1, const std::vector<double>& b2,
    std::vector<std::vector<double>>& intersection_points)
{
    const double eps = GEOMETRY_EPSILON;
    intersection_points.clear();

    // Extract 2D coordinates (first two components).
    const double ax1 = a1[0], ay1 = a1[1];
    const double ax2 = a2[0], ay2 = a2[1];
    const double bx1 = b1[0], by1 = b1[1];
    const double bx2 = b2[0], by2 = b2[1];

    // Direction vectors.
    const double adx = ax2 - ax1, ady = ay2 - ay1;
    const double bdx = bx2 - bx1, bdy = by2 - by1;

    // Cross product for parametric intersection: a1 + s*(a2-a1) = b1 + t*(b2-b1).
    const double cross = adx * bdy - ady * bdx;
    const double abs_cross = (cross >= 0.0) ? cross : -cross;

    if (abs_cross < eps) {
        // Segments are parallel or collinear.
        // Check if they are collinear by testing if b1 lies on line through a1, a2.
        const double dx = bx1 - ax1, dy = by1 - ay1;
        const double test_cross = adx * dy - ady * dx;
        const double abs_test = (test_cross >= 0.0) ? test_cross : -test_cross;
        if (abs_test > eps) {
            // Not collinear -> no intersection.
            return false;
        }

        // Collinear: project onto the line and find overlap.
        const double len2 = adx * adx + ady * ady;
        if (len2 < eps * eps) {
            // Segment a is degenerate (a point).
            // Check if b1 or b2 coincide with a1.
            const double dist_b1 = (bx1 - ax1) * (bx1 - ax1) + (by1 - ay1) * (by1 - ay1);
            const double dist_b2 = (bx2 - ax1) * (bx2 - ax1) + (by2 - ay1) * (by2 - ay1);
            if (dist_b1 < eps * eps || dist_b2 < eps * eps) {
                intersection_points.push_back(a1);
                return true;
            }
            return false;
        }

        // Project b1, b2 onto the line through a1, a2.
        const double t_b1 = ((bx1 - ax1) * adx + (by1 - ay1) * ady) / len2;
        const double t_b2 = ((bx2 - ax1) * adx + (by2 - ay1) * ady) / len2;

        double t_min = (t_b1 < t_b2) ? t_b1 : t_b2;
        double t_max = (t_b1 < t_b2) ? t_b2 : t_b1;

        // Clamp to [0, 1] (segment a).
        if (t_min > 1.0 + eps || t_max < -eps) {
            return false;
        }
        if (t_min < 0.0) t_min = 0.0;
        if (t_max > 1.0) t_max = 1.0;

        // Intersection is the segment from t_min to t_max.
        std::vector<double> p1(a1.size(), 0.0);
        std::vector<double> p2(a1.size(), 0.0);
        for (std::size_t k = 0; k < a1.size(); ++k) {
            p1[k] = a1[k] + t_min * (a2[k] - a1[k]);
            p2[k] = a1[k] + t_max * (a2[k] - a1[k]);
        }
        intersection_points.push_back(p1);
        const double dx_seg = p2[0] - p1[0], dy_seg = p2[1] - p1[1];
        if (dx_seg * dx_seg + dy_seg * dy_seg > eps * eps) {
            intersection_points.push_back(p2);
        }
        return true;
    }

    // Non-parallel: solve for parameters s, t.
    const double dx = bx1 - ax1, dy = by1 - ay1;
    const double s = (dx * bdy - dy * bdx) / cross;
    const double t = (dx * ady - dy * adx) / cross;

    if (s >= -eps && s <= 1.0 + eps && t >= -eps && t <= 1.0 + eps) {
        std::vector<double> p(a1.size(), 0.0);
        for (std::size_t k = 0; k < a1.size(); ++k) {
            p[k] = a1[k] + s * (a2[k] - a1[k]);
        }
        intersection_points.push_back(p);
        return true;
    }

    return false;
}

// Helper: test which side of a line a point is on.
// Returns positive if point is on the left side of the line from p1 to p2.
static double side_of_line_2d(const std::vector<double>& point,
                               const std::vector<double>& p1,
                               const std::vector<double>& p2)
{
    return (p2[0] - p1[0]) * (point[1] - p1[1]) - (p2[1] - p1[1]) * (point[0] - p1[0]);
}

// Helper: compute intersection of line segment s1-s2 with line through p1-p2.
static std::vector<double> line_intersection_2d(
    const std::vector<double>& s1, const std::vector<double>& s2,
    const std::vector<double>& p1, const std::vector<double>& p2,
    double d1, double d2)
{
    const double t = d1 / (d1 - d2);
    std::vector<double> p(2, 0.0);
    p[0] = s1[0] + t * (s2[0] - s1[0]);
    p[1] = s1[1] + t * (s2[1] - s1[1]);
    return p;
}

// Helper: intersect two 2D convex polygons using Sutherland-Hodgman clipping.
// Both polygons are given as ConvexPolyhedron with dimension 2.
// Assumes vertices are in CCW order (as returned by convex_hull).
static bool intersect_convex_polygons_2d(
    const ConvexPolyhedron& poly_a,
    const ConvexPolyhedron& poly_b,
    ConvexPolyhedron& intersection)
{
    const double eps = GEOMETRY_EPSILON;

    if (poly_a.vertices.empty() || poly_b.vertices.empty()) {
        return false;
    }

    // Start with polygon A.
    std::vector<std::vector<double>> subject = poly_a.vertices;

    // Clip against each edge of polygon B.
    const std::vector<std::vector<double>>& clip_verts = poly_b.vertices;
    const std::size_t m = clip_verts.size();

    for (std::size_t i = 0; i < m; ++i) {
        if (subject.empty()) {
            return false;
        }

        const std::vector<double>& edge_p1 = clip_verts[i];
        const std::vector<double>& edge_p2 = clip_verts[(i + 1) % m];

        std::vector<std::vector<double>> clipped;
        clipped.reserve(subject.size() + 1);

        const std::size_t n = subject.size();
        for (std::size_t j = 0; j < n; ++j) {
            const std::vector<double>& current = subject[j];
            const std::vector<double>& next = subject[(j + 1) % n];

            const double d_current = side_of_line_2d(current, edge_p1, edge_p2);
            const double d_next = side_of_line_2d(next, edge_p1, edge_p2);

            const bool inside_current = d_current >= -eps;
            const bool inside_next = d_next >= -eps;

            if (inside_next) {
                if (!inside_current) {
                    // Entering: add intersection point.
                    clipped.push_back(line_intersection_2d(current, next, edge_p1, edge_p2,
                                                           d_current, d_next));
                }
                // Add next vertex (it's inside).
                clipped.push_back(next);
            } else if (inside_current) {
                // Leaving: add intersection point.
                clipped.push_back(line_intersection_2d(current, next, edge_p1, edge_p2,
                                                       d_current, d_next));
            }
        }

        subject = clipped;
    }

    if (subject.empty()) {
        return false;
    }

    intersection.vertices = subject;
    return true;
}

bool intersect_convex_polyhedra(
    const std::vector<ConvexPolyhedron>& polyhedra,
    ConvexPolyhedron& intersection)
{
    if (polyhedra.empty()) {
        return false;
    }

    const std::size_t ambient_dim = polyhedra.front().ambient_dimension();

    // For 2D ambient space, use exact polygon intersection.
    if (ambient_dim == 2u) {
        ConvexPolyhedron result = polyhedra[0];
        for (std::size_t i = 1; i < polyhedra.size(); ++i) {
            if (polyhedra[i].ambient_dimension() != ambient_dim) {
                return false;
            }

            // Special cases for degenerate polyhedra (points, line segments).
            const std::size_t n_result = result.vertices.size();
            const std::size_t n_other = polyhedra[i].vertices.size();

            // Case: both are single points
            if (n_result == 1u && n_other == 1u) {
                const double eps = GEOMETRY_EPSILON;
                const double dist_sq = squared_distance(result.vertices[0], polyhedra[i].vertices[0]);
                if (dist_sq > eps * eps) {
                    // Different points - no intersection
                    return false;
                }
                // Same point - intersection is that point
                result.intrinsic_dim = 0;
                // result.vertices already contains the point
            } else if (n_result == 1u || n_other == 1u) {
                // One is a point, the other is not
                // Check if the point lies in/on the other polyhedron
                const std::vector<double>& pt = (n_result == 1u) ? result.vertices[0] : polyhedra[i].vertices[0];
                const ConvexPolyhedron& poly = (n_result == 1u) ? polyhedra[i] : result;

                // TODO: Implement proper point-in-convex-polygon test.
                // For now, use a simple containment test: check if point is one of the vertices.
                // This is sufficient for many cases but may miss points on edges or in the interior.
                const double eps = GEOMETRY_EPSILON;
                bool found = false;
                for (const std::vector<double>& v : poly.vertices) {
                    if (squared_distance(pt, v) < eps * eps) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    // Point not in polyhedron - no intersection
                    return false;
                }

                // Point is in polyhedron - intersection is the point
                result.vertices = {pt};
                result.intrinsic_dim = 0;
            } else if (n_result == 2u && n_other == 2u) {
                // Both are line segments.
                std::vector<std::vector<double>> inter_pts;
                if (!intersect_line_segments_2d(result.vertices[0], result.vertices[1],
                                                 polyhedra[i].vertices[0], polyhedra[i].vertices[1],
                                                 inter_pts)) {
                    return false;
                }
                result.vertices = inter_pts;
                result.intrinsic_dim = compute_intrinsic_dimension(inter_pts);
            } else if (n_result == 2u || n_other == 2u) {
                // One is a segment, the other is a polygon (or both are segments, already handled above).
                // Clip the segment against the polygon.

                // Make sure we're not in the "both are segments" case (already handled above)
                if (n_result == 2u && n_other == 2u) {
                    // This should have been handled in the previous case
                    return false;
                }

                const ConvexPolyhedron& seg = (n_result == 2u) ? result : polyhedra[i];
                const ConvexPolyhedron& poly = (n_result == 2u) ? polyhedra[i] : result;

                // Debug output
                const bool debug = getGeometryConfig().debug;
                if (debug) {
                    std::cout << "DEBUG: Clipping segment against polygon" << std::endl;
                    std::cout << "  Segment: (" << seg.vertices[0][0] << "," << seg.vertices[0][1]
                              << ") to (" << seg.vertices[1][0] << "," << seg.vertices[1][1] << ")" << std::endl;
                    std::cout << "  Polygon: " << poly.vertices.size() << " vertices" << std::endl;
                    for (std::size_t v = 0; v < poly.vertices.size(); ++v) {
                        std::cout << "    v" << v << ": (" << poly.vertices[v][0] << "," << poly.vertices[v][1] << ")" << std::endl;
                    }
                }

                std::vector<std::vector<double>> inter_pts;
                inter_pts.push_back(seg.vertices[0]);
                inter_pts.push_back(seg.vertices[1]);

                // Clip segment against each edge of the polygon.
                const std::size_t m = poly.vertices.size();
                for (std::size_t j = 0; j < m; ++j) {
                    if (inter_pts.empty()) {
                        if (debug) {
                            std::cout << "  Edge " << j << ": inter_pts is empty, returning false" << std::endl;
                        }
                        return false;
                    }

                    const std::vector<double>& edge_start = poly.vertices[j];
                    const std::vector<double>& edge_end = poly.vertices[(j + 1) % m];

                    const double ex = edge_end[0] - edge_start[0];
                    const double ey = edge_end[1] - edge_start[1];
                    const double nx = -ey;
                    const double ny = ex;

                    if (debug) {
                        std::cout << "  Edge " << j << ": (" << edge_start[0] << "," << edge_start[1]
                                  << ") to (" << edge_end[0] << "," << edge_end[1] << ")" << std::endl;
                        std::cout << "    Normal: (" << nx << "," << ny << ")" << std::endl;
                        std::cout << "    inter_pts before clipping: " << inter_pts.size() << " points" << std::endl;
                        for (std::size_t k = 0; k < inter_pts.size(); ++k) {
                            std::cout << "      p" << k << ": (" << inter_pts[k][0] << "," << inter_pts[k][1] << ")" << std::endl;
                        }
                    }

                    std::vector<std::vector<double>> clipped;
                    for (std::size_t k = 0; k + 1 < inter_pts.size(); ++k) {
                        const std::vector<double>& p1 = inter_pts[k];
                        const std::vector<double>& p2 = inter_pts[k + 1];

                        const double dx1 = p1[0] - edge_start[0];
                        const double dy1 = p1[1] - edge_start[1];
                        const double dist1 = nx * dx1 + ny * dy1;

                        const double dx2 = p2[0] - edge_start[0];
                        const double dy2 = p2[1] - edge_start[1];
                        const double dist2 = nx * dx2 + ny * dy2;

                        const double eps = GEOMETRY_EPSILON;
                        const bool inside1 = dist1 >= -eps;
                        const bool inside2 = dist2 >= -eps;

                        if (inside1 && inside2) {
                            if (clipped.empty()) clipped.push_back(p1);
                            clipped.push_back(p2);
                        } else if (inside1 && !inside2) {
                            if (clipped.empty()) clipped.push_back(p1);
                            const double denom = nx * (p2[0] - p1[0]) + ny * (p2[1] - p1[1]);
                            const double abs_denom = (denom >= 0.0) ? denom : -denom;
                            if (abs_denom > eps) {
                                const double t = -dist1 / denom;
                                std::vector<double> p(2, 0.0);
                                p[0] = p1[0] + t * (p2[0] - p1[0]);
                                p[1] = p1[1] + t * (p2[1] - p1[1]);
                                clipped.push_back(p);
                            }
                        } else if (!inside1 && inside2) {
                            const double denom = nx * (p2[0] - p1[0]) + ny * (p2[1] - p1[1]);
                            const double abs_denom = (denom >= 0.0) ? denom : -denom;
                            if (abs_denom > eps) {
                                const double t = -dist1 / denom;
                                std::vector<double> p(2, 0.0);
                                p[0] = p1[0] + t * (p2[0] - p1[0]);
                                p[1] = p1[1] + t * (p2[1] - p1[1]);
                                clipped.push_back(p);
                            }
                            clipped.push_back(p2);
                        }
                    }
                    inter_pts = clipped;

                    if (debug) {
                        std::cout << "    inter_pts after clipping: " << inter_pts.size() << " points" << std::endl;
                        for (std::size_t k = 0; k < inter_pts.size(); ++k) {
                            std::cout << "      p" << k << ": (" << inter_pts[k][0] << "," << inter_pts[k][1] << ")" << std::endl;
                        }
                    }
                }

                if (inter_pts.empty()) {
                    return false;
                }
                result.vertices = inter_pts;
                result.intrinsic_dim = compute_intrinsic_dimension(inter_pts);
            } else {
                // Both are polygons with 3+ vertices.
                ConvexPolyhedron temp;
                if (!intersect_convex_polygons_2d(result, polyhedra[i], temp)) {
                    return false;
                }
                temp.intrinsic_dim = compute_intrinsic_dimension(temp.vertices);
                result = temp;
            }
        }
        intersection = result;
        intersection.intrinsic_dim = compute_intrinsic_dimension(intersection.vertices);
        return true;
    }

    // Fallback for other dimensions: use axis-aligned bounding boxes.
    std::vector<ConvexPolyhedronBox> boxes;
    boxes.reserve(polyhedra.size());
    for (const ConvexPolyhedron& poly : polyhedra) {
        boxes.push_back(bounding_box(poly));
    }

    ConvexPolyhedronBox inter_box;
    if (!intersect_boxes(boxes, inter_box)) {
        return false;
    }

    // Represent the intersection as the convex hull of the box corners.
    std::vector<std::vector<double>> corners;
    box_corners(inter_box, corners);
    intersection = convex_hull(corners);
    return true;
}

bool intersect_convex_polyhedron_with_last_coordinate_zero(
    const ConvexPolyhedron& polyhedron,
    ConvexPolyhedron& intersection)
{
    if (polyhedron.vertices.empty()) {
        return false;
    }

    const std::size_t ambient_dim = polyhedron.ambient_dimension();
    if (ambient_dim == 0u) {
        return false;
    }

    const std::size_t last = ambient_dim - 1u;
    const double eps = GEOMETRY_EPSILON;

    // For 2D and 3D ambient space, we compute the exact intersection of the convex hull
    // with the hyperplane x_{n+1} = 0 by intersecting segments with that hyperplane
    // and taking the convex hull of the resulting points.
    if (ambient_dim == 2u || ambient_dim == 3u) {
        ConvexPolyhedron hull = convex_hull(polyhedron.vertices);
        const std::vector<std::vector<double>>& verts = hull.vertices;
        const std::size_t n = verts.size();
        if (n == 0u) {
            return false;
        }

        std::vector<std::vector<double>> on_plane_points;
        on_plane_points.reserve(n);

        // Collect vertices that lie on the hyperplane.
        for (std::size_t i = 0; i < n; ++i) {
            const double val = verts[i][last];
            const double abs_val = (val >= 0.0) ? val : -val;
            if (abs_val <= eps) {
                std::vector<double> p = verts[i];
                p[last] = 0.0;
                on_plane_points.push_back(p);
            }
        }

        // For 2D: vertices are in counter-clockwise order, so adjacent vertices form edges.
        // For 3D: we need to check all pairs since we don't have explicit edge information.
        if (ambient_dim == 2u) {
            // Check only adjacent pairs (edges) for 2D polygons.
            // This is O(n) instead of O(n^2).
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t j = (i + 1) % n;  // Next vertex (wraps around)

                const double zi = verts[i][last];
                const double zj = verts[j][last];
                const double abs_zi = (zi >= 0.0) ? zi : -zi;
                const double abs_zj = (zj >= 0.0) ? zj : -zj;

                // Skip if either endpoint is (almost) on the plane; those were
                // already handled above.
                if (abs_zi <= eps || abs_zj <= eps) {
                    continue;
                }

                const bool above_i = zi > 0.0;
                const bool above_j = zj > 0.0;
                if (above_i == above_j) {
                    continue;
                }

                const double t = zi / (zi - zj);
                std::vector<double> p(ambient_dim, 0.0);
                for (std::size_t k = 0; k < ambient_dim; ++k) {
                    p[k] = verts[i][k] + t * (verts[j][k] - verts[i][k]);
                }
                p[last] = 0.0;
                on_plane_points.push_back(p);
            }
        } else {
            // For 3D: check all pairs of vertices (we don't have explicit edges).
            for (std::size_t i = 0; i < n; ++i) {
                const double zi = verts[i][last];
                const double abs_zi = (zi >= 0.0) ? zi : -zi;
                for (std::size_t j = i + 1; j < n; ++j) {
                    const double zj = verts[j][last];
                    const double abs_zj = (zj >= 0.0) ? zj : -zj;

                    // Skip if either endpoint is (almost) on the plane; those were
                    // already handled above.
                    if (abs_zi <= eps || abs_zj <= eps) {
                        continue;
                    }

                    const bool above_i = zi > 0.0;
                    const bool above_j = zj > 0.0;
                    if (above_i == above_j) {
                        continue;
                    }

                    const double t = zi / (zi - zj);
                    std::vector<double> p(ambient_dim, 0.0);
                    for (std::size_t k = 0; k < ambient_dim; ++k) {
                        p[k] = verts[i][k] + t * (verts[j][k] - verts[i][k]);
                    }
                    p[last] = 0.0;
                    on_plane_points.push_back(p);
                }
            }
        }

        if (on_plane_points.empty()) {
            return false;
        }

        intersection = convex_hull(on_plane_points);
        return true;
    }

    // Fallback for other dimensions: use axis-aligned bounding boxes.
    const ConvexPolyhedronBox box = bounding_box(polyhedron);
    ConvexPolyhedronBox sliced_box;
    if (!intersect_box_with_last_coordinate_zero_impl(box, sliced_box)) {
        return false;
    }

    std::vector<std::vector<double>> corners;
    box_corners(sliced_box, corners);
    intersection = convex_hull(corners);
    return true;
}

} // namespace polynomial_solver

