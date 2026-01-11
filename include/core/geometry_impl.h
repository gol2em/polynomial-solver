#ifndef GEOMETRY_IMPL_H
#define GEOMETRY_IMPL_H

/**
 * @file detail/geometry_impl.h
 * @brief Template implementation of geometry algorithms
 *
 * This file contains the templated implementations of geometric algorithms.
 * It is included by geometry_base.h and should not be included directly.
 */

#include <algorithm>
#include <cmath>
#include <limits>

namespace polynomial_solver {

//=============================================================================
// Helper functions
//=============================================================================

namespace detail {

/// Compute squared distance between two points
template<typename Scalar>
Scalar squared_distance(const std::vector<Scalar>& a, const std::vector<Scalar>& b) {
    Scalar sum = Scalar(0);
    const std::size_t dim = a.size();
    for (std::size_t i = 0; i < dim; ++i) {
        const Scalar diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum;
}

/// 2D cross product (z-component of 3D cross product)
template<typename Scalar>
Scalar cross2d(const std::vector<Scalar>& o,
               const std::vector<Scalar>& a,
               const std::vector<Scalar>& b) {
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
}

/// Check if points are collinear
template<typename Scalar>
bool are_collinear(const std::vector<std::vector<Scalar>>& points, Scalar eps) {
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
        return true;
    }

    // Direction vector
    std::vector<Scalar> dir(dim);
    Scalar dir_len_sq = Scalar(0);
    for (std::size_t j = 0; j < dim; ++j) {
        dir[j] = points[i1][j] - points[i0][j];
        dir_len_sq += dir[j] * dir[j];
    }

    // Check if all other points lie on the line
    for (std::size_t i = 0; i < n; ++i) {
        if (i == i0 || i == i1) continue;

        std::vector<Scalar> v(dim);
        for (std::size_t j = 0; j < dim; ++j) {
            v[j] = points[i][j] - points[i0][j];
        }

        Scalar cross_sq = Scalar(0);
        if (dim == 2) {
            Scalar cross_val = v[0] * dir[1] - v[1] * dir[0];
            cross_sq = cross_val * cross_val;
        } else if (dim == 3) {
            Scalar cx = v[1] * dir[2] - v[2] * dir[1];
            Scalar cy = v[2] * dir[0] - v[0] * dir[2];
            Scalar cz = v[0] * dir[1] - v[1] * dir[0];
            cross_sq = cx * cx + cy * cy + cz * cz;
        } else {
            // General case: perpendicular distance
            Scalar dot = Scalar(0);
            for (std::size_t j = 0; j < dim; ++j) {
                dot += v[j] * dir[j];
            }
            Scalar proj_len = dot / dir_len_sq;
            Scalar perp_sq = Scalar(0);
            for (std::size_t j = 0; j < dim; ++j) {
                Scalar perp_j = v[j] - proj_len * dir[j];
                perp_sq += perp_j * perp_j;
            }
            cross_sq = perp_sq;
        }

        if (cross_sq > eps * eps * dir_len_sq) {
            return false;
        }
    }
    return true;
}

/// Check if points are coplanar (3D only)
template<typename Scalar>
bool are_coplanar_3d(const std::vector<std::vector<Scalar>>& points, Scalar eps) {
    if (points.size() <= 3) {
        return true;
    }

    const std::size_t n = points.size();
    std::size_t i0 = 0, i1 = 1, i2 = 2;

    // Find i1 != i0
    while (i1 < n && squared_distance(points[i0], points[i1]) < eps * eps) {
        ++i1;
    }
    if (i1 >= n) return true;

    // Find i2 not collinear
    while (i2 < n) {
        if (i2 == i0 || i2 == i1) { ++i2; continue; }
        std::vector<std::vector<Scalar>> three = {points[i0], points[i1], points[i2]};
        if (!are_collinear(three, eps)) break;
        ++i2;
    }
    if (i2 >= n) return true;

    // Compute plane normal
    std::vector<Scalar> v1(3), v2(3);
    for (std::size_t j = 0; j < 3; ++j) {
        v1[j] = points[i1][j] - points[i0][j];
        v2[j] = points[i2][j] - points[i0][j];
    }

    std::vector<Scalar> normal(3);
    normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
    normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
    normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

    Scalar normal_len_sq = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
    if (normal_len_sq < eps * eps) return true;

    // Check all other points
    for (std::size_t i = 0; i < n; ++i) {
        if (i == i0 || i == i1 || i == i2) continue;
        std::vector<Scalar> v(3);
        for (std::size_t j = 0; j < 3; ++j) {
            v[j] = points[i][j] - points[i0][j];
        }
        Scalar dot = v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2];
        Scalar dist_sq = (dot * dot) / normal_len_sq;
        if (dist_sq > eps * eps) return false;
    }
    return true;
}

} // namespace detail

//=============================================================================
// Intrinsic dimension computation
//=============================================================================

template<typename Scalar>
std::size_t compute_intrinsic_dimension_impl(
    const std::vector<std::vector<Scalar>>& points,
    Scalar eps)
{
    if (points.empty()) return 0;
    if (points.size() == 1) return 0;

    const std::size_t ambient_dim = points[0].size();

    // Check if all points are identical
    bool all_same = true;
    for (std::size_t i = 1; i < points.size() && all_same; ++i) {
        if (detail::squared_distance(points[0], points[i]) > eps * eps) {
            all_same = false;
        }
    }
    if (all_same) return 0;

    // Check if collinear
    if (detail::are_collinear(points, eps)) return 1;

    // For 3D, check coplanar
    if (ambient_dim == 3 && detail::are_coplanar_3d(points, eps)) return 2;

    // Otherwise, full dimension (up to ambient_dim - 1 for affine hull)
    return (ambient_dim >= 2) ? 2 : 1;
}

//=============================================================================
// Bounding box computation
//=============================================================================

template<typename Scalar>
ConvexPolyhedronBoxBase<Scalar>
bounding_box_impl(const ConvexPolyhedronBase<Scalar>& polyhedron)
{
    ConvexPolyhedronBoxBase<Scalar> box;
    if (polyhedron.vertices.empty()) {
        return box;
    }

    const std::size_t dim = polyhedron.ambient_dimension();
    box.min_coords.assign(dim, std::numeric_limits<Scalar>::max());
    box.max_coords.assign(dim, std::numeric_limits<Scalar>::lowest());

    for (const auto& v : polyhedron.vertices) {
        for (std::size_t j = 0; j < dim; ++j) {
            if (v[j] < box.min_coords[j]) box.min_coords[j] = v[j];
            if (v[j] > box.max_coords[j]) box.max_coords[j] = v[j];
        }
    }
    return box;
}

/// Build bounding box from points
template<typename Scalar>
ConvexPolyhedronBoxBase<Scalar>
bounding_box_from_points(const std::vector<std::vector<Scalar>>& points)
{
    ConvexPolyhedronBoxBase<Scalar> box;
    if (points.empty()) return box;

    const std::size_t dim = points[0].size();
    box.min_coords.assign(dim, std::numeric_limits<Scalar>::max());
    box.max_coords.assign(dim, std::numeric_limits<Scalar>::lowest());

    for (const auto& p : points) {
        if (p.size() != dim) continue;
        for (std::size_t j = 0; j < dim; ++j) {
            if (p[j] < box.min_coords[j]) box.min_coords[j] = p[j];
            if (p[j] > box.max_coords[j]) box.max_coords[j] = p[j];
        }
    }
    return box;
}

/// Generate corners of a bounding box
template<typename Scalar>
void box_corners(const ConvexPolyhedronBoxBase<Scalar>& box,
                 std::vector<std::vector<Scalar>>& corners)
{
    const std::size_t dim = box.dimension();
    corners.clear();
    if (dim == 0) return;

    const std::size_t num_corners = std::size_t(1) << dim;
    corners.reserve(num_corners);

    for (std::size_t mask = 0; mask < num_corners; ++mask) {
        std::vector<Scalar> vertex(dim);
        for (std::size_t j = 0; j < dim; ++j) {
            const bool use_max = ((mask >> j) & 1) != 0;
            vertex[j] = use_max ? box.max_coords[j] : box.min_coords[j];
        }
        corners.push_back(vertex);
    }
}

//=============================================================================
// 2D Convex Hull (Monotone Chain Algorithm)
//=============================================================================

template<typename Scalar>
ConvexPolyhedronBase<Scalar>
convex_hull_2d(std::vector<std::vector<Scalar>>& pts, Scalar eps)
{
    ConvexPolyhedronBase<Scalar> poly;

    // Sort lexicographically
    std::sort(pts.begin(), pts.end(),
              [](const std::vector<Scalar>& a, const std::vector<Scalar>& b) {
                  if (a[0] < b[0]) return true;
                  if (a[0] > b[0]) return false;
                  return a[1] < b[1];
              });

    // Remove duplicates
    auto last = std::unique(pts.begin(), pts.end(),
                           [eps](const std::vector<Scalar>& a, const std::vector<Scalar>& b) {
                               Scalar dx = a[0] - b[0];
                               Scalar dy = a[1] - b[1];
                               return (dx*dx + dy*dy) < eps*eps;
                           });
    pts.erase(last, pts.end());

    if (pts.empty()) {
        poly.intrinsic_dim = 0;
        return poly;
    }
    if (pts.size() == 1) {
        poly.vertices = pts;
        poly.intrinsic_dim = 0;
        return poly;
    }
    if (pts.size() == 2) {
        poly.vertices = pts;
        poly.intrinsic_dim = 1;
        return poly;
    }

    // Build lower hull
    std::vector<std::vector<Scalar>> lower;
    for (const auto& p : pts) {
        while (lower.size() >= 2 &&
               detail::cross2d(lower[lower.size()-2], lower[lower.size()-1], p) <= Scalar(0)) {
            lower.pop_back();
        }
        lower.push_back(p);
    }

    // Build upper hull
    std::vector<std::vector<Scalar>> upper;
    for (auto it = pts.rbegin(); it != pts.rend(); ++it) {
        while (upper.size() >= 2 &&
               detail::cross2d(upper[upper.size()-2], upper[upper.size()-1], *it) <= Scalar(0)) {
            upper.pop_back();
        }
        upper.push_back(*it);
    }

    // Concatenate hulls (remove duplicate endpoints)
    poly.vertices.reserve(lower.size() + upper.size() - 2);
    for (std::size_t i = 0; i < lower.size(); ++i) {
        poly.vertices.push_back(lower[i]);
    }
    for (std::size_t i = 1; i + 1 < upper.size(); ++i) {
        poly.vertices.push_back(upper[i]);
    }

    // Determine intrinsic dimension
    if (poly.vertices.size() <= 1) {
        poly.intrinsic_dim = 0;
    } else if (poly.vertices.size() == 2) {
        poly.intrinsic_dim = 1;
    } else {
        poly.intrinsic_dim = detail::are_collinear(poly.vertices, eps) ? 1 : 2;
    }

    return poly;
}

//=============================================================================
// Main convex hull function
//=============================================================================

template<typename Scalar>
ConvexPolyhedronBase<Scalar>
convex_hull_impl(const std::vector<std::vector<Scalar>>& points)
{
    const Scalar eps = GeometryToleranceTraits<Scalar>::epsilon();
    ConvexPolyhedronBase<Scalar> poly;

    if (points.empty()) {
        poly.intrinsic_dim = 0;
        return poly;
    }

    const std::size_t dim = points[0].size();

    // Filter mismatched dimensions
    std::vector<std::vector<Scalar>> pts;
    pts.reserve(points.size());
    for (const auto& p : points) {
        if (p.size() == dim) pts.push_back(p);
    }

    if (pts.empty()) {
        poly.intrinsic_dim = 0;
        return poly;
    }

    // Compute intrinsic dimension
    poly.intrinsic_dim = compute_intrinsic_dimension_impl(pts, eps);

    // For 2D, use monotone chain
    if (dim == 2) {
        return convex_hull_2d(pts, eps);
    }

    // For other dimensions, return all points (conservative approximation)
    // 3D convex hull is complex and left for future work
    poly.vertices = pts;
    return poly;
}

//=============================================================================
// Box intersection with hyperplane x_n = 0
//=============================================================================

template<typename Scalar>
bool intersect_box_with_last_coordinate_zero(
    const ConvexPolyhedronBoxBase<Scalar>& box,
    ConvexPolyhedronBoxBase<Scalar>& intersection)
{
    const std::size_t dim = box.dimension();
    if (dim == 0) return false;

    intersection = box;
    const std::size_t last = dim - 1;
    const Scalar min_v = intersection.min_coords[last];
    const Scalar max_v = intersection.max_coords[last];

    if (min_v > Scalar(0) || max_v < Scalar(0)) {
        return false;
    }

    intersection.min_coords[last] = Scalar(0);
    intersection.max_coords[last] = Scalar(0);
    return true;
}

//=============================================================================
// Polyhedron intersection with hyperplane x_n = 0
//=============================================================================

template<typename Scalar>
bool intersect_convex_polyhedron_with_last_coordinate_zero_impl(
    const ConvexPolyhedronBase<Scalar>& polyhedron,
    ConvexPolyhedronBase<Scalar>& intersection)
{
    if (polyhedron.vertices.empty()) return false;

    const std::size_t ambient_dim = polyhedron.ambient_dimension();
    if (ambient_dim == 0) return false;

    const std::size_t last = ambient_dim - 1;
    const Scalar eps = GeometryToleranceTraits<Scalar>::epsilon();

    std::vector<std::vector<Scalar>> on_plane_points;

    // Collect vertices on the hyperplane
    for (const auto& v : polyhedron.vertices) {
        Scalar z = v[last];
        if (z >= -eps && z <= eps) {
            std::vector<Scalar> p = v;
            p[last] = Scalar(0);
            on_plane_points.push_back(p);
        }
    }

    // Check all pairs of vertices for edge crossings
    const std::size_t n = polyhedron.vertices.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const auto& vi = polyhedron.vertices[i];
            const auto& vj = polyhedron.vertices[j];
            Scalar zi = vi[last];
            Scalar zj = vj[last];

            // Check if edge crosses the hyperplane
            if ((zi < -eps && zj > eps) || (zi > eps && zj < -eps)) {
                // Compute intersection point
                Scalar t = -zi / (zj - zi);
                std::vector<Scalar> p(ambient_dim);
                for (std::size_t k = 0; k < ambient_dim; ++k) {
                    p[k] = vi[k] + t * (vj[k] - vi[k]);
                }
                p[last] = Scalar(0);
                on_plane_points.push_back(p);
            }
        }
    }

    if (on_plane_points.empty()) {
        return false;
    }

    intersection = convex_hull_impl(on_plane_points);
    return true;
}

//=============================================================================
// Intersect collection of convex polyhedra (using bounding boxes for now)
//=============================================================================

template<typename Scalar>
bool intersect_boxes(
    const std::vector<ConvexPolyhedronBoxBase<Scalar>>& boxes,
    ConvexPolyhedronBoxBase<Scalar>& intersection)
{
    if (boxes.empty()) return false;

    const std::size_t dim = boxes[0].dimension();
    intersection = boxes[0];

    for (std::size_t k = 1; k < boxes.size(); ++k) {
        const auto& box = boxes[k];
        if (box.dimension() != dim) return false;

        for (std::size_t j = 0; j < dim; ++j) {
            if (box.min_coords[j] > intersection.min_coords[j])
                intersection.min_coords[j] = box.min_coords[j];
            if (box.max_coords[j] < intersection.max_coords[j])
                intersection.max_coords[j] = box.max_coords[j];
        }
    }

    // Check for empty intersection
    for (std::size_t j = 0; j < dim; ++j) {
        if (intersection.min_coords[j] > intersection.max_coords[j]) {
            return false;
        }
    }
    return true;
}

template<typename Scalar>
bool intersect_convex_polyhedra_impl(
    const std::vector<ConvexPolyhedronBase<Scalar>>& polyhedra,
    ConvexPolyhedronBase<Scalar>& intersection)
{
    if (polyhedra.empty()) return false;

    // Use bounding box intersection as conservative approximation
    std::vector<ConvexPolyhedronBoxBase<Scalar>> boxes;
    boxes.reserve(polyhedra.size());

    for (const auto& poly : polyhedra) {
        boxes.push_back(bounding_box_impl(poly));
    }

    ConvexPolyhedronBoxBase<Scalar> inter_box;
    if (!intersect_boxes(boxes, inter_box)) {
        return false;
    }

    // Convert box to polyhedron
    std::vector<std::vector<Scalar>> corners;
    box_corners(inter_box, corners);
    intersection = convex_hull_impl(corners);
    return true;
}

} // namespace polynomial_solver

#endif // GEOMETRY_IMPL_H

