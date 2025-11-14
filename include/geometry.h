#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cstddef>
#include <vector>

/**
 * @file geometry.h
 * @brief Geometry tools for polynomial solver
 *
 * This module provides geometric utilities for working with polynomial curves
 * and their geometric envelopes. In particular, it provides simple
 * representations of points, intervals, and axis-aligned convex polyhedra
 * (hyper-rectangles) together with basic operations on them.
 */

namespace polynomial_solver {

/**
 * @brief Default geometric tolerance for numerical comparisons.
 *
 * This tolerance is used for:
 * - Detecting collinear/coplanar points
 * - Filtering duplicate points
 * - Comparing distances in geometric algorithms
 */
constexpr double GEOMETRY_EPSILON = 1e-15;

/**
 * @brief Global configuration for geometry operations.
 */
struct GeometryConfig {
    bool debug;        ///< Enable debug output for geometric operations (default: false)

    GeometryConfig() : debug(false) {}
};

/**
 * @brief Get the global geometry configuration.
 */
GeometryConfig& getGeometryConfig();

/**
 * @class Point2D
 * @brief Represents a 2D point
 */
class Point2D {
public:
    Point2D();
    Point2D(double x, double y);
    ~Point2D();

private:
    double x_;
    double y_;
};

/**
 * @class Interval
 * @brief Represents a 1D interval
 */
class Interval {
public:
    Interval();
    Interval(double min, double max);
    ~Interval();

private:
    double min_;
    double max_;
};

/**
 * @struct ConvexPolyhedron
 * @brief Generic convex polyhedron in R^{d}, represented by its vertices.
 *
 * Semantically, this is the convex hull of the stored points.
 *
 * The structure tracks both:
 * - ambient_dim: The dimension of the ambient space (size of each vertex vector)
 * - intrinsic_dim: The intrinsic dimension of the polytope itself
 *   (0 for point, 1 for line segment, 2 for polygon, etc.)
 *
 * For example, a line segment in 3D has ambient_dim=3 and intrinsic_dim=1.
 */
struct ConvexPolyhedron {
    std::vector<std::vector<double>> vertices;
    std::size_t intrinsic_dim; ///< Intrinsic dimension of the polytope

    ConvexPolyhedron() : intrinsic_dim(0) {}

    /// Returns the ambient dimension (dimension of the space)
    std::size_t ambient_dimension() const {
        return vertices.empty() ? 0u : vertices.front().size();
    }

    /// Returns the intrinsic dimension (dimension of the polytope itself)
    std::size_t dimension() const {
        return intrinsic_dim;
    }
};

/**
 * @struct ConvexPolyhedronBox
 * @brief Simple axis-aligned bounding box in R^{d}.
 */
struct ConvexPolyhedronBox {
    std::vector<double> min_coords; ///< Lower bounds in each coordinate.
    std::vector<double> max_coords; ///< Upper bounds in each coordinate.

    std::size_t dimension() const { return min_coords.size(); }
};

/**
 * @brief Compute the intrinsic dimension of a point set.
 *
 * Returns the dimension of the affine hull of the points.
 * For example, collinear points have intrinsic dimension 1,
 * coplanar points have intrinsic dimension 2, etc.
 */
std::size_t compute_intrinsic_dimension(
    const std::vector<std::vector<double>>& points);

/**
 * @brief Construct a convex polyhedron as the convex hull of a set of points.
 *
 * Uses Chan's algorithm for 2D and 3D cases, which is output-sensitive
 * and handles degenerate cases robustly.
 * The intrinsic dimension is automatically computed and stored.
 */
ConvexPolyhedron
convex_hull(const std::vector<std::vector<double>>& points);

/**
 * @brief Intersect a collection of convex polyhedra.
 *
 * On success, @p intersection stores a convex polyhedron that contains the
 * true geometric intersection (it may be a conservative approximation).
 */
bool intersect_convex_polyhedra(
    const std::vector<ConvexPolyhedron>& polyhedra,
    ConvexPolyhedron& intersection);

/**
 * @brief Intersect a convex polyhedron with the hyperplane x_{n+1} = 0.
 *
 * The dimension is taken from the polyhedron. If the last coordinate interval
 * does not cross 0, the function returns false.
 */
bool intersect_convex_polyhedron_with_last_coordinate_zero(
    const ConvexPolyhedron& polyhedron,
    ConvexPolyhedron& intersection);

/**
 * @brief Compute the axis-aligned bounding box of a convex polyhedron.
 */
ConvexPolyhedronBox
bounding_box(const ConvexPolyhedron& polyhedron);

} // namespace polynomial_solver

#endif // GEOMETRY_H

