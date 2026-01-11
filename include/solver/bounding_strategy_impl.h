#ifndef BOUNDING_STRATEGY_IMPL_H
#define BOUNDING_STRATEGY_IMPL_H

/**
 * @file detail/bounding_strategy_impl.h
 * @brief Implementation of bounding strategy classes
 */

#include <algorithm>
#include <cmath>
#include <map>

namespace polynomial_solver {

//=============================================================================
// ProjectedPolyhedralStrategy - Direction-by-direction projection
//=============================================================================

/**
 * @class ProjectedPolyhedralStrategy
 * @brief Projects control points direction-by-direction to compute bounds
 *
 * For each direction i:
 * 1. Project graph control points to 2D (coord i + function value)
 * 2. Compute convex hull of these 2D points
 * 3. Intersect with y=0 axis
 * 4. Project to 1D to get interval bound
 * 5. Intersect all intervals from all equations
 */
template<typename Scalar>
class ProjectedPolyhedralStrategy : public BoundingStrategyBase<Scalar> {
public:
    using Poly = PolynomialBase<Scalar>;
    using Result = BoundingResult<Scalar>;

    Result computeBounds(
        const std::vector<Poly>& polys,
        std::size_t dim,
        const std::vector<Scalar>& current_lower,
        const std::vector<Scalar>& current_upper) const override
    {
        Result result;
        result.lower.resize(dim, Scalar(0));
        result.upper.resize(dim, Scalar(1));
        result.has_root = true;
        
        Scalar total_old_volume = Scalar(1);
        Scalar total_new_volume = Scalar(1);
        
        // Process each direction independently
        for (std::size_t dir = 0; dir < dim; ++dir) {
            Scalar dir_min = Scalar(0);
            Scalar dir_max = Scalar(1);
            
            for (std::size_t eq_idx = 0; eq_idx < polys.size(); ++eq_idx) {
                const Poly& poly = polys[eq_idx];
                
                // Get graph control points
                std::vector<Scalar> control_points;
                poly.graphControlPoints(control_points);
                
                const std::size_t point_dim = dim + 1;
                const std::size_t num_points = control_points.size() / point_dim;
                
                // Group by x-coordinate and find extremes
                std::map<Scalar, std::vector<Scalar>> vertical_groups;
                for (std::size_t i = 0; i < num_points; ++i) {
                    Scalar x = control_points[i * point_dim + dir];
                    Scalar y = control_points[i * point_dim + dim];
                    vertical_groups[x].push_back(y);
                }
                
                // Build projected 2D points
                std::vector<std::vector<Scalar>> projected_2d;
                projected_2d.reserve(2 * vertical_groups.size());
                
                for (auto it = vertical_groups.begin(); it != vertical_groups.end(); ++it) {
                    Scalar x = it->first;
                    const std::vector<Scalar>& y_values = it->second;
                    
                    Scalar y_min = *std::min_element(y_values.begin(), y_values.end());
                    Scalar y_max = *std::max_element(y_values.begin(), y_values.end());
                    
                    projected_2d.push_back({x, y_min});
                    if (y_min != y_max) {
                        projected_2d.push_back({x, y_max});
                    }
                }
                
                // Compute convex hull
                ConvexPolyhedronBase<Scalar> hull = convex_hull_impl(projected_2d);
                
                // Find intersection with y=0 axis
                Scalar eq_min, eq_max;
                if (!findYZeroIntersection(hull, eq_min, eq_max)) {
                    result.has_root = false;
                    return result;
                }
                
                // Intersect with current bounds
                dir_min = std::max(dir_min, eq_min);
                dir_max = std::min(dir_max, eq_max);
                
                if (dir_min > dir_max) {
                    result.has_root = false;
                    return result;
                }
            }
            
            Scalar old_width = current_upper.empty() ? Scalar(1) 
                             : (current_upper[dir] - current_lower[dir]);
            Scalar new_width = dir_max - dir_min;
            
            total_old_volume *= old_width;
            total_new_volume *= new_width;
            
            result.lower[dir] = dir_min;
            result.upper[dir] = dir_max;
        }
        
        result.contraction_ratio = (total_old_volume > Scalar(0)) 
                                 ? (total_new_volume / total_old_volume) 
                                 : Scalar(0);
        return result;
    }

    const char* name() const override { return "ProjectedPolyhedral"; }

private:
    bool findYZeroIntersection(
        const ConvexPolyhedronBase<Scalar>& hull,
        Scalar& x_min, Scalar& x_max) const
    {
        const std::vector<std::vector<Scalar>>& verts = hull.vertices;
        if (verts.size() < 2) {
            return false;
        }
        
        std::vector<Scalar> x_crossings;
        std::size_t n = verts.size();
        
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t j = (i + 1) % n;
            Scalar y1 = verts[i][1];
            Scalar y2 = verts[j][1];
            
            // Check if edge crosses y=0
            if ((y1 <= Scalar(0) && y2 >= Scalar(0)) || 
                (y1 >= Scalar(0) && y2 <= Scalar(0))) {
                if (y1 == y2) {
                    x_crossings.push_back(verts[i][0]);
                    x_crossings.push_back(verts[j][0]);
                } else {
                    Scalar t = -y1 / (y2 - y1);
                    Scalar x_cross = verts[i][0] + t * (verts[j][0] - verts[i][0]);
                    x_crossings.push_back(x_cross);
                }
            }
        }
        
        if (x_crossings.empty()) {
            return false;
        }
        
        x_min = *std::min_element(x_crossings.begin(), x_crossings.end());
        x_max = *std::max_element(x_crossings.begin(), x_crossings.end());
        return true;
    }
};

//=============================================================================
// GraphHullStrategy - Exact convex hull in graph space (1D/2D only)
//=============================================================================

/**
 * @class GraphHullStrategy
 * @brief Uses exact convex hull of graph control points
 *
 * For each equation:
 * 1. Compute convex hull of graph control points in R^{n+1}
 * 2. Intersect with hyperplane x_{n+1} = 0
 * 3. Project to R^n and intersect all results
 * 4. Return bounding box
 *
 * Only implemented for 1D and 2D systems.
 */
template<typename Scalar>
class GraphHullStrategy : public BoundingStrategyBase<Scalar> {
public:
    using Poly = PolynomialBase<Scalar>;
    using Result = BoundingResult<Scalar>;

    Result computeBounds(
        const std::vector<Poly>& polys,
        std::size_t dim,
        const std::vector<Scalar>& current_lower,
        const std::vector<Scalar>& current_upper) const override
    {
        Result result;
        result.lower.resize(dim, Scalar(0));
        result.upper.resize(dim, Scalar(1));
        result.has_root = true;

        // Only implemented for 1D and 2D
        if (dim > 2) {
            // Fall back to input box
            if (!current_lower.empty()) {
                result.lower = current_lower;
                result.upper = current_upper;
            }
            result.contraction_ratio = Scalar(1);
            return result;
        }

        // Compute bounding box from each equation's graph hull
        for (std::size_t eq_idx = 0; eq_idx < polys.size(); ++eq_idx) {
            const Poly& poly = polys[eq_idx];

            // Get graph control points in R^{n+1}
            std::vector<Scalar> control_points;
            poly.graphControlPoints(control_points);

            const std::size_t point_dim = dim + 1;
            const std::size_t num_points = control_points.size() / point_dim;

            // Convert to vector of vectors
            std::vector<std::vector<Scalar>> graph_points;
            graph_points.reserve(num_points);
            for (std::size_t i = 0; i < num_points; ++i) {
                std::vector<Scalar> pt(point_dim);
                for (std::size_t j = 0; j < point_dim; ++j) {
                    pt[j] = control_points[i * point_dim + j];
                }
                graph_points.push_back(pt);
            }

            // Compute convex hull in R^{n+1}
            ConvexPolyhedronBase<Scalar> hull = convex_hull_impl(graph_points);

            // Intersect with hyperplane x_{n+1} = 0
            ConvexPolyhedronBase<Scalar> intersection =
                intersect_convex_polyhedron_with_last_coordinate_zero_impl(hull);

            if (intersection.vertices.empty()) {
                result.has_root = false;
                return result;
            }

            // Project to R^n and compute bounding box
            std::vector<std::vector<Scalar>> projected;
            projected.reserve(intersection.vertices.size());
            for (const auto& v : intersection.vertices) {
                std::vector<Scalar> proj(dim);
                for (std::size_t j = 0; j < dim; ++j) {
                    proj[j] = v[j];
                }
                projected.push_back(proj);
            }

            ConvexPolyhedronBoxBase<Scalar> bbox = bounding_box_impl(projected);

            // Intersect with current bounds
            for (std::size_t j = 0; j < dim; ++j) {
                result.lower[j] = std::max(result.lower[j], bbox.lower[j]);
                result.upper[j] = std::min(result.upper[j], bbox.upper[j]);

                if (result.lower[j] > result.upper[j]) {
                    result.has_root = false;
                    return result;
                }
            }
        }

        // Compute contraction ratio
        Scalar total_old_volume = Scalar(1);
        Scalar total_new_volume = Scalar(1);
        for (std::size_t j = 0; j < dim; ++j) {
            Scalar old_width = current_upper.empty() ? Scalar(1)
                             : (current_upper[j] - current_lower[j]);
            Scalar new_width = result.upper[j] - result.lower[j];
            total_old_volume *= old_width;
            total_new_volume *= new_width;
        }

        result.contraction_ratio = (total_old_volume > Scalar(0))
                                 ? (total_new_volume / total_old_volume)
                                 : Scalar(0);
        return result;
    }

    const char* name() const override { return "GraphHull"; }

    bool supportsDimension(std::size_t dim) const override {
        return dim <= 2;
    }
};

//=============================================================================
// IntervalArithmeticStrategy - Placeholder for interval-based bounds
//=============================================================================

/**
 * @class IntervalArithmeticStrategy
 * @brief Placeholder for interval arithmetic bounding (future implementation)
 */
template<typename Scalar>
class IntervalArithmeticStrategy : public BoundingStrategyBase<Scalar> {
public:
    using Poly = PolynomialBase<Scalar>;
    using Result = BoundingResult<Scalar>;

    Result computeBounds(
        const std::vector<Poly>& /*polys*/,
        std::size_t dim,
        const std::vector<Scalar>& current_lower,
        const std::vector<Scalar>& current_upper) const override
    {
        // TODO: Implement interval arithmetic bounding
        // For now, just return the current box
        Result result;
        result.has_root = true;
        result.lower = current_lower.empty() ? std::vector<Scalar>(dim, Scalar(0)) : current_lower;
        result.upper = current_upper.empty() ? std::vector<Scalar>(dim, Scalar(1)) : current_upper;
        result.contraction_ratio = Scalar(1);
        return result;
    }

    const char* name() const override { return "IntervalArithmetic"; }
};

//=============================================================================
// Factory function to create strategies
//=============================================================================

/**
 * @brief Create a bounding strategy by type
 */
template<typename Scalar>
std::unique_ptr<BoundingStrategyBase<Scalar>> createBoundingStrategy(RootBoundingMethodBase method) {
    switch (method) {
        case RootBoundingMethodBase::ProjectedPolyhedral:
            return std::unique_ptr<BoundingStrategyBase<Scalar>>(new ProjectedPolyhedralStrategy<Scalar>());
        case RootBoundingMethodBase::GraphHull:
            return std::unique_ptr<BoundingStrategyBase<Scalar>>(new GraphHullStrategy<Scalar>());
        case RootBoundingMethodBase::IntervalArithmetic:
            return std::unique_ptr<BoundingStrategyBase<Scalar>>(new IntervalArithmeticStrategy<Scalar>());
        case RootBoundingMethodBase::None:
        default:
            return std::unique_ptr<BoundingStrategyBase<Scalar>>(new NoBoundingStrategy<Scalar>());
    }
}

} // namespace polynomial_solver

#endif // BOUNDING_STRATEGY_IMPL_H

