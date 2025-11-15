#include "solver.h"
#include <queue>
#include <map>


/**
 * @file solver.cpp
 * @brief Implementation of the main solver interface
 */

namespace polynomial_solver {

PolynomialSystem::PolynomialSystem()
    : dimension_(0u)
{
}

PolynomialSystem::PolynomialSystem(const std::vector<Polynomial>& equations)
    : dimension_(0u),
      equations_(equations)
{
    if (!equations_.empty()) {
        dimension_ = equations_[0].dimension();
        // TODO: Optionally verify that all equations share this dimension.
    }
}

std::size_t PolynomialSystem::dimension() const {
    return dimension_;
}

std::size_t PolynomialSystem::equationCount() const {
    return equations_.size();
}

const Polynomial& PolynomialSystem::equation(std::size_t i) const {
    return equations_[i];
}

const std::vector<Polynomial>& PolynomialSystem::equations() const {
    return equations_;
}
std::vector<GraphControlNet> PolynomialSystem::graphControlNets() const {
    std::vector<GraphControlNet> nets;
    nets.reserve(equations_.size());

    for (const Polynomial& poly : equations_) {
        GraphControlNet net;
        net.degrees = poly.degrees();
        poly.graphControlPoints(net.control_points);
        nets.push_back(std::move(net));
    }

    return nets;
}




Solver::Solver() {
    // TODO: Implement constructor
}

Solver::~Solver() {
    // TODO: Implement destructor
}

namespace {

struct SubdivisionNode {
    std::vector<double> box_lower;  ///< Global box lower corner in [0,1]^n.
    std::vector<double> box_upper;  ///< Global box upper corner in [0,1]^n.
    unsigned int depth;             ///< Subdivision depth.

    // Polynomials restricted to this node's box and re-parameterised to [0,1]^n.
    std::vector<polynomial_solver::Polynomial> polys;
};

struct NodeQueueEntry {
    unsigned int depth;
    std::size_t index;  ///< Tie-breaker to keep ordering stable.
    SubdivisionNode node;
};

struct NodeQueueCompare {
    bool operator()(const NodeQueueEntry& lhs, const NodeQueueEntry& rhs) const {
        // std::priority_queue is a max-heap. We want the smallest depth first.
        if (lhs.depth != rhs.depth) {
            return lhs.depth > rhs.depth;
        }
        return lhs.index > rhs.index;
    }
};

/**
 * @brief Compute root bounding box using graph convex hull method.
 *
 * For each equation f_i(x) = 0, we build the graph control net (x, f_i(x)) in R^{n+1}.
 * The roots must lie in the convex hull of these control points. We:
 * 1. Compute the convex hull of each equation's graph control points in R^{n+1}.
 * 2. Intersect each hull with the hyperplane x_{n+1} = 0.
 * 3. Project each result to R^n by dropping the last coordinate (which is 0).
 * 4. Intersect all these projected polyhedra in R^n.
 * 5. Compute the axis-aligned bounding box in R^n (parameter space).
 *
 * Returns true if a non-empty bounding box is found, false otherwise.
 */
bool compute_graph_hull_bounds(
    const std::vector<polynomial_solver::Polynomial>& polys,
    std::size_t dim,
    std::vector<double>& local_bound_lower,
    std::vector<double>& local_bound_upper)
{
    if (polys.empty() || dim == 0u) {
        return false;
    }

    // Only implemented for 1D and 2D systems (graphs in R^2 and R^3).
    if (dim > 2u) {
        return false;
    }

    // Step 1: For each equation, compute convex hull of graph control points in R^{n+1}.
    // Step 2: For each equation, intersect that hull with hyperplane x_{n+1} = 0.
    // This gives us polyhedra in R^n (parameter space with last coord = 0).
    std::vector<polynomial_solver::ConvexPolyhedron> hyperplane_intersections;
    hyperplane_intersections.reserve(polys.size());

    for (const polynomial_solver::Polynomial& poly : polys) {
        std::vector<double> control_points;
        poly.graphControlPoints(control_points);

        // Convert flat array to vector of points in R^{n+1}.
        const std::size_t num_coeffs = poly.coefficientCount();
        const std::size_t point_dim = dim + 1u;

        std::vector<std::vector<double>> points;
        points.reserve(num_coeffs);

        for (std::size_t i = 0; i < num_coeffs; ++i) {
            std::vector<double> pt(point_dim, 0.0);
            for (std::size_t j = 0; j < point_dim; ++j) {
                pt[j] = control_points[i * point_dim + j];
            }
            points.push_back(pt);
        }

        // Compute convex hull of graph control points in R^{n+1}.
        polynomial_solver::ConvexPolyhedron hull = polynomial_solver::convex_hull(points);

        // Intersect this hull with hyperplane x_{n+1} = 0.
        polynomial_solver::ConvexPolyhedron hyperplane_intersection;
        if (!polynomial_solver::intersect_convex_polyhedron_with_last_coordinate_zero(
                hull, hyperplane_intersection)) {
            return false;
        }

        // Project to R^n by dropping the last coordinate (which is 0).
        std::vector<std::vector<double>> projected_points;
        projected_points.reserve(hyperplane_intersection.vertices.size());
        for (const std::vector<double>& v : hyperplane_intersection.vertices) {
            std::vector<double> proj(dim, 0.0);
            for (std::size_t j = 0; j < dim; ++j) {
                proj[j] = v[j];
            }
            projected_points.push_back(proj);
        }

        // Recompute convex hull in R^n to ensure proper vertex ordering.
        polynomial_solver::ConvexPolyhedron projected = polynomial_solver::convex_hull(projected_points);

        hyperplane_intersections.push_back(projected);
    }

    // Step 3: Intersect all the projected polyhedra (now in R^n).
    polynomial_solver::ConvexPolyhedron intersection;
    if (!polynomial_solver::intersect_convex_polyhedra(hyperplane_intersections, intersection)) {
        return false;
    }

    // Compute axis-aligned bounding box (now in R^n).
    polynomial_solver::ConvexPolyhedronBox bbox =
        polynomial_solver::bounding_box(intersection);

    if (bbox.dimension() != dim) {
        return false;
    }

    // Extract bounds (already in parameter space R^n).
    local_bound_lower.resize(dim);
    local_bound_upper.resize(dim);

    for (std::size_t i = 0; i < dim; ++i) {
        local_bound_lower[i] = bbox.min_coords[i];
        local_bound_upper[i] = bbox.max_coords[i];

        // Clamp to [0, 1] (local parameter space).
        if (local_bound_lower[i] < 0.0) local_bound_lower[i] = 0.0;
        if (local_bound_upper[i] > 1.0) local_bound_upper[i] = 1.0;
        if (local_bound_lower[i] > 1.0) local_bound_lower[i] = 1.0;
        if (local_bound_upper[i] < 0.0) local_bound_upper[i] = 0.0;

        // Check for empty interval.
        if (local_bound_lower[i] > local_bound_upper[i]) {
            return false;
        }
    }

    return true;
}

} // namespace

std::vector<SubdivisionBoxResult>
Solver::subdivisionSolve(const PolynomialSystem& system,
                         const SubdivisionConfig& config,
                         RootBoundingMethod method) const
{
    std::vector<SubdivisionBoxResult> results;

    const std::size_t dim = system.dimension();
    if (dim == 0u) {
        return results;
    }

    // Root node: full box [0,1]^dim with the original equations.
    SubdivisionNode root;
    root.box_lower.assign(dim, 0.0);
    root.box_upper.assign(dim, 1.0);
    root.depth = 0u;
    root.polys = system.equations();

    std::priority_queue<NodeQueueEntry, std::vector<NodeQueueEntry>, NodeQueueCompare> queue;
    std::size_t next_index = 0u;
    queue.push(NodeQueueEntry{root.depth, next_index++, root});

    const double tolerance = config.tolerance;

    // Track boxes per depth for degeneracy detection
    std::map<unsigned int, unsigned int> boxes_per_depth;

    while (!queue.empty()) {
        NodeQueueEntry entry = queue.top();
        queue.pop();

        SubdivisionNode node = std::move(entry.node);

        // Track boxes at this depth
        boxes_per_depth[node.depth]++;

        // Check for degeneracy: too many boxes at this depth
        if (boxes_per_depth[node.depth] > config.max_boxes_per_depth) {
            // Degenerate case detected (multiplicity, infinite roots, etc.)
            // Return all results found so far with a warning
            SubdivisionBoxResult warning_box;
            warning_box.lower = node.box_lower;
            warning_box.upper = node.box_upper;
            warning_box.depth = node.depth;
            warning_box.converged = false;
            results.push_back(std::move(warning_box));

            // Add a marker result with depth = max_depth + 1 to signal degeneracy
            SubdivisionBoxResult degeneracy_marker;
            degeneracy_marker.lower.assign(dim, -1.0);  // Invalid box as marker
            degeneracy_marker.upper.assign(dim, -1.0);
            degeneracy_marker.depth = config.max_depth + 1u;
            degeneracy_marker.converged = false;
            results.push_back(std::move(degeneracy_marker));

            return results;
        }

        // Check depth limit
        if (node.depth >= config.max_depth) {
            SubdivisionBoxResult res;
            res.lower = node.box_lower;
            res.upper = node.box_upper;
            res.depth = node.depth;
            res.converged = false;
            results.push_back(std::move(res));
            continue;
        }

        // --- Iterative bounding step ---------------------------------------
        // Workflow: Compute bounding box, if empty quit, if small enough iterate,
        // else subdivide.

        bool converged = false;
        const unsigned int max_iterations = 100u;  // Prevent infinite loops

        for (unsigned int iter = 0; iter < max_iterations; ++iter) {
            // Step 1: Compute bounding box of roots in local [0,1]^n space
            std::vector<double> local_bound_lower(dim, 0.0);
            std::vector<double> local_bound_upper(dim, 1.0);

            bool has_roots = true;
            if (method == RootBoundingMethod::GraphHull) {
                // Use exact convex hull of graph control points
                if (!compute_graph_hull_bounds(node.polys, dim,
                                               local_bound_lower, local_bound_upper)) {
                    // Step 2: If empty, quit (no roots in this box)
                    has_roots = false;
                }
            }
            // For RootBoundingMethod::None, keep the default [0,1]^n bounds

            if (!has_roots) {
                // No roots in this box; discard it
                break;
            }

            // Step 3: Check if bounding box is small enough
            bool all_small = true;
            for (std::size_t i = 0; i < dim; ++i) {
                const double local_width = local_bound_upper[i] - local_bound_lower[i];
                const double old_low = node.box_lower[i];
                const double old_high = node.box_upper[i];
                const double old_width = old_high - old_low;
                const double global_width = local_width * old_width;

                if (global_width > tolerance) {
                    all_small = false;
                    break;
                }
            }

            if (all_small) {
                // Step 4: Box is small enough, converged
                converged = true;
                break;
            }

            // Step 5: Box not small enough, contract and iterate
            // Make bounding box the new region
            for (std::size_t i = 0; i < dim; ++i) {
                const double a = local_bound_lower[i];
                const double b = local_bound_upper[i];

                const double old_low = node.box_lower[i];
                const double old_high = node.box_upper[i];
                const double old_width = old_high - old_low;

                const double new_low = old_low + a * old_width;
                const double new_high = old_low + b * old_width;

                node.box_lower[i] = new_low;
                node.box_upper[i] = new_high;

                // Restrict polynomials to new interval
                for (std::size_t eq = 0; eq < node.polys.size(); ++eq) {
                    node.polys[eq] = node.polys[eq].restrictedToInterval(i, a, b);
                }
            }

            // Continue iteration (go back to step 1)
        }

        if (converged) {
            // Box converged, add to results
            SubdivisionBoxResult res;
            res.lower = node.box_lower;
            res.upper = node.box_upper;
            res.depth = node.depth;
            res.converged = true;
            results.push_back(std::move(res));
            continue;
        }

        // --- Subdivision step ----------------------------------------------
        // Subdivide along axes that are not small enough.
        // Only subdivide dimensions whose width is larger than tolerance.
        std::vector<bool> split_dim(dim, false);
        for (std::size_t i = 0; i < dim; ++i) {
            const double width = node.box_upper[i] - node.box_lower[i];
            if (width > tolerance) {
                split_dim[i] = true;
            }
        }

        // Seed for children: start from the (possibly contracted) node.
        std::vector<SubdivisionNode> children;
        SubdivisionNode seed = node;
        seed.depth = node.depth + 1u;
        children.push_back(seed);

        for (std::size_t axis = 0; axis < dim; ++axis) {
            if (!split_dim[axis]) {
                continue;
            }

            std::vector<SubdivisionNode> next_children;
            next_children.reserve(children.size() * 2u);

            for (const SubdivisionNode& child : children) {
                const double old_low = child.box_lower[axis];
                const double old_high = child.box_upper[axis];
                const double mid = 0.5 * (old_low + old_high);

                SubdivisionNode left_child = child;
                SubdivisionNode right_child = child;

                left_child.box_upper[axis] = mid;
                right_child.box_lower[axis] = mid;

                left_child.polys.clear();
                right_child.polys.clear();
                left_child.polys.reserve(child.polys.size());
                right_child.polys.reserve(child.polys.size());

                for (const Polynomial& poly : child.polys) {
                    left_child.polys.push_back(poly.restrictedToInterval(axis, 0.0, 0.5));
                    right_child.polys.push_back(poly.restrictedToInterval(axis, 0.5, 1.0));
                }

                next_children.push_back(std::move(left_child));
                next_children.push_back(std::move(right_child));
            }

            children.swap(next_children);
        }

        for (SubdivisionNode& child : children) {
            queue.push(NodeQueueEntry{child.depth, next_index++, std::move(child)});
        }
    }

    return results;
}

} // namespace polynomial_solver

