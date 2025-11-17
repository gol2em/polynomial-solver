#include "solver.h"
#include <queue>
#include <map>
#include <iostream>
#include <limits>


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

void PolynomialSystem::evaluate(const std::vector<double>& point, std::vector<double>& values) const {
    values.resize(equations_.size());
    for (std::size_t i = 0; i < equations_.size(); ++i) {
        values[i] = equations_[i].evaluate(point);
    }
}

bool PolynomialSystem::isApproximateRoot(const std::vector<double>& point, double tolerance) const {
    if (point.size() != dimension_) {
        return false;
    }

    for (const Polynomial& eq : equations_) {
        const double value = eq.evaluate(point);
        const double abs_value = (value >= 0.0) ? value : -value;
        if (abs_value > tolerance) {
            return false;
        }
    }

    return true;
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

// Helper: Analyze bounding box to detect degenerate cases
// Returns:
//   0: Empty box (should not happen if has_roots is true)
//   1: Single point (all dimensions < tolerance)
//   2: Line segment (exactly one dimension >= tolerance)
//   3: Higher dimensional box (multiple dimensions >= tolerance)
int analyze_box_dimension(const std::vector<double>& lower,
                          const std::vector<double>& upper,
                          double tolerance,
                          std::size_t& active_axis)
{
    const std::size_t dim = lower.size();
    std::size_t num_active = 0;
    active_axis = 0;

    for (std::size_t i = 0; i < dim; ++i) {
        const double width = upper[i] - lower[i];
        if (width > tolerance) {
            num_active++;
            active_axis = i;  // Remember the last active axis
        }
    }

    if (num_active == 0) {
        return 1;  // Single point
    } else if (num_active == 1) {
        return 2;  // Line segment (1D problem)
    } else {
        return 3;  // Higher dimensional box
    }
}

/**
 * @brief Compute root bounding box using the projected polyhedral method.
 *
 * This function implements the following workflow for each direction i:
 * 1. For each equation, project all graph control points to 2D (coordinate i + function value).
 * 2. For each equation, compute convex hull of these 2D points.
 * 3. For each equation, intersect with horizontal axis (function value = 0).
 * 4. For each equation, project to 1D to get an interval.
 * 5. Intersect all intervals from all equations to get the bound in direction i.
 *
 * Returns true if a non-empty bounding box is found, false otherwise.
 */
bool compute_projected_polyhedral_bounds(
    const std::vector<polynomial_solver::Polynomial>& polys,
    std::size_t dim,
    std::vector<double>& local_bound_lower,
    std::vector<double>& local_bound_upper)
{
    if (polys.empty() || dim == 0u) {
        return false;
    }

    local_bound_lower.resize(dim);
    local_bound_upper.resize(dim);

    // Process each direction independently
    for (std::size_t dir = 0; dir < dim; ++dir) {
        // For this direction, we'll compute the intersection of intervals from all equations
        double dir_min = 0.0;
        double dir_max = 1.0;

        for (const polynomial_solver::Polynomial& poly : polys) {
            // Get graph control points in R^{n+1}
            std::vector<double> control_points;
            poly.graphControlPoints(control_points);

            const std::size_t num_coeffs = poly.coefficientCount();
            const std::size_t point_dim = dim + 1u;

            // Project to 2D: keep coordinate 'dir' and the last coordinate (function value)
            std::vector<std::vector<double>> projected_2d;
            projected_2d.reserve(num_coeffs);

            for (std::size_t i = 0; i < num_coeffs; ++i) {
                std::vector<double> pt_2d(2);
                pt_2d[0] = control_points[i * point_dim + dir];  // coordinate in direction 'dir'
                pt_2d[1] = control_points[i * point_dim + dim];  // function value (last coordinate)
                projected_2d.push_back(pt_2d);
            }

            // Compute convex hull in 2D
            polynomial_solver::ConvexPolyhedron hull_2d = polynomial_solver::convex_hull(projected_2d);

            // Intersect with horizontal axis (y = 0, i.e., last coordinate = 0)
            polynomial_solver::ConvexPolyhedron intersection_1d;
            if (!polynomial_solver::intersect_convex_polyhedron_with_last_coordinate_zero(
                    hull_2d, intersection_1d)) {
                // No intersection with axis for this equation means no roots
                return false;
            }

            // Project to 1D by taking the first coordinate (drop the second which is 0)
            // Find min and max of the first coordinate
            if (intersection_1d.vertices.empty()) {
                return false;
            }

            double eq_min = intersection_1d.vertices[0][0];
            double eq_max = intersection_1d.vertices[0][0];

            for (const std::vector<double>& v : intersection_1d.vertices) {
                if (v[0] < eq_min) eq_min = v[0];
                if (v[0] > eq_max) eq_max = v[0];
            }

            // Intersect with current bounds for this direction
            if (eq_min > dir_min) dir_min = eq_min;
            if (eq_max < dir_max) dir_max = eq_max;

            // Check if intersection is empty
            if (dir_min > dir_max) {
                return false;
            }
        }

        // Store the bounds for this direction
        local_bound_lower[dir] = dir_min;
        local_bound_upper[dir] = dir_max;

        // Clamp to [0, 1] (local parameter space)
        if (local_bound_lower[dir] < 0.0) local_bound_lower[dir] = 0.0;
        if (local_bound_upper[dir] > 1.0) local_bound_upper[dir] = 1.0;
        if (local_bound_lower[dir] > 1.0) local_bound_lower[dir] = 1.0;
        if (local_bound_upper[dir] < 0.0) local_bound_upper[dir] = 0.0;

        // Check for empty interval
        if (local_bound_lower[dir] > local_bound_upper[dir]) {
            return false;
        }
    }

    return true;
}

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

SubdivisionSolverResult
Solver::subdivisionSolve(const PolynomialSystem& system,
                         const SubdivisionConfig& config,
                         RootBoundingMethod method) const
{
    SubdivisionSolverResult result;
    result.num_resolved = 0;
    result.degeneracy_detected = false;

    std::vector<SubdivisionBoxResult> resolved_boxes;
    std::vector<SubdivisionBoxResult> unresolved_boxes;

    const std::size_t dim = system.dimension();
    if (dim == 0u) {
        return result;
    }

    // Compute expected maximum number of roots: product of degrees
    std::size_t expected_max_roots = 1;
    for (const Polynomial& eq : system.equations()) {
        const std::vector<unsigned int>& degrees = eq.degrees();
        for (unsigned int d : degrees) {
            expected_max_roots *= d;
        }
    }
    const std::size_t degeneracy_threshold =
        static_cast<std::size_t>(config.degeneracy_multiplier * static_cast<double>(expected_max_roots));

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
    std::map<unsigned int, std::size_t> boxes_per_depth;

    while (!queue.empty()) {
        NodeQueueEntry entry = queue.top();
        queue.pop();

        SubdivisionNode node = std::move(entry.node);

        // Track boxes at this depth
        boxes_per_depth[node.depth]++;

        // Check for degeneracy: too many boxes at this depth compared to expected root count
        if (boxes_per_depth[node.depth] > degeneracy_threshold) {
            // Degenerate case detected (multiplicity, infinite roots, etc.)
            result.degeneracy_detected = true;

            // Add current node to unresolved
            SubdivisionBoxResult box;
            box.lower = node.box_lower;
            box.upper = node.box_upper;
            box.center.resize(dim);
            box.max_error.resize(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                box.center[i] = 0.5 * (node.box_lower[i] + node.box_upper[i]);
                const double half_width = 0.5 * (node.box_upper[i] - node.box_lower[i]);
                box.max_error[i] = (half_width < std::numeric_limits<double>::epsilon())
                    ? std::numeric_limits<double>::epsilon()
                    : half_width;
            }
            box.depth = node.depth;
            box.converged = false;
            unresolved_boxes.push_back(std::move(box));

            // Collect all remaining boxes in queue as unresolved
            while (!queue.empty()) {
                NodeQueueEntry remaining_entry = queue.top();
                queue.pop();
                SubdivisionNode& remaining_node = remaining_entry.node;

                SubdivisionBoxResult remaining_box;
                remaining_box.lower = remaining_node.box_lower;
                remaining_box.upper = remaining_node.box_upper;
                remaining_box.center.resize(dim);
                remaining_box.max_error.resize(dim);
                for (std::size_t i = 0; i < dim; ++i) {
                    remaining_box.center[i] = 0.5 * (remaining_node.box_lower[i] + remaining_node.box_upper[i]);
                    const double half_width = 0.5 * (remaining_node.box_upper[i] - remaining_node.box_lower[i]);
                    remaining_box.max_error[i] = (half_width < std::numeric_limits<double>::epsilon())
                        ? std::numeric_limits<double>::epsilon()
                        : half_width;
                }
                remaining_box.depth = remaining_node.depth;
                remaining_box.converged = false;
                unresolved_boxes.push_back(std::move(remaining_box));
            }

            // Stop processing
            break;
        }

        // Check depth limit (should rarely happen - indicates problem)
        if (node.depth >= config.max_depth) {
            std::cerr << "Warning: Maximum depth " << config.max_depth
                      << " reached. This may indicate numerical issues or insufficient tolerance."
                      << std::endl;

            SubdivisionBoxResult box;
            box.lower = node.box_lower;
            box.upper = node.box_upper;
            box.center.resize(dim);
            box.max_error.resize(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                box.center[i] = 0.5 * (node.box_lower[i] + node.box_upper[i]);
                const double half_width = 0.5 * (node.box_upper[i] - node.box_lower[i]);
                box.max_error[i] = (half_width < std::numeric_limits<double>::epsilon())
                    ? std::numeric_limits<double>::epsilon()
                    : half_width;
            }
            box.depth = node.depth;
            box.converged = false;
            unresolved_boxes.push_back(std::move(box));
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
            } else if (method == RootBoundingMethod::ProjectedPolyhedral) {
                // Use projected polyhedral method (direction-by-direction)
                if (!compute_projected_polyhedral_bounds(node.polys, dim,
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

            // Step 3: Contract the box to the bounding box
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

            // Step 4: Check if contracted box is small enough
            bool all_small = true;
            for (std::size_t i = 0; i < dim; ++i) {
                const double width = node.box_upper[i] - node.box_lower[i];
                if (width > tolerance) {
                    all_small = false;
                    break;
                }
            }

            if (all_small) {
                // Step 5: Box is small enough, converged
                converged = true;
                break;
            }

            // Continue iteration (go back to step 1)
        }

        if (converged) {
            // Box converged
            // Only do degenerate box handling when using GraphHull or ProjectedPolyhedral method
            if (method == RootBoundingMethod::GraphHull || method == RootBoundingMethod::ProjectedPolyhedral) {
                std::size_t active_axis = 0;
                int box_dim = analyze_box_dimension(node.box_lower, node.box_upper,
                                                    tolerance, active_axis);

                if (box_dim == 1) {
                    // Single point: check if it's actually a root
                    std::vector<double> center(dim);
                    for (std::size_t i = 0; i < dim; ++i) {
                        center[i] = 0.5 * (node.box_lower[i] + node.box_upper[i]);
                    }

                    // Create a PolynomialSystem from the current polynomials
                    PolynomialSystem local_system(node.polys);

                    // Check if center is approximately a root
                    // Use a tolerance based on the box size
                    double max_box_width = 0.0;
                    for (std::size_t i = 0; i < dim; ++i) {
                        const double w = node.box_upper[i] - node.box_lower[i];
                        if (w > max_box_width) {
                            max_box_width = w;
                        }
                    }
                    const double root_tolerance = std::max(1e-6, max_box_width);

                    if (local_system.isApproximateRoot(center, root_tolerance)) {
                        // It's a root, add to results
                        SubdivisionBoxResult box;
                        box.lower = node.box_lower;
                        box.upper = node.box_upper;
                        box.center = center;
                        box.max_error.resize(dim);
                        for (std::size_t i = 0; i < dim; ++i) {
                            const double half_width = 0.5 * (node.box_upper[i] - node.box_lower[i]);
                            box.max_error[i] = (half_width < std::numeric_limits<double>::epsilon())
                                ? std::numeric_limits<double>::epsilon()
                                : half_width;
                        }
                        box.depth = node.depth;
                        box.converged = true;
                        resolved_boxes.push_back(std::move(box));
                    } else {
                        // Not a root, discard (false positive from bounding)
                        // This can happen due to numerical errors or conservative bounding
                    }
                    continue;

                } else if (box_dim == 2) {
                    // Line segment (1D degenerate case)
                    // The problem is now 1D along active_axis
                    // For now, mark as converged and let user handle
                    // TODO: Implement 1D subdivision along the active axis
                    SubdivisionBoxResult box;
                    box.lower = node.box_lower;
                    box.upper = node.box_upper;
                    box.center.resize(dim);
                    box.max_error.resize(dim);
                    for (std::size_t i = 0; i < dim; ++i) {
                        box.center[i] = 0.5 * (node.box_lower[i] + node.box_upper[i]);
                        const double half_width = 0.5 * (node.box_upper[i] - node.box_lower[i]);
                        box.max_error[i] = (half_width < std::numeric_limits<double>::epsilon())
                            ? std::numeric_limits<double>::epsilon()
                            : half_width;
                    }
                    box.depth = node.depth;
                    box.converged = true;
                    resolved_boxes.push_back(std::move(box));
                    continue;
                }
                // Fall through to normal case for box_dim == 3
            }

            // Normal case: add converged box to results
            SubdivisionBoxResult box;
            box.lower = node.box_lower;
            box.upper = node.box_upper;
            box.center.resize(dim);
            box.max_error.resize(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                box.center[i] = 0.5 * (node.box_lower[i] + node.box_upper[i]);
                const double half_width = 0.5 * (node.box_upper[i] - node.box_lower[i]);
                box.max_error[i] = (half_width < std::numeric_limits<double>::epsilon())
                    ? std::numeric_limits<double>::epsilon()
                    : half_width;
            }
            box.depth = node.depth;
            box.converged = true;
            resolved_boxes.push_back(std::move(box));
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

    // Assemble final result: resolved boxes first, then unresolved
    result.boxes.reserve(resolved_boxes.size() + unresolved_boxes.size());
    result.boxes.insert(result.boxes.end(), resolved_boxes.begin(), resolved_boxes.end());
    result.boxes.insert(result.boxes.end(), unresolved_boxes.begin(), unresolved_boxes.end());
    result.num_resolved = resolved_boxes.size();

    return result;
}

} // namespace polynomial_solver

