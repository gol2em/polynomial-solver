#include "solver.h"
#include <queue>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <sys/stat.h>
#include <sys/types.h>


/**
 * @file solver.cpp
 * @brief Implementation of the main solver interface
 */

namespace polynomial_solver {

// Helper function to create directory if it doesn't exist
static void ensure_directory_exists(const std::string& path) {
    std::size_t pos = path.find_last_of("/\\");
    if (pos != std::string::npos) {
        std::string dir = path.substr(0, pos);
        // Try to create directory (ignore errors if it already exists)
        #ifdef _WIN32
            _mkdir(dir.c_str());
        #else
            mkdir(dir.c_str(), 0755);
        #endif
    }
}

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

// Forward declaration for internal helper function
// This function is always available, but dump I/O is only compiled when ENABLE_GEOMETRY_DUMP is defined
bool compute_projected_polyhedral_bounds_with_dump(
    const std::vector<polynomial_solver::Polynomial>& polys,
    std::size_t dim,
    std::vector<double>& local_bound_lower,
    std::vector<double>& local_bound_upper,
    const std::string& dump_file,
    const std::vector<double>& global_box_lower,
    const std::vector<double>& global_box_upper,
    unsigned int depth,
    unsigned int iteration,
    const std::string& decision);

struct SubdivisionNode {
    std::vector<double> box_lower;  ///< Global box lower corner in [0,1]^n.
    std::vector<double> box_upper;  ///< Global box upper corner in [0,1]^n.
    unsigned int depth;             ///< Subdivision depth.

    // Polynomials restricted to this node's box and re-parameterised to [0,1]^n.
    std::vector<polynomial_solver::Polynomial> polys;

    // Original polynomials (on [0,1]^n, never modified) for direct contraction.
    std::vector<polynomial_solver::Polynomial> original_polys;
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
    // Forward to the dump version with empty dump file
    // This ensures we maintain a single implementation
    // When ENABLE_GEOMETRY_DUMP is disabled, the dump I/O code is compiled out
    std::vector<double> dummy_box_lower(dim, 0.0);
    std::vector<double> dummy_box_upper(dim, 1.0);
    return compute_projected_polyhedral_bounds_with_dump(
        polys, dim, local_bound_lower, local_bound_upper,
        "",  // empty dump file means no dumping
        dummy_box_lower, dummy_box_upper,
        0, 0, "CONTRACT");
}

/**
 * @brief Compute root bounding box using projected polyhedral method with optional geometry dump.
 *
 * Same as compute_projected_polyhedral_bounds but with optional dump of geometric data.
 * If dump_file is not empty AND ENABLE_GEOMETRY_DUMP is defined, writes:
 * - Projected 2D points for each direction and equation
 * - Convex hull vertices
 * - Intersection with axis
 * - Final bounding box
 *
 * When ENABLE_GEOMETRY_DUMP is not defined, all dump I/O code is removed at compile time.
 */
bool compute_projected_polyhedral_bounds_with_dump(
    const std::vector<polynomial_solver::Polynomial>& polys,
    std::size_t dim,
    std::vector<double>& local_bound_lower,
    std::vector<double>& local_bound_upper,
    const std::string& dump_file,
    const std::vector<double>& global_box_lower,
    const std::vector<double>& global_box_upper,
    unsigned int depth,
    unsigned int iteration,
    const std::string& decision = "CONTRACT")
{
    if (polys.empty() || dim == 0u) {
        return false;
    }

#ifdef ENABLE_GEOMETRY_DUMP
    // Dump I/O setup (only compiled when dump is enabled)
    std::ofstream dump;
    bool do_dump = !dump_file.empty();
    if (do_dump) {
        dump.open(dump_file.c_str(), std::ios::app);
        if (!dump.is_open()) {
            std::cerr << "Warning: Failed to open dump file: " << dump_file << std::endl;
            do_dump = false;
        } else {
            dump << std::setprecision(16);
            dump << "# Iteration " << iteration << ", Depth " << depth << "\n";
            dump << "# Decision: " << decision << "\n";
            dump << "# Global box: [";
            for (std::size_t i = 0; i < dim; ++i) {
                if (i > 0) dump << ", ";
                dump << global_box_lower[i] << ", " << global_box_upper[i];
            }
            dump << "]\n";
        }
    }
#endif

    local_bound_lower.resize(dim);
    local_bound_upper.resize(dim);

    // Process each direction independently
    for (std::size_t dir = 0; dir < dim; ++dir) {
#ifdef ENABLE_GEOMETRY_DUMP
        if (do_dump) {
            dump << "\n## Direction " << dir << "\n";
        }
#endif

        // For this direction, we'll compute the intersection of intervals from all equations
        double dir_min = 0.0;
        double dir_max = 1.0;

        for (std::size_t eq_idx = 0; eq_idx < polys.size(); ++eq_idx) {
            const polynomial_solver::Polynomial& poly = polys[eq_idx];

#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                dump << "\n### Equation " << eq_idx << "\n";
            }
#endif

            // Get graph control points in R^{n+1}
            std::vector<double> control_points;
            poly.graphControlPoints(control_points);

            const std::size_t num_coeffs = poly.coefficientCount();
            const std::size_t point_dim = dim + 1u;

#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                // Dump 3D control points (for 2D problems: x, y, f(x,y))
                dump << "Control_Points_3D " << num_coeffs << "\n";
                for (std::size_t i = 0; i < num_coeffs; ++i) {
                    for (std::size_t j = 0; j < point_dim; ++j) {
                        if (j > 0) dump << " ";
                        dump << control_points[i * point_dim + j];
                    }
                    dump << "\n";
                }
            }
#endif

            // Project to 2D: keep coordinate 'dir' and the last coordinate (function value)
            // Optimization: Control points are structured in vertical groups when projected.
            // For convex hull, we only need the extreme (min/max) points in each group.
            // This reduces the point set from O(d₀·d₁) to O(d_dir) where d_dir is the degree
            // in the projection direction.

            std::vector<std::vector<double>> projected_2d_all;
            projected_2d_all.reserve(num_coeffs);

#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                dump << "Projected_Points " << num_coeffs << "\n";
            }
#endif

            // First pass: collect all projected points and group by x-coordinate
            std::map<double, std::vector<double>> vertical_groups;

            for (std::size_t i = 0; i < num_coeffs; ++i) {
                double x = control_points[i * point_dim + dir];  // coordinate in direction 'dir'
                double y = control_points[i * point_dim + dim];  // function value (last coordinate)

                vertical_groups[x].push_back(y);

#ifdef ENABLE_GEOMETRY_DUMP
                if (do_dump) {
                    dump << x << " " << y << "\n";
                }
#endif
            }

            // Second pass: extract only extreme points (min/max) from each vertical group
            std::vector<std::vector<double>> projected_2d;
            projected_2d.reserve(2 * vertical_groups.size());

            for (std::map<double, std::vector<double>>::const_iterator it = vertical_groups.begin();
                 it != vertical_groups.end(); ++it) {
                const double x = it->first;
                const std::vector<double>& y_values = it->second;

                double y_min = *std::min_element(y_values.begin(), y_values.end());
                double y_max = *std::max_element(y_values.begin(), y_values.end());

                std::vector<double> pt_min(2);
                pt_min[0] = x;
                pt_min[1] = y_min;
                projected_2d.push_back(pt_min);

                // Only add max if it's different from min
                if (y_min != y_max) {
                    std::vector<double> pt_max(2);
                    pt_max[0] = x;
                    pt_max[1] = y_max;
                    projected_2d.push_back(pt_max);
                }
            }

#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                dump << "Reduced_Points " << projected_2d.size()
                     << " (from " << num_coeffs << " control points)\n";
            }
#endif

            // Compute convex hull in 2D (now with reduced point set)
            polynomial_solver::ConvexPolyhedron hull_2d = polynomial_solver::convex_hull(projected_2d);

#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                dump << "ConvexHull " << hull_2d.vertices.size() << "\n";
                for (const std::vector<double>& v : hull_2d.vertices) {
                    dump << v[0] << " " << v[1] << "\n";
                }
            }
#endif

            // Intersect with horizontal axis (y = 0, i.e., last coordinate = 0)
            polynomial_solver::ConvexPolyhedron intersection_1d;
            if (!polynomial_solver::intersect_convex_polyhedron_with_last_coordinate_zero(
                    hull_2d, intersection_1d)) {
                // No intersection with axis for this equation means no roots
#ifdef ENABLE_GEOMETRY_DUMP
                if (do_dump) {
                    dump << "Intersection EMPTY\n";
                    dump.close();
                }
#endif
                return false;
            }

#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                dump << "Intersection " << intersection_1d.vertices.size() << "\n";
                for (const std::vector<double>& v : intersection_1d.vertices) {
                    dump << v[0] << " " << v[1] << "\n";
                }
            }
#endif

            // Project to 1D by taking the first coordinate (drop the second which is 0)
            // Find min and max of the first coordinate
            if (intersection_1d.vertices.empty()) {
#ifdef ENABLE_GEOMETRY_DUMP
                if (do_dump) {
                    dump.close();
                }
#endif
                return false;
            }

            double eq_min = intersection_1d.vertices[0][0];
            double eq_max = intersection_1d.vertices[0][0];

            for (const std::vector<double>& v : intersection_1d.vertices) {
                if (v[0] < eq_min) eq_min = v[0];
                if (v[0] > eq_max) eq_max = v[0];
            }

#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                dump << "Interval [" << eq_min << ", " << eq_max << "]\n";
            }
#endif

            // Intersect with current bounds for this direction
            if (eq_min > dir_min) dir_min = eq_min;
            if (eq_max < dir_max) dir_max = eq_max;

            // Check if intersection is empty
            if (dir_min > dir_max) {
#ifdef ENABLE_GEOMETRY_DUMP
                if (do_dump) {
                    dump << "Direction_Interval EMPTY\n";
                    dump.close();
                }
#endif
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
#ifdef ENABLE_GEOMETRY_DUMP
            if (do_dump) {
                dump << "Final_Interval EMPTY\n";
                dump.close();
            }
#endif
            return false;
        }

#ifdef ENABLE_GEOMETRY_DUMP
        if (do_dump) {
            dump << "Final_Interval [" << local_bound_lower[dir] << ", " << local_bound_upper[dir] << "]\n";
        }
#endif
    }

#ifdef ENABLE_GEOMETRY_DUMP
    if (do_dump) {
        dump << "\n# Bounding Box (local): [";
        for (std::size_t i = 0; i < dim; ++i) {
            if (i > 0) dump << ", ";
            dump << local_bound_lower[i] << ", " << local_bound_upper[i];
        }
        dump << "]\n";

        // Compute global bounding box
        dump << "# Bounding Box (global): [";
        for (std::size_t i = 0; i < dim; ++i) {
            if (i > 0) dump << ", ";
            double global_low = global_box_lower[i] + local_bound_lower[i] * (global_box_upper[i] - global_box_lower[i]);
            double global_high = global_box_lower[i] + local_bound_upper[i] * (global_box_upper[i] - global_box_lower[i]);
            dump << global_low << ", " << global_high;
        }
        dump << "]\n\n";

        dump.close();
    }
#endif

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
    root.original_polys = system.equations();  // Store original for direct contraction

    std::priority_queue<NodeQueueEntry, std::vector<NodeQueueEntry>, NodeQueueCompare> queue;
    std::size_t next_index = 0u;
    queue.push(NodeQueueEntry{root.depth, next_index++, root});

    const double tolerance = config.tolerance;

    // Setup dump file if requested
#ifdef ENABLE_GEOMETRY_DUMP
    std::string dump_file;
    unsigned int dump_iteration = 0;
    if (config.dump_geometry) {
        dump_file = config.dump_prefix + "_geometry.txt";
        // Ensure directory exists
        ensure_directory_exists(dump_file);
        // Clear the file
        std::ofstream clear_dump(dump_file.c_str());
        if (clear_dump.is_open()) {
            clear_dump << "# Polynomial Solver Geometry Dump\n";
            clear_dump << "# Method: ";
            if (method == RootBoundingMethod::ProjectedPolyhedral) {
                clear_dump << "ProjectedPolyhedral\n";
            } else if (method == RootBoundingMethod::GraphHull) {
                clear_dump << "GraphHull\n";
            } else {
                clear_dump << "None\n";
            }
            clear_dump << "# Dimension: " << dim << "\n";
            clear_dump << "# Number of equations: " << system.equationCount() << "\n";
            clear_dump << "\n";
            clear_dump.close();
        }
    }
#endif

    // Track boxes that need subdivision per depth for degeneracy detection
    // Only count boxes that would be subdivided, not pruned or converged boxes
    std::map<unsigned int, std::size_t> subdivision_boxes_per_depth;
    bool degeneracy_mode = false;  // Flag to indicate we're in degeneracy mode

    while (!queue.empty()) {
        NodeQueueEntry entry = queue.top();
        queue.pop();

        SubdivisionNode node = std::move(entry.node);

        // Check depth limit (should rarely happen - indicates problem)
        if (node.depth >= config.max_depth) {
            std::cerr << "Warning: Maximum depth " << config.max_depth
                      << " reached. This may indicate numerical issues or insufficient tolerance."
                      << std::endl;

#ifdef ENABLE_GEOMETRY_DUMP
            if (config.dump_geometry && method == RootBoundingMethod::ProjectedPolyhedral) {
                std::ofstream dump(dump_file.c_str(), std::ios::app);
                if (dump.is_open()) {
                    dump << "# Iteration " << dump_iteration++ << ", Depth " << node.depth << "\n";
                    dump << "# Decision: MAX_DEPTH_REACHED\n";
                    dump << "# Global box: [";
                    for (std::size_t i = 0; i < dim; ++i) {
                        if (i > 0) dump << ", ";
                        dump << node.box_lower[i] << ", " << node.box_upper[i];
                    }
                    dump << "]\n";
                    dump << "# FINAL_DECISION: UNRESOLVED (max depth reached)\n\n";
                    dump.close();
                }
            }
#endif

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
        std::vector<bool> final_needs_subdivision(dim, false);  // Track which dimensions need subdivision
        bool box_was_pruned = false;  // Track if box was pruned (to skip subdivision)

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
#ifdef ENABLE_GEOMETRY_DUMP
                if (config.dump_geometry) {
                    // Determine decision based on strategy
                    std::ostringstream decision_stream;
                    decision_stream << "BOUNDING (";
                    if (config.strategy == SubdivisionStrategy::ContractFirst) {
                        decision_stream << "ContractFirst";
                    } else if (config.strategy == SubdivisionStrategy::SubdivideFirst) {
                        decision_stream << "SubdivideFirst";
                    } else if (config.strategy == SubdivisionStrategy::Simultaneous) {
                        decision_stream << "Simultaneous";
                    }
                    decision_stream << ", iter " << iter << ")";

                    // Use dump version
                    if (!compute_projected_polyhedral_bounds_with_dump(
                            node.polys, dim,
                            local_bound_lower, local_bound_upper,
                            dump_file,
                            node.box_lower, node.box_upper,
                            node.depth, dump_iteration++,
                            decision_stream.str())) {
                        // Step 2: If empty, quit (no roots in this box)
                        has_roots = false;
                    }
                } else
#endif
                {
                    // Use regular version
                    if (!compute_projected_polyhedral_bounds(node.polys, dim,
                                                             local_bound_lower, local_bound_upper)) {
                        // Step 2: If empty, quit (no roots in this box)
                        has_roots = false;
                    }
                }
            }
            // For RootBoundingMethod::None, keep the default [0,1]^n bounds

            if (!has_roots) {
                // No roots in this box; discard it
#ifdef ENABLE_GEOMETRY_DUMP
                if (config.dump_geometry && method == RootBoundingMethod::ProjectedPolyhedral) {
                    std::ofstream dump(dump_file.c_str(), std::ios::app);
                    if (dump.is_open()) {
                        dump << "# FINAL_DECISION: PRUNED (empty bounding box, no roots)\n\n";
                        dump.close();
                    }
                }
#endif
                box_was_pruned = true;
                break;
            }

            // Step 3: Analyze contraction per direction
            std::vector<double> contraction_ratio(dim);
            std::vector<bool> needs_subdivision(dim, false);
            bool any_needs_subdivision = false;

            for (std::size_t i = 0; i < dim; ++i) {
                const double old_width = node.box_upper[i] - node.box_lower[i];
                const double new_width = (local_bound_upper[i] - local_bound_lower[i]) * old_width;
                contraction_ratio[i] = (old_width > 0.0) ? (new_width / old_width) : 0.0;

                // Check if this direction contracted enough
                if (contraction_ratio[i] > config.contraction_threshold && old_width > tolerance) {
                    needs_subdivision[i] = true;
                    any_needs_subdivision = true;
                }
            }

#ifdef ENABLE_GEOMETRY_DUMP
            // Dump contraction analysis
            if (config.dump_geometry && method == RootBoundingMethod::ProjectedPolyhedral) {
                std::ofstream dump(dump_file.c_str(), std::ios::app);
                if (dump.is_open()) {
                    dump << std::setprecision(6);
                    dump << "# Contraction Analysis:\n";
                    for (std::size_t i = 0; i < dim; ++i) {
                        dump << "#   Direction " << i << ": ratio=" << contraction_ratio[i]
                             << " (threshold=" << config.contraction_threshold << ")";
                        if (needs_subdivision[i]) {
                            dump << " -> NEEDS_SUBDIVISION";
                        } else {
                            dump << " -> OK";
                        }
                        dump << "\n";
                    }
                    dump.close();
                }
            }
#endif

            // Step 4: Apply strategy-specific logic
            bool should_subdivide_now = false;
            bool should_contract = true;

            if (config.strategy == SubdivisionStrategy::SubdivideFirst) {
                // SubdivideFirst: If any direction didn't contract enough, subdivide immediately
                if (any_needs_subdivision) {
                    should_subdivide_now = true;
                    should_contract = false;
                }
            } else if (config.strategy == SubdivisionStrategy::Simultaneous) {
                // Simultaneous: Subdivide in non-contracting directions, contract in others
                if (any_needs_subdivision) {
                    should_subdivide_now = true;
                    should_contract = true;  // Will contract only in directions that don't need subdivision
                }
            } else if (config.strategy == SubdivisionStrategy::ContractFirst) {
                // ContractFirst: Contract first, but subdivide if contraction doesn't make progress
                // This prevents wasting iterations when the box contains multiple roots
                if (any_needs_subdivision) {
                    should_subdivide_now = true;
                    should_contract = false;
                }
                // Otherwise, keep contracting (default behavior)
            }

            // Step 5: Contract the box (if strategy allows)
            if (should_contract) {
                // First, update global box coordinates
                for (std::size_t i = 0; i < dim; ++i) {
                    // For Simultaneous strategy, only contract directions that don't need subdivision
                    if (config.strategy == SubdivisionStrategy::Simultaneous && needs_subdivision[i]) {
                        continue;  // Skip contraction for this direction
                    }

                    const double a = local_bound_lower[i];
                    const double b = local_bound_upper[i];

                    const double old_low = node.box_lower[i];
                    const double old_high = node.box_upper[i];
                    const double old_width = old_high - old_low;

                    const double new_low = old_low + a * old_width;
                    const double new_high = old_low + b * old_width;

                    node.box_lower[i] = new_low;
                    node.box_upper[i] = new_high;
                }

                // Direct contraction: Recompute polynomials from original using global box
                // This eliminates error accumulation from repeated restrictions
                for (std::size_t eq = 0; eq < node.polys.size(); ++eq) {
                    node.polys[eq] = node.original_polys[eq];
                    for (std::size_t i = 0; i < dim; ++i) {
                        // Restrict from original [0,1] to global [box_lower[i], box_upper[i]]
                        node.polys[eq] = node.polys[eq].restrictedToInterval(
                            i, node.box_lower[i], node.box_upper[i]);
                    }
                }
            }

            // Step 6: If strategy says subdivide now, break out of iteration loop
            if (should_subdivide_now) {
                // Save which dimensions need subdivision for later
                final_needs_subdivision = needs_subdivision;
                // Will subdivide after the iteration loop
                break;
            }

            // Step 7: Check if contracted box is small enough
            bool all_small = true;
            for (std::size_t i = 0; i < dim; ++i) {
                const double width = node.box_upper[i] - node.box_lower[i];
                if (width > tolerance) {
                    all_small = false;
                    break;
                }
            }

            if (all_small) {
                // Step 8: Box is small enough, converged
#ifdef ENABLE_GEOMETRY_DUMP
                if (config.dump_geometry && method == RootBoundingMethod::ProjectedPolyhedral) {
                    std::ofstream dump(dump_file.c_str(), std::ios::app);
                    if (dump.is_open()) {
                        dump << "# FINAL_DECISION: CONVERGED (iter=" << iter << ")\n\n";
                        dump.close();
                    }
                }
#endif
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
                    // Single point: trust the algorithm, add to results
                    // The bounding method guarantees that if a box converges to a single point,
                    // it contains a root (or is very close to one within tolerance)
                    std::vector<double> center(dim);
                    for (std::size_t i = 0; i < dim; ++i) {
                        center[i] = 0.5 * (node.box_lower[i] + node.box_upper[i]);
                    }

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

        // Skip subdivision if box was pruned
        if (box_was_pruned) {
            continue;
        }

        // --- Subdivision step ----------------------------------------------
        // Subdivide along axes that are not small enough.
        // Only subdivide dimensions whose width is larger than tolerance.
        std::vector<bool> split_dim(dim, false);

        if (config.strategy == SubdivisionStrategy::Simultaneous) {
            // For Simultaneous strategy, only subdivide dimensions that didn't contract enough
            split_dim = final_needs_subdivision;
        } else {
            // For ContractFirst and SubdivideFirst, subdivide all dimensions not small enough
            for (std::size_t i = 0; i < dim; ++i) {
                const double width = node.box_upper[i] - node.box_lower[i];
                if (width > tolerance) {
                    split_dim[i] = true;
                }
            }
        }

        // Track subdivision boxes per depth for degeneracy detection
        subdivision_boxes_per_depth[node.depth]++;

        // Check for degeneracy: too many boxes needing subdivision at this depth
        if (!degeneracy_mode && subdivision_boxes_per_depth[node.depth] > degeneracy_threshold) {
            degeneracy_mode = true;
            result.degeneracy_detected = true;

            std::cerr << "Warning: Degeneracy detected at depth " << node.depth
                      << " (subdivision boxes: " << subdivision_boxes_per_depth[node.depth]
                      << ", threshold: " << degeneracy_threshold << ")" << std::endl;
        }

        // If in degeneracy mode, add box to unresolved instead of subdividing
        if (degeneracy_mode) {
#ifdef ENABLE_GEOMETRY_DUMP
            if (config.dump_geometry && method == RootBoundingMethod::ProjectedPolyhedral) {
                std::ofstream dump(dump_file.c_str(), std::ios::app);
                if (dump.is_open()) {
                    dump << "# FINAL_DECISION: UNRESOLVED (degeneracy mode, would subdivide in [";
                    bool first = true;
                    for (std::size_t i = 0; i < dim; ++i) {
                        if (split_dim[i]) {
                            if (!first) dump << ", ";
                            dump << "axis " << i;
                            first = false;
                        }
                    }
                    dump << "])\n\n";
                    dump.close();
                }
            }
#endif

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

        // Dump subdivision decision
#ifdef ENABLE_GEOMETRY_DUMP
        if (config.dump_geometry && method == RootBoundingMethod::ProjectedPolyhedral) {
            std::ofstream dump(dump_file.c_str(), std::ios::app);
            if (dump.is_open()) {
                dump << "# FINAL_DECISION: SUBDIVIDE in [";
                bool first = true;
                for (std::size_t i = 0; i < dim; ++i) {
                    if (split_dim[i]) {
                        if (!first) dump << ", ";
                        dump << "axis " << i;
                        first = false;
                    }
                }
                dump << "]\n\n";
                dump.close();
            }
        }
#endif

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

                // Propagate original polynomials to children for direct contraction
                left_child.original_polys = child.original_polys;
                right_child.original_polys = child.original_polys;

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

