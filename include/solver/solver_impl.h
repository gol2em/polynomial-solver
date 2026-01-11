#ifndef SOLVER_IMPL_H
#define SOLVER_IMPL_H

/**
 * @file detail/solver_impl.h
 * @brief Template implementation for SolverBase<Scalar>
 *
 * This file contains the implementation of SolverBase methods.
 * It should only be included from solver_base.h.
 */

#include "core/geometry.h"       // For original convex hull (double fallback)
#include "core/geometry_base.h"  // For templated geometry

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <queue>

// Debug flag for tracing solver execution
// #define SOLVER_DEBUG_TRACE

namespace polynomial_solver {

//=============================================================================
// Helper: Analyze box dimension
//=============================================================================

template<typename Scalar>
int SolverBase<Scalar>::analyzeBoxDimension(
    const std::vector<Scalar>& lower,
    const std::vector<Scalar>& upper,
    const Scalar& tolerance,
    std::size_t& active_axis)
{
    const std::size_t dim = lower.size();
    std::size_t num_active = 0;
    active_axis = 0;

    for (std::size_t i = 0; i < dim; ++i) {
        Scalar width = upper[i] - lower[i];
        if (width > tolerance) {
            num_active++;
            active_axis = i;
        }
    }

    if (num_active == 0) return 1;      // Single point
    else if (num_active == 1) return 2; // Line segment
    else return 3;                       // Higher dimensional
}

//=============================================================================
// Helper: Compute projected polyhedral bounds
//=============================================================================

template<typename Scalar>
bool SolverBase<Scalar>::computeProjectedPolyhedralBounds(
    const std::vector<Poly>& polys,
    std::size_t dim,
    std::vector<Scalar>& local_bound_lower,
    std::vector<Scalar>& local_bound_upper)
{
    if (polys.empty() || dim == 0u) {
        return false;
    }

    local_bound_lower.resize(dim);
    local_bound_upper.resize(dim);

    // Process each direction independently
    for (std::size_t dir = 0; dir < dim; ++dir) {
        Scalar dir_min = Scalar(0);
        Scalar dir_max = Scalar(1);

        for (std::size_t eq_idx = 0; eq_idx < polys.size(); ++eq_idx) {
            const Poly& poly = polys[eq_idx];

            // Get graph control points
            std::vector<Scalar> control_points;
            poly.graphControlPoints(control_points);

            const std::size_t point_dim = dim + 1u;
            const std::size_t num_points = control_points.size() / point_dim;

            // Project to 2D and find extreme points
            std::map<Scalar, std::vector<Scalar>> vertical_groups;

            for (std::size_t i = 0; i < num_points; ++i) {
                Scalar x = control_points[i * point_dim + dir];
                Scalar y = control_points[i * point_dim + dim];
                vertical_groups[x].push_back(y);
            }

            // Extract extreme points and compute 2D convex hull using templated geometry
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

            // Compute convex hull using templated geometry
            ConvexPolyhedronBase<Scalar> hull_2d = convex_hull_impl(projected_2d);

            // Intersect with y=0 axis (matching original solver's tolerance handling)
            Scalar eq_min = Scalar(0);
            Scalar eq_max = Scalar(1);
            bool found_interval = false;

            const std::vector<std::vector<Scalar>>& verts = hull_2d.vertices;
            if (verts.size() >= 2) {
                std::size_t n = verts.size();
                std::vector<Scalar> x_crossings;
                const Scalar eps = GeometryToleranceTraits<Scalar>::epsilon();

                // First, collect vertices that are (almost) on y=0
                for (std::size_t i = 0; i < n; ++i) {
                    Scalar y = verts[i][1];
                    Scalar abs_y = (y >= Scalar(0)) ? y : -y;
                    if (abs_y <= eps) {
                        x_crossings.push_back(verts[i][0]);
                    }
                }

                // Then check edges for crossings (skip edges with endpoints on y=0)
                for (std::size_t i = 0; i < n; ++i) {
                    std::size_t j = (i + 1) % n;
                    Scalar y1 = verts[i][1];
                    Scalar y2 = verts[j][1];
                    Scalar abs_y1 = (y1 >= Scalar(0)) ? y1 : -y1;
                    Scalar abs_y2 = (y2 >= Scalar(0)) ? y2 : -y2;

                    // Skip if either endpoint is (almost) on y=0 - already handled above
                    if (abs_y1 <= eps || abs_y2 <= eps) {
                        continue;
                    }

                    // Check if edge crosses y=0 (endpoints on opposite sides)
                    bool above1 = y1 > Scalar(0);
                    bool above2 = y2 > Scalar(0);
                    if (above1 != above2) {
                        Scalar t = y1 / (y1 - y2);
                        Scalar x_cross = verts[i][0] + t * (verts[j][0] - verts[i][0]);
                        x_crossings.push_back(x_cross);
                    }
                }

                if (!x_crossings.empty()) {
                    eq_min = *std::min_element(x_crossings.begin(), x_crossings.end());
                    eq_max = *std::max_element(x_crossings.begin(), x_crossings.end());
                    found_interval = true;
                }
            }

            if (!found_interval) {
                // No intersection - box can be pruned
                return false;
            }

            // Intersect with current bounds
            dir_min = (eq_min > dir_min) ? eq_min : dir_min;
            dir_max = (eq_max < dir_max) ? eq_max : dir_max;

            if (dir_min > dir_max) {
                return false;  // Empty intersection
            }
        }

        local_bound_lower[dir] = dir_min;
        local_bound_upper[dir] = dir_max;
    }

    return true;
}

//=============================================================================
// Main subdivision solve implementation
//=============================================================================

template<typename Scalar>
typename SolverBase<Scalar>::SolverResult
SolverBase<Scalar>::subdivisionSolve(
    const System& system,
    const Config& config,
    RootBoundingMethodBase method) const
{
    SolverResult result;
    result.num_resolved = 0;
    result.degeneracy_detected = false;

    std::vector<BoxResult> resolved_boxes;
    std::vector<BoxResult> unresolved_boxes;

    const std::size_t dim = system.dimension();
    if (dim == 0u) {
        return result;
    }

    // Compute expected maximum number of roots
    std::size_t expected_max_roots = 1;
    for (const Poly& eq : system.equations()) {
        const std::vector<unsigned int>& degrees = eq.degrees();
        for (unsigned int d : degrees) {
            expected_max_roots *= d;
        }
    }
    const std::size_t degeneracy_threshold =
        static_cast<std::size_t>(static_cast<double>(config.degeneracy_multiplier) *
                                 static_cast<double>(expected_max_roots));

    // Root node: full box [0,1]^dim
    SubdivisionNode root;
    root.box_lower.assign(dim, Scalar(0));
    root.box_upper.assign(dim, Scalar(1));
    root.depth = 0u;

    // Convert polynomials to Bernstein
    root.polys.reserve(system.equations().size());
    for (const Poly& poly : system.equations()) {
        Poly bernstein_poly = poly;
        bernstein_poly.ensureBernsteinPrimary();
        root.polys.push_back(bernstein_poly);
    }
    root.original_polys = root.polys;

    // Priority queue
    std::priority_queue<NodeQueueEntry, std::vector<NodeQueueEntry>, NodeQueueCompare> queue;
    std::size_t next_index = 0u;
    queue.push(NodeQueueEntry{root.depth, next_index++, root});

    const Scalar tolerance = config.tolerance;
    std::map<unsigned int, std::size_t> subdivision_boxes_per_depth;

    while (!queue.empty()) {
        NodeQueueEntry entry = queue.top();
        queue.pop();

        SubdivisionNode node = std::move(entry.node);

        // Check depth limit
        if (node.depth >= config.max_depth) {
            if (config.verbose_warnings) {
                std::cerr << "Warning: Maximum depth " << config.max_depth << " reached." << std::endl;
            }

            BoxResult box;
            box.lower = node.box_lower;
            box.upper = node.box_upper;
            box.center.resize(dim);
            box.max_error.resize(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                box.center[i] = Scalar(0.5) * (node.box_lower[i] + node.box_upper[i]);
                Scalar half_width = Scalar(0.5) * (node.box_upper[i] - node.box_lower[i]);
                const Scalar machine_eps = MachineEpsilonTraits<Scalar>::epsilon();
                box.max_error[i] = (half_width < machine_eps) ? machine_eps : half_width;
            }
            box.depth = node.depth;
            box.converged = false;
            unresolved_boxes.push_back(std::move(box));
            continue;
        }

        // Iterative bounding step
        bool converged = false;
        const unsigned int max_iterations = 100u;
        std::vector<bool> final_needs_subdivision(dim, false);
        bool box_was_pruned = false;

        for (unsigned int iter = 0; iter < max_iterations; ++iter) {
            std::vector<Scalar> local_bound_lower(dim, Scalar(0));
            std::vector<Scalar> local_bound_upper(dim, Scalar(1));

            bool has_roots = true;
            if (method == RootBoundingMethodBase::ProjectedPolyhedral) {
                has_roots = computeProjectedPolyhedralBounds(
                    node.polys, dim, local_bound_lower, local_bound_upper);
            }

#ifdef SOLVER_DEBUG_TRACE
            std::cerr << "  BOUND: box [" << static_cast<double>(node.box_lower[0])
                      << ", " << static_cast<double>(node.box_upper[0]) << "] iter=" << iter
                      << " has_roots=" << has_roots;
            if (has_roots) {
                std::cerr << " local=[" << static_cast<double>(local_bound_lower[0])
                          << ", " << static_cast<double>(local_bound_upper[0]) << "]";
            }
            std::cerr << "\n";
#endif

            if (!has_roots) {
                box_was_pruned = true;
                break;
            }

            // Map local bounds to global coordinates
            std::vector<Scalar> new_global_lower(dim);
            std::vector<Scalar> new_global_upper(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                Scalar width = node.box_upper[i] - node.box_lower[i];
                new_global_lower[i] = node.box_lower[i] + local_bound_lower[i] * width;
                new_global_upper[i] = node.box_lower[i] + local_bound_upper[i] * width;
            }

            // Check contraction and convergence - matches original solver.cpp logic
            std::vector<Scalar> contraction_ratio(dim);
            std::vector<bool> needs_subdivision(dim, false);
            bool any_needs_subdivision = false;

            for (std::size_t i = 0; i < dim; ++i) {
                Scalar old_width = node.box_upper[i] - node.box_lower[i];
                Scalar new_width = (local_bound_upper[i] - local_bound_lower[i]) * old_width;
                contraction_ratio[i] = (old_width > Scalar(0))
                    ? (new_width / old_width) : Scalar(0);

                // Check if this direction contracted enough
                if (contraction_ratio[i] > config.contraction_threshold &&
                    old_width > tolerance) {
                    needs_subdivision[i] = true;
                    any_needs_subdivision = true;
                }
            }

            // Apply strategy-specific logic
            bool should_subdivide_now = false;
            bool should_contract = false;

            if (config.strategy == SubdivisionStrategyBase::SubdivideFirst) {
                // SubdivideFirst: If ANY direction didn't contract enough, subdivide ALL immediately
                if (any_needs_subdivision) {
                    should_subdivide_now = true;
                    should_contract = false;
                } else {
                    // All directions contracted well - contract and continue iteration
                    should_contract = true;
                }
            } else if (config.strategy == SubdivisionStrategyBase::Simultaneous) {
                // Simultaneous: Contract in good directions, subdivide in bad ones - same step
                if (any_needs_subdivision) {
                    should_subdivide_now = true;
                    should_contract = true;  // Contract only in well-contracting directions
                } else {
                    // All directions contracted well - just contract and continue
                    should_contract = true;
                }
            } else if (config.strategy == SubdivisionStrategyBase::ContractFirst) {
                // ContractFirst: If ANY direction has good ratio, contract those and continue iteration
                // Only subdivide when NO direction contracts well anymore
                bool any_good_ratio = false;
                for (std::size_t i = 0; i < dim; ++i) {
                    Scalar old_width = node.box_upper[i] - node.box_lower[i];
                    if (old_width > tolerance && contraction_ratio[i] <= config.contraction_threshold) {
                        any_good_ratio = true;
                        break;
                    }
                }

                if (any_good_ratio) {
                    // Contract in good directions, leave others unchanged, continue iteration
                    should_contract = true;
                    should_subdivide_now = false;  // Keep iterating!
                } else if (any_needs_subdivision) {
                    // No good contractions and box not converged - NOW subdivide
                    should_subdivide_now = true;
                    should_contract = false;
                }
            }

            // Contract the box (if strategy allows)
            if (should_contract) {
                for (std::size_t i = 0; i < dim; ++i) {
                    // Determine which directions to contract based on strategy
                    bool skip_this_direction = false;
                    if (config.strategy == SubdivisionStrategyBase::Simultaneous) {
                        // Skip directions that need subdivision
                        skip_this_direction = needs_subdivision[i];
                    } else if (config.strategy == SubdivisionStrategyBase::ContractFirst) {
                        // Skip directions with bad ratio (ratio > threshold)
                        Scalar old_width = node.box_upper[i] - node.box_lower[i];
                        skip_this_direction = (old_width > tolerance &&
                                               contraction_ratio[i] > config.contraction_threshold);
                    }
                    // SubdivideFirst: contract all directions (no skip)

                    if (skip_this_direction) {
                        continue;
                    }

                    Scalar a = local_bound_lower[i];
                    Scalar b = local_bound_upper[i];
                    Scalar old_low = node.box_lower[i];
                    Scalar old_high = node.box_upper[i];
                    Scalar old_width = old_high - old_low;

                    node.box_lower[i] = old_low + a * old_width;
                    node.box_upper[i] = old_low + b * old_width;
                }

                // Direct contraction: Recompute polynomials from original using global box
                // This eliminates error accumulation from repeated restrictions
                // (matches original solver.cpp behavior)
                for (std::size_t eq_idx = 0; eq_idx < node.polys.size(); ++eq_idx) {
                    node.polys[eq_idx] = node.original_polys[eq_idx];
                    for (std::size_t i = 0; i < dim; ++i) {
                        // Restrict from original [0,1] to global [box_lower[i], box_upper[i]]
                        node.polys[eq_idx] = node.polys[eq_idx].restrictedToInterval(
                            i, node.box_lower[i], node.box_upper[i]);
                    }
                }
            }

            // If strategy says subdivide now, break out of iteration loop
            if (should_subdivide_now) {
                final_needs_subdivision = needs_subdivision;
                break;
            }
            // Otherwise, continue iteration (ContractFirst will keep contracting until no progress)

            // Check if contracted box is small enough
            bool all_small = true;
            for (std::size_t i = 0; i < dim; ++i) {
                Scalar width = node.box_upper[i] - node.box_lower[i];
                if (width > tolerance) {
                    all_small = false;
                    break;
                }
            }

            if (all_small) {
                converged = true;
                break;
            }
        }

        if (box_was_pruned) {
#ifdef SOLVER_DEBUG_TRACE
            std::cerr << "PRUNED: box [" << static_cast<double>(node.box_lower[0])
                      << ", " << static_cast<double>(node.box_upper[0]) << "] depth=" << node.depth << "\n";
#endif
            continue;
        }

        if (converged) {
            BoxResult box;
            box.lower = node.box_lower;
            box.upper = node.box_upper;
            box.center.resize(dim);
            box.max_error.resize(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                box.center[i] = Scalar(0.5) * (node.box_lower[i] + node.box_upper[i]);
                Scalar half_width = Scalar(0.5) * (node.box_upper[i] - node.box_lower[i]);
                const Scalar machine_eps = MachineEpsilonTraits<Scalar>::epsilon();
                box.max_error[i] = (half_width < machine_eps) ? machine_eps : half_width;
            }
            box.depth = node.depth;
            box.converged = true;
#ifdef SOLVER_DEBUG_TRACE
            std::cerr << "CONVERGED: box [" << static_cast<double>(node.box_lower[0])
                      << ", " << static_cast<double>(node.box_upper[0]) << "] depth=" << node.depth << "\n";
#endif
            resolved_boxes.push_back(std::move(box));
            continue;
        }

        // Determine which dimensions to subdivide (matching original solver.cpp logic)
        std::vector<bool> split_dims(dim, false);

        if (config.strategy == SubdivisionStrategyBase::Simultaneous) {
            // Only subdivide dimensions that didn't contract enough
            split_dims = final_needs_subdivision;
        } else {
            // ContractFirst and SubdivideFirst: subdivide all dimensions not small enough
            for (std::size_t i = 0; i < dim; ++i) {
                Scalar width = node.box_upper[i] - node.box_lower[i];
                if (width > tolerance) {
                    split_dims[i] = true;
                }
            }
        }

        // Track subdivision boxes per depth for degeneracy detection
        subdivision_boxes_per_depth[node.depth]++;

        // Check for degeneracy: too many boxes needing subdivision at this depth
        if (subdivision_boxes_per_depth[node.depth] > degeneracy_threshold) {
            if (!result.degeneracy_detected) {
                result.degeneracy_detected = true;
                // Once degeneracy is detected, stop subdividing and add remaining boxes as unresolved
            }
        }

        // If degeneracy detected, add box to unresolved instead of subdividing further
        if (result.degeneracy_detected) {
            BoxResult box;
            box.lower.resize(dim);
            box.upper.resize(dim);
            box.center.resize(dim);
            box.max_error.resize(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                box.lower[i] = static_cast<double>(node.box_lower[i]);
                box.upper[i] = static_cast<double>(node.box_upper[i]);
                box.center[i] = 0.5 * (box.lower[i] + box.upper[i]);
                double half_width = 0.5 * (box.upper[i] - box.lower[i]);
                box.max_error[i] = (half_width < std::numeric_limits<double>::epsilon())
                    ? std::numeric_limits<double>::epsilon()
                    : half_width;
            }
            box.depth = node.depth;
            box.converged = false;
            unresolved_boxes.push_back(std::move(box));
            continue;
        }

        // Create children by subdividing in all marked dimensions
        // This produces 2^k children for k dimensions being split
        std::vector<SubdivisionNode> children;
        SubdivisionNode seed = node;
        seed.depth = node.depth + 1;
        children.push_back(seed);

        for (std::size_t axis = 0; axis < dim; ++axis) {
            if (!split_dims[axis]) continue;

            std::vector<SubdivisionNode> next_children;
            next_children.reserve(children.size() * 2);

            for (const SubdivisionNode& child : children) {
                Scalar old_low = child.box_lower[axis];
                Scalar old_high = child.box_upper[axis];
                Scalar mid = Scalar(0.5) * (old_low + old_high);

                SubdivisionNode left_child = child;
                SubdivisionNode right_child = child;

                left_child.box_upper[axis] = mid;
                right_child.box_lower[axis] = mid;

                // Restrict polynomials using LOCAL [0, 0.5] and [0.5, 1]
                // The polynomial's region tracking handles the global mapping
                left_child.polys.clear();
                right_child.polys.clear();
                for (const Poly& poly : child.polys) {
                    left_child.polys.push_back(
                        poly.restrictedToInterval(axis, Scalar(0), Scalar(0.5)));
                    right_child.polys.push_back(
                        poly.restrictedToInterval(axis, Scalar(0.5), Scalar(1)));
                }

                // Propagate original polynomials to children for direct contraction
                // (matches original solver.cpp behavior)
                left_child.original_polys = child.original_polys;
                right_child.original_polys = child.original_polys;

                next_children.push_back(std::move(left_child));
                next_children.push_back(std::move(right_child));
            }
            children = std::move(next_children);
        }

        // Add all children to queue
#ifdef SOLVER_DEBUG_TRACE
        std::cerr << "SUBDIVIDE: box [" << static_cast<double>(node.box_lower[0])
                  << ", " << static_cast<double>(node.box_upper[0]) << "] depth=" << node.depth
                  << " -> " << children.size() << " children\n";
        for (const auto& child : children) {
            std::cerr << "  child [" << static_cast<double>(child.box_lower[0])
                      << ", " << static_cast<double>(child.box_upper[0]) << "]\n";
        }
#endif
        for (auto& child : children) {
            queue.push(NodeQueueEntry{child.depth, next_index++, std::move(child)});
        }
    }

    // Combine results
    result.boxes.reserve(resolved_boxes.size() + unresolved_boxes.size());
    for (auto& box : resolved_boxes) {
        result.boxes.push_back(std::move(box));
    }
    result.num_resolved = resolved_boxes.size();
    for (auto& box : unresolved_boxes) {
        result.boxes.push_back(std::move(box));
    }

    return result;
}

} // namespace polynomial_solver

#endif // SOLVER_IMPL_H

