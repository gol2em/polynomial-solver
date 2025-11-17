#ifndef SOLVER_H
#define SOLVER_H

/**
 * @file solver.h
 * @brief Main solver interface for polynomial root finding
 *
 * This module provides the main interface for solving polynomial equations
 * using various algorithms including De Casteljau subdivision.
 */

#include <vector>
#include <string>

#include "polynomial.h"
#include "geometry.h"
#include "de_casteljau.h"

namespace polynomial_solver {

/**
 * @brief Control net for the graph of a scalar polynomial f: [0,1]^n -> R.
 *
 * The control points live in R^{n+1}. The underlying Bernstein degrees per
 * variable are stored in degrees, and control_points is a flat array of
 * length prod_i (degrees[i] + 1) * (degrees.size() + 1).
 */
struct GraphControlNet {
    /// Degrees per variable of the underlying polynomial.
    std::vector<unsigned int> degrees;

    /// Control points of the graph hypersurface (x, f(x)) on [0,1]^n.
    ///
    /// Layout: for each coefficient index q in
    /// [0, prod_i (degrees[i] + 1)), the block
    ///   control_points[q * (degrees.size() + 1) + j]   (0 <= j < degrees.size())
    /// stores the j-th coordinate in [0,1], and
    ///   control_points[q * (degrees.size() + 1) + degrees.size()]
    /// stores the function value.
    std::vector<double> control_points;
};

/**
 * @class PolynomialSystem
 * @brief Represents a system of polynomial equations F: [0,1]^n -> R^m.
 */
class PolynomialSystem {
public:
    /// Construct an empty system with dimension 0.
    PolynomialSystem();

    /// Construct a system from a list of equations.
    explicit PolynomialSystem(const std::vector<Polynomial>& equations);

    /// Dimension of the domain (n in [0,1]^n).
    std::size_t dimension() const;

    /// Number of scalar equations (m in R^m).
    std::size_t equationCount() const;

    /// Access a single equation by index.
    const Polynomial& equation(std::size_t i) const;

    /// Access all equations.
    const std::vector<Polynomial>& equations() const;

    /**
     * @brief Build graph control nets (x, f(x)) for all equations.
     *
     * The system is assumed to be defined on [0,1]^dimension(). Each
     * equation's graph is a control net in R^{n+1}, where the first n
     * coordinates encode the normalized domain position and the last
     * coordinate encodes the polynomial value.
     */
    std::vector<GraphControlNet> graphControlNets() const;

    /**
     * @brief Evaluate all equations at a point.
     *
     * @param point Parameter values for each variable (size = dimension()).
     * @param values Output vector of equation values (size = equationCount()).
     */
    void evaluate(const std::vector<double>& point, std::vector<double>& values) const;

    /**
     * @brief Check if a point is approximately a root of the system.
     *
     * A point is considered a root if all equation values are within tolerance.
     *
     * @param point Parameter values for each variable (size = dimension()).
     * @param tolerance Maximum absolute value for each equation (default: 1e-6).
     * @return true if point is approximately a root, false otherwise.
     */
    bool isApproximateRoot(const std::vector<double>& point, double tolerance = 1e-6) const;

private:
    std::size_t dimension_;
    std::vector<Polynomial> equations_;
};

/**
 * @brief Root bounding method selector for the subdivision solver.
 *
 * - None: No contraction; root bounding box equals the current box.
 * - GraphHull: Uses exact convex hull of graph control points (implemented for 1D and 2D systems).
 * - ProjectedPolyhedral: Projects graph control points direction-by-direction to compute bounds.
 * - IntervalArithmetic: Placeholder for interval-arithmetic-based bounding (not yet implemented).
 */
enum class RootBoundingMethod {
    None,                ///< No contraction; root bounding box equals the current box.
    GraphHull,           ///< Exact convex-hull-based bounding in graph space (1D and 2D only).
    ProjectedPolyhedral, ///< Direction-by-direction projection of graph control points.
    IntervalArithmetic   ///< Placeholder for interval-arithmetic-based bounding.
};

/**
 * @brief Subdivision strategy selector for the subdivision solver.
 *
 * - ContractFirst: Always contract first, subdivide only when contraction doesn't shrink enough (default).
 * - SubdivideFirst: When any direction doesn't shrink enough, immediately subdivide in all directions.
 * - Simultaneous: In one step, subdivide in directions that didn't shrink enough, contract in others.
 */
enum class SubdivisionStrategy {
    ContractFirst,   ///< Contract iteratively, subdivide when stuck (default).
    SubdivideFirst,  ///< Subdivide immediately when contraction is insufficient.
    Simultaneous     ///< Subdivide non-contracting directions, contract others simultaneously.
};

/**
 * @brief Configuration parameters for the subdivision-based solver.
 */
struct SubdivisionConfig {
    double tolerance;          ///< Absolute minimum width of a box in each dimension (convergence criterion).
    unsigned int max_depth;    ///< Maximum allowed subdivision depth (for safety/debug, should be large).
    double degeneracy_multiplier; ///< Multiplier for expected root count to detect degeneracy (default: 5.0).
    double contraction_threshold; ///< Minimum contraction ratio to consider progress (default: 0.9).
    SubdivisionStrategy strategy; ///< Subdivision strategy to use (default: ContractFirst).
    bool dump_geometry;        ///< If true, dump detailed geometric information during solving.
    std::string dump_prefix;   ///< Prefix for dump files (default: "dump").

    SubdivisionConfig()
        : tolerance(1e-8),
          max_depth(100u),
          degeneracy_multiplier(5.0),
          contraction_threshold(0.9),
          strategy(SubdivisionStrategy::ContractFirst),
          dump_geometry(false),
          dump_prefix("dump")
    {
    }
};

/**
 * @brief Resulting box from the subdivision solver.
 *
 * Each box represents a region in the normalized domain [0,1]^n.
 * The center of the box is the estimated root, and the half-widths are the max errors.
 */
struct SubdivisionBoxResult {
    std::vector<double> lower;      ///< Lower corner of the box in each dimension.
    std::vector<double> upper;      ///< Upper corner of the box in each dimension.
    std::vector<double> center;     ///< Center of the box (estimated root location).
    std::vector<double> max_error;  ///< Half-width of the box in each dimension (max error, machine epsilon if point).
    unsigned int depth;             ///< Subdivision depth at which this box was produced.
    bool converged;                 ///< True if the box was terminated because it was small enough.
};

/**
 * @brief Complete result from the subdivision solver.
 *
 * Contains resolved roots (converged boxes) and unresolved roots (degenerate or max depth).
 * The first `num_resolved` boxes are resolved roots with small error.
 * The remaining boxes are unresolved (degenerate cases or max depth reached).
 */
struct SubdivisionSolverResult {
    std::vector<SubdivisionBoxResult> boxes;  ///< All boxes (resolved first, then unresolved).
    std::size_t num_resolved;                 ///< Number of resolved roots (first num_resolved boxes).
    bool degeneracy_detected;                 ///< True if degeneracy was detected during solving.
};

/**
 * @class Solver
 * @brief Main interface for polynomial solving
 */
class Solver {
public:
    /**
     * @brief Default constructor
     */
    Solver();

    /**
     * @brief Destructor
     */
    ~Solver();

    /**
     * @brief Run a subdivision-based search for root regions.
     *
     * The system is assumed to be defined on the normalized domain [0,1]^n.
     * The algorithm maintains a priority queue of boxes ordered by subdivision
     * depth, contracts them using a root bounding method, and either refines or
     * subdivides them, until all boxes are small enough or the maximum depth is
     * reached.
     *
     * Returns a result containing:
     * - Resolved roots (converged boxes with error < tolerance)
     * - Unresolved roots (degenerate cases or max depth reached)
     * - A counter separating the two groups
     */
    SubdivisionSolverResult
    subdivisionSolve(const PolynomialSystem& system,
                     const SubdivisionConfig& config,
                     RootBoundingMethod method = RootBoundingMethod::None) const;

private:
    // TODO: Add member variables
};

} // namespace polynomial_solver

#endif // SOLVER_H

