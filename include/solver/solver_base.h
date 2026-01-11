#ifndef SOLVER_BASE_H
#define SOLVER_BASE_H

/**
 * @file solver_base.h
 * @brief Templated solver interface for polynomial root finding
 *
 * This header provides SolverBase<Scalar>, a fully templated solver that works
 * with any scalar type (double, mpreal, etc.). It unifies the functionality of
 * the non-templated Solver class into a flexible template.
 *
 * Key classes:
 * - PolynomialSystemBase<Scalar>: System of polynomial equations
 * - SolverBase<Scalar>: Main solver interface
 *
 * Usage:
 * @code
 * using Poly = PolynomialBase<double>;
 * Poly p1 = Poly::fromPower({2}, {-1.0, 0.0, 1.0});  // x^2 - 1
 * PolynomialSystemBase<double> system({p1});
 *
 * SolverBase<double> solver;
 * auto result = solver.subdivisionSolve(system, config);
 * @endcode
 */

#include "core/polynomial_base.h"
#include "core/geometry.h"
#include "solver/bounding_strategy.h"
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <algorithm>
#include <limits>

namespace polynomial_solver {

//=============================================================================
// Forward declarations
//=============================================================================

template<typename Scalar> class PolynomialSystemBase;
template<typename Scalar> class SolverBase;

//=============================================================================
// GraphControlNetBase - Control net for graph of scalar polynomial
//=============================================================================

/**
 * @brief Control net for the graph of a scalar polynomial f: [0,1]^n -> R.
 *
 * @tparam Scalar The coefficient/evaluation type
 */
template<typename Scalar>
struct GraphControlNetBase {
    std::vector<unsigned int> degrees;   ///< Degrees per variable
    std::vector<Scalar> control_points;  ///< Control points in R^{n+1}
};

//=============================================================================
// PolynomialSystemBase - System of polynomial equations
//=============================================================================

/**
 * @class PolynomialSystemBase
 * @brief Templated system of polynomial equations F: [0,1]^n -> R^m
 *
 * @tparam Scalar The coefficient/evaluation type
 */
template<typename Scalar>
class PolynomialSystemBase {
public:
    using Poly = PolynomialBase<Scalar>;

    /// Construct an empty system
    PolynomialSystemBase() : dimension_(0) {}

    /// Construct from a list of equations
    explicit PolynomialSystemBase(const std::vector<Poly>& equations)
        : dimension_(0), equations_(equations)
    {
        if (!equations_.empty()) {
            dimension_ = equations_[0].dimension();
        }
    }

    /// Dimension of the domain (n in [0,1]^n)
    std::size_t dimension() const { return dimension_; }

    /// Number of scalar equations
    std::size_t equationCount() const { return equations_.size(); }

    /// Access a single equation by index
    const Poly& equation(std::size_t i) const { return equations_[i]; }

    /// Access all equations
    const std::vector<Poly>& equations() const { return equations_; }

    /// Build graph control nets for all equations
    std::vector<GraphControlNetBase<Scalar>> graphControlNets() const {
        std::vector<GraphControlNetBase<Scalar>> nets;
        nets.reserve(equations_.size());

        for (const Poly& poly : equations_) {
            GraphControlNetBase<Scalar> net;
            net.degrees = poly.degrees();
            poly.graphControlPoints(net.control_points);
            nets.push_back(std::move(net));
        }

        return nets;
    }

    /// Evaluate all equations at a point
    void evaluate(const std::vector<Scalar>& point, std::vector<Scalar>& values) const {
        values.resize(equations_.size());
        for (std::size_t i = 0; i < equations_.size(); ++i) {
            values[i] = equations_[i].evaluate(point);
        }
    }

    /// Check if a point is approximately a root
    bool isApproximateRoot(const std::vector<Scalar>& point, 
                          const Scalar& tolerance = Scalar(1e-6)) const {
        if (point.size() != dimension_) {
            return false;
        }

        for (const Poly& eq : equations_) {
            Scalar value = eq.evaluate(point);
            Scalar abs_value = (value >= Scalar(0)) ? value : -value;
            if (abs_value > tolerance) {
                return false;
            }
        }

        return true;
    }

private:
    std::size_t dimension_;
    std::vector<Poly> equations_;
};

//=============================================================================
// Enums
//=============================================================================

// RootBoundingMethodBase is defined in bounding_strategy.h

/**
 * @brief Subdivision strategy selector
 */
enum class SubdivisionStrategyBase {
    ContractFirst,   ///< Contract iteratively, subdivide when stuck
    SubdivideFirst,  ///< Subdivide immediately when contraction is insufficient
    Simultaneous     ///< Mixed strategy
};

//=============================================================================
// SubdivisionConfigBase - Configuration for subdivision solver
//=============================================================================

/**
 * @brief Configuration for subdivision-based solver
 */
template<typename Scalar>
struct SubdivisionConfigBase {
    Scalar tolerance;                    ///< Minimum box width for convergence
    unsigned int max_depth;              ///< Maximum subdivision depth
    Scalar degeneracy_multiplier;        ///< Multiplier for degeneracy detection
    Scalar contraction_threshold;        ///< Minimum contraction ratio
    SubdivisionStrategyBase strategy;    ///< Subdivision strategy
    bool dump_geometry;                  ///< Enable geometry dump
    std::string dump_prefix;             ///< Prefix for dump files
    bool verbose_warnings;               ///< Print warnings

    SubdivisionConfigBase()
        : tolerance(Scalar(1e-8)),
          max_depth(100u),
          degeneracy_multiplier(Scalar(5)),
          contraction_threshold(Scalar(0.9)),
          strategy(SubdivisionStrategyBase::ContractFirst),
          dump_geometry(false),
          dump_prefix("dump"),
          verbose_warnings(false)
    {}
};

//=============================================================================
// SubdivisionBoxResultBase - Result box from solver
//=============================================================================

/**
 * @brief Result box from subdivision solver
 */
template<typename Scalar>
struct SubdivisionBoxResultBase {
    std::vector<Scalar> lower;       ///< Lower corner
    std::vector<Scalar> upper;       ///< Upper corner
    std::vector<Scalar> center;      ///< Center (estimated root)
    std::vector<Scalar> max_error;   ///< Half-width (max error)
    unsigned int depth;              ///< Subdivision depth
    bool converged;                  ///< True if converged
};

//=============================================================================
// SubdivisionSolverResultBase - Complete solver result
//=============================================================================

/**
 * @brief Complete result from subdivision solver
 */
template<typename Scalar>
struct SubdivisionSolverResultBase {
    std::vector<SubdivisionBoxResultBase<Scalar>> boxes;  ///< All result boxes
    std::size_t num_resolved;                             ///< Number of resolved roots
    bool degeneracy_detected;                             ///< Degeneracy flag
};

//=============================================================================
// SolverBase - Main solver interface
//=============================================================================

/**
 * @class SolverBase
 * @brief Templated polynomial root solver
 *
 * @tparam Scalar The coefficient/evaluation type
 */
template<typename Scalar>
class SolverBase {
public:
    using Poly = PolynomialBase<Scalar>;
    using System = PolynomialSystemBase<Scalar>;
    using Config = SubdivisionConfigBase<Scalar>;
    using BoxResult = SubdivisionBoxResultBase<Scalar>;
    using SolverResult = SubdivisionSolverResultBase<Scalar>;

    SolverBase() {}
    ~SolverBase() {}

    /**
     * @brief Run subdivision-based search for root regions
     *
     * The system is assumed to be defined on [0,1]^n.
     */
    SolverResult subdivisionSolve(
        const System& system,
        const Config& config,
        RootBoundingMethodBase method = RootBoundingMethodBase::None) const;

private:
    // Internal node structure for subdivision
    struct SubdivisionNode {
        std::vector<Scalar> box_lower;
        std::vector<Scalar> box_upper;
        unsigned int depth;
        std::vector<Poly> polys;
        std::vector<Poly> original_polys;
    };

    // Priority queue entry
    struct NodeQueueEntry {
        unsigned int depth;
        std::size_t index;
        SubdivisionNode node;
    };

    // Comparator for min-heap by depth
    struct NodeQueueCompare {
        bool operator()(const NodeQueueEntry& lhs, const NodeQueueEntry& rhs) const {
            if (lhs.depth != rhs.depth) return lhs.depth > rhs.depth;
            return lhs.index > rhs.index;
        }
    };

    // Helper: analyze box dimension
    static int analyzeBoxDimension(
        const std::vector<Scalar>& lower,
        const std::vector<Scalar>& upper,
        const Scalar& tolerance,
        std::size_t& active_axis);

    // Helper: compute projected polyhedral bounds
    static bool computeProjectedPolyhedralBounds(
        const std::vector<Poly>& polys,
        std::size_t dim,
        std::vector<Scalar>& local_bound_lower,
        std::vector<Scalar>& local_bound_upper);
};

} // namespace polynomial_solver

// Include template implementation
#include "solver/solver_impl.h"

#endif // SOLVER_BASE_H

