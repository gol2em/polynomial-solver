#ifndef BOUNDING_STRATEGY_H
#define BOUNDING_STRATEGY_H

/**
 * @file bounding_strategy.h
 * @brief Modular bounding strategy interface for polynomial root solvers
 *
 * This header provides a pluggable interface for root bounding methods.
 * Different strategies can be swapped to determine if a root exists in a box
 * and to contract the bounding box around potential roots.
 *
 * Strategies include:
 * - NoBoundingStrategy: No contraction (baseline)
 * - ProjectedPolyhedralStrategy: Direction-by-direction projection
 * - GraphHullStrategy: Exact convex hull (1D/2D only)
 * - IntervalArithmeticStrategy: Interval-based bounds (future)
 */

#include "core/polynomial_base.h"
#include "core/geometry_base.h"
#include <vector>
#include <memory>

namespace polynomial_solver {

//=============================================================================
// RootBoundingMethodBase enum
//=============================================================================

/**
 * @brief Root bounding method selector
 */
enum class RootBoundingMethodBase {
    None,                ///< No contraction
    GraphHull,           ///< Exact convex-hull-based bounding (1D and 2D only)
    ProjectedPolyhedral, ///< Direction-by-direction projection
    IntervalArithmetic   ///< Placeholder for interval-based bounding
};

//=============================================================================
// BoundingResult - Result from a bounding computation
//=============================================================================

/**
 * @brief Result of a bounding strategy computation
 * @tparam Scalar The coefficient/evaluation type
 */
template<typename Scalar>
struct BoundingResult {
    bool has_root;                    ///< True if box may contain a root
    std::vector<Scalar> lower;        ///< Contracted lower bounds
    std::vector<Scalar> upper;        ///< Contracted upper bounds
    Scalar contraction_ratio;         ///< Ratio of volume reduction
    
    BoundingResult() : has_root(true), contraction_ratio(Scalar(1)) {}
};

//=============================================================================
// BoundingStrategyBase - Abstract interface for bounding strategies
//=============================================================================

/**
 * @class BoundingStrategyBase
 * @brief Abstract interface for root bounding strategies
 *
 * Implementations determine:
 * 1. Whether a box may contain roots (exclusion test)
 * 2. A tighter bounding box if roots are possible (contraction)
 *
 * @tparam Scalar The coefficient/evaluation type
 */
template<typename Scalar>
class BoundingStrategyBase {
public:
    using Poly = PolynomialBase<Scalar>;
    using Result = BoundingResult<Scalar>;

    virtual ~BoundingStrategyBase() = default;

    /**
     * @brief Compute bounds for roots of polynomials in a box
     *
     * @param polys Polynomials restricted to [0,1]^n
     * @param dim Dimension of the domain
     * @param current_lower Current box lower bounds (in [0,1]^n)
     * @param current_upper Current box upper bounds (in [0,1]^n)
     * @return BoundingResult with has_root=false if box provably empty
     */
    virtual Result computeBounds(
        const std::vector<Poly>& polys,
        std::size_t dim,
        const std::vector<Scalar>& current_lower,
        const std::vector<Scalar>& current_upper) const = 0;

    /**
     * @brief Get the strategy name for logging/debugging
     */
    virtual const char* name() const = 0;

    /**
     * @brief Check if strategy supports the given dimension
     * @param dim Dimension to check
     * @return True if strategy works for this dimension
     */
    virtual bool supportsDimension(std::size_t dim) const { 
        (void)dim;
        return true; 
    }
};

//=============================================================================
// NoBoundingStrategy - No contraction (baseline)
//=============================================================================

/**
 * @class NoBoundingStrategy
 * @brief No-op bounding strategy that returns the input box unchanged
 */
template<typename Scalar>
class NoBoundingStrategy : public BoundingStrategyBase<Scalar> {
public:
    using Poly = PolynomialBase<Scalar>;
    using Result = BoundingResult<Scalar>;

    Result computeBounds(
        const std::vector<Poly>& /*polys*/,
        std::size_t dim,
        const std::vector<Scalar>& current_lower,
        const std::vector<Scalar>& current_upper) const override
    {
        Result result;
        result.has_root = true;
        result.lower = current_lower;
        result.upper = current_upper;
        result.contraction_ratio = Scalar(1);
        
        // Initialize if empty
        if (result.lower.empty()) {
            result.lower.resize(dim, Scalar(0));
            result.upper.resize(dim, Scalar(1));
        }
        
        return result;
    }

    const char* name() const override { return "None"; }
};

//=============================================================================
// Forward declarations for strategy implementations
//=============================================================================

template<typename Scalar> class ProjectedPolyhedralStrategy;
template<typename Scalar> class GraphHullStrategy;
template<typename Scalar> class IntervalArithmeticStrategy;

} // namespace polynomial_solver

// Include strategy implementations
#include "solver/bounding_strategy_impl.h"

#endif // BOUNDING_STRATEGY_H

