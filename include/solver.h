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

private:
    std::size_t dimension_;
    std::vector<Polynomial> equations_;
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

    // TODO: Add solver operations
    // - Find all roots in an interval
    // - Find roots with specified precision
    // - Isolate roots
    // - Refine root approximations

private:
    // TODO: Add member variables
};

} // namespace polynomial_solver

#endif // SOLVER_H

