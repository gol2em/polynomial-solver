#ifndef SOLVER_H
#define SOLVER_H

/**
 * @file solver.h
 * @brief Main solver interface for polynomial root finding
 * 
 * This module provides the main interface for solving polynomial equations
 * using various algorithms including De Casteljau subdivision.
 */

#include "polynomial.h"
#include "geometry.h"
#include "de_casteljau.h"

namespace polynomial_solver {

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

