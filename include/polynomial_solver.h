/**
 * @file polynomial_solver.h
 * @brief Unified header for the Polynomial Solver library
 * 
 * This is the main header file that users should include to access all
 * functionality of the Polynomial Solver library.
 * 
 * @section usage Basic Usage
 * 
 * @code{.cpp}
 * #include <polynomial_solver.h>
 * 
 * int main() {
 *     // 1. Define polynomial system
 *     std::vector<Polynomial> system = { ... };
 *     
 *     // 2. Solve (fast, double precision)
 *     SubdivisionSolver solver;
 *     auto result = solver.solve(system);
 *     
 *     // 3. Refine (high precision, 1e-15)
 *     ResultRefiner refiner;
 *     auto refined = refiner.refine1D(result, system[0]);
 *     
 *     // 4. Check results
 *     for (const auto& root : refined.roots) {
 *         if (root.needs_higher_precision) {
 *             // Use higher precision arithmetic (future feature)
 *         }
 *     }
 *     
 *     return 0;
 * }
 * @endcode
 * 
 * @section parameters Parameters
 * 
 * See docs/PARAMETERS.md for detailed parameter documentation.
 * 
 * Key parameters:
 * - tolerance: Box size threshold (default: 1e-8)
 * - max_depth: Maximum subdivision depth (default: 100)
 * - degeneracy_multiplier: Degeneracy detection threshold (default: 5.0)
 * 
 * @section conditioning Conditioning
 * 
 * See docs/CONDITIONING_AND_PRECISION.md for conditioning strategy.
 * 
 * The solver automatically detects ill-conditioned problems and flags roots
 * that need higher precision arithmetic.
 * 
 * @author gol2em
 * @date 2024
 */

#ifndef POLYNOMIAL_SOLVER_H
#define POLYNOMIAL_SOLVER_H

// Core polynomial functionality
#include "polynomial.h"

// Solver functionality
#include "solver.h"

// Result refinement
#include "result_refiner.h"

// Differentiation (for advanced users)
#include "differentiation.h"

// Geometry algorithms (for advanced users)
#include "geometry.h"

// De Casteljau algorithm (for advanced users)
#include "de_casteljau.h"

/**
 * @namespace polynomial_solver
 * @brief Main namespace for the Polynomial Solver library
 * 
 * All classes and functions are in the global namespace for simplicity.
 * Future versions may use a namespace.
 */

/**
 * @mainpage Polynomial Solver Library
 * 
 * @section intro Introduction
 * 
 * A high-performance C++11 library for solving systems of multivariate
 * polynomial equations using subdivision methods with Bernstein basis
 * representation.
 * 
 * @section features Key Features
 * 
 * - Multivariate polynomial support (arbitrary dimensions and degrees)
 * - Bernstein basis representation (numerically stable)
 * - Multiple root bounding methods (GraphHull, ProjectedPolyhedral)
 * - Automatic degeneracy detection
 * - High-precision result refinement (1e-15 with double precision)
 * - Condition number estimation (detects ill-conditioned problems)
 * - Header-only option (future feature)
 * - No external dependencies (pure C++11)
 * 
 * @section workflow Typical Workflow
 * 
 * 1. **Define polynomial system**
 *    - Create Polynomial objects from control points
 *    - Or convert from power basis
 * 
 * 2. **Solve** (fast, double precision)
 *    - Use SubdivisionSolver to find boxes containing roots
 *    - Configure tolerance, max_depth, etc.
 * 
 * 3. **Refine** (high precision)
 *    - Use ResultRefiner to refine roots to 1e-15 precision
 *    - Automatic multiplicity detection
 *    - Condition number estimation
 * 
 * 4. **Check conditioning**
 *    - If needs_higher_precision flag is set, use higher precision arithmetic
 *    - See docs/CONDITIONING_AND_PRECISION.md
 * 
 * @section examples Examples
 * 
 * See examples/ directory for complete examples:
 * - examples/cubic_1d_roots.cpp: 1D cubic polynomial
 * - examples/multiplicity_1d_roots.cpp: Multiple roots
 * - examples/wilkinson_1d_roots.cpp: Ill-conditioned problem
 * - examples/circle_ellipse_intersection.cpp: 2D system
 * 
 * @section docs Documentation
 * 
 * - docs/PARAMETERS.md: Parameter reference
 * - docs/CONDITIONING_AND_PRECISION.md: Conditioning strategy
 * - docs/QUICKSTART.md: Quick start guide
 * - docs/GEOMETRY_ALGORITHMS.md: Geometric algorithms
 * - docs/result_refinement_design.md: Refinement design
 * 
 * @section building Building
 * 
 * @code{.sh}
 * # Build library and examples
 * ./build.sh
 * 
 * # Build with tests
 * ./build.sh --test
 * 
 * # Link against library
 * g++ -std=c++11 my_code.cpp -I include -L build/lib -lpolynomial_solver
 * @endcode
 * 
 * @section license License
 * 
 * TBD
 */

#endif // POLYNOMIAL_SOLVER_H

