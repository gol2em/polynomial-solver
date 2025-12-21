/**
 * @file polynomial_solver.h
 * @brief Unified header for the Polynomial Solver library
 * 
 * This is the main header file that users should include to access all
 * functionality of the Polynomial Solver library.
 * 
 * @section usage Basic Usage
 *
 * @subsection usage_simple Simple Usage (with defaults)
 *
 * @code{.cpp}
 * #include <polynomial_solver.h>
 * using namespace polynomial_solver;
 *
 * int main() {
 *     // 1. Define polynomial system
 *     std::vector<unsigned int> degrees = {3};
 *     std::vector<double> power_coeffs = {-0.08, 0.66, -1.5, 1.0};
 *     Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
 *     PolynomialSystem system(std::vector<Polynomial>{poly});
 *
 *     // 2. Solve (fast, double precision) - uses default config
 *     Solver solver;
 *     auto result = solver.subdivisionSolve(system, defaultSolverConfig(),
 *                                           RootBoundingMethod::ProjectedPolyhedral);
 *
 *     // 3. Refine (high precision, 1e-15) - uses default config
 *     ResultRefiner refiner;
 *     auto refined = refiner.refine(result, system, defaultRefinementConfig());
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
 * @subsection usage_custom Custom Configuration
 *
 * @code{.cpp}
 * // Customize solver parameters
 * SubdivisionConfig solver_config = defaultSolverConfig();
 * solver_config.tolerance = 1e-10;  // Tighter tolerance
 * solver_config.max_depth = 150;    // Allow deeper subdivision
 *
 * // Customize refinement parameters
 * RefinementConfig refine_config = defaultRefinementConfig();
 * refine_config.target_tolerance = 1e-12;   // Target error for refined roots
 * refine_config.residual_tolerance = 1e-12; // Maximum residual |f(x)|
 *
 * auto result = solver.subdivisionSolve(system, solver_config,
 *                                       RootBoundingMethod::ProjectedPolyhedral);
 * auto refined = refiner.refine(result, system, refine_config);
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

// Polynomial interpolation from function samples
#include "interpolation.h"

/**
 * @namespace polynomial_solver
 * @brief Main namespace for the Polynomial Solver library
 */

namespace polynomial_solver {

/**
 * @brief Create default solver configuration
 *
 * Returns a SubdivisionConfig with sensible defaults:
 * - tolerance: 1e-8 (box size threshold for convergence)
 * - max_depth: 100 (maximum subdivision depth)
 * - degeneracy_multiplier: 5.0 (degeneracy detection threshold)
 * - contraction_threshold: 0.8 (minimum contraction ratio)
 * - strategy: SubdivideFirst (balanced approach)
 *
 * @return Default solver configuration
 */
inline SubdivisionConfig defaultSolverConfig() {
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    config.max_depth = 100;
    config.degeneracy_multiplier = 5.0;
    config.contraction_threshold = 0.8;
    config.strategy = SubdivisionStrategy::SubdivideFirst;
    return config;
}

/**
 * @brief Create default refinement configuration
 *
 * Returns a RefinementConfig with sensible defaults:
 *
 * - target_tolerance: 1e-15 (target error tolerance)
 *   Used for TWO purposes:
 *   1. Condition-aware convergence: Root is accepted only if estimated_error < target_tolerance
 *   2. Exclusion radius computation: radius ≈ target_tolerance / |f'(x)| for simple roots
 *
 *   The refiner uses CONDITION-AWARE CONVERGENCE to prevent accepting inaccurate roots:
 *   - When |f(x)| < residual_tolerance, estimate condition number κ ≈ |f''| / |f'|²
 *   - Estimate actual error: error ≈ κ × |f(x)| / |f'(x)|
 *   - Accept root only if estimated_error < target_tolerance
 *
 *   This ensures ill-conditioned roots are rejected even when residual is small.
 *
 * - residual_tolerance: 1e-15 (residual threshold for convergence check)
 *   Newton's method checks convergence when |f(x)| < residual_tolerance.
 *   However, the root is NOT automatically accepted - the condition-aware criterion
 *   also checks if the estimated error is within target_tolerance.
 *
 *   For well-conditioned problems: small residual → small error (accepted)
 *   For ill-conditioned problems: small residual ≠ small error (rejected)
 *
 *   See docs/CONVERGENCE_CRITERIA_ANALYSIS.md for details.
 *
 * - max_newton_iters: 50 (maximum Newton iterations)
 * - max_multiplicity: 10 (maximum multiplicity to check)
 * - exclusion_multiplier: 3.0 (multiplier for exclusion radius)
 *
 * @return Default refinement configuration
 */
inline RefinementConfig defaultRefinementConfig() {
    RefinementConfig config;
    config.target_tolerance = 1e-15;
    config.residual_tolerance = 1e-15;
    config.max_newton_iters = 50;
    config.max_multiplicity = 10;
    config.exclusion_multiplier = 3.0;
    return config;
}

} // namespace polynomial_solver

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

