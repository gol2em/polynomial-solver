#ifndef RESULT_REFINER_H
#define RESULT_REFINER_H

/**
 * @file result_refiner.h
 * @brief Post-processing tool for refining and consolidating solver results
 *
 * This module provides tools to:
 * - Verify roots at high precision (1e-15)
 * - Estimate root multiplicity from derivatives
 * - Eliminate duplicate/nearby boxes
 * - Consolidate results into unique roots
 */

#include <vector>
#include <cstddef>
#include "polynomial.h"
#include "solver.h"

namespace polynomial_solver {

/**
 * @brief Configuration for result refinement
 */
struct RefinementConfig {
    double target_tolerance;        ///< Target precision for refined roots (default: 1e-15)
    double residual_tolerance;      ///< Max residual |f(x)| to accept as root (default: 1e-12)
    unsigned int max_newton_iters;  ///< Maximum Newton iterations (default: 50)
    unsigned int max_multiplicity;  ///< Maximum multiplicity to check (default: 10)
    double exclusion_multiplier;    ///< Multiplier for exclusion radius (default: 3.0)

    RefinementConfig()
        : target_tolerance(1e-15),
          residual_tolerance(1e-12),
          max_newton_iters(50u),
          max_multiplicity(10u),
          exclusion_multiplier(3.0)
    {}
};

/**
 * @brief A refined root with multiplicity and verification info
 */
struct RefinedRoot {
    std::vector<double> location;      ///< Root location in [0,1]^n
    std::vector<double> residual;      ///< f(root) for each equation
    std::vector<double> max_error;     ///< Error bounds per dimension
    unsigned int multiplicity;         ///< Estimated multiplicity (1 = simple root)
    std::vector<std::size_t> source_boxes;  ///< Indices of original boxes merged into this root
    bool verified;                     ///< True if passed high-precision verification
    unsigned int depth;                ///< Subdivision depth of primary source box
};

/**
 * @brief Result of refinement process
 */
struct RefinementResult {
    std::vector<RefinedRoot> roots;           ///< Verified and consolidated roots
    std::vector<std::size_t> cancelled_boxes; ///< Indices of boxes eliminated as duplicates
    std::vector<std::size_t> unverified_boxes; ///< Indices of boxes that didn't pass verification
};

/**
 * @brief Tool for refining and consolidating solver results
 */
class ResultRefiner {
public:
    /**
     * @brief Constructor
     */
    ResultRefiner();
    
    /**
     * @brief Destructor
     */
    ~ResultRefiner();
    
    /**
     * @brief Refine solver results using Newton's method with sign checking
     *
     * This method:
     * 1. For each resolved box, use Newton's method to refine to target_tolerance
     * 2. Use subdivision with sign checking to ensure convergence
     * 3. Verify convergence by checking residual
     * 4. Estimate multiplicity from derivatives
     * 5. Cancel nearby boxes within exclusion radius
     * 6. Returns consolidated unique roots
     *
     * Currently only supports 1D systems. 2D systems may have infinite roots
     * in degenerate regions.
     *
     * @param solver_result Raw result from subdivision solver
     * @param original_system Original polynomial system (on [0,1]^n domain)
     * @param config Refinement configuration
     * @return Refined and consolidated roots
     */
    RefinementResult refine(
        const SubdivisionSolverResult& solver_result,
        const PolynomialSystem& original_system,
        const RefinementConfig& config) const;

private:
    /**
     * @brief Refine a 1D root using Newton's method with sign checking
     *
     * Uses Newton iteration: x_{n+1} = x_n - f(x_n)/f'(x_n)
     * Checks sign changes to ensure convergence within the box.
     * Subdivides the box if Newton step goes outside.
     *
     * @param box Initial box containing root
     * @param poly Polynomial to refine (1D)
     * @param config Refinement configuration
     * @param refined_location Output: refined root location
     * @param residual Output: residual at refined location
     * @return True if refinement succeeded
     */
    bool refineRoot1D(
        const SubdivisionBoxResult& box,
        const Polynomial& poly,
        const RefinementConfig& config,
        double& refined_location,
        double& residual) const;

    /**
     * @brief Check if Newton step is valid (stays within box and reduces residual)
     *
     * @param x_old Old position
     * @param x_new New position from Newton step
     * @param lower Box lower bound
     * @param upper Box upper bound
     * @param f_old Old residual
     * @param f_new New residual
     * @return True if step is valid
     */
    bool isValidNewtonStep(
        double x_old, double x_new,
        double lower, double upper,
        double f_old, double f_new) const;
    
    /**
     * @brief Estimate root multiplicity from derivatives
     *
     * Checks derivatives up to max_order and returns the order of the
     * first non-zero derivative. For simple roots, returns 1.
     * For multiple roots, returns the multiplicity estimate.
     *
     * @param point Point to check (in [0,1]^n)
     * @param system Original polynomial system
     * @param max_order Maximum order to check
     * @return Estimated multiplicity
     */
    unsigned int estimateMultiplicity(
        const std::vector<double>& point,
        const PolynomialSystem& system,
        unsigned int max_order) const;
    
    /**
     * @brief Compute exclusion radius based on multiplicity
     *
     * For simple roots: radius = multiplier * tolerance
     * For multiple roots: radius = multiplier * tolerance^(1/m)
     *
     * @param multiplicity Root multiplicity
     * @param tolerance Verification tolerance
     * @param multiplier Scaling factor
     * @return Exclusion radius
     */
    double computeExclusionRadius(
        unsigned int multiplicity,
        double tolerance,
        double multiplier) const;
    
    /**
     * @brief Compute Euclidean distance between two points
     *
     * @param p1 First point
     * @param p2 Second point
     * @return Euclidean distance
     */
    double computeDistance(
        const std::vector<double>& p1,
        const std::vector<double>& p2) const;
};

} // namespace polynomial_solver

#endif // RESULT_REFINER_H

