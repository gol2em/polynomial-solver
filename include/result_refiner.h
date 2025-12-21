#ifndef RESULT_REFINER_H
#define RESULT_REFINER_H

/**
 * @file result_refiner.h
 * @brief Post-processing tool for refining and consolidating solver results
 *
 * This module provides tools to:
 * - Verify roots at high precision (1e-15) using condition-aware convergence
 * - Estimate root multiplicity from derivatives
 * - Eliminate duplicate/nearby boxes
 * - Consolidate results into unique roots
 * - Detect ill-conditioned problems requiring higher precision
 *
 * ## Condition-Aware Convergence
 *
 * The refiner uses a robust convergence criterion that checks BOTH residual and
 * estimated error. This prevents accepting inaccurate roots for ill-conditioned problems.
 *
 * Traditional convergence: |f(x)| < residual_tolerance
 * Problem: For ill-conditioned problems, small residual ≠ small error
 *
 * Condition-aware convergence:
 *   1. Check if |f(x)| < residual_tolerance
 *   2. Estimate condition number: κ ≈ |f''(x)| / |f'(x)|²
 *   3. Estimate actual error: error ≈ κ × |f(x)| / |f'(x)|
 *   4. Accept root only if estimated_error < target_tolerance
 *
 * This ensures roots are rejected when residual is small but error is large,
 * which occurs for ill-conditioned problems requiring higher precision arithmetic.
 */

#include <vector>
#include <cstddef>
#include <functional>
#include "polynomial.h"
#include "solver.h"

namespace polynomial_solver {

/**
 * @brief Configuration for result refinement
 */
struct RefinementConfig {
    double target_tolerance;        ///< Target precision for both error estimation and exclusion radius (default: 1e-15)
    double residual_tolerance;      ///< Max residual |f(x)| threshold for convergence check (default: 1e-15)
    unsigned int max_newton_iters;  ///< Maximum Newton iterations (default: 50)
    unsigned int max_multiplicity;  ///< Maximum multiplicity to check (default: 10)
    double exclusion_multiplier;    ///< Multiplier for exclusion radius (default: 3.0)

    RefinementConfig()
        : target_tolerance(1e-15),
          residual_tolerance(1e-15),
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
    double first_nonzero_derivative;   ///< Value of first non-zero derivative
    double exclusion_radius;           ///< Computed exclusion radius based on multiplicity and derivative
    std::vector<std::size_t> source_boxes;  ///< Indices of original boxes merged into this root
    bool verified;                     ///< True if passed high-precision verification
    unsigned int depth;                ///< Subdivision depth of primary source box
    double condition_estimate;         ///< Estimated condition number (error/residual ratio)
    bool needs_higher_precision;       ///< True if condition number suggests double precision is insufficient
};

/**
 * @brief A problematic region formed by merging nearby unverified boxes
 *
 * When refinement fails for a cluster of adjacent boxes (e.g., due to ill-conditioning
 * or multiple roots), they are merged into a problematic region. The refiner attempts
 * to resolve the region using modified Newton with multiplicity estimation.
 */
struct ProblematicRegion {
    double lower;                      ///< Lower bound of the region (1D only)
    double upper;                      ///< Upper bound of the region (1D only)
    std::vector<std::size_t> box_indices;  ///< Indices of boxes merged into this region

    // Refinement results (if attempted)
    bool refinement_attempted;         ///< True if refinement was attempted
    bool refinement_succeeded;         ///< True if refinement succeeded with condition-aware convergence
    double refined_root;               ///< Refined root location (if successful)
    unsigned int multiplicity;         ///< Estimated multiplicity at refined root
    double residual;                   ///< Residual |f(x)| at refined root
    double condition_estimate;         ///< Estimated condition number
    bool needs_higher_precision;       ///< True if condition suggests double precision insufficient

    ProblematicRegion()
        : lower(0.0), upper(0.0),
          refinement_attempted(false), refinement_succeeded(false),
          refined_root(0.0), multiplicity(0), residual(0.0),
          condition_estimate(1.0), needs_higher_precision(false)
    {}
};

/**
 * @brief Result of refinement process
 */
struct RefinementResult {
    std::vector<RefinedRoot> roots;           ///< Verified and consolidated roots
    std::vector<std::size_t> cancelled_boxes; ///< Indices of boxes eliminated as duplicates
    std::vector<std::size_t> unverified_boxes; ///< Indices of boxes that didn't pass verification (kept for backward compatibility)
    std::vector<ProblematicRegion> problematic_regions; ///< Merged regions from unverified boxes (1D only)
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
     * 1. For each resolved box, use Newton's method to refine until |f(x)| < residual_tolerance
     * 2. Use subdivision with sign checking to ensure convergence
     * 3. Verify convergence by checking residual |f(x)| < residual_tolerance
     * 4. Estimate multiplicity from derivatives
     * 5. Compute exclusion radius using target_tolerance and derivatives
     * 6. Cancel nearby boxes within exclusion radius
     * 7. Returns consolidated unique roots
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

    /**
     * @brief Verify root and determine multiplicity from derivatives (Taylor method)
     *
     * Checks derivatives up to max_order and returns the order of the
     * first non-zero derivative. For simple roots, returns 1.
     * For multiple roots, returns the multiplicity estimate.
     *
     * For a root at x with multiplicity m:
     * - f(x) = 0 (already verified by Newton refinement)
     * - f'(x) = 0, f''(x) = 0, ..., f^(m-1)(x) = 0
     * - f^(m)(x) ≠ 0 (first non-zero derivative)
     *
     * This provides a rigorous verification of the root's multiplicity.
     *
     * @param point Point to check (in [0,1]^n)
     * @param system Original polynomial system
     * @param max_order Maximum order to check
     * @param derivative_threshold Threshold for considering derivative as zero (default: 1e-10)
     * @param first_nonzero_deriv Output: value of first non-zero derivative
     * @return Estimated multiplicity (1 = simple root, 2+ = multiple root)
     */
    unsigned int estimateMultiplicity(
        const std::vector<double>& point,
        const PolynomialSystem& system,
        unsigned int max_order,
        double derivative_threshold,
        double& first_nonzero_deriv) const;

    /**
     * @brief Estimate multiplicity using Ostrowski's method from Newton iterates
     *
     * Uses 3 consecutive Newton iterates to estimate multiplicity via:
     *   p = 1/2 + (x₁ - x₂) / (x₃ - 2x₂ + x₁)
     *   m = floor(p)
     *
     * This method is more robust to distance from root and multiple roots
     * compared to the Taylor derivative method.
     *
     * @param x1 First Newton iterate
     * @param x2 Second Newton iterate
     * @param x3 Third Newton iterate
     * @return Estimated multiplicity (1 = simple root, 2+ = multiple root)
     */
    unsigned int estimateMultiplicityOstrowski(
        double x1, double x2, double x3) const;

    /**
     * @brief Estimate multiplicity using Ostrowski's method from a starting point
     *
     * Performs 3 Newton iterations from x0, then uses Ostrowski's formula.
     * This is the recommended method for multiplicity detection in double precision.
     *
     * @param x0 Starting point
     * @param poly Polynomial (1D)
     * @return Estimated multiplicity (1 = simple root, 2+ = multiple root)
     */
    unsigned int estimateMultiplicityOstrowskiFromPoint(
        double x0, const Polynomial& poly) const;

    /**
     * @brief Refine a 1D root from a starting point using Newton's method
     *
     * Uses the new Ostrowski-based workflow:
     * 1. Do 3 standard Newton iterations
     * 2. Estimate multiplicity using Ostrowski
     * 3. Use modified Newton with detected multiplicity
     *
     * @param x0 Initial guess
     * @param poly Polynomial (1D)
     * @param config Refinement configuration
     * @param refined_location Output: refined root location
     * @param residual Output: residual at refined location
     * @return True if refinement succeeded
     */
    bool refineRoot1D_fromPoint(
        double x0,
        const Polynomial& poly,
        const RefinementConfig& config,
        double& refined_location,
        double& residual) const;

    /**
     * @brief Refine a 1D root with known multiplicity using modified Newton method
     *
     * For a root with multiplicity m, uses the modified Newton iteration:
     *   x_{n+1} = x_n - m * f(x_n) / f'(x_n)
     *
     * This converges quadratically even for multiple roots, unlike standard Newton
     * which only converges linearly for multiple roots.
     *
     * @param initial_guess Initial guess for root location
     * @param lower Lower bound of search interval
     * @param upper Upper bound of search interval
     * @param poly Polynomial to refine (1D)
     * @param multiplicity Known or estimated multiplicity
     * @param config Refinement configuration
     * @param refined_location Output: refined root location
     * @param residual Output: residual at refined location
     * @return True if refinement succeeded
     */
    bool refineRoot1DWithMultiplicity(
        double initial_guess,
        double lower,
        double upper,
        const Polynomial& poly,
        unsigned int multiplicity,
        const RefinementConfig& config,
        double& refined_location,
        double& residual) const;

    /**
     * @brief Estimate condition number for a 1D root
     *
     * The condition number κ relates the error in the root to the residual:
     *   |x - x_true| ≈ κ * |f(x)| / |f'(x_true)|
     *
     * For well-conditioned problems: κ ≈ 1
     * For ill-conditioned problems: κ >> 1
     *
     * We estimate κ by computing:
     *   κ ≈ ||f|| / (|f'(x)| * spacing)
     * where ||f|| is a norm of the polynomial coefficients and spacing is
     * the typical distance between roots.
     *
     * A simpler practical estimate uses the ratio of higher derivatives:
     *   κ ≈ |f''(x)| / |f'(x)|^2 * typical_scale
     *
     * @param location Root location
     * @param poly Polynomial (1D)
     * @param derivative_value Value of f'(x) at the root
     * @return Estimated condition number
     */
    double estimateConditionNumber1D(
        double location,
        const Polynomial& poly,
        double derivative_value) const;

#ifdef ENABLE_HIGH_PRECISION
    /**
     * @brief Refine a 1D root with automatic precision escalation
     *
     * This method first attempts refinement with double precision.
     * If the condition number is too large (estimated error > 1e-10),
     * it automatically switches to high precision with appropriate bit depth.
     *
     * Precision selection based on condition number:
     * - κ < 1e5:  256 bits (~77 decimal digits)
     * - κ < 1e10: 512 bits (~154 decimal digits)
     * - κ ≥ 1e10: 1024 bits (~308 decimal digits)
     *
     * @param initial_guess Initial guess for root location
     * @param poly Polynomial to refine (1D, double precision)
     * @param config Double precision refinement configuration
     * @param refined_root Output: refined root information
     * @return True if refinement succeeded (either in double or high precision)
     */
    bool refineRoot1DWithPrecisionEscalation(
        double initial_guess,
        const Polynomial& poly,
        const RefinementConfig& config,
        RefinedRoot& refined_root) const;
#endif

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
     * @brief Compute exclusion radius based on multiplicity and derivative value
     *
     * For a root with multiplicity m and first non-zero derivative D:
     * The exclusion radius is estimated as: r ≈ (tolerance * m! / |D|)^(1/m)
     *
     * For simple roots (m=1): r = multiplier * tolerance / |D|
     * For multiple roots: r = multiplier * (tolerance / |D|)^(1/m)
     *
     * @param multiplicity Root multiplicity
     * @param first_nonzero_deriv Value of first non-zero derivative
     * @param tolerance Verification tolerance
     * @param multiplier Scaling factor
     * @return Exclusion radius
     */
    double computeExclusionRadiusFromDerivative(
        unsigned int multiplicity,
        double first_nonzero_deriv,
        double tolerance,
        double multiplier) const;

    /**
     * @brief Compute exclusion radius based on multiplicity (legacy)
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

    /**
     * @brief Merge nearby unverified boxes into problematic regions (1D only)
     *
     * Adjacent boxes from subdivision that failed refinement are merged into
     * continuous regions. Each region is then refined using modified Newton
     * with multiplicity estimation.
     *
     * @param unverified_indices Indices of unverified boxes
     * @param solver_result Original solver result containing boxes
     * @param poly Polynomial (1D)
     * @param config Refinement configuration
     * @return Vector of problematic regions with refinement results
     */
    std::vector<ProblematicRegion> mergeUnverifiedBoxes1D(
        const std::vector<std::size_t>& unverified_indices,
        const SubdivisionSolverResult& solver_result,
        const Polynomial& poly,
        const RefinementConfig& config) const;

    /**
     * @brief Attempt to refine a problematic region
     *
     * Tries Newton refinement from the center of the region with multiplicity
     * estimation and condition-aware convergence checking.
     *
     * @param region Region to refine (modified in place with results)
     * @param poly Polynomial (1D)
     * @param config Refinement configuration
     */
    void refineProblematicRegion1D(
        ProblematicRegion& region,
        const Polynomial& poly,
        const RefinementConfig& config) const;

};

//=============================================================================
// CurveRefiner - Refinement for curves (1 equation in 2D)
//=============================================================================

/**
 * @brief Result of curve point refinement
 */
struct CurveRefinedPoint {
    double x;                    ///< Refined x coordinate
    double y;                    ///< Refined y coordinate
    double residual;             ///< |g(x,y)| at refined point
    bool converged;              ///< True if refinement converged
    unsigned int iterations;     ///< Number of iterations used
};

/**
 * @brief Configuration for curve refinement
 */
struct CurveRefinementConfig {
    double residual_tolerance;      ///< Target |g(x,y)| tolerance (default: 1e-14)
    unsigned int max_iterations;    ///< Maximum iterations (default: 50)
    double min_gradient_norm;       ///< Minimum |∇g| to avoid singular points (default: 1e-20)

    CurveRefinementConfig()
        : residual_tolerance(1e-14)
        , max_iterations(50)
        , min_gradient_norm(1e-20)
    {}
};

/**
 * @brief Refiner for projecting points onto curves defined by g(x,y) = 0
 *
 * For a curve defined by a bivariate polynomial g(x,y) = 0, this class
 * provides methods to project initial guesses onto the curve using
 * gradient projection (Newton's method along the gradient direction).
 *
 * The algorithm:
 *   1. Given initial point (x₀, y₀) near the curve
 *   2. Evaluate g and ∇g = (gₓ, gᵧ) at current point
 *   3. Move along gradient: (x,y) ← (x,y) - g/|∇g|² × ∇g
 *   4. Repeat until |g| < tolerance
 *
 * This finds the point on the curve closest to the initial guess
 * (in the gradient direction), which is the orthogonal projection
 * onto the curve for small displacements.
 *
 * Usage:
 * @code
 * Polynomial g = ...;  // 2D polynomial defining the curve g(x,y) = 0
 * CurveRefiner refiner(g);
 *
 * double x = 0.5, y = 0.5;  // Initial guess
 * CurveRefinedPoint result = refiner.refine(x, y);
 *
 * if (result.converged) {
 *     std::cout << "Point on curve: (" << result.x << ", " << result.y << ")\n";
 * }
 * @endcode
 */
class CurveRefiner {
public:
    /**
     * @brief Construct a curve refiner for a given polynomial
     *
     * @param curve_poly Bivariate polynomial g(x,y) defining the curve g = 0
     * @throws std::invalid_argument if polynomial is not 2-dimensional
     */
    explicit CurveRefiner(const Polynomial& curve_poly);

    /**
     * @brief Refine a point onto the curve
     *
     * Projects the initial guess onto the curve g(x,y) = 0 using
     * gradient projection.
     *
     * @param x0 Initial x coordinate
     * @param y0 Initial y coordinate
     * @param config Refinement configuration
     * @return Refinement result with converged point
     */
    CurveRefinedPoint refine(
        double x0, double y0,
        const CurveRefinementConfig& config = CurveRefinementConfig()) const;

    /**
     * @brief Refine multiple points onto the curve
     *
     * @param points Vector of (x, y) pairs to refine
     * @param config Refinement configuration
     * @return Vector of refinement results
     */
    std::vector<CurveRefinedPoint> refineMultiple(
        const std::vector<std::pair<double, double>>& points,
        const CurveRefinementConfig& config = CurveRefinementConfig()) const;

    /**
     * @brief Get the curve polynomial
     */
    const Polynomial& polynomial() const { return curve_poly_; }

private:
    Polynomial curve_poly_;      ///< The curve polynomial g(x,y)
    Polynomial dg_dx_;           ///< Partial derivative ∂g/∂x
    Polynomial dg_dy_;           ///< Partial derivative ∂g/∂y
};

//=============================================================================
// Function-based curve refinement (uses numerical gradient)
//=============================================================================

/**
 * @brief Refine a point onto a curve using only function evaluations
 *
 * Uses Newton's method with numerical gradient (central differences) to
 * project a point onto the zero set of g(x,y) = 0.
 *
 * This is useful when only function evaluations are available (no analytic
 * derivatives or polynomial representation).
 *
 * @param g Function g(x,y) defining the curve g = 0
 * @param x0 Initial x coordinate
 * @param y0 Initial y coordinate
 * @param config Refinement configuration
 * @return Refinement result with converged point
 *
 * Usage:
 * @code
 * auto g = [](double x, double y) { return x*x + y*y - 1.0; };  // Unit circle
 * CurveRefinedPoint result = refineCurveNumerical(g, 0.8, 0.7);
 * if (result.converged) {
 *     // result.x, result.y is on the circle
 * }
 * @endcode
 */
CurveRefinedPoint refineCurveNumerical(
    const std::function<double(double, double)>& g,
    double x0, double y0,
    const CurveRefinementConfig& config = CurveRefinementConfig());

/**
 * @brief Refine multiple points onto a curve using only function evaluations
 *
 * @param g Function g(x,y) defining the curve g = 0
 * @param points Vector of (x, y) pairs to refine
 * @param config Refinement configuration
 * @return Vector of refinement results
 */
std::vector<CurveRefinedPoint> refineCurveNumericalMultiple(
    const std::function<double(double, double)>& g,
    const std::vector<std::pair<double, double>>& points,
    const CurveRefinementConfig& config = CurveRefinementConfig());

} // namespace polynomial_solver

#endif // RESULT_REFINER_H

