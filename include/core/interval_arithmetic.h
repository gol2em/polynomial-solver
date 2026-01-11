#ifndef INTERVAL_ARITHMETIC_H
#define INTERVAL_ARITHMETIC_H

/**
 * @file interval_arithmetic.h
 * @brief Interval arithmetic for rigorous polynomial evaluation bounds
 *
 * This module provides interval arithmetic for high-precision polynomial
 * evaluation. Key features:
 * - Rigorous enclosure of polynomial values over an interval
 * - Power basis evaluation via monomial range analysis
 * - Used by interval Newton method for guaranteed error bounds
 *
 * ## Algorithm for Monomial Evaluation
 *
 * For monomial c*x^n on interval [a,b]:
 * - Odd n: monotonic, range is [c*a^n, c*b^n] (swapped if c < 0)
 * - Even n: 
 *   - If 0 ∈ [a,b]: minimum at 0, range includes 0
 *   - If 0 ∉ [a,b]: monotonic on interval
 *
 * Only available when ENABLE_HIGH_PRECISION is defined.
 */

#ifdef ENABLE_HIGH_PRECISION

#include "hp/high_precision_types.h"
#include "hp/polynomial_hp.h"
#include <vector>
#include <algorithm>

namespace polynomial_solver {

/**
 * @class Interval
 * @brief Represents a closed interval [lower, upper] with rigorous arithmetic
 */
class Interval {
public:
    mpreal lower;  ///< Lower bound of interval
    mpreal upper;  ///< Upper bound of interval

    /**
     * @brief Default constructor: [0, 0]
     */
    Interval() : lower(0), upper(0) {}

    /**
     * @brief Construct interval from bounds
     * @param lo Lower bound
     * @param hi Upper bound (must be >= lo)
     */
    Interval(const mpreal& lo, const mpreal& hi) : lower(lo), upper(hi) {}

    /**
     * @brief Construct point interval [x, x]
     */
    explicit Interval(const mpreal& x) : lower(x), upper(x) {}

    /**
     * @brief Interval width (upper - lower)
     */
    mpreal width() const { return upper - lower; }

    /**
     * @brief Interval midpoint
     */
    mpreal midpoint() const { return (lower + upper) / mpreal(2); }

    /**
     * @brief Check if interval contains a point
     */
    bool contains(const mpreal& x) const { return x >= lower && x <= upper; }

    /**
     * @brief Check if interval contains zero
     */
    bool containsZero() const { return lower <= mpreal(0) && upper >= mpreal(0); }

    /**
     * @brief Interval addition: [a,b] + [c,d] = [a+c, b+d]
     */
    Interval operator+(const Interval& other) const {
        return Interval(lower + other.lower, upper + other.upper);
    }

    /**
     * @brief Interval subtraction: [a,b] - [c,d] = [a-d, b-c]
     */
    Interval operator-(const Interval& other) const {
        return Interval(lower - other.upper, upper - other.lower);
    }

    /**
     * @brief Interval multiplication
     * [a,b] * [c,d] = [min(ac,ad,bc,bd), max(ac,ad,bc,bd)]
     */
    Interval operator*(const Interval& other) const;

    /**
     * @brief Scalar multiplication: c * [a,b]
     */
    Interval operator*(const mpreal& scalar) const {
        if (scalar >= mpreal(0)) {
            return Interval(lower * scalar, upper * scalar);
        } else {
            return Interval(upper * scalar, lower * scalar);
        }
    }

    /**
     * @brief Interval division (other must not contain zero)
     */
    Interval operator/(const Interval& other) const;

    /**
     * @brief Compound addition
     */
    Interval& operator+=(const Interval& other) {
        lower += other.lower;
        upper += other.upper;
        return *this;
    }
};

// Scalar * Interval
inline Interval operator*(const mpreal& scalar, const Interval& interval) {
    return interval * scalar;
}

/**
 * @brief Compute interval power x^n for interval x
 *
 * For x^n where x = [a,b]:
 * - Odd n: monotonic, range = [a^n, b^n]
 * - Even n: check if 0 ∈ [a,b]
 *
 * @param base Interval base
 * @param n Non-negative integer exponent
 * @return Interval enclosing all x^n for x in base
 */
Interval intervalPow(const Interval& base, unsigned int n);

/**
 * @brief Evaluate univariate polynomial on interval (power basis)
 *
 * Computes rigorous enclosure of {p(x) : x ∈ interval} by evaluating
 * each monomial c_k * x^k on the interval and summing.
 *
 * @param poly Polynomial (must be 1D, uses power coefficients)
 * @param interval Input interval
 * @return Interval enclosing all polynomial values on input interval
 */
Interval evaluateOnInterval(const PolynomialHP& poly, const Interval& interval);

/**
 * @brief Compute minimum absolute value of polynomial on interval
 *
 * Returns min{|p(x)| : x ∈ interval}. Used for bounding f' away from zero
 * in interval Newton method.
 *
 * @param poly Polynomial
 * @param interval Input interval
 * @return Lower bound on |p(x)| over interval (may be 0)
 */
mpreal minAbsOnInterval(const PolynomialHP& poly, const Interval& interval);

/**
 * @brief Result of interval Newton iteration
 */
struct IntervalNewtonResult {
    Interval enclosure;         ///< Final interval enclosing the root
    mpreal root;                ///< Best estimate of root (midpoint of final interval)
    mpreal error_bound;         ///< Guaranteed error bound (half-width of final interval)
    bool unique_root;           ///< True if unique root existence is verified
    bool converged;             ///< True if iteration converged
    unsigned int iterations;    ///< Number of iterations performed
    std::string status;         ///< Status message
};

/**
 * @brief Perform interval Newton iteration for root refinement
 *
 * The interval Newton operator is:
 *   N(X) = x* - f(x*) / F'(X)
 * where x* is the midpoint of X and F'(X) is the interval evaluation of f'.
 *
 * Properties:
 * - If N(X) ∩ X is empty, there is no root in X
 * - If N(X) ⊂ X (strictly), there is exactly one root in X
 * - The intersection X ∩ N(X) contains all roots in X
 *
 * @param poly Polynomial (1D, uses power coefficients)
 * @param initial_interval Starting interval (should contain a root)
 * @param derivative Derivative polynomial
 * @param max_iterations Maximum number of iterations
 * @param tolerance Stop when interval width < tolerance
 * @return IntervalNewtonResult with enclosure, root estimate, and error bound
 */
IntervalNewtonResult intervalNewton(
    const PolynomialHP& poly,
    const Interval& initial_interval,
    const PolynomialHP& derivative,
    unsigned int max_iterations = 100,
    const mpreal& tolerance = mpreal("1e-100"));

/**
 * @brief Compute intersection of two intervals
 * @return Intersection interval, or empty interval if no intersection
 */
Interval intersect(const Interval& a, const Interval& b);

/**
 * @brief Check if interval a is strictly contained in interval b
 */
bool strictlyContainedIn(const Interval& a, const Interval& b);

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

#endif // INTERVAL_ARITHMETIC_H

