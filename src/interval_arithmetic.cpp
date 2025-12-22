#include "interval_arithmetic.h"

#ifdef ENABLE_HIGH_PRECISION

namespace polynomial_solver {

Interval Interval::operator*(const Interval& other) const {
    // [a,b] * [c,d] = [min(ac,ad,bc,bd), max(ac,ad,bc,bd)]
    mpreal ac = lower * other.lower;
    mpreal ad = lower * other.upper;
    mpreal bc = upper * other.lower;
    mpreal bd = upper * other.upper;

    mpreal lo = min(min(ac, ad), min(bc, bd));
    mpreal hi = max(max(ac, ad), max(bc, bd));

    return Interval(lo, hi);
}

Interval Interval::operator/(const Interval& other) const {
    // Division requires other not containing zero
    // [a,b] / [c,d] where c,d have same sign
    // = [a,b] * [1/d, 1/c] if c,d > 0
    // = [a,b] * [1/d, 1/c] if c,d < 0 (note: 1/d > 1/c when d < c < 0)

    if (other.containsZero()) {
        // Return infinite interval if dividing by interval containing zero
        // This is a conservative approach
        mpreal inf = mpreal("1e308");
        return Interval(-inf, inf);
    }

    mpreal inv_lo, inv_hi;
    if (other.lower > mpreal(0)) {
        // Both positive: 1/d <= 1/x <= 1/c
        inv_lo = mpreal(1) / other.upper;
        inv_hi = mpreal(1) / other.lower;
    } else {
        // Both negative: 1/c <= 1/x <= 1/d
        inv_lo = mpreal(1) / other.upper;
        inv_hi = mpreal(1) / other.lower;
    }

    return (*this) * Interval(inv_lo, inv_hi);
}

Interval intervalPow(const Interval& base, unsigned int n) {
    if (n == 0) {
        // x^0 = 1 for all x
        return Interval(mpreal(1), mpreal(1));
    }

    if (n == 1) {
        return base;
    }

    mpreal a = base.lower;
    mpreal b = base.upper;

    // Compute a^n and b^n
    mpreal a_n = pow(a, static_cast<unsigned long>(n));
    mpreal b_n = pow(b, static_cast<unsigned long>(n));

    if (n % 2 == 1) {
        // Odd power: monotonically increasing
        // Range is [a^n, b^n]
        return Interval(a_n, b_n);
    } else {
        // Even power: x^n is U-shaped with minimum at x=0
        if (base.containsZero()) {
            // 0 is in [a,b], so minimum is 0
            // Maximum is max(|a|^n, |b|^n) = max(a^n, b^n) since n is even
            mpreal hi = max(a_n, b_n);
            return Interval(mpreal(0), hi);
        } else if (a >= mpreal(0)) {
            // Interval is entirely non-negative: monotonically increasing
            return Interval(a_n, b_n);
        } else {
            // Interval is entirely negative: monotonically decreasing
            // a < b < 0, so a^n > b^n
            return Interval(b_n, a_n);
        }
    }
}

Interval evaluateOnInterval(const PolynomialHP& poly, const Interval& interval) {
    // Evaluate polynomial p(x) = sum_{k=0}^n c_k * x^k on interval
    // by computing interval for each monomial and summing

    if (poly.dimension() != 1) {
        // Only support univariate polynomials for now
        return Interval(mpreal(0), mpreal(0));
    }

    const std::vector<mpreal>& coeffs = poly.powerCoefficients();
    if (coeffs.empty()) {
        return Interval(mpreal(0), mpreal(0));
    }

    // Start with constant term (c_0 * x^0 = c_0)
    Interval result(coeffs[0], coeffs[0]);

    // Add each monomial c_k * x^k
    for (size_t k = 1; k < coeffs.size(); ++k) {
        Interval x_k = intervalPow(interval, static_cast<unsigned int>(k));
        Interval term = coeffs[k] * x_k;
        result += term;
    }

    return result;
}

mpreal minAbsOnInterval(const PolynomialHP& poly, const Interval& interval) {
    // Compute min{|p(x)| : x ∈ interval}
    //
    // Strategy:
    // 1. Evaluate polynomial on interval to get [lo, hi]
    // 2. If 0 ∈ [lo, hi], then min|p(x)| = 0
    // 3. Otherwise, min|p(x)| = min(|lo|, |hi|)

    Interval p_interval = evaluateOnInterval(poly, interval);

    if (p_interval.containsZero()) {
        return mpreal(0);
    }

    // Both bounds have same sign, return the one closer to zero
    return min(abs(p_interval.lower), abs(p_interval.upper));
}

Interval intersect(const Interval& a, const Interval& b) {
    mpreal lo = max(a.lower, b.lower);
    mpreal hi = min(a.upper, b.upper);

    if (lo > hi) {
        // Empty intersection - return a degenerate interval
        // (caller should check for this)
        return Interval(mpreal(1), mpreal(0));  // lower > upper indicates empty
    }
    return Interval(lo, hi);
}

bool strictlyContainedIn(const Interval& a, const Interval& b) {
    // a ⊂ b strictly means b.lower < a.lower and a.upper < b.upper
    return b.lower < a.lower && a.upper < b.upper;
}

IntervalNewtonResult intervalNewton(
    const PolynomialHP& poly,
    const Interval& initial_interval,
    const PolynomialHP& derivative,
    unsigned int max_iterations,
    const mpreal& tolerance)
{
    IntervalNewtonResult result;
    result.enclosure = initial_interval;
    result.unique_root = false;
    result.converged = false;
    result.iterations = 0;

    Interval X = initial_interval;

    for (unsigned int iter = 0; iter < max_iterations; ++iter) {
        result.iterations = iter + 1;

        // Compute midpoint
        mpreal x_mid = X.midpoint();

        // Evaluate f at midpoint (point evaluation)
        mpreal f_mid = poly.evaluate(x_mid);

        // Evaluate f' on the interval
        Interval df_interval = evaluateOnInterval(derivative, X);

        // Check if derivative interval contains zero
        if (df_interval.containsZero()) {
            // Cannot proceed - derivative might be zero
            // This could indicate multiple root or need to subdivide
            result.status = "Derivative interval contains zero";
            result.enclosure = X;
            result.root = x_mid;
            result.error_bound = X.width() / mpreal(2);
            return result;
        }

        // Compute Newton operator: N(X) = x_mid - f(x_mid) / F'(X)
        Interval f_mid_interval(f_mid, f_mid);
        Interval newton_correction = f_mid_interval / df_interval;

        // N(X) = x_mid - correction
        // Since x_mid is a point and correction is an interval:
        // N(X) = [x_mid - correction.upper, x_mid - correction.lower]
        Interval N_X(x_mid - newton_correction.upper, x_mid - newton_correction.lower);

        // Intersect with current interval
        Interval X_new = intersect(X, N_X);

        // Check for empty intersection (no root in X)
        if (X_new.lower > X_new.upper) {
            result.status = "Empty intersection - no root in interval";
            result.enclosure = X;
            result.root = x_mid;
            result.error_bound = X.width() / mpreal(2);
            result.converged = false;
            return result;
        }

        // Check for strict containment (unique root verified)
        if (!result.unique_root && strictlyContainedIn(N_X, X)) {
            result.unique_root = true;
        }

        // Check convergence
        mpreal new_width = X_new.width();
        if (new_width < tolerance) {
            result.converged = true;
            result.enclosure = X_new;
            result.root = X_new.midpoint();
            result.error_bound = new_width / mpreal(2);
            result.status = "Converged";
            return result;
        }

        // Check for stagnation (no progress)
        mpreal old_width = X.width();
        if (new_width >= old_width * mpreal("0.99")) {
            // Very slow convergence - might be multiple root
            result.status = "Slow convergence - possible multiple root";
            result.enclosure = X_new;
            result.root = X_new.midpoint();
            result.error_bound = new_width / mpreal(2);
            return result;
        }

        X = X_new;
    }

    // Max iterations reached
    result.status = "Max iterations reached";
    result.enclosure = X;
    result.root = X.midpoint();
    result.error_bound = X.width() / mpreal(2);
    return result;
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

