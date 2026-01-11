# Subdivision Design for Both Representations

## Core Principle

**Bernstein basis**: Modify coefficients (De Casteljau)
**Power basis**: Transform coordinates (keep coefficients unchanged)

---

## Bernstein Subdivision (Coefficient-Based)

### Algorithm
For polynomial in Bernstein basis on [0,1]^n, restrict to [a,b] along axis i:
- Use De Casteljau algorithm to compute new Bernstein coefficients
- Result is polynomial on [0,1]^n representing restricted region

### Properties
- ✅ Numerically stable (convex combinations only)
- ✅ Natural for Bernstein basis
- ✅ Error grows slowly: O(ε · depth)
- ✅ Coefficients have geometric meaning (control points)

### Implementation
```cpp
PolynomialBase restrictedToInterval(size_t axis, Scalar a, Scalar b) const {
    assert(basis_ == Basis::BERNSTEIN);
    std::vector<Scalar> new_coeffs = de_casteljau_restrict(degrees_, coeffs_, axis, a, b);
    return PolynomialBase::fromBernstein(degrees_, new_coeffs);
}
```

---

## Power Subdivision (Coordinate-Based)

### Algorithm
For polynomial in power basis, DON'T modify coefficients. Instead:
- Store region bounds [lower, upper] ∈ R^n alongside polynomial
- Transform evaluation points: x = lower + (upper - lower) ⊙ u
- Evaluate using Horner's method at transformed point

### Properties
- ✅ **Minimal rounding error**: Only coordinate transformation
- ✅ **Preserves exact coefficients**: No binomial expansion
- ✅ **Composable**: Subdivisions compose as affine transformations
- ✅ **Error grows linearly**: O(ε · depth) vs O(ε · degree² · 2^depth)

### Data Structure
```cpp
template<typename Scalar>
struct PowerPolynomialWithRegion {
    std::vector<Scalar> coeffs;           // Power basis coefficients
    std::vector<unsigned int> degrees;    // Degrees per dimension
    std::vector<Scalar> region_lower;     // Lower bounds per dimension
    std::vector<Scalar> region_upper;     // Upper bounds per dimension
    
    // Evaluate at u ∈ [0,1]^n
    Scalar evaluate(const std::vector<Scalar>& u) const {
        // Transform to actual region
        std::vector<Scalar> x(u.size());
        for (size_t i = 0; i < u.size(); ++i) {
            x[i] = region_lower[i] + (region_upper[i] - region_lower[i]) * u[i];
        }
        // Evaluate using Horner
        return horner_eval(degrees, coeffs, x);
    }
    
    // Subdivide along axis
    PowerPolynomialWithRegion restrictedToInterval(size_t axis, Scalar a, Scalar b) const {
        PowerPolynomialWithRegion result = *this;
        // Update region bounds (coefficients unchanged!)
        Scalar width = region_upper[axis] - region_lower[axis];
        result.region_lower[axis] = region_lower[axis] + a * width;
        result.region_upper[axis] = region_lower[axis] + b * width;
        return result;
    }
};
```

---

## Unified Interface

Both representations should support the same operations:

```cpp
enum class Basis { POWER, BERNSTEIN };

template<typename Scalar, typename ErrorScalar = double>
class PolynomialBase {
private:
    std::vector<Scalar> coeffs_;
    Basis basis_;
    std::vector<unsigned int> degrees_;
    
    // For POWER basis: region bounds (default [0,1]^n)
    std::vector<Scalar> region_lower_;
    std::vector<Scalar> region_upper_;
    
    ErrorBound<ErrorScalar> error_bound_;
    
public:
    // Subdivision (works for both bases)
    PolynomialBase restrictedToInterval(size_t axis, Scalar a, Scalar b) const {
        if (basis_ == Basis::BERNSTEIN) {
            // Modify coefficients using De Casteljau
            auto new_coeffs = de_casteljau_restrict(degrees_, coeffs_, axis, a, b);
            return PolynomialBase::fromBernstein(degrees_, new_coeffs);
        } else {
            // Modify region bounds, keep coefficients
            PolynomialBase result = *this;
            Scalar width = region_upper_[axis] - region_lower_[axis];
            result.region_lower_[axis] = region_lower_[axis] + a * width;
            result.region_upper_[axis] = region_lower_[axis] + b * width;
            return result;
        }
    }
    
    // Evaluation (works for both bases)
    Scalar evaluate(const std::vector<Scalar>& u) const {
        if (basis_ == Basis::BERNSTEIN) {
            // u is already in [0,1]^n
            return de_casteljau_eval(degrees_, coeffs_, u);
        } else {
            // Transform u from [0,1]^n to region
            std::vector<Scalar> x(u.size());
            for (size_t i = 0; i < u.size(); ++i) {
                x[i] = region_lower_[i] + (region_upper_[i] - region_lower_[i]) * u[i];
            }
            return horner_eval(degrees_, coeffs_, x);
        }
    }
    
    // Conversion
    PolynomialBase convertToBernstein() const {
        if (basis_ == Basis::BERNSTEIN) return *this;
        
        // Must first normalize region to [0,1]^n by evaluating at transformed points
        // This is the ONE place where we modify power coefficients
        // Use high-precision arithmetic here!
        auto normalized = normalizeRegionToUnit<HigherPrecision>();
        auto bern_coeffs = power_to_bernstein(normalized.degrees_, normalized.coeffs_);
        return PolynomialBase::fromBernstein(degrees_, bern_coeffs);
    }
};
```

---

## High-Precision Conversion Strategy

### Key Insight
The ONLY numerically sensitive operation is **power-to-Bernstein conversion**.
All other operations preserve precision:
- Bernstein subdivision: stable (convex combinations)
- Power subdivision: stable (coordinate transform only)
- Horner evaluation: stable
- De Casteljau evaluation: stable

### Strategy
**Always perform basis conversion in higher precision:**

```cpp
template<typename Scalar, typename ErrorScalar>
PolynomialBase<Scalar, ErrorScalar> 
PolynomialBase<Scalar, ErrorScalar>::convertToBernstein() const {
    using HighPrecision = /* mpreal or higher precision type */;
    
    // 1. Lift to high precision
    auto hp_poly = liftPrecision<HighPrecision>();
    
    // 2. Convert in high precision (numerically stable)
    auto hp_coeffs = power_to_bernstein_hp(hp_poly.degrees_, hp_poly.coeffs_);
    
    // 3. Round back to target precision
    std::vector<Scalar> bern_coeffs;
    for (const auto& c : hp_coeffs) {
        bern_coeffs.push_back(static_cast<Scalar>(c));
    }
    
    // 4. Track error from rounding
    ErrorScalar conversion_error = estimate_conversion_error(...);
    
    return PolynomialBase::fromBernstein(degrees_, bern_coeffs)
        .withError(error_bound_.absolute_error + conversion_error);
}
```

This ensures conversion error stays near machine epsilon of the target type.

---

## Summary

| Operation | Bernstein | Power | Precision Strategy |
|-----------|-----------|-------|-------------------|
| Subdivision | Modify coeffs (De Casteljau) | Modify region bounds | Both stable |
| Evaluation | De Casteljau on [0,1]^n | Horner on transformed point | Both stable |
| Conversion | N/A | Power→Bernstein in HP | **Use higher precision** |
| Error growth | O(ε·depth) | O(ε·depth) | Linear in both |

**Key takeaway**: Power basis subdivision via coordinate transformation is vastly superior to coefficient transformation.

