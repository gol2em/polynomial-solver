# Unified API Design: Consistent Operations Across Both Bases

## The Challenge

Different operations are natural in different bases:
- **Bernstein basis**: Subdivision (De Casteljau), evaluation, convex hull
- **Power basis**: Differentiation, division with remainder, Sturm sequences, arithmetic

We need a **unified API** where all operations work regardless of current basis.

---

## Core Design Principle: Implicit Normalization

**Key Insight**: Power basis polynomials with region bounds [a,b]^n are NOT in standard form for algebraic operations.

### The Problem

```cpp
// Power polynomial on region [0.3, 0.7]
PolynomialBase p_power = /* ... */;
p_power.region_lower = {0.3};
p_power.region_upper = {0.7};

// What does differentiation mean here?
auto dp = differentiate(p_power, 0);  // ??? 
```

**Issue**: d/dx p(x) where x ∈ [0.3, 0.7] is NOT the same as d/du p(0.3 + 0.4u) where u ∈ [0,1]!

### The Solution: Normalize Before Algebraic Operations

For power basis polynomials on non-standard regions, we must:
1. **Normalize to [0,1]^n** (transform coefficients)
2. **Perform operation** (differentiation, division, etc.)
3. **Result is on [0,1]^n**

---

## Unified API Design

```cpp
enum class Basis { POWER, BERNSTEIN };

template<typename Scalar, typename ErrorScalar = double>
class PolynomialBase {
private:
    std::vector<Scalar> coeffs_;
    Basis basis_;
    std::vector<unsigned int> degrees_;
    std::size_t dimension_;
    
    // Region bounds (only meaningful for POWER basis)
    // For BERNSTEIN: always [0,1]^n (not stored)
    // For POWER: current region (default [0,1]^n)
    std::vector<Scalar> region_lower_;  // Empty if Bernstein or standard region
    std::vector<Scalar> region_upper_;  // Empty if Bernstein or standard region
    
    ErrorBound<ErrorScalar> error_bound_;
    
    // Check if region is standard [0,1]^n
    bool hasStandardRegion() const {
        if (basis_ == Basis::BERNSTEIN) return true;
        if (region_lower_.empty()) return true;  // Default is [0,1]^n
        
        for (size_t i = 0; i < dimension_; ++i) {
            if (region_lower_[i] != Scalar(0) || region_upper_[i] != Scalar(1)) {
                return false;
            }
        }
        return true;
    }
    
public:
    //=========================================================================
    // SUBDIVISION (works for both bases, always returns polynomial on [0,1]^n)
    //=========================================================================
    
    PolynomialBase restrictedToInterval(size_t axis, Scalar a, Scalar b) const {
        if (basis_ == Basis::BERNSTEIN) {
            // Bernstein: modify coefficients using De Casteljau
            auto new_coeffs = de_casteljau_restrict(degrees_, coeffs_, axis, a, b);
            ErrorScalar sub_err = errorFromSubdivision(degrees_, new_coeffs);
            
            PolynomialBase result;
            result.coeffs_ = new_coeffs;
            result.basis_ = Basis::BERNSTEIN;
            result.degrees_ = degrees_;
            result.dimension_ = dimension_;
            result.error_bound_.absolute_error = error_bound_.absolute_error + sub_err;
            return result;
            
        } else {
            // Power: modify region bounds (coefficients unchanged)
            PolynomialBase result = *this;
            
            // Initialize region if not set
            if (result.region_lower_.empty()) {
                result.region_lower_.assign(dimension_, Scalar(0));
                result.region_upper_.assign(dimension_, Scalar(1));
            }
            
            // Update region bounds
            Scalar width = result.region_upper_[axis] - result.region_lower_[axis];
            result.region_lower_[axis] = result.region_lower_[axis] + a * width;
            result.region_upper_[axis] = result.region_lower_[axis] + b * width;
            
            // Error from coordinate transformation (minimal)
            ErrorScalar transform_err = errorFromCoordinateTransform(a, b);
            result.error_bound_.absolute_error = error_bound_.absolute_error + transform_err;
            
            return result;
        }
    }
    
    //=========================================================================
    // EVALUATION (works for both bases, u ∈ [0,1]^n)
    //=========================================================================
    
    Scalar evaluate(const std::vector<Scalar>& u) const {
        if (basis_ == Basis::BERNSTEIN) {
            // u is already in [0,1]^n
            return de_casteljau_eval(degrees_, coeffs_, u);
            
        } else {
            // Transform u from [0,1]^n to actual region
            std::vector<Scalar> x(u.size());
            
            if (region_lower_.empty()) {
                // Standard region [0,1]^n
                x = u;
            } else {
                // Non-standard region
                for (size_t i = 0; i < u.size(); ++i) {
                    x[i] = region_lower_[i] + (region_upper_[i] - region_lower_[i]) * u[i];
                }
            }
            
            return horner_eval(degrees_, coeffs_, x);
        }
    }
    
    //=========================================================================
    // DIFFERENTIATION (works for both bases!)
    //=========================================================================

    PolynomialBase differentiate(size_t axis, unsigned int order = 1) const {
        if (basis_ == Basis::BERNSTEIN) {
            // Bernstein differentiation: d_i = n * (b_{i+1} - b_i)
            // Very simple and numerically stable!
            auto deriv_coeffs = differentiate_bernstein(degrees_, coeffs_, axis, order);

            ErrorScalar deriv_err = errorFromDifferentiation(degrees_, coeffs_);

            PolynomialBase result;
            result.coeffs_ = deriv_coeffs;
            result.basis_ = Basis::BERNSTEIN;
            result.degrees_ = degrees_;  // Degree reduced by order
            result.dimension_ = dimension_;
            result.error_bound_.absolute_error = error_bound_.absolute_error + deriv_err;
            return result;

        } else {
            // Power differentiation: d/dx(a_i * x^i) = i * a_i * x^(i-1)
            // Must be on standard region [0,1]^n
            PolynomialBase p_normalized = ensurePowerStandardRegion();

            auto deriv_coeffs = differentiate_power(p_normalized.degrees_,
                                                     p_normalized.coeffs_,
                                                     axis, order);

            ErrorScalar deriv_err = errorFromDifferentiation(p_normalized.degrees_,
                                                              p_normalized.coeffs_);

            PolynomialBase result;
            result.coeffs_ = deriv_coeffs;
            result.basis_ = Basis::POWER;
            result.degrees_ = p_normalized.degrees_;  // Degree reduced
            result.dimension_ = dimension_;
            result.error_bound_.absolute_error = error_bound_.absolute_error + deriv_err;
            return result;
        }
    }
    
    //=========================================================================
    // DIVISION WITH REMAINDER (requires power basis on [0,1]^n)
    //=========================================================================

    std::pair<PolynomialBase, PolynomialBase>
    divideWithRemainder(const PolynomialBase& divisor) const {
        // Both must be in power basis on standard region
        PolynomialBase dividend_norm = ensurePowerStandardRegion();
        PolynomialBase divisor_norm = divisor.ensurePowerStandardRegion();

        auto [quot_coeffs, rem_coeffs] = divide_power(
            dividend_norm.degrees_, dividend_norm.coeffs_,
            divisor_norm.degrees_, divisor_norm.coeffs_);

        return {
            PolynomialBase::fromPower(quot_degrees, quot_coeffs),
            PolynomialBase::fromPower(rem_degrees, rem_coeffs)
        };
    }

    //=========================================================================
    // STURM SEQUENCE (requires power basis on [0,1]^n)
    //=========================================================================

    std::vector<PolynomialBase> sturmSequence() const {
        // Must be univariate
        if (dimension_ != 1) {
            throw std::invalid_argument("Sturm sequence only for univariate polynomials");
        }

        // Ensure power basis on standard region
        PolynomialBase p0 = ensurePowerStandardRegion();
        PolynomialBase p1 = p0.differentiate(0, 1);

        std::vector<PolynomialBase> sequence;
        sequence.push_back(p0);
        sequence.push_back(p1);

        // Build sequence: p_{i+1} = -rem(p_{i-1}, p_i)
        while (!sequence.back().isZero()) {
            auto [quot, rem] = sequence[sequence.size() - 2].divideWithRemainder(sequence.back());
            if (rem.isZero()) break;
            sequence.push_back(-rem);  // Negate remainder
        }

        return sequence;
    }

    //=========================================================================
    // CONVERSION (always returns polynomial on [0,1]^n)
    //=========================================================================

    PolynomialBase convertToBernstein() const {
        if (basis_ == Basis::BERNSTEIN) return *this;

        // Must normalize to [0,1]^n first
        PolynomialBase p_norm = ensurePowerStandardRegion();

        // Convert in high precision to minimize error
        auto bern_coeffs = convertPowerToBernsteinHP(p_norm.degrees_, p_norm.coeffs_);

        ErrorScalar conversion_err = errorFromConversion(p_norm.degrees_, bern_coeffs);

        PolynomialBase result;
        result.coeffs_ = bern_coeffs;
        result.basis_ = Basis::BERNSTEIN;
        result.degrees_ = degrees_;
        result.dimension_ = dimension_;
        result.error_bound_.absolute_error = error_bound_.absolute_error + conversion_err;
        return result;
    }

    PolynomialBase convertToPower() const {
        if (basis_ == Basis::POWER) return *this;

        // Bernstein is always on [0,1]^n, so result is too
        auto power_coeffs = convertBernsteinToPowerHP(degrees_, coeffs_);

        ErrorScalar conversion_err = errorFromConversion(degrees_, power_coeffs);

        PolynomialBase result;
        result.coeffs_ = power_coeffs;
        result.basis_ = Basis::POWER;
        result.degrees_ = degrees_;
        result.dimension_ = dimension_;
        // Region is standard [0,1]^n (empty vectors)
        result.error_bound_.absolute_error = error_bound_.absolute_error + conversion_err;
        return result;
    }

private:
    //=========================================================================
    // INTERNAL: Ensure power basis on standard region [0,1]^n
    //=========================================================================

    PolynomialBase ensurePowerStandardRegion() const {
        // Already power on standard region?
        if (basis_ == Basis::POWER && hasStandardRegion()) {
            return *this;
        }

        // Bernstein basis?
        if (basis_ == Basis::BERNSTEIN) {
            // Convert to power (result is on [0,1]^n)
            return convertToPower();
        }

        // Power basis on non-standard region - must normalize
        // This is the expensive operation!
        return normalizePowerRegion();
    }

    //=========================================================================
    // INTERNAL: Normalize power polynomial from [a,b]^n to [0,1]^n
    //=========================================================================

    PolynomialBase normalizePowerRegion() const {
        assert(basis_ == Basis::POWER);
        assert(!hasStandardRegion());

        // For p(x) on [a,b], compute q(u) = p(a + (b-a)u) on [0,1]
        // This requires expanding and collecting terms - numerically sensitive!

        // Use high precision for this operation
        auto normalized_coeffs = normalizePowerRegionHP(
            degrees_, coeffs_, region_lower_, region_upper_);

        ErrorScalar normalization_err = errorFromRegionNormalization(
            degrees_, region_lower_, region_upper_);

        PolynomialBase result;
        result.coeffs_ = normalized_coeffs;
        result.basis_ = Basis::POWER;
        result.degrees_ = degrees_;
        result.dimension_ = dimension_;
        // Region is now standard [0,1]^n (empty vectors)
        result.error_bound_.absolute_error = error_bound_.absolute_error + normalization_err;
        return result;
    }
};
```

---

## Key Design Decisions

### 1. **Subdivision: Basis-Specific, Always Returns [0,1]^n**

- **Bernstein**: Modify coefficients (De Casteljau) → result on [0,1]^n
- **Power**: Modify region bounds → result on sub-region, but still evaluates on [0,1]^n

Both approaches give **consistent API**: `restrictedToInterval(axis, a, b)` always returns a polynomial that evaluates on [0,1]^n.

### 2. **Algebraic Operations: Require Normalization**

Operations like differentiation, division, Sturm sequences **require power basis on [0,1]^n**.

The API handles this transparently:
```cpp
// User doesn't care about current representation
auto dp = poly.differentiate(0);  // Works for any basis/region
```

Internally:
1. Check if already power on [0,1]^n → fast path
2. If Bernstein → convert to power
3. If power on [a,b]^n → normalize region (expensive!)

### 3. **Normalization is Expensive - Cache When Possible**

```cpp
// BAD: Normalizes 3 times!
auto dp = poly.differentiate(0);
auto ddp = poly.differentiate(0, 2);
auto sturm = poly.sturmSequence();

// GOOD: Normalize once, reuse
auto poly_power = poly.convertToPower();  // Normalize once
auto dp = poly_power.differentiate(0);
auto ddp = poly_power.differentiate(0, 2);
auto sturm = poly_power.sturmSequence();
```

### 4. **Error Tracking Through Operations**

Each operation adds error:
- **Subdivision (Bernstein)**: Small error from convex combinations
- **Subdivision (Power)**: Tiny error from coordinate transform
- **Conversion**: Moderate error (use high precision!)
- **Normalization**: Large error (use high precision!)
- **Differentiation**: No additional error (exact operation)
- **Division**: Moderate error (use high precision!)

---

## Usage Examples

### Example 1: Subdivision Workflow (Solver)

```cpp
// Start with power polynomial
auto poly = PolynomialBase<double>::fromPower(degrees, coeffs);

// Convert to Bernstein for subdivision (one-time cost)
auto bern = poly.convertToBernstein();

// Subdivide many times (cheap - no normalization needed!)
auto sub1 = bern.restrictedToInterval(0, 0.0, 0.5);
auto sub2 = sub1.restrictedToInterval(1, 0.3, 0.7);
auto sub3 = sub2.restrictedToInterval(0, 0.2, 0.8);

// All results are on [0,1]^n, can evaluate consistently
double val = sub3.evaluate({0.5, 0.5});  // u ∈ [0,1]^2
```

### Example 2: Sturm Sequence (Root Counting)

```cpp
// Polynomial in any basis/region
PolynomialBase<mpreal> poly = /* ... */;

// Sturm sequence (automatically normalizes if needed)
auto sturm = poly.sturmSequence();

// Count roots in [0.3, 0.7]
int count = countRootsInInterval(sturm, mpreal("0.3"), mpreal("0.7"));
```

### Example 3: Mixed Operations

```cpp
// Start in Bernstein (from subdivision)
auto bern = /* ... */;

// Need derivative (converts to power on [0,1]^n)
auto dbern = bern.differentiate(0);

// Subdivide derivative (stays in power, uses region bounds)
auto sub_deriv = dbern.restrictedToInterval(0, 0.2, 0.8);

// Evaluate (works regardless of basis)
double val = sub_deriv.evaluate({0.5});
```

---

## Performance Characteristics

| Operation | Bernstein | Power [0,1]^n | Power [a,b]^n |
|-----------|-----------|---------------|---------------|
| Subdivision | O(n) | O(1) | O(1) |
| Evaluation | O(n²) | O(n) | O(n) |
| Differentiation | O(n³) normalize | O(n) | O(n³) normalize |
| Division | O(n³) normalize | O(n²) | O(n³) normalize |
| Conversion | O(n²) | N/A | O(n³) normalize |

**Key insight**: Power basis on non-standard regions is **fast for subdivision/evaluation** but **slow for algebraic operations**.

---

## Summary

✅ **Unified API**: All operations work on any polynomial regardless of basis/region
✅ **Consistent semantics**: All polynomials evaluate on [0,1]^n
✅ **Automatic normalization**: Algebraic operations normalize transparently
✅ **Performance hints**: User can pre-normalize to avoid repeated cost
✅ **Error tracking**: Every operation updates error bound
✅ **High precision**: Sensitive operations use HP internally

---

## Key Design Decisions

### 1. Differentiation Works in BOTH Bases!

**Important discovery**: Differentiation has simple, stable formulas in both bases:

**Bernstein basis** (degree n → n-1):
```
d_i = n · (b_{i+1} - b_i)    for i = 0, ..., n-1
```

**Power basis** (degree n → n-1):
```
d_i = (i+1) · a_{i+1}        for i = 0, ..., n-1
```

**Implication**: We DON'T need to convert to power for differentiation!

```cpp
PolynomialBase differentiate(size_t axis, unsigned int order = 1) const {
    if (basis_ == Basis::BERNSTEIN) {
        // Use Bernstein formula (simple and stable!)
        auto deriv_coeffs = differentiate_bernstein(degrees_, coeffs_, axis, order);
        return PolynomialBase::fromBernstein(degrees_, deriv_coeffs);
    } else {
        // Use power formula (must be on [0,1]^n)
        PolynomialBase p_normalized = ensurePowerStandardRegion();
        auto deriv_coeffs = differentiate_power(p_normalized.degrees_,
                                                 p_normalized.coeffs_, axis, order);
        return PolynomialBase::fromPower(degrees_, deriv_coeffs);
    }
}
```

### 2. Division Requires Power Basis

Polynomial division in Bernstein basis is complex and not commonly implemented.

**Recommendation**: Convert to power basis for division operations.

```cpp
PolynomialBase divideWithRemainder(const PolynomialBase& divisor) const {
    // Convert both to power on [0,1]^n
    PolynomialBase dividend_power = ensurePowerStandardRegion();
    PolynomialBase divisor_power = divisor.ensurePowerStandardRegion();

    // Divide in power basis
    auto [quot, rem] = divide_power(dividend_power, divisor_power);
    return {quot, rem};
}
```

### 3. Smart Basis Selection

| Operation | Best Basis | Auto-Convert? | Reason |
|-----------|-----------|---------------|--------|
| Differentiation | **Either!** | No | Simple formula in both |
| Subdivision | Bernstein | Yes | De Casteljau algorithm |
| Division | Power | Yes | Standard algorithm |
| Evaluation | Either | No | Use current basis |
| Sturm sequence | Power | Yes | Needs division |

**Performance win**: Differentiating subdivided polynomials no longer requires conversion!

```cpp
// Old design (required 2 conversions):
auto bern = poly.convertToBernstein();  // Conversion 1
auto sub = bern.restrictedToInterval(0, 0.3, 0.7);
auto deriv = sub.differentiate(0);  // Conversion 2 (Bernstein → Power)

// New design (only 1 conversion):
auto bern = poly.convertToBernstein();  // Conversion 1
auto sub = bern.restrictedToInterval(0, 0.3, 0.7);
auto deriv = sub.differentiate(0);  // No conversion! (Bernstein formula)
```

**Savings**: 50% fewer conversions in subdivision-heavy solvers!

