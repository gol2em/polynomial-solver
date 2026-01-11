# Key Insight: Differentiation Works in Both Bases!

## The Discovery

While researching polynomial division in Bernstein basis, I discovered that **differentiation has simple, efficient formulas in BOTH Bernstein and power bases**.

This is a **major design simplification** because:
1. We don't need to convert to power for differentiation
2. Subdivision-heavy solvers can stay in Bernstein throughout
3. Performance improves significantly (fewer conversions)

---

## Differentiation Formulas

### Bernstein Basis (Univariate)

For polynomial of degree n with Bernstein coefficients b₀, b₁, ..., bₙ:

```
d/dx B(x) = n · Σ(b_{i+1} - b_i) · B_{i,n-1}(x)
```

**Derivative coefficients**:
```
d_i = n · (b_{i+1} - b_i)    for i = 0, ..., n-1
```

**Properties**:
- ✅ Simple: Just differences and scaling
- ✅ Stable: No binomial coefficients
- ✅ Efficient: O(n) operations
- ✅ Degree reduction: Automatic (n → n-1)

**Example**:
```
B(x) = [1, 2, 4, 5]  (degree 3)
dB/dx = 3·[(2-1), (4-2), (5-4)] = [3, 6, 3]  (degree 2)
```

### Power Basis (Univariate)

For polynomial of degree n with power coefficients a₀, a₁, ..., aₙ:

```
d/dx P(x) = Σ i·a_i·x^(i-1)
```

**Derivative coefficients**:
```
d_i = (i+1) · a_{i+1}    for i = 0, ..., n-1
```

**Properties**:
- ✅ Simple: Just scaling
- ✅ Stable: No cancellation
- ✅ Efficient: O(n) operations
- ✅ Degree reduction: Automatic (n → n-1)

**Example**:
```
P(x) = [1, 2, 4, 5]  (1 + 2x + 4x² + 5x³)
dP/dx = [2, 8, 15]   (2 + 8x + 15x²)
```

---

## Implementation

### Unified Interface

```cpp
template<typename Scalar, typename ErrorScalar>
class PolynomialBase {
public:
    PolynomialBase differentiate(size_t axis, unsigned int order = 1) const {
        if (basis_ == Basis::BERNSTEIN) {
            // Use Bernstein formula
            auto deriv_coeffs = differentiate_bernstein(degrees_, coeffs_, axis, order);
            
            PolynomialBase result;
            result.coeffs_ = deriv_coeffs;
            result.basis_ = Basis::BERNSTEIN;
            result.degrees_ = degrees_;  // Degree reduced by order
            result.dimension_ = dimension_;
            result.error_bound_ = /* updated */;
            return result;
            
        } else {
            // Use power formula (must be on [0,1]^n)
            PolynomialBase p_normalized = ensurePowerStandardRegion();
            auto deriv_coeffs = differentiate_power(p_normalized.degrees_, 
                                                     p_normalized.coeffs_, axis, order);
            
            PolynomialBase result;
            result.coeffs_ = deriv_coeffs;
            result.basis_ = Basis::POWER;
            result.degrees_ = p_normalized.degrees_;
            result.dimension_ = dimension_;
            result.error_bound_ = /* updated */;
            return result;
        }
    }
};
```

### Multivariate Case

For tensor-product polynomials, differentiate along each axis:

**Bernstein**:
```cpp
// Differentiate along axis k
for (size_t i = 0; i < num_coeffs; ++i) {
    auto idx = multiIndexFromLinear(i, degrees_);
    if (idx[k] < degrees_[k]) {
        auto idx_next = idx;
        idx_next[k]++;
        
        size_t i_next = linearIndexFromMulti(idx_next, degrees_);
        deriv_coeffs[i] = degrees_[k] * (coeffs_[i_next] - coeffs_[i]);
    }
}
```

**Power**:
```cpp
// Differentiate along axis k
for (size_t i = 0; i < num_coeffs; ++i) {
    auto idx = multiIndexFromLinear(i, degrees_);
    if (idx[k] > 0) {
        auto idx_prev = idx;
        idx_prev[k]--;
        
        size_t i_prev = linearIndexFromMulti(idx_prev, degrees_);
        deriv_coeffs[i_prev] = idx[k] * coeffs_[i];
    }
}
```

---

## Performance Impact

### Before (Old Design)

Differentiation required conversion to power basis:

```cpp
// Workflow: Subdivide → Differentiate → Subdivide → ...
auto bern = poly.convertToBernstein();  // Conversion 1
auto sub1 = bern.restrictedToInterval(0, 0.3, 0.7);

auto power = sub1.convertToPower();     // Conversion 2 ❌
auto deriv = power.differentiate(0);

auto bern2 = deriv.convertToBernstein(); // Conversion 3 ❌
auto sub2 = bern2.restrictedToInterval(0, 0.4, 0.6);
```

**Total conversions**: 3

### After (New Design)

Differentiation works in Bernstein:

```cpp
// Workflow: Subdivide → Differentiate → Subdivide → ...
auto bern = poly.convertToBernstein();  // Conversion 1
auto sub1 = bern.restrictedToInterval(0, 0.3, 0.7);

auto deriv = sub1.differentiate(0);     // No conversion! ✅

auto sub2 = deriv.restrictedToInterval(0, 0.4, 0.6);
```

**Total conversions**: 1

**Improvement**: 66% fewer conversions!

---

## Numerical Stability

Both formulas are numerically stable:

### Bernstein
- Only uses **differences** and **scaling**
- No binomial coefficients (which can be large)
- Coefficients remain bounded by original polynomial

### Power
- Only uses **scaling**
- No cancellation
- Straightforward formula

---

## References

### Bernstein Differentiation

From Farin's "Curves and Surfaces for CAGD" (5th edition):

> "The derivative of a Bernstein polynomial is given by:
> d/dx B(x) = n · Σ Δb_i · B_{i,n-1}(x)
> where Δb_i = b_{i+1} - b_i are the forward differences."

### Implementation

Already implemented in `src/polynomial/differentiation.cpp`:

```cpp
template<typename Scalar>
std::vector<Scalar> differentiate_bernstein(
    const std::vector<unsigned int>& degrees,
    const std::vector<Scalar>& coeffs,
    size_t axis,
    unsigned int order) {
    // Implementation using forward differences
    // ...
}
```

---

## Summary

✅ **Differentiation works in BOTH bases**
✅ **Simple formulas** (differences for Bernstein, scaling for power)
✅ **Numerically stable** (no cancellation)
✅ **Performance win** (fewer conversions in subdivision-heavy solvers)
✅ **Already implemented** (in `differentiation.cpp`)

**Key takeaway**: We can stay in Bernstein basis for subdivision + differentiation workflows!

