# Operations in Bernstein Basis

## Summary

| Operation | Bernstein Basis | Recommendation |
|-----------|----------------|----------------|
| **Differentiation** | ✅ **Easy and stable** | Use Bernstein formula directly |
| **Subdivision** | ✅ **Easy and stable** (De Casteljau) | Use Bernstein |
| **Evaluation** | ✅ **Stable** (De Casteljau) | Use Bernstein |
| **Division** | ❌ **Complex** | Convert to power first |
| **GCD** | ❌ **Complex** | Convert to power first |

---

## Differentiation in Bernstein Basis

### Formula (Univariate)

For Bernstein polynomial of degree n with coefficients b₀, b₁, ..., bₙ:

```
d/dx B(x) = n · Σ(b_{i+1} - b_i) · B_{i,n-1}(x)
```

Where the derivative has degree n-1 with coefficients:
```
d_i = n · (b_{i+1} - b_i)    for i = 0, ..., n-1
```

### Properties

✅ **Numerically stable**: Only differences and scaling
✅ **Simple**: No binomial coefficients needed
✅ **Degree reduction**: Automatic (degree n → n-1)
✅ **Already implemented**: See `differentiation.cpp`

### Example

```cpp
// Bernstein polynomial of degree 3: b = [1, 2, 4, 5]
// Derivative has degree 2 with coefficients:
// d[0] = 3 * (2 - 1) = 3
// d[1] = 3 * (4 - 2) = 6
// d[2] = 3 * (5 - 4) = 3
// Result: d = [3, 6, 3]
```

### Multivariate Case

For tensor-product Bernstein polynomials, differentiate along each axis separately:

```cpp
∂/∂x_k B(x₁,...,xₙ) = n_k · Σ(b_{...,i+1,...} - b_{...,i,...}) · B_{i,n_k-1}(x_k) · ∏_{j≠k} B_{·,n_j}(x_j)
```

---

## Division in Bernstein Basis

### UPDATE: Division Algorithm Exists!

**Important discovery**: There IS a division algorithm for Bernstein polynomials!

From Busé et al. (2008), the algorithm uses:
1. **Homogenization**: Convert to homogeneous polynomials
2. **Division formula**: `g(t) = q(t)·f(t) + (1-t)^(e-d+1)·r(t)`
3. **Key operations**: Multiply/divide by powers of `t` and `(1-t)`

See [BERNSTEIN_DIVISION_ALGORITHM.md](BERNSTEIN_DIVISION_ALGORITHM.md) for full details.

### Key Formulas

**Multiply by `t^d`**:
```
c^(n+d)_k = ∏(i=0 to d-1) (k-i)/(n+d-i) · c^n_(k-d)
```

**Multiply by `(1-t)^d`**:
```
c^(n+d)_k = ∏(i=0 to d-1) (n+d-k-i)/(n+d-i) · c^n_k
```

**Divide by `(1-t)^j`**:
```
c^(n-j)_k = ∏(i=0 to j-1) (n-i)/(n-k-i) · c^n_k
```

### Implementation Strategy

**Option 1: Implement Bernstein division** (recommended for pure Bernstein workflows)

```cpp
PolynomialBase divideWithRemainder(const PolynomialBase& divisor) const {
    if (basis_ == Basis::BERNSTEIN && divisor.basis_ == Basis::BERNSTEIN) {
        // Use Bernstein division algorithm
        auto [quot_coeffs, rem_coeffs] = divide_bernstein(
            degrees_, coeffs_,
            divisor.degrees_, divisor.coeffs_);

        return {
            PolynomialBase::fromBernstein(quot_degrees, quot_coeffs),
            PolynomialBase::fromBernstein(rem_degrees, rem_coeffs)
        };
    } else {
        // Convert to power basis
        PolynomialBase dividend_power = ensurePowerStandardRegion();
        PolynomialBase divisor_power = divisor.ensurePowerStandardRegion();

        auto [quot_coeffs, rem_coeffs] = divide_power(
            dividend_power.degrees_, dividend_power.coeffs_,
            divisor_power.degrees_, divisor_power.coeffs_);

        return {
            PolynomialBase::fromPower(quot_degrees, quot_coeffs),
            PolynomialBase::fromPower(rem_degrees, rem_coeffs)
        };
    }
}
```

**Option 2: Convert to power** (simpler, but requires conversion)

```cpp
PolynomialBase divideWithRemainder(const PolynomialBase& divisor) const {
    // Convert both to power basis on [0,1]^n
    PolynomialBase dividend_power = ensurePowerStandardRegion();
    PolynomialBase divisor_power = divisor.ensurePowerStandardRegion();

    // Divide in power basis (standard algorithm)
    auto [quot_coeffs, rem_coeffs] = divide_power(
        dividend_power.degrees_, dividend_power.coeffs_,
        divisor_power.degrees_, divisor_power.coeffs_);

    return {
        PolynomialBase::fromPower(quot_degrees, quot_coeffs),
        PolynomialBase::fromPower(rem_degrees, rem_coeffs)
    };
}
```

### Recommendation

**For now**: Use Option 2 (convert to power) for simplicity.

**Future**: Implement Option 1 (Bernstein division) if:
- Sturm sequences are computed frequently in Bernstein basis
- Conversion overhead becomes significant
- Need to maintain Bernstein representation throughout

---

## Sturm Sequence

### Requirements

Sturm sequence requires:
1. **Differentiation**: ✅ Works in Bernstein!
2. **Division with remainder**: ❌ Needs power basis

### Implementation Strategy

```cpp
std::vector<PolynomialBase> sturmSequence() const {
    // Convert to power basis once
    PolynomialBase p0 = ensurePowerStandardRegion();
    
    // Differentiate (can stay in power or convert to Bernstein)
    PolynomialBase p1 = p0.differentiate(0, 1);
    
    std::vector<PolynomialBase> sequence;
    sequence.push_back(p0);
    sequence.push_back(p1);
    
    // Build sequence using division (requires power basis)
    while (!sequence.back().isZero()) {
        auto [quot, rem] = sequence[sequence.size() - 2].divideWithRemainder(sequence.back());
        if (rem.isZero()) break;
        sequence.push_back(-rem);
    }
    
    return sequence;
}
```

**Note**: Since division requires power basis anyway, the entire Sturm sequence is computed in power basis.

---

## Updated Design Decisions

### When to Use Each Basis

**Use Bernstein for:**
- ✅ Subdivision (De Casteljau)
- ✅ Evaluation
- ✅ Differentiation (simple formula!)
- ✅ Convex hull properties
- ✅ Geometric operations

**Use Power for:**
- ✅ Division with remainder
- ✅ GCD computation
- ✅ Sturm sequences
- ✅ Polynomial arithmetic (when needed)

### Conversion Strategy

1. **Start in power basis** (input from user)
2. **Convert to Bernstein** for subdivision-heavy solver
3. **Differentiate in Bernstein** (no conversion needed!)
4. **Convert to power** only when division is needed (Sturm sequence)

---

## Performance Implications

### Before (Old Design)

```cpp
// Differentiation required conversion to power
auto bern = poly.convertToBernstein();  // Conversion 1
auto sub = bern.restrictedToInterval(0, 0.3, 0.7);
auto deriv = sub.differentiate(0);  // Conversion 2 (Bernstein → Power)
```

**Cost**: 2 conversions

### After (New Design)

```cpp
// Differentiation works in Bernstein!
auto bern = poly.convertToBernstein();  // Conversion 1
auto sub = bern.restrictedToInterval(0, 0.3, 0.7);
auto deriv = sub.differentiate(0);  // No conversion needed!
```

**Cost**: 1 conversion

**Savings**: 50% fewer conversions when differentiating subdivided polynomials!

---

## Summary

✅ **Differentiation**: Use Bernstein formula (simple, stable, already implemented)
❌ **Division**: Convert to power basis (complex in Bernstein, not worth implementing)
✅ **Overall**: Bernstein is great for subdivision + differentiation, power for division

