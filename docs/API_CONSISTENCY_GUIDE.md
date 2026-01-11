# API Consistency Guide: How Subdivision Works Across Both Bases

## The Core Question

**How do we keep subdivision consistent when:**
- Bernstein basis modifies coefficients (De Casteljau)
- Power basis should NOT modify coefficients (numerical instability)
- Both need to support differentiation, division, Sturm sequences

## The Solution: Region Bounds + Auto-Normalization

---

## Subdivision: Two Strategies, One API

### Bernstein Basis: Coefficient Modification

```cpp
// Input: Bernstein polynomial on [0,1]^n
PolynomialBase<double> bern = PolynomialBase::fromBernstein(degrees, coeffs);

// Subdivide: Modify coefficients using De Casteljau
auto sub = bern.restrictedToInterval(0, 0.3, 0.7);

// Result: New Bernstein polynomial on [0,1]^n
// Coefficients are different, but polynomial evaluates on [0,1]^n
```

**What happens:**
1. De Casteljau algorithm computes new Bernstein coefficients
2. New coefficients represent the restricted polynomial
3. Result is on [0,1]^n (always!)

### Power Basis: Region Modification

```cpp
// Input: Power polynomial on [0,1]^n
PolynomialBase<double> power = PolynomialBase::fromPower(degrees, coeffs);

// Subdivide: Modify region bounds (coefficients unchanged!)
auto sub = power.restrictedToInterval(0, 0.3, 0.7);

// Result: Same coefficients, but region is now [0.3, 0.7] × [0,1]^{n-1}
// Still evaluates on [0,1]^n via coordinate transformation
```

**What happens:**
1. Region bounds updated: `region_lower[0] = 0.3`, `region_upper[0] = 0.7`
2. Coefficients unchanged (preserves exactness!)
3. Evaluation transforms: `x[0] = 0.3 + 0.4 * u[0]` where u ∈ [0,1]

### Unified Evaluation

**Both bases evaluate at u ∈ [0,1]^n:**

```cpp
// Bernstein: Direct evaluation
double val_bern = bern_sub.evaluate({0.5, 0.5});  // u ∈ [0,1]^2

// Power: Transform then evaluate
double val_power = power_sub.evaluate({0.5, 0.5});  // u ∈ [0,1]^2
// Internally: x = [0.3 + 0.4*0.5, 0.5] = [0.5, 0.5]
//             val = horner_eval(coeffs, x)
```

**Result: Same API, same semantics, different implementation!**

---

## Algebraic Operations: The Normalization Challenge

### The Problem

Differentiation, division, and Sturm sequences require **power basis on [0,1]^n**.

But after subdivision, power polynomials are on **[a,b]^n**!

```cpp
// Power polynomial on [0.3, 0.7]
auto power_sub = power.restrictedToInterval(0, 0.3, 0.7);

// What does differentiation mean here?
auto dpoly = power_sub.differentiate(0);  // ???
```

**Issue**: d/dx p(x) where x ∈ [0.3, 0.7] is NOT what we want!

We want: d/du p(0.3 + 0.4u) where u ∈ [0,1]

### The Solution: Auto-Normalization

```cpp
PolynomialBase differentiate(size_t axis, unsigned int order) const {
    // Step 1: Ensure power basis on [0,1]^n
    PolynomialBase p_norm = ensurePowerStandardRegion();
    
    // Step 2: Differentiate (standard formula)
    auto deriv_coeffs = differentiate_power(p_norm.degrees_, 
                                             p_norm.coeffs_, 
                                             axis, order);
    
    // Step 3: Return result on [0,1]^n
    return PolynomialBase::fromPower(deriv_degrees, deriv_coeffs);
}
```

**What `ensurePowerStandardRegion()` does:**

1. **Already power on [0,1]^n?** → Return as-is (fast path)
2. **Bernstein?** → Convert to power (result on [0,1]^n)
3. **Power on [a,b]^n?** → Normalize region (expensive!)

### Region Normalization (The Expensive Operation)

For p(x) on [a,b], compute q(u) = p(a + (b-a)u) on [0,1]:

```cpp
// Example: p(x) = c₀ + c₁x + c₂x² on [0.3, 0.7]
// Want: q(u) = p(0.3 + 0.4u) on [0,1]

// Substitute x = 0.3 + 0.4u:
// q(u) = c₀ + c₁(0.3 + 0.4u) + c₂(0.3 + 0.4u)²
//      = c₀ + c₁·0.3 + c₁·0.4u + c₂·0.09 + c₂·0.24u + c₂·0.16u²
//      = (c₀ + 0.3c₁ + 0.09c₂) + (0.4c₁ + 0.24c₂)u + (0.16c₂)u²

// New coefficients: [c₀ + 0.3c₁ + 0.09c₂, 0.4c₁ + 0.24c₂, 0.16c₂]
```

**This is numerically sensitive!** → Use high precision

---

## Complete API Flow

### Scenario 1: Subdivision-Heavy (Solver)

```cpp
// Start with power
auto poly = PolynomialBase<double>::fromPower(degrees, coeffs);

// Convert to Bernstein ONCE
auto bern = poly.convertToBernstein();

// Subdivide many times (cheap - no normalization!)
for (int i = 0; i < 100; ++i) {
    bern = bern.restrictedToInterval(axis, a, b);
}

// Evaluate (cheap)
double val = bern.evaluate({0.5, 0.5});
```

**Cost**: 1 conversion + 100 cheap subdivisions

### Scenario 2: Algebra-Heavy (Sturm)

```cpp
// Start with any polynomial
auto poly = /* ... */;

// Ensure power on [0,1]^n ONCE
auto poly_norm = poly.convertToPower();  // Or ensurePowerStandardRegion()

// Many algebraic operations (cheap - no normalization!)
auto dp = poly_norm.differentiate(0);
auto ddp = poly_norm.differentiate(0, 2);
auto sturm = poly_norm.sturmSequence();
```

**Cost**: 1 normalization + many cheap operations

### Scenario 3: Mixed Operations

```cpp
// Subdivided Bernstein polynomial
auto bern_sub = /* ... */;

// Differentiate (converts to power on [0,1]^n)
auto dbern = bern_sub.differentiate(0);

// Subdivide derivative (uses region bounds - cheap!)
auto sub_deriv = dbern.restrictedToInterval(0, 0.2, 0.8);

// Evaluate (transforms coordinates)
double val = sub_deriv.evaluate({0.5});
```

**Cost**: 1 conversion + cheap operations

---

## Performance Summary

| Operation | Bernstein | Power [0,1]^n | Power [a,b]^n |
|-----------|-----------|---------------|---------------|
| `restrictedToInterval` | O(n) modify coeffs | O(1) update bounds | O(1) update bounds |
| `evaluate` | O(n²) De Casteljau | O(n) Horner | O(n) transform + Horner |
| `differentiate` | O(n³) normalize | O(n) direct | O(n³) normalize |
| `divideWithRemainder` | O(n³) normalize | O(n²) direct | O(n³) normalize |
| `sturmSequence` | O(n³) normalize | O(n²) direct | O(n³) normalize |

**Key insight**: Power on [a,b]^n is great for subdivision, but expensive for algebra!

---

## Best Practices

### ✅ DO: Pre-normalize for algebra

```cpp
auto poly_norm = poly.convertToPower();  // Normalize once
auto dp = poly_norm.differentiate(0);
auto sturm = poly_norm.sturmSequence();
```

### ✅ DO: Use Bernstein for subdivision

```cpp
auto bern = poly.convertToBernstein();  // Convert once
for (int i = 0; i < 1000; ++i) {
    bern = bern.restrictedToInterval(axis, a, b);
}
```

### ❌ DON'T: Repeatedly normalize

```cpp
// BAD: Normalizes 3 times!
auto dp = poly.differentiate(0);
auto ddp = poly.differentiate(0, 2);
auto sturm = poly.sturmSequence();
```

### ✅ DO: Let API handle it automatically

```cpp
// OK: API handles normalization transparently
auto dp = poly.differentiate(0);  // Works for any basis/region
```

---

## Summary

✅ **Unified API**: All operations work on any polynomial
✅ **Consistent semantics**: All polynomials evaluate on [0,1]^n
✅ **Optimal performance**: Choose right basis for your workload
✅ **Automatic normalization**: Algebraic ops normalize transparently
✅ **Error tracking**: Every operation updates error bound
✅ **Numerical stability**: Power subdivision via coordinates, not coefficients

