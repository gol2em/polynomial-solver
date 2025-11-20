# Complexity Analysis: Derivatives vs. Subdivision

## Problem Statement

When working with polynomial derivatives on subregions, we have two approaches:

**Approach A**: Restrict polynomial to subregion first, then differentiate
```
p → restrict([a,b]) → p_restricted → differentiate → dp_restricted
```

**Approach B**: Differentiate first, then restrict to subregion
```
p → differentiate → dp → restrict([a,b]) → dp_restricted
```

This document analyzes the computational complexity of both approaches to determine the optimal strategy.

---

## 1. Complexity of Basic Operations

### 1.1 De Casteljau Subdivision (1D)

**Operation**: `subdivide1D(coeffs, t, left, right)`

**Algorithm** (from `src/de_casteljau.cpp:41-76`):
```cpp
// Build De Casteljau triangle
for (r = 1; r <= degree; ++r) {
    for (i = 0; i < degree - r + 1; ++i) {
        b[r][i] = (1-t) * b[r-1][i] + t * b[r-1][i+1];
    }
}
```

**Complexity**:
- Outer loop: `degree` iterations
- Inner loop: `(degree - r + 1)` iterations
- Total operations: `Σ(r=1 to n) (n - r + 1) = n + (n-1) + ... + 1 = n(n+1)/2`

**Result**: **O(n²)** where n = degree

**Memory**: O(n²) for the triangle (can be optimized to O(n))

### 1.2 Restriction to Interval [a,b] (1D)

**Operation**: `restrictedToInterval(axis, a, b)` for 1D polynomial

**Algorithm** (from `src/polynomial.cpp:360-435`):
```cpp
// Two subdivisions:
if (a > 0.0) {
    subdivide1D(segment, a, left, right);  // O(n²)
    segment = right;
}
subdivide1D(segment, t_rel, left, right);  // O(n²)
```

**Complexity**: **O(n²)** (two subdivisions, both O(n²))

### 1.3 Differentiation (1D, Bernstein Basis)

**Operation**: Compute derivative of Bernstein polynomial

**Algorithm**:
```cpp
// For degree n with coefficients b[0..n]:
// Derivative has degree n-1 with coefficients:
// new_coeffs[i] = n * (b[i+1] - b[i]) for i = 0..n-1
for (i = 0; i < n; ++i) {
    new_coeffs[i] = n * (coeffs[i+1] - coeffs[i]);
}
```

**Complexity**: **O(n)** where n = degree

**Memory**: O(n) for new coefficients

---

## 2. Multivariate Complexity Analysis

For tensor-product polynomials with dimension d and degrees (n₀, n₁, ..., n_{d-1}):

**Total coefficients**: N = ∏(nᵢ + 1)

### 2.1 Restriction Along One Axis (Multivariate)

**Algorithm** (from `src/polynomial.cpp:399-432`):
```cpp
// Iterate over all 1D slices perpendicular to the axis
num_slices = N / (n_axis + 1)
for each slice:
    gather coefficients along axis      // O(n_axis)
    subdivide1D twice                   // O(n_axis²)
    scatter coefficients back           // O(n_axis)
```

**Complexity**: **O(N · n_axis)** where:
- N = total coefficients
- n_axis = degree along the restricted axis

**Example**: For 2D polynomial with degrees (10, 20):
- N = 11 × 21 = 231 coefficients
- Restrict along axis 0: O(231 × 10) = O(2,310) operations
- Restrict along axis 1: O(231 × 20) = O(4,620) operations

### 2.2 Differentiation Along One Axis (Multivariate)

**Algorithm**:
```cpp
// Iterate over all 1D slices perpendicular to the axis
num_slices = N / (n_axis + 1)
for each slice:
    gather coefficients along axis      // O(n_axis)
    differentiate 1D                    // O(n_axis)
    scatter coefficients back           // O(n_axis - 1)
```

**Complexity**: **O(N)** (linear in total coefficients)

**Note**: After differentiation, the new polynomial has fewer coefficients:
- N_new = N × n_axis / (n_axis + 1)

**Example**: For 2D polynomial with degrees (10, 20):
- N = 231 coefficients
- Differentiate along axis 0: O(231) operations, result has 10 × 21 = 210 coefficients
- Differentiate along axis 1: O(231) operations, result has 11 × 20 = 220 coefficients

---

## 3. Comparison: Approach A vs. Approach B

### 3.1 Single Derivative on Subregion (1D)

**Approach A**: Restrict then differentiate
```
restrict: O(n²)
differentiate: O(n-1) ≈ O(n)
Total: O(n²)
```

**Approach B**: Differentiate then restrict
```
differentiate: O(n)
restrict: O((n-1)²) ≈ O(n²)
Total: O(n²)
```

**Winner**: **Approach B** (slightly better constant factors)
- Differentiation reduces degree first, making restriction cheaper
- Differentiation is O(n) vs. restriction's O(n²)

### 3.2 Multiple Derivatives on Same Subregion (1D)

**Scenario**: Compute f, f', f'', f''' on subregion [a,b]

**Approach A**: Restrict once, then differentiate multiple times
```
restrict: O(n²)
diff 1st: O(n)
diff 2nd: O(n-1)
diff 3rd: O(n-2)
Total: O(n²) + O(n) + O(n) + O(n) = O(n²)
```

**Approach B**: Differentiate multiple times, then restrict each
```
diff 1st: O(n)
restrict: O((n-1)²)
diff 2nd: O(n-1)
restrict: O((n-2)²)
diff 3rd: O(n-2)
restrict: O((n-3)²)
Total: O(n²) + O(n²) + O(n²) = O(3n²)
```

**Winner**: **Approach A** (3× better)
- Restrict once, differentiate many times
- Each differentiation is cheap O(n)

### 3.3 Multivariate Case (2D Example)

**Setup**: Polynomial with degrees (n, m), total N = (n+1)(m+1) coefficients
**Task**: Compute ∂f/∂x on subregion [a,b] × [c,d]

**Approach A**: Restrict both axes, then differentiate
```
restrict axis 0: O(N · n)
restrict axis 1: O(N · m)  [N unchanged]
differentiate axis 0: O(N)
Total: O(N · (n + m + 1))
```

**Approach B**: Differentiate, then restrict both axes
```
differentiate axis 0: O(N)
restrict axis 0: O(N' · n)  [N' = N · n/(n+1)]
restrict axis 1: O(N' · m)
Total: O(N + N' · (n + m)) ≈ O(N · (n + m))
```

**Winner**: **Approach B** (slightly better)
- Differentiation reduces coefficient count before restriction

### 3.4 Gradient on Subregion (Multivariate)

**Setup**: d-dimensional polynomial, compute all d partial derivatives on subregion

**Approach A**: Restrict to subregion, then compute gradient
```
restrict d axes: O(d · N · n_avg)  [n_avg = average degree]
compute d derivatives: O(d · N)
Total: O(d · N · (n_avg + 1))
```

**Approach B**: Compute gradient, then restrict each derivative
```
compute d derivatives: O(d · N)
restrict d derivatives, d axes each: O(d² · N' · n_avg)  [N' ≈ N]
Total: O(d² · N · n_avg)
```

**Winner**: **Approach A** when d is large
- Approach A: O(d · N · n_avg)
- Approach B: O(d² · N · n_avg)
- Approach A is d times better

### 3.5 Hessian on Subregion (Multivariate)

**Setup**: d-dimensional polynomial, compute all d² second derivatives on subregion

**Approach A**: Restrict to subregion, then compute Hessian
```
restrict d axes: O(d · N · n_avg)
compute d² second derivatives: O(d² · N)
Total: O(d · N · (d · n_avg + n_avg))
```

**Approach B**: Compute Hessian, then restrict each
```
compute d² second derivatives: O(d² · N)
restrict d² derivatives, d axes each: O(d³ · N' · n_avg)
Total: O(d³ · N · n_avg)
```

**Winner**: **Approach A** (d² times better)
- Approach A: O(d² · N · n_avg)
- Approach B: O(d³ · N · n_avg)

---

## 4. Summary and Recommendations

### 4.1 Complexity Summary Table

| Operation | Approach A (Restrict → Diff) | Approach B (Diff → Restrict) | Winner |
|-----------|------------------------------|------------------------------|--------|
| Single derivative (1D) | O(n²) | O(n²) | **B** (better constants) |
| k derivatives (1D) | O(n²) | O(k · n²) | **A** (k× better) |
| Single derivative (d-D) | O(N · n_avg) | O(N · n_avg) | **B** (better constants) |
| Gradient (d-D) | O(d · N · n_avg) | O(d² · N · n_avg) | **A** (d× better) |
| Hessian (d-D) | O(d² · N · n_avg) | O(d³ · N · n_avg) | **A** (d²× better) |

Where:
- n = degree (1D)
- d = dimension
- N = total coefficients = ∏(nᵢ + 1)
- n_avg = average degree across dimensions

### 4.2 Practical Recommendations

#### **Use Approach A (Restrict → Differentiate) when:**

1. **Multiple derivatives needed** on the same subregion
   - Computing f, f', f'', ... on [a,b]
   - Computing gradient (all ∂f/∂xᵢ) on a box
   - Computing Hessian (all ∂²f/∂xᵢ∂xⱼ) on a box

2. **High-dimensional problems** (d ≥ 3)
   - Gradient computation: d× speedup
   - Hessian computation: d²× speedup

3. **Iterative refinement** scenarios
   - Newton's method: need gradient/Hessian at each iteration
   - Multiplicity detection: need f, f', f'', ... at suspected root

#### **Use Approach B (Differentiate → Restrict) when:**

1. **Single derivative** on a subregion
   - Only need one derivative, not multiple

2. **Different subregions** for each derivative
   - Need f' on [a₁,b₁] and f'' on [a₂,b₂]

3. **Caching derivatives** for reuse across many subregions
   - Compute f' once, restrict to many different boxes
   - Useful in subdivision algorithms

### 4.3 Hybrid Approach for Subdivision Algorithms

In the context of your subdivision solver, the optimal strategy is:

**Strategy**: Cache derivatives globally, restrict on-demand per box

```cpp
// At solver initialization:
DerivativeCache global_cache(polynomial);
global_cache.precomputeUpToOrder(max_order);  // e.g., order 3

// For each subdivision box:
for (auto& box : boxes) {
    // Get cached derivatives (already computed)
    const Polynomial& f = global_cache.get({0, 0});
    const Polynomial& fx = global_cache.get({1, 0});
    const Polynomial& fy = global_cache.get({0, 1});

    // Restrict to current box (cheap, done once per box)
    Polynomial f_box = f.restrictedToInterval(0, box.lower[0], box.upper[0])
                        .restrictedToInterval(1, box.lower[1], box.upper[1]);
    Polynomial fx_box = fx.restrictedToInterval(0, box.lower[0], box.upper[0])
                          .restrictedToInterval(1, box.lower[1], box.upper[1]);
    // ... use for degeneracy detection
}
```

**Benefits**:
- Derivatives computed once: O(d² · N) for Hessian
- Restriction per box: O(d · N · n_avg) per box
- Total: O(d² · N + B · d · N · n_avg) where B = number of boxes
- Much better than computing derivatives per box: O(B · d² · N)

---

## 5. Numerical Considerations

### 5.1 Numerical Stability

**Differentiation**:
- Amplifies high-frequency errors by factor of n (degree)
- Higher derivatives amplify more: k-th derivative amplifies by n!/(n-k)!

**Subdivision**:
- Generally stable (convex combinations)
- Error grows slowly with subdivision depth

**Recommendation**: For high-order derivatives (k ≥ 3), consider:
- Using higher precision arithmetic
- Monitoring condition numbers
- Switching to alternative methods (e.g., finite differences) for verification

### 5.2 Memory Considerations

**Caching all derivatives up to order k**:
- 1D: k+1 polynomials
- 2D: (k+1)(k+2)/2 polynomials (triangular)
- d-D: C(d+k, k) polynomials (binomial coefficient)

**Example** (d=3, k=3):
- Number of derivatives: C(6,3) = 20 polynomials
- If each has ~1000 coefficients: ~20KB per polynomial = 400KB total

**Recommendation**: For high dimensions or high orders, use sparse caching (only compute what's needed).

---

## 6. Conclusion

**For your degeneracy detection use case:**

✅ **Use Approach B with global caching** (Differentiate → Restrict):

1. **At solver initialization**: Compute and cache all needed derivatives globally
   ```cpp
   DerivativeCache cache(polynomial);
   cache.precomputeUpToOrder(3);  // For multiplicity detection
   ```

2. **Per subdivision box**: Restrict cached derivatives to box
   ```cpp
   Polynomial f_box = cache.get({0}).restrictedToInterval(...);
   Polynomial df_box = cache.get({1}).restrictedToInterval(...);
   ```

3. **For degeneracy detection**: Evaluate restricted derivatives at box center
   ```cpp
   double f_val = f_box.evaluate(box.center);
   double df_val = df_box.evaluate(box.center);
   // Check if |f_val|, |df_val|, ... are all small
   ```

**Complexity**: O(d^k · N) for initial caching + O(B · d · N · n_avg) for B boxes
**Memory**: O(C(d+k, k) · N) for cached derivatives
**Numerical accuracy**: 10-5500× better than Approach A (see numerical_stability_analysis.md)

This is optimal for your use case where you need multiple derivatives on many boxes.

**Why this strategy wins**:
1. ✅ **Computational**: k× faster than computing derivatives per box
2. ✅ **Numerical**: 10-5500× more accurate than Approach A (experimental results)
3. ✅ **Memory**: Sparse caching, only compute what's needed
4. ✅ **Simplicity**: Differentiate once, restrict many times

