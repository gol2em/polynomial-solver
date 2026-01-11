# Power vs Bernstein: Detailed Comparison

## Subdivision Cost Analysis

### Power Basis (Region Tracking)

**Algorithm**:
```cpp
// Subdivide [a,b] to [c,d] where c,d ∈ [0,1]
region_lower_new = region_lower + c * (region_upper - region_lower)
region_upper_new = region_lower + d * (region_upper - region_lower)
```

**Cost**:
- Time: **O(1)** - 4 arithmetic operations per axis
- Space: **O(1)** - update 2 numbers
- Numerical error: **0** - coefficients unchanged

**Example** (degree 100 polynomial):
- Subdivide 1000 times: 4000 operations total
- Coefficients: unchanged (101 values)

### Bernstein Basis (De Casteljau)

**Algorithm**:
```cpp
// De Casteljau subdivision
for (int r = 1; r <= n; r++) {
    for (int i = 0; i <= n - r; i++) {
        b[i] = (1-t) * b[i] + t * b[i+1];
    }
}
```

**Cost**:
- Time: **O(n²)** - n(n+1)/2 operations
- Space: **O(n)** - store n+1 new coefficients
- Numerical error: **Accumulates** - each subdivision adds error

**Example** (degree 100 polynomial):
- Subdivide once: ~5050 operations
- Subdivide 1000 times: ~5,050,000 operations
- Coefficients: recomputed 1000 times (101 values each)

### Comparison

| Metric | Power | Bernstein | Ratio |
|--------|-------|-----------|-------|
| Time (1 subdivision) | O(1) | O(n²) | **1000x faster** (n=100) |
| Time (1000 subdivisions) | O(1) | O(1000n²) | **1000x faster** |
| Space | O(1) | O(n) | **100x less** |
| Numerical error | 0 | Accumulates | **Infinite** |

**Winner**: Power basis by a huge margin!

---

## Evaluation Cost Analysis

### Power Basis (with Region Rescaling)

**Algorithm**:
```cpp
// Rescale point from [0,1] to [a,b]
t_rescaled = a + t * (b - a)  // 1 multiply, 1 add per axis

// Horner evaluation
result = coeffs[n]
for (i = n-1; i >= 0; i--) {
    result = result * t_rescaled + coeffs[i]
}
```

**Cost**:
- Rescaling: **O(d)** - d multiplies + d adds (d = dimension)
- Evaluation: **O(n)** - n multiplies + n adds
- Total: **O(n + d)**

**Example** (degree 100, 2D):
- Rescaling: 2 multiplies + 2 adds = 4 ops
- Horner: 100 multiplies + 100 adds = 200 ops
- Total: 204 ops

### Bernstein Basis (De Casteljau)

**Algorithm**:
```cpp
// De Casteljau evaluation
for (int r = 1; r <= n; r++) {
    for (int i = 0; i <= n - r; i++) {
        b[i] = (1-t) * b[i] + t * b[i+1];
    }
}
return b[0];
```

**Cost**:
- Time: **O(n²)** - n(n+1)/2 multiplies + n(n+1)/2 adds
- Space: **O(n)** - temporary array

**Example** (degree 100):
- Operations: ~5050 multiplies + ~5050 adds = 10,100 ops

### Comparison

| Metric | Power | Bernstein | Ratio |
|--------|-------|-----------|-------|
| Time | O(n) | O(n²) | **50x faster** (n=100) |
| Operations (n=100) | ~200 | ~10,000 | **50x fewer** |
| Space | O(1) | O(n) | **100x less** |

**Winner**: Power basis (much faster)

---

## Differentiation Cost Analysis

### Power Basis

**Algorithm**:
```cpp
// d/dx (a₀ + a₁x + a₂x² + ... + aₙxⁿ) = a₁ + 2a₂x + ... + naₙx^(n-1)
for (int i = 0; i < n; i++) {
    deriv_coeffs[i] = (i+1) * coeffs[i+1];
}
```

**Cost**:
- Time: **O(n)** - n multiplies
- Space: **O(n)** - n new coefficients
- Numerical error: **Minimal** (one multiply per coefficient)

### Bernstein Basis

**Algorithm**:
```cpp
// d/dx B(x) = n * Σ(b_{i+1} - b_i) * B_{i,n-1}(x)
for (int i = 0; i < n; i++) {
    deriv_coeffs[i] = n * (coeffs[i+1] - coeffs[i]);
}
```

**Cost**:
- Time: **O(n)** - n subtracts + n multiplies
- Space: **O(n)** - n new coefficients
- Numerical error: **Minimal** (one subtract + one multiply per coefficient)

### Comparison

| Metric | Power | Bernstein | Ratio |
|--------|-------|-----------|-------|
| Time | O(n) | O(n) | **Tie** |
| Operations | n multiplies | n subtracts + n multiplies | **Bernstein slightly more** |
| Numerical stability | Excellent | Excellent | **Tie** |

**Winner**: Tie (both excellent)

---

## Division Cost Analysis

### Power Basis (Standard Algorithm)

**Algorithm**: Euclidean division

**Cost**:
- Time: **O(d·e)** where d = deg(divisor), e = deg(dividend)
- Space: **O(e)**
- Complexity: **Simple** (well-known algorithm)

### Bernstein Basis (Homogenization)

**Algorithm**: Busé et al. (2008)

**Cost**:
- Time: **O(d·e)** (similar to power)
- Space: **O(e)**
- Complexity: **Complex** (homogenization, binomial coefficients)

### Comparison

| Metric | Power | Bernstein | Ratio |
|--------|-------|-----------|-------|
| Time | O(d·e) | O(d·e) | **Tie** |
| Implementation complexity | Simple | Complex | **Power simpler** |
| Numerical stability | Excellent | Good | **Power better** |

**Winner**: Power (simpler and better understood)

---

## Convex Hull Property

### Power Basis

**Not available** - coefficients don't bound values

### Bernstein Basis

**Available** - coefficients form convex hull of polynomial values

**Algorithm**:
```cpp
min_value = *std::min_element(coeffs.begin(), coeffs.end());
max_value = *std::max_element(coeffs.begin(), coeffs.end());
// Polynomial values are guaranteed to be in [min_value, max_value]
```

**Cost**: O(n)

**Winner**: Bernstein (only option)

---

## Overall Comparison

| Operation | Power | Bernstein | Winner | Reason |
|-----------|-------|-----------|--------|--------|
| **Subdivision** | O(1) | O(n²) | **Power** | 1000x faster |
| **Evaluation** | O(n) | O(n²) | **Power** | 50x faster |
| **Differentiation** | O(n) | O(n) | **Tie** | Both excellent |
| **Division** | O(d·e) | O(d·e) | **Power** | Simpler |
| **Convex hull** | ❌ | ✅ | **Bernstein** | Only option |

---

## Recommendation

### For Subdivision-Based Solver

**Use Power Basis**:
1. ✅ Subdivision is O(1) vs O(n²)
2. ✅ Evaluation is O(n) vs O(n²)
3. ✅ No numerical error accumulation
4. ✅ Simpler implementation

**Use Bernstein Only When**:
- Need convex hull property for bounds
- Specific algorithm requires Bernstein representation

### Typical Workflow

```cpp
// 1. Start in power basis
auto poly = PolynomialBase<mpreal>::fromPower(degrees, coeffs);

// 2. Subdivide many times (O(1) each!)
for (int i = 0; i < 1000; i++) {
    auto sub = poly.restrictedToInterval(0, a[i], b[i]);
    auto val = sub.evaluate({0.5});
    // ...
}

// 3. Convert to Bernstein only if needed
if (need_convex_hull) {
    auto bern = poly.convertToBernstein();
    auto [min, max] = bern.getConvexHullBounds();
}
```

---

## Summary

**Power basis is superior for subdivision-based solvers**:
- ✅ 1000x faster subdivision
- ✅ 50x faster evaluation
- ✅ No error accumulation
- ✅ Simpler implementation

**Bernstein basis is only needed for**:
- Convex hull property
- Specific algorithms that require it

**Design decision**: Use power basis as primary representation, convert to Bernstein only when explicitly needed.

