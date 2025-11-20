# Differentiation Tools: Design Summary

## Executive Summary

This document summarizes the complexity analysis and design decisions for implementing differentiation tools for Bernstein basis polynomials in the polynomial solver library.

---

## Key Findings

### 1. Complexity Analysis Results

| Scenario | Approach A<br/>(Restrict → Diff) | Approach B<br/>(Diff → Restrict) | Winner | Speedup |
|----------|----------------------------------|----------------------------------|--------|---------|
| **Single derivative (1D)** | O(n²) | O(n²) | B | 1.2× (constants) |
| **k derivatives (1D)** | O(n²) | O(k·n²) | **A** | **k×** |
| **Gradient (d-D)** | O(d·N·n) | O(d²·N·n) | **A** | **d×** |
| **Hessian (d-D)** | O(d²·N·n) | O(d³·N·n) | **A** | **d²×** |

**Legend**:
- n = polynomial degree
- d = dimension
- N = total coefficients
- k = number of derivatives

### 2. Recommended Strategy

✅ **For subdivision solver with degeneracy detection:**

**Use Approach B with global caching (Differentiate → Restrict):**
1. Cache all derivatives globally (once) - **numerically stable**
2. Restrict cached derivatives per box (cheap)
3. Evaluate restricted derivatives at box center

**Complexity**:
- Initial: O(C(d+k,k) · N) - compute all derivatives once
- Per box: O(d · N · n_avg) - restrict to box
- **Total for B boxes**: O(C(d+k,k) · N + B · d · N · n_avg)

**Numerical accuracy**:
- **10-5500× more accurate** than Approach A (restrict→diff)
- Experimental results show Approach B wins for all test cases
- See `docs/numerical_stability_analysis.md` for details

**Alternative (BAD)**: Compute derivatives per box
- **Total for B boxes**: O(B · C(d+k,k) · N)
- **B times slower!**

---

## API Design

### Two-Tier Architecture

```cpp
// Tier 1: Stateless (one-off operations)
class Differentiation {
    static Polynomial derivative(const Polynomial& p, size_t axis, unsigned int order);
    static std::vector<Polynomial> gradient(const Polynomial& p);
    static std::vector<std::vector<Polynomial>> hessian(const Polynomial& p);
};

// Tier 2: Stateful (caching)
class DerivativeCache {
    explicit DerivativeCache(const Polynomial& p);
    const Polynomial& get(const std::vector<unsigned int>& orders);
    const Polynomial& getPartial(size_t axis);
    void precomputeUpToOrder(unsigned int maxOrder);
};
```

### Design Principles

1. **Iterative computation**: Higher-order derivatives computed by iterating first-order
2. **Lazy evaluation**: Compute on-demand, cache results
3. **Optimal for subdivision**: Cache globally, restrict per box
4. **Memory efficient**: Sparse storage (std::map)

---

## Usage Examples

### Example 1: Multiplicity Detection

```cpp
Polynomial p = ...;  // (x-0.6)^6
DerivativeCache cache(p);
std::vector<double> point = {0.6};

int multiplicity = 0;
for (int k = 0; k <= 6; ++k) {
    double dk_val = cache.get(0, k).evaluate(point);
    if (std::abs(dk_val) < 1e-10) {
        multiplicity = k + 1;
    } else {
        break;
    }
}
// Result: multiplicity = 6
```

### Example 2: Subdivision with Degeneracy Detection

```cpp
// At solver initialization
DerivativeCache global_cache(polynomial);
global_cache.precomputeUpToOrder(3);

// For each subdivision box
for (auto& box : boxes) {
    // Restrict cached derivatives (cheap)
    Polynomial f_box = global_cache.get({0,0})
        .restrictedToInterval(0, box.lower[0], box.upper[0])
        .restrictedToInterval(1, box.lower[1], box.upper[1]);
    
    Polynomial fx_box = global_cache.get({1,0})
        .restrictedToInterval(0, box.lower[0], box.upper[0])
        .restrictedToInterval(1, box.lower[1], box.upper[1]);
    
    // Evaluate at box center
    double f_val = f_box.evaluate(box.center);
    double fx_val = fx_box.evaluate(box.center);
    
    // Check degeneracy
    if (std::abs(f_val) < tol && std::abs(fx_val) < tol) {
        // High-multiplicity root detected
    }
}
```

### Example 3: Newton's Method

```cpp
DerivativeCache cache(f);
std::vector<double> x = {0.5, 0.5};

for (int iter = 0; iter < 20; ++iter) {
    std::vector<double> grad = Differentiation::evaluateGradient(cache, x);
    std::vector<std::vector<double>> hess = Differentiation::evaluateHessian(cache, x);
    
    // Solve: H * delta = -grad
    // Update: x = x + delta
}
```

---

## Implementation Plan

### Phase 1: Core Differentiation (Week 1)
- [ ] Implement `Differentiation::derivative()` for 1D
- [ ] Extend to multivariate (tensor-product)
- [ ] Write unit tests for basic differentiation

### Phase 2: Caching (Week 2)
- [ ] Implement `DerivativeCache` class
- [ ] Implement iterative computation logic
- [ ] Add precompute methods
- [ ] Write caching tests

### Phase 3: High-Level API (Week 3)
- [ ] Implement `gradient()` and `hessian()`
- [ ] Implement evaluation helpers
- [ ] Add Jacobian support
- [ ] Write integration tests

### Phase 4: Integration (Week 4)
- [ ] Integrate with subdivision solver
- [ ] Add degeneracy detection using derivatives
- [ ] Update multiplicity example
- [ ] Performance benchmarks

---

## Files to Create

```
include/differentiation.h          - API declarations
src/differentiation.cpp            - Core implementation
src/derivative_cache.cpp           - Cache implementation
tests/test_differentiation.cpp     - Unit tests
tests/test_derivative_cache.cpp    - Cache tests
examples/multiplicity_detection.cpp - Example usage
docs/derivative_complexity_analysis.md - Detailed analysis (✓ created)
docs/derivative_api_design.md      - API documentation (✓ created)
```

---

## Performance Expectations

### Memory Usage

For d-dimensional polynomial with degrees ~n, derivatives up to order k:

| Dimension | Order | # Derivatives | Memory (n=20) |
|-----------|-------|---------------|---------------|
| 1D | 3 | 4 | ~4 KB |
| 2D | 3 | 10 | ~100 KB |
| 3D | 3 | 20 | ~2 MB |
| 2D | 5 | 21 | ~200 KB |

### Computation Time

For Wilkinson polynomial (degree 19, 1D):
- Compute all derivatives up to order 3: ~0.1 ms
- Restrict to 100 boxes: ~10 ms
- **Total**: ~10 ms for 100 boxes with 4 derivatives each

For 2D polynomial (degrees 10×10):
- Compute gradient + Hessian: ~1 ms
- Restrict to 100 boxes: ~50 ms
- **Total**: ~50 ms for 100 boxes with 6 derivatives each

---

## Next Steps

1. **Review this design** with the team
2. **Implement Phase 1** (core differentiation)
3. **Test on extreme cases** (Wilkinson, multiplicity)
4. **Integrate with solver** for degeneracy detection
5. **Benchmark performance** on real problems

---

## References

- `docs/derivative_complexity_analysis.md` - Detailed complexity analysis
- `docs/derivative_api_design.md` - Complete API specification
- Farin, G. "Curves and Surfaces for CAGD" - Bernstein basis theory
- De Casteljau algorithm - Subdivision and evaluation

