# Differentiation Implementation Summary

## Overview

Implemented a comprehensive differentiation API for Bernstein polynomials with caching support. The implementation computes and caches derivative coefficients without applying restriction to regions, following the optimal **Approach B (Differentiate → Restrict)** strategy validated by numerical stability analysis.

---

## Files Created

### 1. `include/differentiation.h` (150 lines)

**Two main classes:**

#### `Differentiation` (Static Utility Class)
- `derivative(p, axis, order)` - Compute ∂^order f / ∂x_axis^order
- `gradient(p)` - Compute [∂f/∂x_0, ..., ∂f/∂x_{n-1}]
- `hessian(p)` - Compute matrix H[i][j] = ∂²f/∂x_i∂x_j
- `differentiateAxis(p, axis)` - Core function: first derivative along one axis

#### `DerivativeCache` (Caching Class)
- `DerivativeCache(p)` - Constructor: stores original polynomial
- `get(orders)` - Get derivative by multi-index {k_0, k_1, ..., k_{n-1}}
- `getPartial(axis, order)` - Convenience: ∂^order f / ∂x_axis^order
- `precomputeUpToOrder(maxOrder)` - Precompute all derivatives up to total order k
- `dimension()` - Get polynomial dimension

**Key design features:**
- Sparse caching using `std::map<std::vector<unsigned int>, Polynomial>`
- Lazy evaluation: derivatives computed only when requested
- Iterative computation: higher-order derivatives via repeated first-order differentiation

---

### 2. `src/differentiation.cpp` (286 lines)

**Core algorithm: Bernstein derivative formula**

For a Bernstein polynomial of degree n with coefficients b_0, ..., b_n:
```
Derivative has degree n-1 with coefficients: d_i = n * (b_{i+1} - b_i)
```

**Implementation details:**

1. **`differentiateAxis(p, axis)`** - O(N) where N = total coefficients
   - Iterates over all 1D slices along the specified axis
   - Applies Bernstein formula: `d[i] = n * (b[i+1] - b[i])`
   - Reduces degree by 1 along the differentiation axis
   - Handles tensor-product layout (last dimension varies fastest)

2. **`derivative(p, axis, order)`** - Iterative application
   - Calls `differentiateAxis()` repeatedly `order` times
   - Efficient for moderate orders (1-5)

3. **`gradient(p)`** - Computes all first-order partials
   - Returns vector of size `dimension`
   - Each element is `differentiateAxis(p, i)`

4. **`hessian(p)`** - Computes all second-order partials
   - Returns `dimension × dimension` matrix
   - H[i][j] = `differentiateAxis(differentiateAxis(p, i), j)`

5. **`DerivativeCache::computeIfNeeded(orders)`** - Recursive caching
   - Checks if derivative already cached
   - If not, finds previous derivative (reduce one order)
   - Differentiates once more and caches result
   - Ensures all derivatives computed via iterative path

---

### 3. `tests/test_differentiation.cpp` (377 lines)

**Comprehensive test suite with 7 tests:**

| Test | Description | Validates |
|------|-------------|-----------|
| 1 | Derivative of x³ | 1st, 2nd, 3rd derivatives |
| 2 | Derivative of x⁵ | Higher-degree polynomials |
| 3 | Gradient of x² + y² | Multivariate first derivatives |
| 4 | DerivativeCache for x⁴ | Caching and reference stability |
| 5 | Mixed partials of x²y² | ∂²f/∂x∂y computation |
| 6 | Hessian of x² + y² | Full Hessian matrix |
| 7 | Precompute up to order 3 | Batch precomputation |

**All tests PASS ✅**

---

## Mathematical Background

### Bernstein Derivative Formula

For univariate Bernstein polynomial:
```
B(t) = Σ_{i=0}^n b_i * B_i^n(t)
```

The derivative is:
```
B'(t) = n * Σ_{i=0}^{n-1} (b_{i+1} - b_i) * B_i^{n-1}(t)
```

### Multivariate Case

For tensor-product polynomial f(x_0, ..., x_{d-1}):
- Partial derivative ∂f/∂x_k: differentiate along axis k, keep others unchanged
- Mixed partial ∂²f/∂x_i∂x_j: differentiate along i, then along j
- Gradient: [∂f/∂x_0, ..., ∂f/∂x_{d-1}]
- Hessian: H[i][j] = ∂²f/∂x_i∂x_j (symmetric matrix)

---

## Numerical Stability

**Why Approach B (Differentiate → Restrict) is optimal:**

1. **Differentiation in Bernstein basis is stable**
   - Formula `d_i = n * (b_{i+1} - b_i)` involves only one subtraction
   - Bernstein coefficients are well-spaced (convex hull property)
   - Error per operation: O(ε · ||b||)

2. **Experimental validation**
   - Test 1 (x⁵): Approach B is 4-46× more accurate
   - Test 2 (degree 20): Approach B is 1.6-5500× more accurate
   - Test 3 (ill-conditioned): Approach B is 2.4-11× more accurate

3. **Computational efficiency**
   - Differentiation: O(n) per operation
   - Subdivision: O(n²) per operation
   - Caching derivatives globally: k× faster for k derivatives

---

## Usage Examples

### Example 1: One-off differentiation
```cpp
Polynomial p = Polynomial::fromPower({3}, {0, 0, 0, 1});  // x^3
Polynomial dp = Differentiation::derivative(p, 0, 1);     // 3x^2
Polynomial d2p = Differentiation::derivative(p, 0, 2);    // 6x
```

### Example 2: Gradient computation
```cpp
Polynomial f = /* f(x,y) = x^2 + y^2 */;
std::vector<Polynomial> grad = Differentiation::gradient(f);
// grad[0] = ∂f/∂x = 2x
// grad[1] = ∂f/∂y = 2y
```

### Example 3: Caching for multiple derivatives
```cpp
Polynomial p = /* some polynomial */;
DerivativeCache cache(p);

// Compute derivatives as needed (cached automatically)
const Polynomial& dp = cache.getPartial(0, 1);
const Polynomial& d2p = cache.getPartial(0, 2);
const Polynomial& d3p = cache.getPartial(0, 3);

// Getting again returns same reference (no recomputation)
const Polynomial& dp_again = cache.getPartial(0, 1);
assert(&dp == &dp_again);  // Same object!
```

### Example 4: Precompute for subdivision solver
```cpp
// At solver initialization
DerivativeCache global_cache(polynomial);
global_cache.precomputeUpToOrder(3);  // For multiplicity detection

// For each subdivision box (many times)
for (auto& box : boxes) {
    // Restrict cached derivatives to box (Approach B)
    Polynomial f_box = global_cache.get({0})
        .restrictedToInterval(0, box.lower[0], box.upper[0]);
    Polynomial df_box = global_cache.get({1})
        .restrictedToInterval(0, box.lower[0], box.upper[0]);
    
    // Evaluate for degeneracy detection
    double f_val = f_box.evaluate(box.center);
    double df_val = df_box.evaluate(box.center);
    
    if (std::abs(f_val) < tol && std::abs(df_val) < tol) {
        // Potential high-multiplicity root!
    }
}
```

---

## Next Steps

1. **Apply to degeneracy detection** in subdivision solver
2. **Implement multiplicity detection** using derivative vanishing
3. **Newton refinement** using gradient/Hessian
4. **Update README** with differentiation API documentation

---

## Build and Test

```bash
cd build
cmake ..
make differentiation test_differentiation
./bin/test_differentiation
ctest  # Run all tests
```

**Result:** All 13 tests pass ✅

