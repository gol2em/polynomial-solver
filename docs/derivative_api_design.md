# Differentiation API Design for Bernstein Polynomials

## Overview

This document describes the design of differentiation tools for Bernstein basis polynomials, optimized for iterative computation and caching.

---

## Design Principles

1. **Iterative-friendly**: Compute higher-order derivatives by iterating first-order derivatives
2. **Caching**: Store computed derivatives to avoid recomputation
3. **Lazy evaluation**: Only compute what's needed
4. **Optimal for subdivision**: Cache globally, restrict per box

---

## API Structure

### Core Classes

```
Differentiation (static utility)
├── derivative()           - Stateless differentiation
├── gradient()            - All first-order partials
├── hessian()             - All second-order partials
├── evaluateGradient()    - Evaluate gradient at point
├── evaluateHessian()     - Evaluate Hessian at point
└── createCache()         - Factory for DerivativeCache

DerivativeCache (stateful caching)
├── get(orders)           - Get derivative by multi-index
├── getPartial(axis)      - Get ∂f/∂xᵢ
├── isCached()            - Check if cached
├── precompute()          - Precompute derivatives
└── clearCache()          - Clear cache
```

---

## Usage Patterns

### Pattern 1: One-off Differentiation

**Use case**: Need a single derivative, no caching needed

```cpp
Polynomial p = ...;

// First derivative
Polynomial dp = Differentiation::derivative(p, 0);

// Second derivative (iterative)
Polynomial d2p = Differentiation::derivative(p, 0, 2);

// Gradient (2D)
std::vector<Polynomial> grad = Differentiation::gradient(p);
```

**Complexity**: O(N) per derivative, where N = total coefficients

---

### Pattern 2: Multiple Derivatives with Caching

**Use case**: Need multiple derivatives of same polynomial

```cpp
Polynomial p = ...;  // degree 20 polynomial
DerivativeCache cache(p);

// Compute as needed (cached automatically)
const Polynomial& dp = cache.getPartial(0);      // ∂f/∂x
const Polynomial& d2p = cache.get(0, 2);         // ∂²f/∂x²
const Polynomial& d3p = cache.get(0, 3);         // ∂³f/∂x³

// Mixed derivatives (2D)
const Polynomial& dxy = cache.get({1, 1});       // ∂²f/∂x∂y
const Polynomial& d2xy = cache.get({2, 1});      // ∂³f/∂x²∂y
```

**Complexity**: O(N) per unique derivative, O(1) for cached access

---

### Pattern 3: Precompute All Derivatives

**Use case**: Know you'll need many derivatives, compute upfront

```cpp
Polynomial p = ...;
DerivativeCache cache(p);

// Precompute all derivatives up to order 3
cache.precomputeUpToOrder(3);

// Now all derivatives are cached:
// Order 0: f
// Order 1: ∂f/∂x, ∂f/∂y
// Order 2: ∂²f/∂x², ∂²f/∂x∂y, ∂²f/∂y²
// Order 3: ∂³f/∂x³, ∂³f/∂x²∂y, ∂³f/∂x∂y², ∂³f/∂y³

// Access is O(1)
const Polynomial& d3x = cache.get({3, 0});
```

**Complexity**: O(C(d+k, k) · N) where d=dimension, k=max order

---

### Pattern 4: Subdivision with Derivatives (RECOMMENDED)

**Use case**: Subdivision solver with degeneracy detection

```cpp
// At solver initialization
Polynomial p = ...;
DerivativeCache global_cache(p);
global_cache.precomputeUpToOrder(3);  // For multiplicity up to 3

// For each subdivision box
for (auto& box : subdivision_boxes) {
    // Restrict cached derivatives to box (cheap)
    Polynomial f_box = global_cache.get({0, 0})
        .restrictedToInterval(0, box.lower[0], box.upper[0])
        .restrictedToInterval(1, box.lower[1], box.upper[1]);
    
    Polynomial fx_box = global_cache.get({1, 0})
        .restrictedToInterval(0, box.lower[0], box.upper[0])
        .restrictedToInterval(1, box.lower[1], box.upper[1]);
    
    Polynomial fy_box = global_cache.get({0, 1})
        .restrictedToInterval(0, box.lower[0], box.upper[0])
        .restrictedToInterval(1, box.lower[1], box.upper[1]);
    
    // Evaluate at box center for degeneracy check
    std::vector<double> center = box.center;
    double f_val = f_box.evaluate(center);
    double fx_val = fx_box.evaluate(center);
    double fy_val = fy_box.evaluate(center);
    
    // Check degeneracy
    if (std::abs(f_val) < tol && 
        std::abs(fx_val) < tol && 
        std::abs(fy_val) < tol) {
        // Potential high-multiplicity root
        // Check higher derivatives...
    }
}
```

**Complexity**:
- Initial caching: O(C(d+k, k) · N) - done once
- Per box: O(d · N · n_avg) - restriction cost
- Total for B boxes: O(C(d+k, k) · N + B · d · N · n_avg)

**Alternative (BAD)**: Computing derivatives per box
- Per box: O(C(d+k, k) · N) - recompute derivatives
- Total for B boxes: O(B · C(d+k, k) · N)
- **B times slower!**

---

### Pattern 5: Multiplicity Detection

**Use case**: Determine multiplicity of a root

```cpp
Polynomial p = ...;  // Suspected root at x = 0.6
DerivativeCache cache(p);
std::vector<double> point = {0.6};

const double tol = 1e-10;
int multiplicity = 0;

// Check derivatives until one is non-zero
for (int k = 0; k <= 6; ++k) {
    const Polynomial& dk = cache.get(0, k);
    double dk_val = dk.evaluate(point);
    
    if (std::abs(dk_val) < tol) {
        multiplicity = k + 1;
    } else {
        break;  // Found first non-zero derivative
    }
}

std::cout << "Multiplicity: " << multiplicity << std::endl;
```

---

### Pattern 6: Newton's Method with Hessian

**Use case**: Newton's method for root refinement

```cpp
Polynomial f = ...;  // 2D polynomial
DerivativeCache cache(f);
std::vector<double> x = {0.5, 0.5};  // Initial guess

const int max_iter = 20;
const double tol = 1e-12;

for (int iter = 0; iter < max_iter; ++iter) {
    // Evaluate function and derivatives at current point
    double f_val = cache.get({0, 0}).evaluate(x);
    
    std::vector<double> grad = 
        Differentiation::evaluateGradient(cache, x);
    
    std::vector<std::vector<double>> hess = 
        Differentiation::evaluateHessian(cache, x);
    
    // Check convergence
    double grad_norm = std::sqrt(grad[0]*grad[0] + grad[1]*grad[1]);
    if (grad_norm < tol) {
        break;
    }
    
    // Solve: H * delta = -grad
    // ... (linear system solver)
    
    // Update: x = x + delta
    // ...
}
```

---

## Implementation Details

### Iterative Computation

Higher-order derivatives computed iteratively:

```cpp
// To compute ∂³f/∂x²∂y:
// Option 1: ∂/∂y(∂²f/∂x²)
// Option 2: ∂/∂x(∂²f/∂x∂y)
// Option 3: ∂²/∂x²(∂f/∂y)

// Choose the path with most cached intermediates
```

### Cache Storage

```cpp
class DerivativeCache {
private:
    Polynomial poly_;
    std::map<std::vector<unsigned int>, Polynomial> cache_;
    
    // Key: multi-index [k₀, k₁, ..., k_{d-1}]
    // Value: ∂^(k₀+k₁+...)/∂x₀^k₀∂x₁^k₁...
};
```

### Memory Usage

For d-dimensional polynomial, derivatives up to order k:

- **Number of derivatives**: C(d+k, k) = (d+k)! / (d! · k!)
- **Examples**:
  - 1D, order 3: 4 polynomials
  - 2D, order 3: 10 polynomials
  - 3D, order 3: 20 polynomials
  - 2D, order 5: 21 polynomials

---

## Files to Create

1. **include/differentiation.h** - API declarations
2. **src/differentiation.cpp** - Implementation
3. **tests/test_differentiation.cpp** - Unit tests
4. **examples/multiplicity_detection.cpp** - Example usage

---

## Next Steps

1. Implement core `Differentiation::derivative()` for 1D
2. Extend to multivariate (tensor-product)
3. Implement `DerivativeCache` with iterative computation
4. Add gradient/Hessian helpers
5. Write comprehensive tests
6. Integrate with subdivision solver for degeneracy detection

