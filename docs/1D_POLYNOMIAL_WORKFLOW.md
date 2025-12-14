# 1D Polynomial Root Refinement Workflow

## Overview

This document describes the recommended workflow for refining roots of 1D polynomials with high precision, including handling of multiple roots.

## Key Components

### 1. Multiplicity Detection

**Recommended Method: Taylor Series Ratio Test**

The Taylor series method with ratio test is the gold standard for multiplicity detection:

```cpp
unsigned int m = ResultRefinerHP::estimateMultiplicity(
    x, poly, max_order, threshold, first_nonzero_deriv);
```

**Why it works:**
- For f(x) = c_m(x-r)^m + ..., the ratio |f^(k+1)| / |f^(k)| reveals if f^(k) vanishes
- If ratio > 100, then f^(k) vanishes at the root
- If ratio < 100, then multiplicity = k
- **Scale-independent**: works at all distances from root
- **Robust**: tested on multiplicities 1-5 with excellent results

**Alternative Methods (NOT recommended):**
- Simple threshold: Fails universally - no single threshold works for all cases
- Ostrowski: Frequently gives wrong estimates (off by ±1), causing slow/no convergence
- Sturm sequence: Works when precision is high enough, but less practical than Taylor

### 2. Modified Newton Iteration

For a root of multiplicity m, use modified Newton:

```cpp
x_new = x - m * f(x) / f'(x)
```

**Convergence:**
- With correct multiplicity: quadratic convergence (error² → 0)
- With wrong multiplicity: linear or no convergence

**Example results from convergence tests:**
- m=2, Taylor estimate: 4 iterations to convergence ✓
- m=2, Ostrowski estimate (m=3): 50 iterations, no convergence ✗
- m=3, Taylor estimate: 4 iterations to convergence ✓
- m=4, Taylor estimate: 4 iterations to convergence ✓

### 3. Error Bounds

Use interval Newton method for rigorous error bounds:

```cpp
mpreal error_bound = abs(f / df) * (1 + condition_factor);
```

**Convergence criterion:**
- Exit when error_bound < tolerance (NOT when |f| < tolerance)
- For multiple roots, |f| can be small even when far from root
- Error bounds provide rigorous guarantees

### 4. Recommended Workflow

```cpp
// 1. Start with initial guess
mpreal x = initial_guess;

// 2. Iterate until convergence
for (unsigned int iter = 0; iter < max_iters; ++iter) {
    // Evaluate polynomial and derivative
    mpreal f = poly.evaluate(x);
    mpreal df = dpoly.evaluate(x);
    
    // Estimate multiplicity using Taylor ratio test
    mpreal first_nonzero_deriv;
    unsigned int m = ResultRefinerHP::estimateMultiplicity(
        x, poly, 10, mpreal("1e-50"), first_nonzero_deriv);
    
    if (m == 0) m = 1;  // Safety check
    
    // Compute error bound
    mpreal error_bound = abs(f / df) * (1 + condition_factor);
    
    // Check convergence
    if (error_bound < tolerance) {
        break;  // Converged!
    }
    
    // Modified Newton step
    mpreal step = mpreal(m) * f / df;
    x = x - step;
}
```

### 5. Precision Requirements

**For multiplicity detection:**
- Minimum 256 bits (~77 decimal digits) recommended
- Higher multiplicity requires higher precision
- Taylor ratio test works well with 256-bit precision for m ≤ 5

**For final refinement:**
- 512 bits (~154 decimal digits) for high-accuracy results
- Adjust based on application requirements

## Test Results Summary

### Multiplicity Detection Accuracy

| Method | m=1 | m=2 | m=3 | m=4 | m=5 | Verdict |
|--------|-----|-----|-----|-----|-----|---------|
| **Taylor (ratio test)** | ✓ | ✓ | ✓ | ✓ | ✓ | **BEST** |
| Simple threshold (all) | ✓ | ✗ | ✗ | ✗ | ✗ | FAILS |
| Ostrowski | ~ | ✗ | ~ | ✗ | ~ | UNRELIABLE |
| Sturm sequence | ✓ | ✓ | ✓ | ✓ | ✓ | Works at root only |

### Convergence Performance

Using estimated multiplicity in modified Newton:

| True m | Taylor Iters | Ostrowski Iters | Winner |
|--------|--------------|-----------------|--------|
| 2 | 4 | 50+ (no conv) | **Taylor** |
| 3 | 4 | 6 | **Taylor** |
| 4 | 4 | 50+ (no conv) | **Taylor** |

## Conclusion

**For 1D polynomial root refinement:**
1. Use **Taylor ratio test** for multiplicity detection
2. Use **modified Newton** with estimated multiplicity
3. Use **error bounds** for convergence criterion
4. Use **256-bit precision** minimum for multiplicity detection
5. Avoid simple thresholds and Ostrowski method

This workflow provides robust, efficient root refinement for 1D polynomials with arbitrary multiplicity.

## References

- `tests/test_multiplicity_methods.cpp` - Comprehensive multiplicity detection comparison
- `tests/test_convergence_comparison.cpp` - Convergence rate comparison
- `src/result_refiner_hp.cpp` - Implementation of all methods
- `docs/MULTIPLICITY_DETECTION_COMPARISON.md` - Detailed method comparison

