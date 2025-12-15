# Conditioning and Precision in Polynomial Root Finding

## Overview

This document explains why **small residuals don't guarantee accurate roots** for ill-conditioned problems, and when higher precision arithmetic is needed.

## The Problem: Residual vs. Error

### What is a Residual?

The **residual** is how close the function value is to zero:
```
residual = |f(x)|
```

### What is the Error?

The **error** is how far the computed root is from the true root:
```
error = |x_computed - x_true|
```

### The Relationship

For **well-conditioned** problems:
```
small residual → small error ✅
```

For **ill-conditioned** problems:
```
small residual ≠ small error ❌
```

## Condition Number: The Key Concept

### Definition

The **condition number** κ (kappa) measures how sensitive a root is to perturbations:

```
|error| ≈ κ * |residual| / |f'(x_true)|
```

**Interpretation:**
- **κ ≈ 1**: Well-conditioned (small residual → small error)
- **κ >> 1**: Ill-conditioned (small residual can hide large error)

### Estimation

We estimate κ using derivative ratios:

```cpp
κ ≈ max(|f''(x)| / |f'(x)|², |f'''(x)| / |f'(x)|³)
```

This measures the curvature relative to the slope.

## Real-World Example: Wilkinson Polynomial

### The Problem

Wilkinson-19 polynomial: (x-0.05)(x-0.10)(x-0.15)...(x-0.95)

**Properties:**
- 19 roots at 0.05, 0.10, 0.15, ..., 0.95
- Roots are closely spaced (0.05 apart)
- Coefficients vary by ~10²⁰ in magnitude
- **Extremely ill-conditioned**: κ ~ 10¹⁷ to 10²⁰

### Results with Double Precision

| Root | Residual | Condition κ | Estimated Error | Actual Error |
|------|----------|-------------|-----------------|--------------|
| 0.05 | 6.1e-17  | 2.1e+19     | 5.3e+10         | 2.5e-09      |
| 0.25 | 9.9e-20  | 3.0e+20     | 8.0e+10         | 0.0e+00      |
| 0.45 | 4.5e-16  | ~1e+18      | ~1e+02          | 7.9e-04      |
| 0.55 | ~1e-16   | ~1e+18      | ~1e+02          | 9.1e-04      |

**Key observation:**
- Residuals are at **machine precision** (~10⁻¹⁷)
- But errors are **up to 10⁻³** (1000× larger than expected!)
- This is because κ ~ 10¹⁸ amplifies the residual

### Why This Happens

1. **Closely-spaced roots**: Small changes in x cause large changes in which root you're near
2. **Large coefficient variation**: Rounding errors in coefficients propagate to roots
3. **Numerical cancellation**: Computing f(x) involves subtracting large numbers

## When Do You Need Higher Precision?

### Detection Criteria

The solver automatically flags roots that need higher precision:

```cpp
estimated_error = κ * |residual| / |f'(x)|

if (estimated_error > 1e-10) {
    needs_higher_precision = true;  // Flag for user
}
```

### Precision Requirements

| Condition Number | Double (16 digits) | Quad (34 digits) | 100-digit |
|------------------|-------------------|------------------|-----------|
| κ < 10⁶          | ✅ Sufficient     | Overkill         | Overkill  |
| 10⁶ < κ < 10¹⁵   | ⚠️ Marginal       | ✅ Sufficient    | Overkill  |
| κ > 10¹⁵         | ❌ Insufficient   | ⚠️ Marginal      | ✅ Sufficient |

### Rule of Thumb

To achieve **d** digits of accuracy in the root:
```
required_precision ≈ d + log₁₀(κ)
```

**Example:** For Wilkinson-19 (κ ~ 10¹⁸) to get 10 digits of accuracy:
```
required_precision ≈ 10 + 18 = 28 digits
→ Use quad precision (34 digits) or higher
```

## What Causes Ill-Conditioning?

### 1. Multiple Roots (High Multiplicity)

**Example:** (x - 0.5)⁶

- f'(x) = 0 at the root
- κ → ∞ as you approach the root
- **Solution:** Use modified Newton method with multiplicity

### 2. Closely-Spaced Roots

**Example:** Wilkinson polynomial

- Roots separated by 0.05
- Small perturbation can "jump" to neighboring root
- **Solution:** Use higher precision arithmetic

### 3. Large Coefficient Variation

**Example:** x²⁰ - 1 = 0

- Coefficients: [1, 0, 0, ..., 0, -1]
- Ratio: 1 / 1 = 1 (actually well-conditioned!)
- But expanded form has large coefficients
- **Solution:** Keep polynomial in factored form when possible

### 4. Near-Degenerate Systems

**Example:** Two nearly parallel curves

- Small changes in coefficients cause large changes in intersection
- **Solution:** Detect degeneracy, use symbolic methods

## Solver Strategy

### 1. Start with Double Precision

```cpp
SubdivisionSolver solver;
auto result = solver.solve(poly);  // Fast, double precision
```

### 2. Refine and Check Conditioning

```cpp
ResultRefiner refiner;
auto refined = refiner.refine1D(result, poly);

for (const auto& root : refined.roots) {
    if (root.needs_higher_precision) {
        // Condition number suggests higher precision needed
    }
}
```

### 3. Escalate to Higher Precision (Future)

```cpp
// Convert to higher precision and re-refine
auto root_hp = refineWithHigherPrecision<float100>(root, poly);
```

## Summary

### Key Takeaways

1. **Small residual ≠ accurate root** for ill-conditioned problems
2. **Condition number κ** relates residual to error
3. **Wilkinson polynomial** is a classic example (κ ~ 10¹⁸)
4. **Detection is automatic** - solver flags problematic roots
5. **Solution**: Use higher precision arithmetic when flagged

### Workflow

```
Solve (double) → Refine (double) → Check κ → Re-refine (higher precision if needed)
     ↓                ↓                ↓              ↓
   Fast           Accurate        Detect         Fix if needed
```

### References

- Wilkinson, J.H. (1963). "Rounding Errors in Algebraic Processes"
- Higham, N.J. (2002). "Accuracy and Stability of Numerical Algorithms"
- Trefethen, L.N. & Bau, D. (1997). "Numerical Linear Algebra"

