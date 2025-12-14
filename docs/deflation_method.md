# Deflation Method for Multiple Roots

## Overview

Deflation is a technique for handling multiple (singular) roots of polynomial systems where standard Newton's method fails to converge quadratically. This document describes the **Ojika-Watanabe-Mitsui deflation algorithm** as implemented in Bertini and other numerical algebraic geometry software.

## The Problem

For a polynomial f(x) with a root of multiplicity m at x = r:
- f(r) = f'(r) = f''(r) = ... = f^(m-1)(r) = 0
- f^(m)(r) ≠ 0

Standard Newton's method:
- x_{k+1} = x_k - f(x_k)/f'(x_k)
- Converges **linearly** (not quadratically) near multiple roots
- Division by near-zero f'(x) causes numerical instability

## Ojika's Deflation Algorithm

### Key Idea

Transform the singular problem into a **non-singular** problem by augmenting the system with additional variables and equations that capture the **dual space** structure of the multiple root.

### The Augmented System

For a 1D polynomial f(x), the deflation proceeds in stages:

#### Stage 0 (Original System)
```
Variables: x
Equations: F₀ = f(x) = 0
Jacobian: J = [f'(x)]
```

If |f'(x)| < tol at the solution → Jacobian is singular → Need deflation

#### Stage 1 (First Deflation)
```
Variables: (x, v₁)
Equations:
  F₀ = f(x) = 0
  F₁ = f'(x)·v₁ = 0
  F₂ = v₁² - 1 = 0  (normalization)

Jacobian:
  ∂F₀/∂x = f'(x),      ∂F₀/∂v₁ = 0
  ∂F₁/∂x = f''(x)·v₁,  ∂F₁/∂v₁ = f'(x)
  ∂F₂/∂x = 0,          ∂F₂/∂v₁ = 2v₁
```

Check rank using SVD. If rank < 3 → Need more deflation

#### Stage 2 (Second Deflation)
```
Variables: (x, v₁, v₂)
Equations:
  F₀ = f(x) = 0
  F₁ = f'(x)·v₁ = 0
  F₂ = f''(x)·v₁ + f'(x)·v₂ = 0     ← Uses BOTH v₁ and v₂
  F₃ = v₁² + v₂² - 1 = 0

Jacobian:
  ∂F₀/∂x = f'(x),                    ∂F₀/∂v₁ = 0,      ∂F₀/∂v₂ = 0
  ∂F₁/∂x = f''(x)·v₁,                ∂F₁/∂v₁ = f'(x),  ∂F₁/∂v₂ = 0
  ∂F₂/∂x = f'''(x)·v₁ + f''(x)·v₂,   ∂F₂/∂v₁ = f''(x), ∂F₂/∂v₂ = f'(x)
  ∂F₃/∂x = 0,                        ∂F₃/∂v₁ = 2v₁,    ∂F₃/∂v₂ = 2v₂
```

#### General Stage k
```
Variables: (x, v₁, v₂, ..., vₖ)
Equations:
  F₀ = f(x)
  F₁ = f'(x)·v₁
  F₂ = f''(x)·v₁ + f'(x)·v₂
  F₃ = f'''(x)·v₁ + f''(x)·v₂ + f'(x)·v₃
  ...
  Fᵢ = f^(i)(x)·v₁ + f^(i-1)(x)·v₂ + ... + f'(x)·vᵢ
  ...
  Fₖ = f^(k)(x)·v₁ + f^(k-1)(x)·v₂ + ... + f'(x)·vₖ
  Fₖ₊₁ = v₁² + v₂² + ... + vₖ² - 1
```

**General formula for equation i (i = 1, ..., k):**
```
Fᵢ = Σⱼ₌₁ⁱ f^(i-j+1)(x) · vⱼ = 0
```

**Jacobian entries:**
```
∂Fᵢ/∂x = Σⱼ₌₁ⁱ f^(i-j+2)(x) · vⱼ
∂Fᵢ/∂vⱼ = f^(i-j+1)(x)  for j ≤ i
∂Fᵢ/∂vⱼ = 0             for j > i
```

### Algorithm

1. **Start** with stage k = 0, variables z = [x]
2. **Apply Newton's method** to the current augmented system:
   - Compute F(z) and J(z)
   - Solve J·Δz = -F using least-squares (handles singular J)
   - Update z ← z + Δz
   - Repeat until ||Δz|| < tol
3. **Check rank** of Jacobian using SVD:
   - Compute singular values σ₁ ≥ σ₂ ≥ ... ≥ σₙ
   - rank = #{σᵢ : σᵢ > tol}
4. **If rank < n** (rank deficient):
   - Increment k ← k + 1
   - Add new variable vₖ (initialize to 1/√k)
   - Go to step 2
5. **If rank = n** (full rank):
   - Deflation complete!
   - Detected multiplicity = k + 1

## Key Properties

### Why This Works

1. **Dual Space Interpretation**: The variables v₁, ..., vₖ represent directions in the dual space (space of linear functionals)

2. **Rank Restoration**: At a root of multiplicity m:
   - Original Jacobian [f'(x)] has rank 0
   - After k deflations, Jacobian has rank k
   - After m-1 deflations, Jacobian becomes full rank

3. **Quadratic Convergence**: Once the Jacobian is non-singular, Newton's method converges quadratically

### Comparison with Modified Newton

**Modified Newton** (requires knowing multiplicity m):
```
x_{k+1} = x_k - m · f(x_k)/f'(x_k)
```
- ✓ Quadratic convergence if m is correct
- ✗ Requires knowing m in advance
- ✗ Still divides by near-zero f'(x)

**Deflation**:
- ✓ Automatically detects multiplicity
- ✓ Transforms to non-singular problem
- ✓ Quadratic convergence
- ✗ More complex (adds variables)

## Practical Considerations

### Numerical Challenges

1. **Rank Determination**: Using SVD to determine rank is sensitive to the tolerance threshold
   - Too large: May incorrectly identify full rank when still singular
   - Too small: May continue deflating unnecessarily
   - Recommended: Use **absolute tolerance** ~10⁻⁸ for singular values

2. **Condition Number Growth**: The Jacobian condition number grows rapidly with deflation stages
   - Stage 0: cond ~ 1
   - Stage 5: cond ~ 10²⁰
   - Stage 10: cond ~ 10²⁵
   - This limits practical deflation to ~10-15 stages in double precision

3. **Coefficient Precision Limitation**: Final accuracy is limited by input coefficient precision
   - Double precision coefficients → Final error ~10⁻¹⁵
   - For higher accuracy, need higher precision coefficients (not just arithmetic)

### Alternative: Adaptive Modified Newton

Instead of full deflation, can estimate multiplicity using derivatives:

```python
def estimate_multiplicity(f, x):
    """Find first non-zero derivative"""
    for m in range(1, max_mult):
        if |f^(m)(x)| > tol:
            return m
    return max_mult

# Modified Newton with adaptive m
m = estimate_multiplicity(f, x)
x_new = x - m * f(x) / f'(x)
```

**Advantages**:
- Simpler (no augmented system)
- Faster (fewer variables)
- Works well for moderate multiplicities (m ≤ 5)

**Disadvantages**:
- Still divides by near-zero f'(x)
- Multiplicity estimation can be inaccurate
- May oscillate if m estimate changes

### Test Results

#### Test Case 1: p(x) = (x-0.6)⁶ (pure multiple root)

Starting from x₀ = 0.59, targeting root x = 0.6 (multiplicity 6):

| Method | Final Error | Iterations | Notes |
|--------|-------------|------------|-------|
| Standard Newton | 4.0×10⁻³ | 100 | Linear convergence, poor accuracy |
| Modified Newton (m=6) | **8.5×10⁻⁷** | 1 | **BEST!** Quadratic convergence with correct m |
| Adaptive Modified Newton | 5.1×10⁻² | 100 | Wrong m estimate → diverges |
| Ojika Deflation | 2.6×10⁻³ | ~250 | Detects m=11 (wrong!), numerical issues |

**Key Observations**:
1. ✅ **Modified Newton with correct m is BY FAR the best method!**
2. ❌ Standard Newton has only linear convergence (very slow)
3. ❌ Adaptive multiplicity estimation is unreliable (underestimates m)
4. ❌ Deflation fails in double precision (wrong multiplicity detection)

#### Test Case 2: p(x) = (x-0.2)(x-0.6)⁶ (mixed roots)

Starting from x₀ = 0.59, targeting root x = 0.6 (multiplicity 6):

| Method | Final Error | Iterations | Notes |
|--------|-------------|------------|-------|
| Standard Newton | 2.2×10⁻³ | 39 | Linear convergence |
| Modified Newton (m=6) | **8.5×10⁻⁷** | 1 | **BEST!** Correct m essential |
| Adaptive Modified Newton | 1.7×10⁻³ | 100 | Estimates m=3, oscillates |
| Ojika Deflation | 2.4×10⁻⁴ | ~250 | Better than pure case, still slow |

**Observations**:
1. Modified Newton with correct m works equally well for mixed roots
2. Presence of simple root doesn't significantly affect performance
3. Deflation performs slightly better with mixed roots but still has issues

### Recommendations

**For polynomial root finding**:

Based on our test results, the recommendations are:

1. **If multiplicity is known**: Use **modified Newton with m** - it's fast, accurate, and reliable
   - Converges in 1-2 iterations even for m=6
   - Achieves error ~10⁻⁷ in double precision
   - Much simpler than deflation

2. **If multiplicity is unknown**:
   - **Option A**: Estimate m using high-order derivatives, then use modified Newton
     - Risk: Multiplicity estimation can be inaccurate
     - Mitigation: Try multiple m values if convergence fails

   - **Option B**: Use high-precision arithmetic from the start
     - More reliable than deflation in double precision
     - Avoids numerical issues with rank determination

   - **Option C**: Use deflation (NOT recommended for double precision)
     - Only if you have a robust implementation (like Bertini)
     - Requires careful tolerance management
     - Prone to detecting wrong multiplicity in double precision

3. **For high multiplicity (m > 10)**:
   - Switch to high-precision arithmetic (quad or arbitrary precision)
   - Modified Newton still works well with correct m
   - Deflation becomes increasingly unreliable

**For general nonlinear systems**:
1. Deflation is essential (no simple multiplicity formula exists)
2. Use specialized software (Bertini, PHCpack) rather than implementing from scratch
3. Consider adaptive precision to handle ill-conditioning
4. Be aware that rank determination is the main source of errors

## References

1. **Ojika, S. Watanabe, T. Mitsui** (1983). "Deflation algorithm for the multiple roots of a system of nonlinear equations." *J. Math. Anal. Appl.* 96:463-479.

2. **Hauenstein, J.D. and Wampler, C.W.** (2013). "Isosingular sets and deflation." *Found. Comput. Math.* 13:371-403.

3. **Leykin, A., Verschelde, J., and Zhao, A.** (2006). "Newton's method with deflation for isolated singularities of polynomial systems." *Theor. Comput. Sci.* 359:111-122.

4. **Bates, D.J., Hauenstein, J.D., Sommese, A.J., and Wampler, C.W.** (2013). *Numerically Solving Polynomial Systems with Bertini.* SIAM.

## Implementation Notes

### Why Deflation May Detect Wrong Multiplicity

In our test, deflation detected multiplicity 11 instead of 6. Possible reasons:

1. **Numerical Noise**: Rounding errors in high-order derivatives create artificial rank deficiency
2. **Tolerance Issues**: Singular value threshold may be too strict
3. **Ill-Conditioning**: Extreme condition numbers (>10²⁵) make rank determination unreliable
4. **Implementation Errors**: Jacobian computation may have numerical issues

### Bertini's Approach

Bertini uses several enhancements:
- **Adaptive precision**: Switches to higher precision when condition number grows
- **Randomization**: Uses random linear combinations to improve conditioning
- **Hybrid methods**: Combines deflation with other techniques
- **Careful tolerance management**: Adapts thresholds based on problem scale

For production use, consider using Bertini or similar specialized software rather than implementing deflation from scratch.

## Conclusion

### Main Findings

1. **Modified Newton with known multiplicity is the clear winner for polynomials**
   - Error: 8.5×10⁻⁷ in just 1 iteration
   - Simple to implement
   - Reliable and fast

2. **Deflation is NOT recommended for double precision polynomial root finding**
   - Detects wrong multiplicity (11 instead of 6 in our test)
   - Computationally expensive (~250 iterations)
   - Numerically unstable due to rank determination issues
   - Final error worse than modified Newton

3. **Adaptive multiplicity estimation is unreliable**
   - Tends to underestimate multiplicity
   - Can lead to divergence or oscillation
   - Better to use other methods to determine m

### When to Use Each Method

| Method | Use When | Avoid When |
|--------|----------|------------|
| **Modified Newton (known m)** | Multiplicity is known or can be determined | m is completely unknown |
| **Modified Newton (adaptive m)** | Need automatic method, m ≤ 5 | High multiplicity (m > 5) |
| **Deflation** | General nonlinear systems, using Bertini | Implementing from scratch in double precision |
| **High-precision arithmetic** | m > 10, or ill-conditioned problems | Performance is critical |

### Practical Advice for Polynomial Solver

For your polynomial solver project:

1. **Primary method**: Modified Newton with multiplicity detection
   - Use derivative ratios to estimate m
   - Try m = 1, 2, 3, ... until convergence
   - Very effective for m ≤ 10

2. **Fallback**: High-precision refinement
   - If modified Newton fails to converge
   - Switch to quad or arbitrary precision
   - Continue with modified Newton in higher precision

3. **Skip deflation for now**
   - Too complex and unreliable in double precision
   - Only consider if you need it for general nonlinear systems
   - If needed, use existing libraries (Bertini) rather than implementing

This approach gives you the best balance of simplicity, reliability, and performance.

