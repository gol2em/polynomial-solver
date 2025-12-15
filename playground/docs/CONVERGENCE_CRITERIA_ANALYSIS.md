# Convergence Criteria Analysis: Residual-Only vs. Condition-Aware

## Problem Statement

The current refiner uses **residual-only convergence**: Newton's method stops when `|f(x)| < residual_tolerance`.

**This is insufficient for ill-conditioned problems.**

## Demonstration: (x - 0.5)³

True root: `x = 0.5` (multiplicity 3)

### Current Behavior

```
Root 1:
  Location: x = 0.4999914745418798
  Actual Error: |x - 0.5| = 8.5255e-06  ← LARGE ERROR
  Residual: |f(x)| = 6.1966e-16         ← TINY RESIDUAL (machine epsilon)
  Condition estimate: 5.79e+29           ← ASTRONOMICAL
  Ratio (error/residual): 1.38e+10       ← RED FLAG!
  
Root 2:
  Location: x = 0.5000057753103395
  Actual Error: |x - 0.5| = 5.7753e-06  ← LARGE ERROR
  Residual: |f(x)| = 1.9263e-16         ← TINY RESIDUAL (machine epsilon)
  Condition estimate: 5.99e+30           ← ASTRONOMICAL
  Ratio (error/residual): 3.00e+10       ← RED FLAG!
```

### The Issue

✅ **Good**: System detects ill-conditioning via `needs_higher_precision = true`  
❌ **Bad**: Convergence criterion accepts roots with large errors  
❌ **Bad**: Warning comes AFTER convergence  
❌ **Bad**: Users might miss the warning

## Why Residual-Only Fails

For ill-conditioned problems:
```
|error| ≈ κ × |residual| / |f'(x)|
```

When κ is large (e.g., 10²⁹):
- Residual can be at machine epsilon (10⁻¹⁶)
- But error can still be large (10⁻⁶)
- Newton's method "converges" to wrong answer

## Current Workflow Detection

The system DOES detect ill-conditioning:

1. **After refinement**: Estimates condition number κ
2. **Calculates estimated error**: `κ × |residual| / |f'(x)|`
3. **Flags if needed**: `needs_higher_precision = (estimated_error > 1e-10)`

**But this happens AFTER accepting the root!**

## Proposed Solutions

### Option 1: Condition-Aware Convergence (Recommended)

Check BOTH residual AND estimated error during Newton iteration:

```cpp
// Current (residual-only)
if (std::abs(f) < config.residual_tolerance) {
    return true;  // Accept root
}

// Proposed (condition-aware)
if (std::abs(f) < config.residual_tolerance) {
    // Estimate condition number
    double kappa = estimateConditionNumber1D(x, poly, df);
    double est_error = kappa * std::abs(f) / std::max(std::abs(df), 1e-14);
    
    if (est_error < config.target_tolerance) {
        return true;  // Accept root
    } else {
        // Residual is small but estimated error is large
        // Flag as needing higher precision
        return false;  // Reject root, mark as unverified
    }
}
```

**Benefits:**
- Rejects roots with small residual but large estimated error
- Forces user to use higher precision for ill-conditioned problems
- More robust convergence criterion

**Drawbacks:**
- Slightly more expensive (condition number estimation per iteration)
- May reject more roots (but correctly!)

### Option 2: Stricter Residual Tolerance

Use much stricter residual tolerance for ill-conditioned problems:

```cpp
// Adaptive residual tolerance based on condition number
double adaptive_tol = config.residual_tolerance / std::max(1.0, kappa / 1e10);
if (std::abs(f) < adaptive_tol) {
    return true;
}
```

**Benefits:**
- Simple to implement
- Automatically tightens tolerance for ill-conditioned problems

**Drawbacks:**
- May not converge at all for very ill-conditioned problems
- Doesn't solve fundamental issue (residual ≠ error)

### Option 3: Early Warning System (Current + Enhancement)

Keep current behavior but add prominent warnings:

```cpp
// During refinement, track ill-conditioned roots
if (needs_higher_precision) {
    std::cerr << "WARNING: Root at x=" << x << " needs higher precision!\n";
    std::cerr << "  Residual: " << residual << "\n";
    std::cerr << "  Estimated error: " << est_error << "\n";
    std::cerr << "  Condition number: " << kappa << "\n";
}
```

**Benefits:**
- No change to convergence behavior
- Users get immediate feedback

**Drawbacks:**
- Doesn't prevent accepting inaccurate roots
- Users might ignore warnings

## Recommendation

**Implement Option 1: Condition-Aware Convergence**

This is the most robust solution:
1. Check residual for fast convergence
2. If residual is small, estimate condition number
3. Check if estimated error is acceptable
4. Reject root if error is too large
5. Mark as unverified, requiring higher precision

This ensures that **residual-only convergence never accepts inaccurate roots**.

## Cost Analysis

### Current Condition Number Estimation

`estimateConditionNumber1D()` computes:
- **f'(x)**: Already available from Newton iteration (FREE)
- **f''(x)**: Requires 1 derivative + 1 evaluation
- **f'''(x)**: Requires 1 derivative + 1 evaluation

**Total cost**: 2 derivatives + 2 evaluations per check

### Performance Impact

For **well-conditioned simple roots**:
- Newton converges in 1-3 iterations
- Condition check happens 1-3 times
- Total overhead: ~6 derivative operations
- **Negligible** compared to solver cost

For **ill-conditioned problems**:
- Newton may take many iterations
- But these are the cases we NEED to detect
- Cost is justified to prevent accepting wrong roots

### Cheaper Alternatives

#### Option A: Check Only f''(x) (Recommended)

```cpp
// Simplified condition check using only second derivative
double f_double_prime = ddpoly.evaluate(x);
double kappa_simple = std::abs(f_double_prime) / (std::abs(df) * std::abs(df));
double est_error = kappa_simple * std::abs(f) / std::abs(df);
```

**Cost**: 1 derivative + 1 evaluation (50% reduction)
**Accuracy**: Good for most cases (f''' check is conservative)

#### Option B: Adaptive Checking

```cpp
// Only check condition number if residual is suspiciously small
if (std::abs(f) < 1e-14) {
    // Residual at machine epsilon - check condition
    double kappa = estimateConditionNumber1D(x, poly, df);
    double est_error = kappa * std::abs(f) / std::abs(df);
    if (est_error > config.target_tolerance) {
        return false;  // Reject
    }
}
```

**Cost**: Only pay for condition check when residual is very small
**Benefit**: Zero overhead for well-conditioned problems

#### Option C: Use Newton Step Size as Proxy

```cpp
// Track Newton step size
double step = std::abs(x_new - x);

// If residual is small but step is still large, problem is ill-conditioned
if (std::abs(f) < config.residual_tolerance) {
    if (step > 100.0 * std::abs(f)) {
        // Step >> residual indicates ill-conditioning
        return false;  // Reject, needs higher precision
    }
    return true;  // Accept
}
```

**Cost**: FREE (step size already computed)
**Drawback**: Less accurate, may miss some ill-conditioned cases

### Recommended Approach

**Hybrid: Option B (Adaptive) + Option A (Simplified)**

```cpp
if (std::abs(f) < config.residual_tolerance) {
    // Only check condition if residual is very small (< 1e-12)
    if (std::abs(f) < 1e-12) {
        // Simplified condition check using only f''
        Polynomial ddpoly = Differentiation::derivative(dpoly, 0, 1);
        double f_double_prime = ddpoly.evaluate(x);
        double kappa = std::abs(f_double_prime) / (std::abs(df) * std::abs(df));
        double est_error = kappa * std::abs(f) / std::abs(df);

        if (est_error > config.target_tolerance) {
            return false;  // Reject
        }
    }
    return true;  // Accept
}
```

**Benefits**:
- Zero overhead for well-conditioned problems (residual > 1e-12)
- 50% cheaper condition check (only f'', not f''')
- Still catches ill-conditioned problems
- Optimal cost/benefit tradeoff

## Implementation Impact

- Modify `ResultRefiner::refineRoot1D()` convergence check
- Add adaptive condition number estimation during iteration
- Add `estimated_error_tolerance` parameter to `RefinementConfig`
- Update documentation to explain condition-aware convergence

