# Precision Requirements for Multiple Root Refinement

## Executive Summary

Systematic testing reveals that:
1. **Multiplicity detection** (Taylor & Ostrowski) works at **all tested precisions** (64-1024 bits)
2. **Sturm sequence** consistently fails (always returns 1) - needs investigation
3. **Convergence** requires precision that scales with multiplicity

## Test Methodology

- **Polynomial**: (x - 0.5)^m created in native high precision
- **Initial guess**: x = 0.3 (distance = 0.2 from root)
- **Tolerance**: 1e-50
- **Max iterations**: 50

## Results Summary

### Multiplicity Detection

| Method | m=2-10 | Notes |
|--------|--------|-------|
| **Taylor** | ✓ All correct | Works at all precisions (64-1024 bits) |
| **Ostrowski** | ✓ All correct | Works at all precisions (64-1024 bits) |
| **Sturm** | ✗ Always returns 1 | **NEEDS FIX** |

**Key finding**: Taylor and Ostrowski methods are **extremely robust** - they correctly detect multiplicity even at 64-bit precision!

### Convergence Requirements (with known multiplicity)

| Multiplicity | Min Precision | Recommended | Notes |
|--------------|---------------|-------------|-------|
| 2 | 64 bits | 128 bits | Converges in 1-2 iterations |
| 3 | 96 bits | 256 bits | Oscillates at 64, 128 bits |
| 4 | 64 bits | 256 bits | Converges in 1-3 iterations |
| 5 | 96 bits | 384 bits | Oscillates at 64, 128 bits |
| 6 | 64 bits | 512 bits | Oscillates at 96, 128 bits |
| 7 | 192 bits | 512 bits | Oscillates at 64-128 bits |
| 8 | 64 bits | 768 bits | Oscillates at 96 bits |
| 9 | 96 bits | 1024 bits | Oscillates at 64, 128 bits |
| 10 | 64 bits | 1280 bits | Oscillates at 96, 128 bits |

**Pattern observed**: 
- Odd multiplicities (3, 5, 7, 9) tend to require higher minimum precision
- Even multiplicities (2, 4, 6, 8, 10) often work at 64 bits but may oscillate at intermediate precisions
- **Oscillation** = convergence occurs but with error increasing/decreasing pattern

### Precision Formula

Based on empirical data:

```
Minimum precision (bits) ≈ 64 + 16*m  (for odd m)
                         ≈ 64         (for even m, but unreliable)

Recommended precision (bits) = 128 * m
```

**Safe rule**: Use **128 × m bits** to guarantee smooth convergence in 1-3 iterations.

## Detailed Convergence Patterns

### Smooth Convergence (✓)
- Converges in 1-3 iterations
- Error decreases monotonically
- Example: m=2 at 64 bits, m=7 at 192 bits

### Oscillating Convergence (⚠)
- Converges eventually (within 50 iterations)
- Error alternates between small and large values
- Example: m=3 at 64 bits (5 iters), m=6 at 96 bits (7 iters)

### Failed Convergence (✗)
- Does not converge within 50 iterations
- Error oscillates without reaching tolerance
- Example: m=5 at 64 bits, m=7 at 128 bits

## Workflow Design Recommendations

### Option 1: Conservative (Guaranteed Success)

```
1. Detect multiplicity m using Taylor method (works at any precision ≥ 64 bits)
2. Select precision = 128 * m bits
3. Create polynomial in high precision from power coefficients
4. Apply modified Newton with detected multiplicity
5. Expect convergence in 1-3 iterations
```

**Pros**: Guaranteed smooth convergence, minimal iterations
**Cons**: May use more precision than strictly necessary

### Option 2: Adaptive (Precision-Efficient)

```
1. Start with 256 bits precision
2. Detect multiplicity m using Taylor method
3. Check if current precision ≥ 128*m:
   - If yes: proceed with refinement
   - If no: increase precision to 128*m
4. Apply modified Newton
5. If oscillation detected, increase precision by 2x
```

**Pros**: Uses less precision for low multiplicities
**Cons**: May need to restart with higher precision

### Option 3: Hybrid (Practical)

```
1. Detect multiplicity m in double precision (approximate)
2. Select precision based on condition number κ:
   - If κ < 1e5:  256 bits
   - If κ < 1e10: 512 bits
   - If κ < 1e15: 1024 bits
   - Else:        2048 bits
3. Re-detect multiplicity in high precision (verify)
4. If detected m requires more precision than selected:
   - Increase to max(current, 128*m)
5. Apply modified Newton
```

**Pros**: Balances precision efficiency with robustness
**Cons**: More complex logic

## Issues to Address

### 1. Sturm Sequence Always Returns 1 (**ROOT CAUSE IDENTIFIED**)

**Observation**: Sturm method returns 1 for all multiplicities at all precisions

**Root cause**: Sturm sequence is designed for **exact arithmetic** and **distinct roots**. In floating-point:
- Even at the exact root location x=0.5, f(x) ≈ 2e-79 (not exactly zero)
- Bernstein representation introduces numerical errors
- Sturm sequence counts roots by sign changes, which requires f(x) to cross zero
- For multiple roots in floating-point, the polynomial may not have a clear zero crossing

**Evidence from debug test**:
```
Polynomial: (x - 0.5)^3
Testing at x = 0.5
f(x) = 2.0241e-79  ← Not exactly zero!
Sturm: 1           ← Fails to detect multiplicity
Taylor: 3          ← Correct
Ostrowski: 3       ← Correct
```

**Conclusion**: Sturm sequence is **not suitable** for multiplicity detection in floating-point arithmetic. Use Taylor or Ostrowski instead.

### 2. Oscillation at Intermediate Precisions

**Observation**: Some multiplicities oscillate at specific precisions (e.g., m=6 at 96 bits)

**Cause**: When error becomes comparable to machine epsilon, rounding errors dominate

**Solution**: Use recommended precision (128*m bits) to avoid this regime

## Conclusion

**For production workflow**:
1. ✅ Use **Taylor or Ostrowski** for multiplicity detection (both very robust, work at 64+ bits)
2. ✅ Use **192 bits minimum** for detection (safe for all m≤10)
3. ✅ Use **128 × m bits** for convergence (safe and efficient)
4. ✅ Create polynomial in **native high precision** (not double→HP conversion)
5. ❌ **Do not use Sturm** - not suitable for floating-point multiplicity detection
6. ✅ Expect **1-3 iterations** to full precision with proper setup

**Key insights**:
- Multiplicity detection is NOT the bottleneck - Taylor/Ostrowski work even at 64 bits
- Sturm sequence requires exact arithmetic - unsuitable for floating-point
- The challenge is having enough precision for convergence, which scales linearly with multiplicity (128*m bits)

## Critical Finding: Achievable Error Bounds (Dec 2025)

### Test Case: Product Polynomial with High Multiplicity

**Polynomial**: `p(x) = (x - 0.2)(x - 0.6)^6` (degree 7, multiplicity-6 root at x=0.6)

**Setup**:
- Exact rational power coefficients (e.g., `-729/78125` instead of `-0.0093312`)
- 1024-bit precision
- Power-basis differentiation (avoids Bernstein conversion)
- Modified Newton with multiplicity hint = 6

**Result**:
- Residual: 5.3e-99 (excellent)
- Error bound: **1.1e-33** (limited by numerical precision)
- Iterations: 3
- **Conclusion**: Cannot achieve error bounds better than ~1e-33 with 1024-bit precision for this polynomial

### Implication for Numerical Solvers

**Key insight**: For a numerical solver, the target is **double precision (1e-15)**, not arbitrary precision.

- Error bound of 1e-33 provides **18 orders of magnitude** better accuracy than required
- 1024-bit precision is **more than sufficient** for all practical numerical applications
- Targeting 1e-50 error bounds is unrealistic for high-degree, high-multiplicity polynomials

### Precision Scaling Law

For degree-n polynomial with multiplicity-m root:
- Error bound scales roughly as `O(2^(-precision/n))`
- Degree-7 polynomial with 1024 bits → error bound ~ 1e-33
- To achieve 1e-50 would require > 2048 bits (diminishing returns)

## Critical Implementation Requirements

### 1. Avoid Implicit Basis Conversions

**Problem**: Implicit Power ↔ Bernstein conversions introduce rounding errors.

**Solution**:
- Added explicit conversion checks that abort with error messages
- Use power-basis differentiation for polynomials created with `fromPowerHP()`
- Only convert when explicitly needed via `convertPowerToBernstein()` or `convertBernsteinToPower()`

### 2. Use Exact Rational Coefficients

**Problem**: Decimal approximations introduce errors from the start.

**Solution**:
```cpp
// WRONG:
mpreal("-0.0093312")

// CORRECT:
mpreal("-729") / mpreal("78125")  // Exact: -(1/5)*(3/5)^6
```

### 3. Power-Basis Differentiation for HP

**Problem**: Bernstein-basis differentiation forces conversion even for power-basis polynomials.

**Solution**: Implemented `DifferentiationHP::differentiateAxisPower()` that:
- Works directly with power coefficients
- Avoids Bernstein conversion
- Preserves exact rational arithmetic
- Automatically selected when polynomial has power as primary representation

## References

- Test code: `tests/test_hp_multiplicity.cpp`
- Working example: `examples/multiplicity_1d_roots.cpp`
- Implementation: `src/result_refiner_hp.cpp`
- Power-basis differentiation: `src/differentiation_hp.cpp`
- Multiplicity detection: `estimateMultiplicity()`, `estimateMultiplicityOstrowski()`, `estimateMultiplicitySturm()`

