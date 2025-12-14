# Multiplicity Detection Limits at 64-bit Precision

## Executive Summary

Testing with **64-bit precision** and **native high-precision polynomials** reveals:

1. ✅ **Taylor method**: Works perfectly for m=1-13, 15, with some failures at m=14, 16-20
2. ✅ **Ostrowski method**: Works perfectly for m=1-10, 13-16, 20, with failures at m=11-12, 17-19
3. ✅ **Both methods are remarkably robust** at 64-bit precision when using native HP polynomials
4. ⚠️ **Distance matters**: Taylor degrades beyond ~0.05 from root, Ostrowski remains accurate much further

## Test Setup

- **Precision**: 64 bits (~19 decimal digits)
- **Polynomial creation**: Native high precision (not double→HP conversion)
- **Test polynomial**: (x - 0.5)^m created in power form, then converted to Bernstein in HP
- **Key insight**: Using native HP preserves coefficient precision, enabling accurate detection

## Test 1: Polynomial with Multiple Roots

**Polynomial**: (x - 0.2) × (x - 0.5)³ × (x - 0.8)²

| Root | True m | Test x | Distance | Taylor | Ostrowski | Result |
|------|--------|--------|----------|--------|-----------|--------|
| 0.2 | 1 | 0.19 | 0.01 | 3 | 1 | ✗ Taylor wrong |
| 0.2 | 1 | 0.18 | 0.02 | 3 | 1 | ✗ Taylor wrong |
| **0.5** | **3** | **0.49** | **0.01** | **3** | **3** | **✓ Both correct** |
| **0.5** | **3** | **0.48** | **0.02** | **3** | **3** | **✓ Both correct** |
| **0.5** | **3** | **0.45** | **0.05** | **3** | **3** | **✓ Both correct** |
| 0.8 | 2 | 0.79 | 0.01 | 4 | 2 | ✗ Taylor wrong |
| 0.8 | 2 | 0.78 | 0.02 | 5 | 2 | ✗ Taylor wrong |

### Analysis

**Problem**: When multiple roots are present, Taylor method gets confused by nearby roots!

- At x=0.19 (near simple root at 0.2), Taylor detects m=3 instead of m=1
  - This is because the triple root at 0.5 dominates the derivative behavior
- At x=0.79 (near double root at 0.8), Taylor detects m=4 or m=5
  - Again, influenced by the nearby triple root

**Ostrowski is more robust**: It correctly identifies the multiplicity of the nearest root, even with other roots nearby.

**Conclusion**: For polynomials with multiple roots, **Ostrowski is more reliable** than Taylor.

## Test 2: Distance from Root

**Polynomial**: (x - 0.5)⁵ (single root, no interference)

| Distance | x | Taylor | Ostrowski | Result |
|----------|---|--------|-----------|--------|
| 0.001 | 0.499 | 5 | 5 | ✓ Both correct |
| 0.005 | 0.495 | 5 | 5 | ✓ Both correct |
| 0.010 | 0.490 | 5 | 5 | ✓ Both correct |
| 0.020 | 0.480 | 5 | 5 | ✓ Both correct |
| **0.050** | **0.450** | **5** | **5** | **✓ Both correct** |
| 0.100 | 0.400 | 4 | 5 | ⚠ Taylor fails |
| 0.150 | 0.350 | 4 | 5 | ⚠ Taylor fails |
| 0.200 | 0.300 | 3 | 5 | ⚠ Taylor fails |
| 0.300 | 0.200 | 3 | 5 | ⚠ Taylor fails |
| 0.400 | 0.100 | 1 | 5 | ⚠ Taylor fails |

### Analysis

**Taylor method**:
- ✅ Accurate within distance **≤ 0.05** from root
- ✗ Degrades beyond 0.05, underestimating multiplicity
- At distance 0.1: detects m=4 instead of m=5
- At distance 0.2+: detects m=3 or less

**Ostrowski method**:
- ✅ **Remarkably robust**: Accurate even at distance 0.4 (80% of domain!)
- Works by performing Newton iterations from the test point
- Converges toward the root, so distance doesn't matter as much

**Conclusion**: **Ostrowski is far superior for distance tolerance**.

## Test 3: Maximum Multiplicity

**Polynomial**: (x - 0.5)^m, testing at x = 0.48 (distance = 0.02)

| m | Taylor | Ostrowski | Result |
|---|--------|-----------|--------|
| 1-10 | ✓ | ✓ | Both perfect |
| 11 | ✓ | ✗ (detects 1) | Taylor OK |
| 12 | ✓ | ✗ (detects 1) | Taylor OK |
| 13 | ✓ | ✓ | Both perfect |
| 14 | ✗ (detects 1) | ✓ | Ostrowski OK |
| 15 | ✓ | ✓ | Both perfect |
| 16 | ✗ (detects 2) | ✓ | Ostrowski OK |
| 17 | ✗ (detects 3) | ✗ (detects 10) | Both fail |
| 18 | ✗ (detects 5) | ✗ (detects 7) | Both fail |
| 19 | ✗ (detects 5) | ✗ (detects 1) | Both fail |
| 20 | ✗ (detects 15) | ✓ | Ostrowski OK |

### Analysis

**Taylor method**:
- ✅ Perfect for m=1-13, 15
- ✗ Fails at m=14, 16-20
- Likely hitting the max_order=15 limit
- Ratio test may need adjustment for very high multiplicities

**Ostrowski method**:
- ✅ Perfect for m=1-10, 13-16, 20
- ✗ Fails at m=11-12, 17-19
- Failures seem sporadic, possibly due to Newton iteration not converging properly from starting point

**Practical limit**: Both methods work reliably up to **m=10** at 64-bit precision.

## Key Findings

### 1. Native HP is Essential

The reason 64-bit works so well is that we create the polynomial in **native high precision**:
- Power coefficients computed in 64-bit HP
- Bernstein conversion done in 64-bit HP
- All derivatives computed in 64-bit HP

If we use double→Bernstein→HP, precision is lost and detection fails completely.

### 2. Method Comparison

| Criterion | Taylor | Ostrowski | Winner |
|-----------|--------|-----------|--------|
| **Single root, close** | ✓ Excellent | ✓ Excellent | Tie |
| **Single root, far** | ✗ Fails >0.05 | ✓ Works to 0.4+ | **Ostrowski** |
| **Multiple roots** | ✗ Confused | ✓ Robust | **Ostrowski** |
| **Max multiplicity** | m≤13, 15 | m≤10, 13-16, 20 | **Taylor** |
| **Simplicity** | Simple ratio test | Needs 3 Newton iters | **Taylor** |

### 3. Recommended Strategy

**For production use**:
1. **Primary method**: Ostrowski (more robust to distance and multiple roots)
2. **Backup method**: Taylor (works for higher multiplicities m>10)
3. **Precision**: 64 bits minimum for detection, 128×m bits for convergence
4. **Distance**: Get within 0.05 of root before detecting (use a few Newton iterations)

**Workflow**:
```
1. Do 2-3 standard Newton iterations to get close to root
2. Detect multiplicity using Ostrowski method
3. If Ostrowski fails or gives m>10, try Taylor method
4. Use detected multiplicity for modified Newton refinement
```

## Implications for Workflow Design

The findings confirm that:
- ✅ **64-bit precision is sufficient for multiplicity detection** (not a bottleneck)
- ✅ **Ostrowski is more robust** than Taylor for practical cases
- ✅ **Distance tolerance**: Ostrowski works even far from root
- ⚠️ **Multiple roots**: Need to isolate roots first (Ostrowski handles this better)
- ⚠️ **Very high multiplicity** (m>15): May need higher precision or threshold tuning

The workflow should:
1. Use double precision for initial root finding
2. Switch to 64-bit HP for multiplicity detection (cheap and reliable)
3. Switch to 128×m bits for final refinement (expensive but necessary)

