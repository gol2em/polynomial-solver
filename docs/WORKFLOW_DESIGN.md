# High-Precision Refinement Workflow Design

## Overview

This document describes the recommended workflow for refining polynomial roots with automatic precision escalation, based on systematic testing of multiplicity detection and convergence requirements.

## Key Findings from Testing

### 1. Multiplicity Detection
- **Taylor method**: ✅ Works at all precisions (64-1024 bits), all multiplicities (m=2-10)
- **Ostrowski method**: ✅ Works at all precisions (64-1024 bits), all multiplicities (m=2-10)
- **Sturm sequence**: ❌ Not suitable for floating-point (requires exact arithmetic)

### 2. Convergence Requirements
- **Minimum precision**: Varies by multiplicity (64-192 bits)
- **Recommended precision**: 128 × m bits for smooth convergence
- **Convergence speed**: 1-3 iterations with adequate precision
- **Oscillation risk**: Occurs when precision is insufficient (error alternates between small/large)

### 3. Precision vs Multiplicity

| Multiplicity | Min Precision | Recommended | Iterations |
|--------------|---------------|-------------|------------|
| 2 | 64 bits | 256 bits | 1-2 |
| 3 | 96 bits | 384 bits | 1 |
| 4 | 64 bits | 512 bits | 1-3 |
| 5 | 96 bits | 640 bits | 1 |
| 6 | 64 bits | 768 bits | 1 |
| 7 | 192 bits | 896 bits | 1 |
| 8 | 64 bits | 1024 bits | 1-3 |
| 9 | 96 bits | 1152 bits | 1 |
| 10 | 64 bits | 1280 bits | 1 |

## Recommended Workflow

### Phase 1: Double-Precision Refinement

```
Input: Initial guess x₀, polynomial coefficients (double precision)

1. Attempt standard Newton refinement in double precision
2. Estimate condition number κ = |f''| / |f'|²
3. If converged with κ < 1e4:
   → Return refined root (simple root, well-conditioned)
4. If not converged or κ ≥ 1e4:
   → Proceed to Phase 2 (precision escalation)
```

### Phase 2: Precision Selection

```
Input: Condition number κ from Phase 1

1. Select initial precision based on condition number:
   - If κ < 1e5:  precision = 256 bits
   - If κ < 1e10: precision = 512 bits
   - If κ < 1e15: precision = 1024 bits
   - Else:        precision = 2048 bits

2. This provides a conservative starting point
```

### Phase 3: Multiplicity Detection

```
Input: Initial guess x₀, selected precision

1. Set precision context to selected precision
2. Convert polynomial coefficients to high precision:
   - If power coefficients available: use fromPowerHP()
   - If only Bernstein available: convert carefully
3. Detect multiplicity using Taylor method:
   - Use ratio threshold = 10.0 (works for m ≤ 11)
   - For extreme cases (m > 11), increase threshold to 50
4. Store detected multiplicity m
```

### Phase 4: Precision Adjustment

```
Input: Detected multiplicity m, current precision

1. Calculate required precision: P_req = 128 × m
2. If current precision < P_req:
   - Increase precision to P_req
   - Re-create polynomial in new precision
3. Set convergence tolerance based on precision:
   - For 256 bits:  tolerance = 1e-25
   - For 512 bits:  tolerance = 1e-50
   - For 1024 bits: tolerance = 1e-100
   - For 2048 bits: tolerance = 1e-200
```

### Phase 5: Modified Newton Refinement

```
Input: Initial guess x₀, multiplicity m, polynomial in HP

1. Apply modified Newton iteration:
   x_{n+1} = x_n - m × f(x_n) / f'(x_n)

2. Monitor convergence:
   - Track error: e_n = |x_n - x_{n-1}|
   - Check for oscillation: e_n > 2 × e_{n-1}
   - If oscillation detected after 5 iterations:
     → Double precision and restart from Phase 3

3. Convergence criteria:
   - |f(x)| < tolerance AND
   - |x_n - x_{n-1}| < tolerance
   - Maximum 50 iterations

4. If converged:
   - Compute error bounds using interval arithmetic
   - Return refined root with bounds
```

## Implementation Pseudocode

```cpp
struct RefinedRoot {
    double location;
    double error_bound;
    unsigned int multiplicity;
    bool converged;
    unsigned int precision_used;
};

RefinedRoot refineRootWithPrecisionEscalation(
    double x0,
    const std::vector<double>& power_coeffs,
    double tolerance = 1e-10)
{
    // Phase 1: Try double precision
    auto result_double = refineInDouble(x0, power_coeffs);
    if (result_double.converged && result_double.condition < 1e4) {
        return result_double;
    }

    // Phase 2: Select precision
    unsigned int precision = selectPrecision(result_double.condition);

    // Phase 3: Detect multiplicity
    PrecisionContext ctx(precision);
    PolynomialHP poly = createPolynomialHP(power_coeffs, precision);
    unsigned int m = detectMultiplicity(x0, poly);

    // Phase 4: Adjust precision
    unsigned int required_precision = 128 * m;
    if (precision < required_precision) {
        precision = required_precision;
        ctx.setPrecision(precision);
        poly = createPolynomialHP(power_coeffs, precision);
    }

    // Phase 5: Refine with modified Newton
    mpreal x = mpreal(x0);
    mpreal tol = selectTolerance(precision);
    
    for (unsigned int iter = 0; iter < 50; ++iter) {
        mpreal f = poly.evaluate(x);
        mpreal df = poly.derivative().evaluate(x);
        
        if (abs(f) < tol && abs(step) < tol) {
            return {toDouble(x), computeErrorBound(x, poly), m, true, precision};
        }
        
        mpreal step = mpreal(m) * f / df;
        x = x - step;
    }

    return {toDouble(x), 0.0, m, false, precision};
}
```

## Performance Characteristics

### Time Complexity
- **Multiplicity detection**: O(m × n) where n = polynomial degree
- **Modified Newton**: O(m_iter × n) where m_iter ≈ 1-3 iterations
- **Total**: O(m × n) dominated by multiplicity detection

### Space Complexity
- **Polynomial storage**: O(n) coefficients in high precision
- **Derivatives**: O(m × n) for m derivatives
- **Total**: O(m × n)

### Precision Overhead
- **64 bits → 256 bits**: ~4× slower
- **256 bits → 512 bits**: ~2× slower
- **512 bits → 1024 bits**: ~2× slower

## Edge Cases and Limitations

### 1. Very High Multiplicity (m > 10)
- **Challenge**: Requires very high precision (> 1280 bits)
- **Solution**: Increase Taylor ratio threshold to 50-100
- **Limitation**: May become impractical for m > 15

### 2. Clustered Roots
- **Challenge**: Multiple distinct roots very close together
- **Solution**: Use interval arithmetic to separate roots
- **Limitation**: May require symbolic preprocessing

### 3. Ill-Conditioned Polynomials
- **Challenge**: Condition number > 1e20
- **Solution**: Use even higher precision (4096+ bits)
- **Limitation**: Computational cost increases significantly

## Testing and Validation

All recommendations are based on systematic testing:
- **Test suite**: `test_precision_requirements.cpp`
- **Multiplicities tested**: m = 2-10
- **Precisions tested**: 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048 bits
- **Test polynomial**: (x - 0.5)^m in native high precision
- **Starting point**: x = 0.3 (distance 0.2 from root)

See `docs/PRECISION_REQUIREMENTS.md` for detailed test results.

