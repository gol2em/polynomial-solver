# The Complete Picture: Bernstein vs Power Basis

## Summary of Discoveries

### Initial Understanding
- **Bernstein**: Good for subdivision (De Casteljau)
- **Power**: Good for division and algebraic operations
- **Assumption**: Need to convert between bases for different operations

### Key Discovery #1: Differentiation Works in Both Bases!

**Bernstein formula**:
```
d_i = n · (b_{i+1} - b_i)
```

**Power formula**:
```
d_i = (i+1) · a_{i+1}
```

**Impact**: No conversion needed for differentiation!

### Key Discovery #2: Division Works in Bernstein Too!

**Algorithm** (Busé et al. 2008):
- Uses homogenization: `f^h(t,s) = s^deg(f) · f(t/s)`
- Division formula: `g(t) = q(t)·f(t) + (1-t)^(e-d+1)·r(t)`
- Key operations: Multiply/divide by powers of `t` and `(1-t)`

**Impact**: Division is possible in Bernstein, though more complex than power!

---

## Complete Operation Matrix

| Operation | Bernstein | Power | Recommendation |
|-----------|-----------|-------|----------------|
| **Differentiation** | ✅ Simple (`d_i = n(b_{i+1}-b_i)`) | ✅ Simple (`d_i = (i+1)a_{i+1}`) | **Use current basis** |
| **Subdivision** | ✅ De Casteljau (O(n)) | ⚠️ Region bounds (O(1)) | **Use Bernstein** |
| **Division** | ✅ Homogenization (complex) | ✅ Standard (simple) | **Use Power** |
| **Evaluation** | ✅ De Casteljau (stable) | ✅ Horner (fast) | **Use current basis** |
| **Multiplication** | ✅ Possible | ✅ Standard | **Either** |
| **GCD** | ✅ Possible (uses division) | ✅ Standard | **Use Power** |
| **Sturm sequence** | ✅ Possible (uses division) | ✅ Standard | **Use Power** |

**Key insight**: Almost ALL operations work in BOTH bases!

---

## Design Implications

### Original Design (Assumed)

```
Bernstein: Subdivision only
Power: Everything else
→ Frequent conversions needed
```

### Actual Reality

```
Bernstein: Subdivision, differentiation, evaluation, division (complex)
Power: Division (simple), differentiation, evaluation, arithmetic
→ Minimal conversions needed!
```

### Optimal Strategy

**For subdivision-heavy solvers** (typical use case):

1. **Start in power** (user input)
2. **Convert to Bernstein once**
3. **Stay in Bernstein** for:
   - Subdivision (De Casteljau)
   - Differentiation (simple formula)
   - Evaluation (stable)
4. **Convert to power only if needed** for:
   - Division (simpler algorithm)
   - Sturm sequences (uses division)

**Performance**: Only 1-2 conversions total, instead of many!

---

## Implementation Phases

### Phase 1: Core Operations (Current)

**Implement in BOTH bases**:
- ✅ Differentiation (Bernstein and Power)
- ✅ Evaluation (Bernstein and Power)
- ✅ Subdivision (Bernstein preferred)

**Implement in Power only**:
- ✅ Division (simpler algorithm)
- ✅ Sturm sequences (uses division)

**Rationale**: Get solver working quickly with minimal complexity.

### Phase 2: Optimization (Future)

**Implement Bernstein division** if:
- Profiling shows conversion overhead is significant
- Sturm sequences are computed frequently
- Need to maintain Bernstein throughout

**Algorithm**: Use homogenization method from Busé et al. (2008)

**Benefit**: Eliminate conversions for division-heavy workflows

---

## Numerical Stability Comparison

### Differentiation

| Basis | Formula | Stability |
|-------|---------|-----------|
| Bernstein | `d_i = n(b_{i+1} - b_i)` | ✅ Excellent (only differences) |
| Power | `d_i = (i+1)a_{i+1}` | ✅ Excellent (only scaling) |

**Winner**: Tie (both excellent)

### Subdivision

| Basis | Method | Stability |
|-------|--------|-----------|
| Bernstein | De Casteljau | ✅ Excellent (convex combinations) |
| Power | Region bounds | ⚠️ Good (no coefficient changes) |

**Winner**: Bernstein (numerically superior)

### Division

| Basis | Method | Stability |
|-------|--------|-----------|
| Bernstein | Homogenization | ⚠️ Good (binomial coefficients) |
| Power | Standard algorithm | ✅ Excellent (well-studied) |

**Winner**: Power (simpler and well-understood)

### Evaluation

| Basis | Method | Stability |
|-------|--------|-----------|
| Bernstein | De Casteljau | ✅ Excellent (convex combinations) |
| Power | Horner | ✅ Excellent (minimal operations) |

**Winner**: Tie (both excellent)

---

## Complexity Comparison

### Differentiation

| Basis | Time | Space |
|-------|------|-------|
| Bernstein | O(n) | O(n) |
| Power | O(n) | O(n) |

**Winner**: Tie

### Subdivision

| Basis | Time | Space |
|-------|------|-------|
| Bernstein | O(n²) | O(n) |
| Power | O(1) | O(1) |

**Winner**: Power (but Bernstein is more stable)

### Division

| Basis | Time | Space |
|-------|------|-------|
| Bernstein | O(d·e) | O(e) |
| Power | O(d·e) | O(e) |

**Winner**: Tie (but Power is simpler to implement)

---

## Final Recommendations

### For Polynomial Solver

**Use Bernstein for**:
1. ✅ Subdivision (De Casteljau - numerically superior)
2. ✅ Differentiation (simple formula - no conversion needed)
3. ✅ Evaluation (stable)

**Use Power for**:
1. ✅ Division (simpler algorithm)
2. ✅ Sturm sequences (uses division)
3. ✅ Initial input (user provides power coefficients)

**Conversion strategy**:
- Convert to Bernstein at start of solver
- Stay in Bernstein for subdivision + differentiation
- Convert to Power only when division is needed
- **Total conversions**: 1-2 per solve (minimal!)

### Implementation Priority

**Phase 1** (now):
1. Differentiation in both bases ✅
2. Subdivision in Bernstein ✅
3. Division in Power ✅
4. Evaluation in both bases ✅

**Phase 2** (later, if needed):
1. Bernstein division (for optimization)
2. Bernstein multiplication (for completeness)
3. Bernstein GCD (for completeness)

---

## Summary

✅ **Differentiation**: Works in BOTH bases (simple formulas)
✅ **Division**: Works in BOTH bases (Power is simpler)
✅ **Subdivision**: Best in Bernstein (De Casteljau)
✅ **Evaluation**: Works in BOTH bases (both stable)

**Key insight**: We can do almost EVERYTHING in Bernstein! Only division is significantly simpler in power basis.

**Performance win**: Minimal conversions needed (1-2 per solve instead of many).

**Numerical win**: Stay in Bernstein for subdivision-heavy workflows (more stable).

**Implementation win**: Start simple (power division), optimize later (Bernstein division) if needed.

