# Polynomial Design Summary: Complete Architecture

This document summarizes the complete design for the templated polynomial class with error tracking, dual-basis support, and unified API.

---

## Design Documents

1. **EXACT_COEFFICIENT_DESIGN.md** - Error tracking architecture
2. **SUBDIVISION_DESIGN.md** - Subdivision strategies for both bases
3. **UNIFIED_API_DESIGN.md** - Consistent API across representations

---

## Core Architecture

### Data Structure

```cpp
template<typename Scalar, typename ErrorScalar = double>
class PolynomialBase {
private:
    // Coefficient storage (single representation)
    std::vector<Scalar> coeffs_;
    Basis basis_;  // POWER or BERNSTEIN
    
    // Metadata
    std::vector<unsigned int> degrees_;
    std::size_t dimension_;
    
    // Region bounds (only for POWER basis on non-standard regions)
    std::vector<Scalar> region_lower_;  // Empty if Bernstein or [0,1]^n
    std::vector<Scalar> region_upper_;  // Empty if Bernstein or [0,1]^n
    
    // Error tracking
    ErrorBound<ErrorScalar> error_bound_;
};
```

### Key Principles

1. **Single representation** - Either power OR Bernstein, not both
2. **Error tracking** - Every operation updates error bound
3. **Region bounds for power** - Enables cheap subdivision
4. **Unified API** - All operations work on any polynomial
5. **Automatic normalization** - Algebraic ops normalize transparently

---

## Operation Summary

### Subdivision

| Basis | Method | Cost | Error Growth |
|-------|--------|------|--------------|
| Bernstein | De Casteljau (modify coeffs) | O(n) | O(ε·depth) |
| Power | Update region bounds | O(1) | O(ε·depth) |

**Both return polynomials that evaluate on [0,1]^n**

### Evaluation

| Basis | Method | Notes |
|-------|--------|-------|
| Bernstein | De Casteljau on u ∈ [0,1]^n | Direct |
| Power | Horner on x = a + (b-a)·u | Transform coordinates |

**Unified interface**: `evaluate(u)` where u ∈ [0,1]^n

### Algebraic Operations

| Operation | Bernstein | Power | Auto-Normalize |
|-----------|-----------|-------|----------------|
| Differentiation | ✅ Easy! d_i = n(b_{i+1}-b_i) | ✅ d/dx(x^i) = i·x^(i-1) | Power only |
| Division | ❌ Hard | ✅ Standard algorithm | Yes (convert to power) |
| Sturm sequence | ❌ Hard | ✅ Standard algorithm | Yes (convert to power) |
| Arithmetic (+,-,×) | ✅ Possible | ✅ Standard | Either basis |

**Differentiation works in BOTH bases! Division requires power basis.**

### Conversion

| From | To | Method | Cost |
|------|-----|--------|------|
| Power [0,1]^n | Bernstein | Power→Bern in HP | O(n²) |
| Power [a,b]^n | Bernstein | Normalize + convert | O(n³) |
| Bernstein | Power [0,1]^n | Bern→Power in HP | O(n²) |

**Always use high precision for conversions to minimize error**

---

## Error Tracking

### Error Sources

1. **Conversion** (Power ↔ Bernstein): Moderate - use HP
2. **Region normalization** (Power [a,b]^n → [0,1]^n): Large - use HP
3. **Subdivision** (Bernstein): Small - convex combinations
4. **Coordinate transform** (Power): Tiny - 2 ops per dimension
5. **Algebraic operations**: Varies - use HP for division

### Error Accumulation

```cpp
ErrorBound<ErrorScalar> {
    ErrorScalar absolute_error;  // Conservative upper bound
    
    bool isExact() const { return absolute_error == ErrorScalar(0); }
};
```

**Error grows additively through operations**

### Precision Lifting

When error exceeds tolerance:
```cpp
// Lift to higher precision
auto hp_poly = poly.liftPrecision<mpreal>();

// Continue with higher precision
auto hp_result = hp_poly.restrictedToInterval(0, a, b);
```

---

## Usage Patterns

### Pattern 1: Subdivision Solver (Bernstein)

```cpp
// Start with power polynomial
auto poly = PolynomialBase<double>::fromPower(degrees, coeffs);

// Convert to Bernstein once
auto bern = poly.convertToBernstein();

// Subdivide many times (cheap!)
while (!converged) {
    auto sub = bern.restrictedToInterval(axis, a, b);
    
    // Check error
    if (sub.errorBound() > tolerance) {
        // Lift precision
        auto hp_poly = poly.liftPrecision<mpreal>();
        auto hp_bern = hp_poly.convertToBernstein();
        sub = hp_bern.restrictedToInterval(axis, a, b);
    }
    
    bern = sub;
}
```

### Pattern 2: Root Isolation (Sturm)

```cpp
// Polynomial in any basis
auto poly = PolynomialBase<mpreal>::fromPower(degrees, coeffs);

// Sturm sequence (auto-normalizes)
auto sturm = poly.sturmSequence();

// Count roots in interval
int count = countRootsInInterval(sturm, a, b);
```

### Pattern 3: Mixed Operations

```cpp
// Subdivided polynomial (Bernstein)
auto bern = /* ... */;

// Differentiate (converts to power)
auto dbern = bern.differentiate(0);

// Subdivide derivative (uses region bounds)
auto sub_deriv = dbern.restrictedToInterval(0, 0.3, 0.7);

// Evaluate
double val = sub_deriv.evaluate({0.5});
```

---

## Performance Optimization

### Cache Normalized Form

```cpp
// BAD: Normalizes 3 times
auto dp = poly.differentiate(0);
auto ddp = poly.differentiate(0, 2);
auto sturm = poly.sturmSequence();

// GOOD: Normalize once
auto poly_norm = poly.convertToPower();  // Ensures [0,1]^n
auto dp = poly_norm.differentiate(0);
auto ddp = poly_norm.differentiate(0, 2);
auto sturm = poly_norm.sturmSequence();
```

### Choose Right Basis

- **Subdivision-heavy**: Use Bernstein
- **Algebra-heavy**: Use Power on [0,1]^n
- **Mixed**: Start Power, convert to Bernstein for subdivision

---

## Implementation Priorities

### Phase 1: Core Infrastructure
- [ ] ErrorBound struct and type traits
- [ ] PolynomialBase with single representation
- [ ] Region bounds for power basis
- [ ] Basic evaluation (both bases)

### Phase 2: Subdivision
- [ ] De Casteljau subdivision (Bernstein)
- [ ] Region-based subdivision (Power)
- [ ] Error estimation for subdivision

### Phase 3: Conversion
- [ ] Power → Bernstein in high precision
- [ ] Bernstein → Power in high precision
- [ ] Region normalization in high precision
- [ ] Error estimation for conversion

### Phase 4: Algebraic Operations
- [ ] Differentiation (power basis)
- [ ] Division with remainder (power basis)
- [ ] Polynomial arithmetic (+, -, ×)
- [ ] Auto-normalization wrapper

### Phase 5: Advanced Features
- [ ] Sturm sequence
- [ ] Precision lifting
- [ ] Error-based precision selection
- [ ] Performance optimizations

---

## Key Advantages

✅ **Numerical stability**: Power subdivision via coordinates (not coefficients)
✅ **Error awareness**: Always know accumulated error
✅ **Exact arithmetic support**: Works with GMP rationals
✅ **Unified API**: Same interface for all operations
✅ **Automatic precision**: Lift when error grows
✅ **Performance**: O(1) subdivision for power basis
✅ **Simplicity**: Single representation, no dual storage

---

## Migration from Current Design

1. Remove dual representation (power_coeffs_ + bernstein_coeffs_)
2. Remove validity flags (power_valid_, bernstein_valid_)
3. Add region bounds (region_lower_, region_upper_)
4. Add error tracking (error_bound_)
5. Update subdivision to use region bounds for power
6. Add auto-normalization for algebraic operations
7. Update conversion to use high precision

