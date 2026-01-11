# Polynomial Operations Quick Reference

## When to Use Each Basis

### Use Bernstein Basis For:

✅ **Subdivision** (De Casteljau algorithm)
```cpp
auto bern = poly.convertToBernstein();
auto sub = bern.restrictedToInterval(0, 0.3, 0.7);
```

✅ **Differentiation** (simple formula: `d_i = n(b_{i+1} - b_i)`)
```cpp
auto deriv = bern.differentiate(0);  // No conversion needed!
```

✅ **Evaluation** (numerically stable)
```cpp
double val = bern.evaluate({0.5, 0.3});
```

✅ **Convex hull properties**
```cpp
auto bounds = bern.getBounds();  // Coefficients bound values
```

### Use Power Basis For:

✅ **Division with remainder** (simpler than Bernstein)
```cpp
auto power = poly.ensurePowerStandardRegion();
auto [quot, rem] = power.divideWithRemainder(divisor);
```

**Note**: Division also works in Bernstein basis (see [BERNSTEIN_DIVISION_ALGORITHM.md](../BERNSTEIN_DIVISION_ALGORITHM.md)), but power basis is simpler.

✅ **Sturm sequences** (uses division)
```cpp
auto sturm = power.sturmSequence();
```

✅ **Polynomial arithmetic** (when needed)
```cpp
auto sum = p1.add(p2);  // Works in power basis
```

---

## Operation Cheat Sheet

| I want to... | Best Basis | Auto-Convert? | Example |
|--------------|-----------|---------------|---------|
| Subdivide | Bernstein | ✅ Yes | `poly.restrictedToInterval(0, a, b)` |
| Differentiate | **Either!** | ❌ No | `poly.differentiate(0)` |
| Divide | Power | ✅ Yes | `poly.divideWithRemainder(divisor)` |
| Evaluate | Either | ❌ No | `poly.evaluate({x, y})` |
| Get bounds | Bernstein | ✅ Yes | `poly.getBounds()` |
| Sturm sequence | Power | ✅ Yes | `poly.sturmSequence()` |

---

## Common Workflows

### Workflow 1: Subdivision-Based Solver

```cpp
// 1. Start with power polynomial
auto poly = PolynomialBase<mpreal>::fromPower(degrees, coeffs);

// 2. Convert to Bernstein (once)
auto bern = poly.convertToBernstein();

// 3. Subdivide repeatedly (stays in Bernstein)
auto left = bern.restrictedToInterval(0, 0.0, 0.5);
auto right = bern.restrictedToInterval(0, 0.5, 1.0);

// 4. Differentiate (stays in Bernstein - no conversion!)
auto dleft = left.differentiate(0);
auto dright = right.differentiate(0);

// 5. Evaluate (Bernstein)
double val_left = dleft.evaluate({0.25});
double val_right = dright.evaluate({0.75});
```

**Performance**: Only 1 conversion (Power → Bernstein at start)

### Workflow 2: Root Counting with Sturm

```cpp
// 1. Start with any polynomial
auto poly = /* ... */;

// 2. Build Sturm sequence (auto-converts to power if needed)
auto sturm = poly.sturmSequence();

// 3. Count roots in interval
int num_roots = countRootsInInterval(sturm, a, b);
```

**Performance**: Converts to power once, then stays in power

### Workflow 3: Mixed Operations

```cpp
// 1. Subdivide in Bernstein
auto bern = poly.convertToBernstein();
auto sub = bern.restrictedToInterval(0, 0.3, 0.7);

// 2. Differentiate (stays in Bernstein!)
auto deriv = sub.differentiate(0);

// 3. Need division? Convert to power
auto power = deriv.ensurePowerStandardRegion();
auto [quot, rem] = power.divideWithRemainder(divisor);

// 4. Back to Bernstein for more subdivision
auto bern2 = quot.convertToBernstein();
```

---

## Performance Tips

### ✅ DO: Pre-convert for repeated operations

```cpp
// Good: Convert once, subdivide many times
auto bern = poly.convertToBernstein();
for (int i = 0; i < 100; ++i) {
    auto sub = bern.restrictedToInterval(0, a[i], b[i]);
    // ...
}
```

### ❌ DON'T: Convert repeatedly

```cpp
// Bad: Converts 100 times!
for (int i = 0; i < 100; ++i) {
    auto sub = poly.restrictedToInterval(0, a[i], b[i]);  // Converts each time!
    // ...
}
```

### ✅ DO: Use appropriate basis for operation

```cpp
// Good: Differentiate in current basis
auto deriv = poly.differentiate(0);  // Works in both bases!
```

### ❌ DON'T: Force unnecessary conversions

```cpp
// Bad: Unnecessary conversion
auto power = poly.convertToPower();  // Not needed for differentiation!
auto deriv = power.differentiate(0);
```

---

## Basis Conversion

### When Conversion Happens Automatically

1. **Subdivision** → Converts to Bernstein if needed
2. **Division** → Converts to Power if needed
3. **Sturm sequence** → Converts to Power if needed
4. **Power on [a,b]^n + algebraic op** → Normalizes to [0,1]^n

### When Conversion Does NOT Happen

1. **Differentiation** → Works in current basis
2. **Evaluation** → Works in current basis
3. **Arithmetic** → Works in current basis (usually power)

### Manual Conversion

```cpp
// Convert to Bernstein
auto bern = poly.convertToBernstein();

// Convert to Power on [0,1]^n
auto power = poly.ensurePowerStandardRegion();

// Check current basis
if (poly.getBasis() == Basis::BERNSTEIN) {
    // ...
}
```

---

## Error Tracking

Every operation updates error bounds:

```cpp
auto poly = PolynomialBase<double>::fromPower(degrees, coeffs);
std::cout << "Initial error: " << poly.getAbsoluteError() << std::endl;

auto deriv = poly.differentiate(0);
std::cout << "After differentiation: " << deriv.getAbsoluteError() << std::endl;

auto sub = deriv.restrictedToInterval(0, 0.3, 0.7);
std::cout << "After subdivision: " << sub.getAbsoluteError() << std::endl;
```

---

## Summary

**Key Takeaways**:
1. ✅ Differentiation works in BOTH bases (no conversion needed!)
2. ✅ Subdivision best in Bernstein (auto-converts)
3. ✅ Division requires Power (auto-converts)
4. ✅ Pre-convert for repeated operations
5. ✅ Error bounds tracked automatically

