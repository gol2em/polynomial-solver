# Polynomial Solver Documentation

## Overview

This directory contains comprehensive documentation for the polynomial solver's design and implementation.

---

## Quick Start

**New to the project?** Start here:

1. 🎯 [**REVISED_DESIGN.md**](REVISED_DESIGN.md) - **START HERE!** Complete API design with region tracking
2. 📍 [**REGION_TRACKING_DESIGN.md**](REGION_TRACKING_DESIGN.md) - How region tracking works
3. 📊 [**POWER_VS_BERNSTEIN_COMPARISON.md**](POWER_VS_BERNSTEIN_COMPARISON.md) - Performance comparison
4. 📋 [**QUICK_REFERENCE.md**](QUICK_REFERENCE.md) - Common operations and workflows

---

## Documentation Index

### Core Design Documents

| Document | Purpose | Audience |
|----------|---------|----------|
| [**REVISED_DESIGN.md**](REVISED_DESIGN.md) | **Complete API design** | **All developers** |
| [**REGION_TRACKING_DESIGN.md**](REGION_TRACKING_DESIGN.md) | **Region tracking explained** | **All developers** |
| [**POWER_VS_BERNSTEIN_COMPARISON.md**](POWER_VS_BERNSTEIN_COMPARISON.md) | **Performance analysis** | **All developers** |
| [**DESIGN_SUMMARY.md**](DESIGN_SUMMARY.md) | Architecture overview | All developers |
| [**POLYNOMIAL_DESIGN_SUMMARY.md**](POLYNOMIAL_DESIGN_SUMMARY.md) | Original design notes | Historical reference |

### Operation-Specific Guides

| Document | Purpose | Audience |
|----------|---------|----------|
| [**BERNSTEIN_OPERATIONS.md**](BERNSTEIN_OPERATIONS.md) | Operations in Bernstein basis | Implementers |
| [**BERNSTEIN_DIVISION_ALGORITHM.md**](BERNSTEIN_DIVISION_ALGORITHM.md) | Division algorithm for Bernstein | Implementers |
| [**DIFFERENTIATION_INSIGHT.md**](DIFFERENTIATION_INSIGHT.md) | Why differentiation works in both bases | All developers |

### Quick References

| Document | Purpose | Audience |
|----------|---------|----------|
| [**QUICK_REFERENCE.md**](QUICK_REFERENCE.md) | Common workflows and tips | Users |

---

## Key Concepts

### Two Bases as Solver Tools

The solver supports two polynomial representations:

1. **Power Basis** (on [a,b]^n)
   - Coefficients: Represent polynomial on original region (unchanged after subdivision)
   - Subdivision: O(1) - update region bounds only
   - Best for: Efficient subdivision via region tracking

2. **Bernstein Basis** (always on [0,1]^n)
   - Coefficients: Always represent polynomial on [0,1]^n
   - Subdivision: O(n²) - De Casteljau algorithm
   - Best for: Convex hull property, numerically stable subdivision

**Design principle**: Both bases are tools. Choose based on solver needs. Region tracking enables reconstruction.

### Key Insights

**1. Region Tracking Enables Complete Workflow**

User input on [a,b]^n → Convert to Bernstein → Subdivide → Reconstruct box in original coordinates

**2. Bernstein Coefficients Always on [0,1]^n**

- Simplifies implementation (no "coefficient region" to track)
- `original_box` tracks mapping to original coordinates
- Enables reconstruction after subdivision

**3. Power Basis Supports Arbitrary Regions**

- Coefficients represent polynomial on original region (unchanged)
- `current_region` tracks current region for evaluation
- O(1) subdivision via region tracking

**4. Both Bases Are Tools**

- Power: O(1) subdivision, arbitrary regions
- Bernstein: Convex hull property, numerically stable
- Choose based on solver needs

---

## Reading Guide

### For Users

**Goal**: Use the polynomial solver effectively

1. Start with [REVISED_DESIGN.md](REVISED_DESIGN.md) for corrected design
2. Read [POWER_VS_BERNSTEIN_COMPARISON.md](POWER_VS_BERNSTEIN_COMPARISON.md) for performance analysis
3. Use [QUICK_REFERENCE.md](QUICK_REFERENCE.md) for common operations

### For Implementers

**Goal**: Implement new operations or modify existing ones

1. Read [REVISED_DESIGN.md](REVISED_DESIGN.md) for API design
2. Study [POWER_VS_BERNSTEIN_COMPARISON.md](POWER_VS_BERNSTEIN_COMPARISON.md) for performance
3. Review [BERNSTEIN_DIVISION_ALGORITHM.md](BERNSTEIN_DIVISION_ALGORITHM.md) if implementing Bernstein division

### For Maintainers

**Goal**: Understand design decisions and evolution

1. Start with [POLYNOMIAL_DESIGN_SUMMARY.md](POLYNOMIAL_DESIGN_SUMMARY.md) (original design)
2. Read [DIFFERENTIATION_INSIGHT.md](DIFFERENTIATION_INSIGHT.md) (discovery #1: differentiation)
3. Read [BERNSTEIN_DIVISION_ALGORITHM.md](BERNSTEIN_DIVISION_ALGORITHM.md) (discovery #2: division)
4. Review [COMPLETE_PICTURE.md](COMPLETE_PICTURE.md) (complete understanding)
5. See [UNIFIED_API_DESIGN.md](UNIFIED_API_DESIGN.md) (current design)

---

## Design Evolution

### Phase 1: Original Design
- Single basis (power)
- Manual conversion for subdivision
- See: [POLYNOMIAL_DESIGN_SUMMARY.md](POLYNOMIAL_DESIGN_SUMMARY.md)

### Phase 2: Dual Basis Design
- Added Bernstein basis
- Automatic conversion for operations
- See: [UNIFIED_API_DESIGN.md](UNIFIED_API_DESIGN.md)

### Phase 3: Region Tracking Design (Current)
- **Key insight**: Bernstein coefficients always on [0,1]^n
- **Region tracking**: Enables reconstruction after subdivision
- **Both bases as tools**: No "winner", choose based on needs
- **Complete workflow**: User input → conversion → subdivision → reconstruction
- See: [REVISED_DESIGN.md](REVISED_DESIGN.md) and [REGION_TRACKING_DESIGN.md](REGION_TRACKING_DESIGN.md)

---

## Common Questions

### Q: Which basis should I use?

**A**: Choose based on your needs:

**Power basis** - for efficient subdivision:
```cpp
auto poly = PolynomialBase::fromPower(degrees, coeffs, {2,3}, {5,7});
auto sub = poly.restrictedToInterval(0, 0.3, 0.7);  // O(1)
auto val = sub.evaluate({0.5});  // Automatic rescaling
```

**Bernstein basis** - for convex hull property:
```cpp
auto bern = power.convertToBernstein();
auto [min, max] = bern.getConvexHullBounds();
```

### Q: How does reconstruction work?

**A**: Region tracking enables reconstruction:

```cpp
// User polynomial on [2,5]×[3,7]
auto power = PolynomialBase::fromPower(degrees, coeffs, {2,3}, {5,7});
auto bern = power.convertToBernstein();

// Subdivide
auto sub = bern.restrictedToInterval(0, 0.0, 0.5);

// Reconstruct box in original coordinates
auto [box_lower, box_upper] = sub.getOriginalBox();
// box_lower = [2, 5], box_upper = [3.5, 7]  ✓
```

See [REGION_TRACKING_DESIGN.md](REGION_TRACKING_DESIGN.md) for details.

### Q: How do I optimize performance?

**A**: Stay in power basis!

```cpp
// Excellent: Power basis subdivision is O(1)
auto poly = PolynomialBase::fromPower(degrees, coeffs);
for (int i = 0; i < 1000; ++i) {
    auto sub = poly.restrictedToInterval(0, a[i], b[i]);  // O(1) each!
    auto val = sub.evaluate({0.5});  // O(n) with rescaling
}

// Bad: Bernstein subdivision is O(n²)
auto bern = poly.convertToBernstein();
for (int i = 0; i < 1000; ++i) {
    auto sub = bern.restrictedToInterval(0, a[i], b[i]);  // O(n²) each!
}
```

See [POWER_VS_BERNSTEIN_COMPARISON.md](POWER_VS_BERNSTEIN_COMPARISON.md) for details.

---

## Contributing

When adding new operations:

1. Determine which basis is most natural
2. Implement in that basis
3. Add auto-conversion if needed
4. Update documentation
5. Add tests

See [UNIFIED_API_DESIGN.md](UNIFIED_API_DESIGN.md) for API design guidelines.

---

## Summary

🎯 **START HERE**: [REVISED_DESIGN.md](REVISED_DESIGN.md) - Complete API design
📍 **Region tracking**: [REGION_TRACKING_DESIGN.md](REGION_TRACKING_DESIGN.md) - How it works
📊 **Performance**: [POWER_VS_BERNSTEIN_COMPARISON.md](POWER_VS_BERNSTEIN_COMPARISON.md) - Detailed comparison
📋 **Quick reference**: [QUICK_REFERENCE.md](QUICK_REFERENCE.md)
💡 **Key insights**:
  - Bernstein coefficients always on [0,1]^n
  - Region tracking enables reconstruction
  - Both bases are tools (no "winner")
  - Complete workflow: input → conversion → subdivision → reconstruction
🔧 **Advanced**:
  - [BERNSTEIN_DIVISION_ALGORITHM.md](BERNSTEIN_DIVISION_ALGORITHM.md) - Division in Bernstein
  - [DIFFERENTIATION_INSIGHT.md](DIFFERENTIATION_INSIGHT.md) - Differentiation in both bases

