# Tier 1 API Baseline Documentation (Phase 3)

## Overview

This document provides a comprehensive baseline of the current double-precision API (Tier 1) before implementing high-precision support (Tiers 2 and 3). All APIs documented here use `double` precision and will need high-precision equivalents in future phases.

**Date**: 2025-11-28  
**Phase**: Phase 3 - Verify Tier 1 (Minimum Version)  
**Purpose**: Document current API for reference during high-precision refactorization

---

## Core Classes

### 1. Polynomial Class

**File**: `include/polynomial.h`, `src/polynomial.cpp`

**Purpose**: Represents multivariate polynomials in Bernstein basis

**Key Methods**:

| Method | Signature | Description | HP Needed? |
|--------|-----------|-------------|------------|
| Constructor | `Polynomial(degrees, bernstein_coeffs)` | Construct from Bernstein coefficients | ✅ Yes |
| fromBernstein | `static Polynomial fromBernstein(degrees, coeffs)` | Factory from Bernstein | ✅ Yes |
| fromPower | `static Polynomial fromPower(degrees, coeffs)` | Factory from power basis | ✅ Yes |
| evaluate | `double evaluate(parameters) const` | Evaluate at point | ✅ Yes |
| evaluate | `double evaluate(t) const` | Evaluate univariate at t | ✅ Yes |
| graphControlPoints | `void graphControlPoints(vector<double>&) const` | Get graph control points | ✅ Yes |
| graphCorners | `void graphCorners(vector<vector<double>>&) const` | Get graph corners | ✅ Yes |
| restrictedToInterval | `Polynomial restrictedToInterval(axis, a, b) const` | Restrict to interval | ✅ Yes |
| dimension | `size_t dimension() const` | Get dimension | ⚪ No |
| degrees | `const vector<unsigned int>& degrees() const` | Get degrees | ⚪ No |
| coefficientCount | `size_t coefficientCount() const` | Get coefficient count | ⚪ No |
| bernsteinCoefficients | `const vector<double>& bernsteinCoefficients() const` | Get coefficients | ✅ Yes |

**Data Members**:
- `dimension_`: `size_t` - Number of variables
- `degrees_`: `vector<unsigned int>` - Degrees per variable
- `bernstein_coeffs_`: `vector<double>` - Bernstein coefficients (**needs HP version**)

---

### 2. DeCasteljau Class

**File**: `include/de_casteljau.h`, `src/de_casteljau.cpp`

**Purpose**: De Casteljau algorithm for polynomial evaluation and subdivision

**Key Methods**:

| Method | Signature | Description | HP Needed? |
|--------|-----------|-------------|------------|
| evaluate1D | `static double evaluate1D(coeffs, t)` | Evaluate univariate polynomial | ✅ Yes |
| subdivide1D | `static void subdivide1D(coeffs, t, left, right)` | Subdivide at t | ✅ Yes |
| evaluateTensorProduct | `static double evaluateTensorProduct(degrees, coeffs, params)` | Evaluate multivariate | ✅ Yes |

**Notes**:
- All methods are static utilities
- Core algorithm for polynomial evaluation
- Critical for high-precision: needs template or HP version

---

### 3. Differentiation Class

**File**: `include/differentiation.h`, `src/differentiation.cpp`

**Purpose**: Compute derivatives of Bernstein polynomials

**Key Methods**:

| Method | Signature | Description | HP Needed? |
|--------|-----------|-------------|------------|
| derivative | `static Polynomial derivative(p, axis, order)` | Compute derivative | ✅ Yes |
| gradient | `static vector<Polynomial> gradient(p)` | Compute gradient | ✅ Yes |
| hessian | `static vector<vector<Polynomial>> hessian(p)` | Compute Hessian | ✅ Yes |
| differentiateAxis | `static Polynomial differentiateAxis(p, axis)` | Core differentiation | ✅ Yes |

**Notes**:
- All methods are static utilities
- Returns new Polynomial objects
- Needs HP version for high-precision derivatives

---

### 4. DerivativeCache Class

**File**: `include/differentiation.h`, `src/differentiation.cpp`

**Purpose**: Cache computed derivatives for efficient reuse

**Key Methods**:

| Method | Signature | Description | HP Needed? |
|--------|-----------|-------------|------------|
| Constructor | `DerivativeCache(p)` | Construct cache for polynomial | ✅ Yes |
| get | `const Polynomial& get(orders)` | Get derivative by multi-index | ✅ Yes |
| getPartial | `const Polynomial& getPartial(axis, order)` | Get partial derivative | ✅ Yes |
| precomputeUpToOrder | `void precomputeUpToOrder(maxOrder)` | Precompute derivatives | ✅ Yes |

**Data Members**:
- `dimension_`: `size_t` - Dimension
- `cache_`: `map<vector<unsigned int>, Polynomial>` - Derivative cache (**needs HP version**)

---

### 5. ResultRefiner Class

**File**: `include/result_refiner.h`, `src/result_refiner.cpp`

**Purpose**: Refine solver results using Newton's method with high precision

**Key Methods**:

| Method | Signature | Description | HP Needed? |
|--------|-----------|-------------|------------|
| refineResults | `static RefinementResult refineResults(system, boxes, config)` | Refine all results | ✅ Yes |
| refineRoot1D | `static RefinedRoot refineRoot1D(p, initial, config)` | Refine 1D root | ✅ Yes |
| estimateMultiplicity1D | `static unsigned int estimateMultiplicity1D(p, x, config)` | Estimate multiplicity | ✅ Yes |

**Structs**:
- `RefinementConfig`: Configuration parameters
- `RefinedRoot`: Refined root with multiplicity info
- `ProblematicRegion`: Unresolved region info
- `RefinementResult`: Complete refinement results

**Notes**:
- Uses Newton's method with condition-aware convergence
- Critical for achieving 1e-15 precision
- **High priority for HP implementation**

---

## Supporting Classes

### 6. Solver Class

**File**: `include/solver.h`, `src/solver.cpp`

**Purpose**: Main solver using subdivision methods

**Key Methods**:

| Method | Signature | Description | HP Needed? |
|--------|-----------|-------------|------------|
| subdivisionSolve | `SubdivisionResult subdivisionSolve(system, config, method)` | Solve system | ⚪ No (uses HP polynomials) |

**Notes**:
- Solver itself doesn't need HP version
- Will work with HP polynomials once they're implemented

---

## Summary: Methods Needing High-Precision Versions

### High Priority (Core Algorithms)
1. ✅ **Polynomial**: All evaluation and manipulation methods
2. ✅ **DeCasteljau**: All evaluation methods
3. ✅ **Differentiation**: All derivative computation methods
4. ✅ **DerivativeCache**: All caching methods
5. ✅ **ResultRefiner**: All refinement methods

### Medium Priority (Utilities)
6. ✅ **Conversion utilities**: double ↔ high-precision

### Low Priority (Metadata)
7. ⚪ **Dimension/degree queries**: No HP needed (return size_t/unsigned int)

---

## Next Steps (Phase 4)

1. Create `PolynomialHP` class (fixed high-precision, no templates)
2. Create `DeCasteljauHP` class
3. Create `DifferentiationHP` class
4. Create `DerivativeCache HP` class
5. Create `ResultRefinerHP` class
6. Add conversion utilities between double and HP types

---

## Testing Baseline (Phase 3 Continued)

See `docs/TIER1_TEST_BASELINE.md` for:
- List of all existing tests
- Baseline performance measurements
- Test coverage analysis

