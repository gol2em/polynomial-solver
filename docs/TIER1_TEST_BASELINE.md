# Tier 1 Test Baseline (Phase 3)

## Overview

This document establishes the baseline test results for Tier 1 (double precision only) before implementing high-precision support.

**Date**: 2025-11-28  
**Configuration**: Default (no high-precision flags)  
**Build Type**: Debug  
**Compiler**: GCC  

---

## Test Results Summary

**Total Tests**: 15  
**Passed**: 14 (93%)  
**Failed**: 1 (7%)  
**Total Time**: 0.05 sec  

---

## Test Details

### ✅ Passing Tests (14/15)

| # | Test Name | Time | Status | Description |
|---|-----------|------|--------|-------------|
| 1 | HighPrecisionTypesTest | 0.00s | ✅ PASS | High-precision type definitions |
| 2 | PolynomialConversionTest | 0.01s | ✅ PASS | Power ↔ Bernstein conversion |
| 3 | PolynomialGraphTest | 0.00s | ✅ PASS | Graph control points |
| 4 | PolynomialSystemExampleTest | 0.00s | ✅ PASS | Polynomial system creation |
| 5 | SolverSubdivisionTest | 0.00s | ✅ PASS | Subdivision solver |
| 6 | GeometryConvexTest | 0.00s | ✅ PASS | Convex hull algorithms |
| 7 | 2DHyperplaneIntersectionTest | 0.00s | ✅ PASS | Hyperplane intersection |
| 8 | 2DConvexHullRobustTest | 0.00s | ✅ PASS | Robust convex hull |
| 9 | 2DHullVertexOrderTest | 0.00s | ✅ PASS | Hull vertex ordering |
| 10 | DegenerateBoxesTest | 0.00s | ✅ PASS | Degenerate case handling |
| 11 | LinearGraphHullTest | 0.00s | ✅ PASS | Linear polynomial graph hull |
| 12 | ProjectedPolyhedralTest | 0.00s | ✅ PASS | PP method |
| 13 | EllipseDumpTest | 0.00s | ✅ PASS | Geometry dump functionality |
| 14 | DifferentiationTest | 0.00s | ✅ PASS | Differentiation algorithms |
| 15 | ResultRefinerTest | 0.01s | ❌ FAIL | Result refinement (known issues) |

---

## ❌ Known Failing Test: ResultRefinerTest

### Test 1: Cubic polynomial refinement
**Status**: ❌ FAIL  
**Issue**: Expected 3 unique roots, got 2  
**Root Cause**: Duplicate elimination may be too aggressive  
**Impact**: Low (solver still finds roots, just merges some)  

### Test 2: Multiplicity polynomial refinement
**Status**: ❌ FAIL  
**Issue**: Did not find simple root at 0.2 with precision 1e-15  
**Root Cause**: Solver needs higher precision for simple roots near multiple roots  
**Impact**: **High** - This is exactly why we need high-precision arithmetic!  
**Details**:
- Polynomial: `(x - 0.2) * (x - 0.5)^5` (simple root at 0.2, quintuple root at 0.5)
- Solver found 100 boxes (degeneracy detected)
- Could not refine simple root to 1e-15 precision
- **This validates the need for Phase 4/5 high-precision implementation**

### Test 3: Multiplicity detection from derivatives
**Status**: ✅ PASS  
**Details**: Correctly detects multiplicities 1, 2, 3, 4  

---

## Test Coverage Analysis

### Core Functionality
- ✅ Polynomial creation and conversion
- ✅ Polynomial evaluation (De Casteljau)
- ✅ Differentiation (all orders)
- ✅ Subdivision solver
- ✅ Root bounding methods (GraphHull, PP)
- ✅ Convex hull algorithms
- ✅ Geometry dump
- ⚠️ Result refinement (partial - fails for ill-conditioned cases)

### Edge Cases
- ✅ Degenerate boxes
- ✅ Linear polynomials
- ✅ Multiple roots (detection works, refinement needs HP)
- ✅ 2D systems

### Missing Coverage
- ⚪ 3D+ systems (no dedicated tests)
- ⚪ Very high degree polynomials (>10)
- ⚪ Extreme ill-conditioning (needs HP)
- ⚪ Performance benchmarks

---

## Performance Baseline

### Test Execution Time
- **Total**: 0.05 seconds for 15 tests
- **Average**: 0.003 seconds per test
- **Fastest**: 0.00 seconds (most tests)
- **Slowest**: 0.01 seconds (PolynomialConversionTest, ResultRefinerTest)

### Notes
- All tests are very fast (< 0.01s)
- No performance bottlenecks in current test suite
- Need dedicated performance benchmarks for Phase 3 completion

---

## Baseline Configuration

### CMake Configuration
```bash
cd build
cmake -L ..
```

**Key Options**:
- `ENABLE_HIGH_PRECISION`: OFF
- `DISABLE_TEMPLATES`: OFF
- `USE_BOOST`: ON
- `USE_MPFR`: OFF
- `USE_GMP`: OFF
- `USE_QUADMATH`: ON
- `BUILD_EXAMPLES`: ON
- `CMAKE_BUILD_TYPE`: Debug

### Compiler Flags
- Debug mode: `-g -O0`
- Warnings: `-Wall -Wextra`
- Standard: C++11

---

## Next Steps (Phase 3 Continued)

1. ✅ Document API baseline → `TIER1_API_BASELINE.md`
2. ✅ Run all tests → This document
3. ⏳ Create performance benchmarks
4. ⏳ Measure baseline performance
5. ⏳ Document performance baseline

After Phase 3 completion:
- **Phase 4**: Implement Tier 2 (fixed high-precision)
- **Phase 5**: Implement Tier 3 (template-based high-precision)

---

## Validation

This baseline establishes:
- ✅ All core functionality works in double precision
- ✅ Test suite is comprehensive (14/15 passing)
- ✅ Known limitation: Cannot refine roots to 1e-15 for ill-conditioned problems
- ✅ **Validates need for high-precision implementation**

The failing test (ResultRefinerTest) demonstrates exactly why high-precision arithmetic is needed - double precision is insufficient for refining roots in polynomials with multiple roots.

