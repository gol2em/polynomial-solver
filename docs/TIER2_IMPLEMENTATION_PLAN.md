# Tier 2 Implementation Plan

## Analysis Summary

### Current State

**High-precision is ONLY used for root refinement**, not for solving. The workflow is:

1. **Solve** (double precision) → Locate roots approximately
2. **Refine** (double/high precision) → Refine roots to 1e-15 precision

**Key Finding**: Only `ResultRefiner` needs high-precision support!

### What Needs High-Precision

From `src/result_refiner.cpp` analysis - **CORRECTED**:

#### Core Refinement Operations (Lines 177-860)
1. **Newton iteration** (lines 193-301)
   - `poly.evaluate(x)` - evaluate f(x) at HP point
   - `dpoly.evaluate(x)` - evaluate f'(x) at HP point
   - `ddpoly.evaluate(x)` - evaluate f''(x) at HP point
   - Arithmetic: `x - f/df` in HP

2. **Differentiation** - **NEEDS HP VERSION!**
   - `Differentiation::derivative(poly, axis, order)` - used 16 times in refiner
   - Creates derivative polynomials from double-precision coefficients
   - **Must support HP coefficient conversion**

3. **Condition number estimation** (lines 429-475)
   - Second/third derivative evaluation
   - Arithmetic: `|f''| / |f'|²` in HP

4. **Multiplicity estimation** (lines 477-558)
   - Higher-order derivative evaluation (up to order 10)
   - Uses `Differentiation::derivative()` and `Differentiation::gradient()`
   - Threshold comparison in HP

5. **Modified Newton for multiple roots** (lines 341-427)
   - m-th derivative evaluation
   - Factorial computation
   - Arithmetic: `x - m*f/f'` in HP

### What Does NOT Need High-Precision

- `Solver` class - stays double precision (subdivision works fine)
- `DeCasteljau` - stays double precision (only used in solver, not refiner)
- `Geometry` - stays double precision (only used in solver, not refiner)

### What DOES Need High-Precision (Corrected)

- `Polynomial::evaluate()` - **NEEDS HP VERSION** (evaluate at HP point)
- `Differentiation::derivative()` - **NEEDS HP VERSION** (create HP polynomial from double coefficients)
- `Differentiation::gradient()` - **NEEDS HP VERSION** (for 2D multiplicity estimation)
- `ResultRefiner` - **NEEDS HP VERSION** (all refinement logic)

### Tier 2 Implementation Strategy

**Option 1: Minimal Approach** (Recommended)
- Keep all existing code as-is (double precision)
- Add `ResultRefinerHP` class with high-precision Newton refinement
- User workflow:
  1. Solve with double precision (existing code)
  2. Refine with double precision (existing code)
  3. For roots with `needs_higher_precision=true`, use `ResultRefinerHP`

**Option 2: Full Approach** (Overkill)
- Create HP versions of all classes
- Much more work, but not needed since solving works fine with double

### Required Components for Tier 2 (CORRECTED)

#### 1. `PolynomialHP` class
**Purpose**: Store and evaluate polynomials with high-precision coefficients

**Methods needed**:
- `evaluate(mpreal x)` - evaluate at high-precision point using HP coefficients
- `evaluate(const std::vector<mpreal>& x)` - multivariate evaluation
- Constructor from double-precision `Polynomial` (converts coefficients to HP)
- Accessors: `degrees()`, `coefficients()`, `dimension()`

**Implementation**: Store Bernstein coefficients as `std::vector<mpreal>`

**NOT needed**:
- De Casteljau (not used in refinement)
- Subdivision (not used in refinement)
- Power basis conversion (not used in refinement)

#### 2. `DifferentiationHP` class
**Purpose**: Compute derivatives of HP polynomials

**Methods needed**:
- `derivative(const PolynomialHP& poly, axis, order) -> PolynomialHP`
- `gradient(const PolynomialHP& poly) -> std::vector<PolynomialHP>`

**Implementation**: Copy logic from `Differentiation`, use HP arithmetic for binomial coefficients

#### 3. `ResultRefinerHP` class
**Purpose**: Refine roots using high-precision Newton iteration

**Methods needed**:
- `refineRoot1D(double initial_guess, const PolynomialHP& poly, config) -> mpreal`
- `estimateMultiplicity()` - high-precision derivative evaluation
- `estimateConditionNumber()` - high-precision condition estimation

**Implementation**: Copy logic from `ResultRefiner`, replace `double` with `mpreal`

#### 4. Conversion utilities (Already exist!)
**Purpose**: Convert between double and high-precision

**Existing in `precision_conversion.h`**:
- `toHighPrecision(double x) -> mpreal`
- `toDouble(mpreal x) -> double`
- `toString(mpreal x, int digits) -> std::string`

**New functions needed**:
- `convertToHighPrecision(const Polynomial& poly) -> PolynomialHP`

### Include Organization

#### Current Structure (Good!)
```
include/
├── polynomial_solver.h          # Main header (users include this)
├── polynomial.h                 # Core classes
├── solver.h
├── result_refiner.h
├── differentiation.h
├── geometry.h
├── de_casteljau.h
├── high_precision_types.h       # Type definitions (mpreal, quad)
├── precision_conversion.h       # Conversion utilities
├── high_precision_refiner.h     # HP refinement (Tier 2/3)
└── precision_context.h          # RAII precision management
```

#### Proposed Tier 2 Structure
```
include/
├── polynomial_solver.h          # Main header (Tier 1 only)
├── polynomial_solver_hp.h       # NEW: Tier 2 header (includes HP support)
│
├── polynomial.h                 # Tier 1 (double precision)
├── solver.h                     # Tier 1 (double precision)
├── result_refiner.h             # Tier 1 (double precision)
│
├── polynomial_hp.h              # NEW: Tier 2 (high-precision polynomial)
├── result_refiner_hp.h          # NEW: Tier 2 (high-precision refinement)
│
├── high_precision_types.h       # Type definitions
├── precision_conversion.h       # Conversion utilities
└── precision_context.h          # RAII precision management
```

### User Workflow

#### Tier 1 (Double Precision Only)
```cpp
#include <polynomial_solver.h>
using namespace polynomial_solver;

// Solve and refine with double precision
Solver solver;
auto result = solver.subdivisionSolve(system, config);

ResultRefiner refiner;
auto refined = refiner.refine(result, system, refine_config);

// Check results
for (const auto& root : refined.roots) {
    if (root.needs_higher_precision) {
        std::cout << "Warning: Root may need higher precision\n";
    }
}
```

#### Tier 2 (Fixed High-Precision)
```cpp
#include <polynomial_solver_hp.h>  // Includes Tier 1 + HP support
using namespace polynomial_solver;

// Step 1: Solve with double precision (fast)
Solver solver;
auto result = solver.subdivisionSolve(system, config);

// Step 2: Refine with double precision (fast)
ResultRefiner refiner;
auto refined = refiner.refine(result, system, refine_config);

// Step 3: Re-refine problematic roots with high precision
for (const auto& root : refined.roots) {
    if (root.needs_higher_precision) {
        // Convert to high precision
        setPrecision(256);  // 256 bits = ~77 decimal digits
        PolynomialHP poly_hp = convertToHighPrecision(system.equation(0));
        
        // Refine with high precision
        ResultRefinerHP refiner_hp;
        mpreal root_hp = refiner_hp.refineRoot1D(
            root.location[0], poly_hp, refine_config);
        
        std::cout << "HP root: " << toString(root_hp, 50) << "\n";
    }
}
```

### Implementation Phases

**Phase 4a: Core Types** (2 hours)
- Create `polynomial_hp.h` and `polynomial_hp.cpp`
- Implement `PolynomialHP::evaluate(mpreal)`
- Implement conversion utilities

**Phase 4b: HP Refinement** (4 hours)
- Create `result_refiner_hp.h` and `result_refiner_hp.cpp`
- Implement `ResultRefinerHP::refineRoot1D()`
- Implement multiplicity/condition estimation

**Phase 4c: Testing** (3 hours)
- Test with multiplicity polynomial from Phase 3
- Verify 1e-15 precision achieved
- Test all three backends (MPFR, cpp_dec_float, quadmath)

**Phase 4d: Documentation** (3 hours)
- Update `polynomial_solver_hp.h` header
- Add examples
- Update README

**Total: 12 hours** (as estimated)

### Next Steps

1. ✅ Analyze current code (DONE)
2. Create `polynomial_hp.h/cpp` with evaluate method
3. Create `result_refiner_hp.h/cpp` with Newton refinement
4. Test with ill-conditioned polynomial
5. Document usage

