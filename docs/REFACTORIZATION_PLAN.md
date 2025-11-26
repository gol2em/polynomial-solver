# Step-by-Step Refactorization Plan

## Overview: Three-Tier System

### Tier 1: Minimum Version (Default)
- **No multiprecision dependencies**
- **Only double precision**
- **Always available**
- **Zero overhead**

### Tier 2: Fallback Version (Partial Multiprecision)
- **Fixed high-precision support**
- **No templates** (separate implementation)
- **Enabled with**: `-DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON`
- **Supports**: One or two fixed precision types (e.g., mpreal with 256-bit)

### Tier 3: Full Version (Template-Based)
- **Flexible multiprecision support**
- **Template-based** (supports any numeric type)
- **Enabled with**: `-DENABLE_HIGH_PRECISION=ON` (default when high precision is on)
- **Supports**: double, __float128, mpreal at any precision

## Configuration Logic

```
Default: Tier 1 (minimum version, double only)
    ↓
ENABLE_HIGH_PRECISION=ON
    ↓
    ├─ Check for Boost headers → Not found → ERROR
    ├─ Check for MPFR/GMP → Not found → Use cpp_dec_float
    ↓
    DISABLE_TEMPLATES=OFF (default)
    ├─ YES → Tier 3 (full version, templates)
    └─ NO  → Tier 2 (fallback version, no templates)
```

## Step-by-Step Process

---

## Phase 1: CMake Configuration and Dependency Checking

### Step 1.1: Create CMake Module for Dependency Detection

**File**: `cmake/FindHighPrecision.cmake`

**Purpose**: Detect available multiprecision libraries and set configuration flags

**Tasks**:
1. Create `cmake/` directory
2. Create `FindHighPrecision.cmake` module
3. Implement detection logic:
   - Check for Boost headers
   - Check for MPFR library
   - Check for GMP library
   - Check for quadmath library
   - Determine which backend to use (MPFR vs cpp_dec_float)
4. Set CMake variables:
   - `HIGH_PRECISION_AVAILABLE` (bool)
   - `HIGH_PRECISION_BACKEND` (MPFR or CPP_DEC_FLOAT)
   - `QUADMATH_AVAILABLE` (bool)
   - `HIGH_PRECISION_INCLUDE_DIRS` (list)
   - `HIGH_PRECISION_LIBRARIES` (list)

**Estimated time**: 1 hour

### Step 1.2: Update Root CMakeLists.txt

**File**: `CMakeLists.txt`

**Purpose**: Add configuration options and include dependency detection

**Tasks**:
1. Add configuration options:
   ```cmake
   option(ENABLE_HIGH_PRECISION "Enable high-precision arithmetic" OFF)
   option(DISABLE_TEMPLATES "Disable template-based implementation" OFF)
   option(ENABLE_QUADMATH "Enable quadmath support" OFF)
   ```
2. Include `FindHighPrecision.cmake` module
3. Add validation logic:
   - If `ENABLE_HIGH_PRECISION=ON` but dependencies not found → ERROR with helpful message
   - If `DISABLE_TEMPLATES=ON` but `ENABLE_HIGH_PRECISION=OFF` → WARNING (ignored)
4. Set preprocessor definitions based on configuration:
   - `ENABLE_HIGH_PRECISION`
   - `USE_TEMPLATES` (if templates enabled)
   - `USE_MPFR_BACKEND` or `USE_CPP_DEC_FLOAT_BACKEND`
   - `ENABLE_QUADMATH` (if available)
5. Print configuration summary

**Estimated time**: 1 hour

---

## Phase 2: Create Type Definition Headers

### Step 2.1: Create High-Precision Type Definitions

**File**: `include/high_precision_types.h`

**Purpose**: Define multiprecision types based on configuration

**Tasks**:
1. Create header with conditional compilation:
   ```cpp
   #ifdef ENABLE_HIGH_PRECISION
       #ifdef USE_MPFR_BACKEND
           // MPFR backend
       #elif defined(USE_CPP_DEC_FLOAT_BACKEND)
           // cpp_dec_float backend
       #endif
       
       #ifdef ENABLE_QUADMATH
           // quadmath types
       #endif
   #endif
   ```
2. Define type aliases:
   - `mpreal` (for Boost.Multiprecision)
   - `quad` (for __float128, if available)
3. Define precision constants:
   - `DEFAULT_HP_PRECISION` (e.g., 256 bits)
   - `DEFAULT_HP_DECIMAL_DIGITS` (e.g., 77 digits)

**Estimated time**: 30 minutes

### Step 2.2: Create Configuration Header Template

**File**: `include/config.h.in` (CMake template)

**Purpose**: Generate compile-time configuration header

**Tasks**:
1. Create CMake template with placeholders:
   ```cpp
   #cmakedefine ENABLE_HIGH_PRECISION
   #cmakedefine USE_TEMPLATES
   #cmakedefine USE_MPFR_BACKEND
   #cmakedefine USE_CPP_DEC_FLOAT_BACKEND
   #cmakedefine ENABLE_QUADMATH
   ```
2. Update CMakeLists.txt to generate `config.h`:
   ```cmake
   configure_file(include/config.h.in include/config.h)
   ```
3. Add generated header to include path

**Estimated time**: 30 minutes

---

## Phase 3: Refactor Core Classes (Tier 1 - Minimum Version)

### Step 3.1: Review Current Implementation

**Files**: 
- `include/polynomial.h`
- `include/de_casteljau.h`
- `include/differentiation.h`
- `src/result_refiner.cpp`

**Purpose**: Understand current double-precision implementation

**Tasks**:
1. Document current API signatures
2. Identify all methods that need high-precision versions
3. List all dependencies between classes
4. Create API compatibility checklist

**Estimated time**: 1 hour

### Step 3.2: Ensure Tier 1 Works Without Changes

**Purpose**: Verify minimum version (double-only) still works

**Tasks**:
1. Build with default configuration (no flags)
2. Run all existing tests
3. Verify no performance regression
4. Document baseline performance

**Estimated time**: 30 minutes

---

## Phase 4: Implement Tier 2 (Fallback Version - No Templates)

### Step 4.1: Create High-Precision Polynomial Class

**File**: `include/polynomial_hp.h`

**Purpose**: Fixed high-precision polynomial class (no templates)

**Tasks**:
1. Create `PolynomialHP` class (copy from `Polynomial`)
2. Replace `double` with `mpreal` throughout
3. Keep same API structure as `Polynomial`
4. Wrap entire file in `#ifdef ENABLE_HIGH_PRECISION`
5. Add conversion methods:
   - `static PolynomialHP fromPolynomial(const Polynomial& p)`
   - `Polynomial toPolynomial() const` (with precision loss warning)

**Estimated time**: 2 hours

### Step 4.2: Create High-Precision De Casteljau

**File**: `include/de_casteljau_hp.h` and `src/de_casteljau_hp.cpp`

**Purpose**: Fixed high-precision De Casteljau evaluation

**Tasks**:
1. Create `DeCasteljauHP` class
2. Copy algorithms from `DeCasteljau`
3. Replace `double` with `mpreal`
4. Implement methods:
   - `static mpreal evaluate1D(const std::vector<mpreal>&, const mpreal&)`
   - `static mpreal evaluateTensorProduct(...)`
5. Wrap in `#ifdef ENABLE_HIGH_PRECISION`

**Estimated time**: 2 hours

### Step 4.3: Create High-Precision Differentiation

**File**: `include/differentiation_hp.h` and `src/differentiation_hp.cpp`

**Purpose**: Fixed high-precision differentiation

**Tasks**:
1. Create `DifferentiationHP` class
2. Copy algorithms from `Differentiation`
3. Replace `double` with `mpreal`, `Polynomial` with `PolynomialHP`
4. Implement methods:
   - `static PolynomialHP derivative(const PolynomialHP&, size_t axis, unsigned int order)`
   - `static std::vector<PolynomialHP> gradient(const PolynomialHP&)`
5. Wrap in `#ifdef ENABLE_HIGH_PRECISION`

**Estimated time**: 2 hours

### Step 4.4: Create High-Precision Result Refiner

**File**: `include/result_refiner_hp.h` and `src/result_refiner_hp.cpp`

**Purpose**: Fixed high-precision Newton refinement

**Tasks**:
1. Create `ResultRefinerHP` class
2. Define `RefinedRootHP` struct:
   ```cpp
   struct RefinedRootHP {
       std::vector<mpreal> root;
       mpreal residual;
       int iterations;
       bool converged;
   };
   ```
3. Implement methods:
   - `static RefinedRootHP refineRoot(const PolynomialHP&, const std::vector<mpreal>&)`
   - `static RefinedRootHP refineRoot1D(const PolynomialHP&, const mpreal&)`
4. Copy Newton iteration logic from `ResultRefiner`
5. Replace all double operations with mpreal
6. Wrap in `#ifdef ENABLE_HIGH_PRECISION`

**Estimated time**: 3 hours

### Step 4.5: Add Conversion Utilities

**File**: `include/precision_conversion.h` and `src/precision_conversion.cpp`

**Purpose**: Convert between double and high-precision types

**Tasks**:
1. Create utility functions:
   ```cpp
   #ifdef ENABLE_HIGH_PRECISION
       mpreal toHighPrecision(double x);
       std::vector<mpreal> toHighPrecision(const std::vector<double>& v);
       double toDouble(const mpreal& x);
       std::vector<double> toDouble(const std::vector<mpreal>& v);
   #endif
   ```
2. Add precision loss warnings (optional logging)
3. Wrap in `#ifdef ENABLE_HIGH_PRECISION`

**Estimated time**: 1 hour

### Step 4.6: Test Tier 2 Implementation

**File**: `tests/test_high_precision_fallback.cpp`

**Purpose**: Verify fallback version works correctly

**Tasks**:
1. Create test file (only compiled when `ENABLE_HIGH_PRECISION=ON` and `DISABLE_TEMPLATES=ON`)
2. Test polynomial evaluation
3. Test differentiation
4. Test Newton refinement
5. Test conversion utilities
6. Compare results with double precision
7. Verify precision improvement

**Estimated time**: 2 hours

**Total for Phase 4**: ~12 hours

---

## Phase 5: Implement Tier 3 (Full Version - Templates)

### Step 5.1: Create Template Implementation Headers

**File**: `include/detail/polynomial_impl.h`

**Purpose**: Private template implementation for Polynomial

**Tasks**:
1. Create `include/detail/` directory
2. Create template implementation:
   ```cpp
   namespace detail {
       template<typename T>
       class PolynomialImpl {
           // Template implementation
       };
   }
   ```
3. Move algorithm logic from `Polynomial` to `PolynomialImpl<T>`
4. Keep as header-only (in `detail` namespace)

**Estimated time**: 3 hours

### Step 5.2: Create Template De Casteljau

**File**: `include/detail/de_casteljau_impl.h`

**Purpose**: Private template implementation for De Casteljau

**Tasks**:
1. Create template functions:
   ```cpp
   namespace detail {
       template<typename T>
       T evaluate1D_impl(const std::vector<T>&, T);
       
       template<typename T>
       T evaluateTensorProduct_impl(...);
   }
   ```
2. Move algorithm logic to templates
3. Keep as header-only

**Estimated time**: 2 hours

### Step 5.3: Create Template Differentiation

**File**: `include/detail/differentiation_impl.h`

**Purpose**: Private template implementation for Differentiation

**Tasks**:
1. Create template functions in `detail` namespace
2. Move algorithm logic to templates
3. Handle template polynomial type

**Estimated time**: 2 hours

### Step 5.4: Update Public APIs to Use Templates

**Files**: 
- `include/polynomial.h`
- `include/de_casteljau.h`
- `include/differentiation.h`

**Purpose**: Make public APIs delegate to template implementations

**Tasks**:
1. Update `Polynomial` class:
   ```cpp
   class Polynomial {
   public:
       double evaluate(...) const {
           #ifdef USE_TEMPLATES
               return detail::PolynomialImpl<double>::evaluate(...);
           #else
               // Original implementation
           #endif
       }
   };
   ```
2. Similarly update `DeCasteljau` and `Differentiation`
3. Ensure API remains unchanged
4. Wrap template delegation in `#ifdef USE_TEMPLATES`

**Estimated time**: 3 hours

### Step 5.5: Create High-Precision Template APIs

**Files**:
- `include/polynomial_mp.h`
- `include/de_casteljau_mp.h`
- `include/differentiation_mp.h`

**Purpose**: High-precision APIs that delegate to same templates

**Tasks**:
1. Create `PolynomialMP` class:
   ```cpp
   #ifdef ENABLE_HIGH_PRECISION
       #ifdef USE_TEMPLATES
           class PolynomialMP {
           public:
               mpreal evaluate(...) const {
                   return detail::PolynomialImpl<mpreal>::evaluate(...);
               }
           };
       #endif
   #endif
   ```
2. Similarly create `DeCasteljauMP`, `DifferentiationMP`
3. All delegate to same template implementation
4. Wrap in `#ifdef ENABLE_HIGH_PRECISION && USE_TEMPLATES`

**Estimated time**: 2 hours

### Step 5.6: Add Quadmath Support (Optional)

**Files**:
- `include/polynomial_quad.h`
- `include/de_casteljau_quad.h`
- `include/differentiation_quad.h`

**Purpose**: Quadmath APIs that delegate to same templates

**Tasks**:
1. Create classes for __float128 type
2. Delegate to same template implementation
3. Wrap in `#ifdef ENABLE_QUADMATH && USE_TEMPLATES`

**Estimated time**: 2 hours

### Step 5.7: Test Tier 3 Implementation

**File**: `tests/test_high_precision_templates.cpp`

**Purpose**: Verify template version works correctly

**Tasks**:
1. Create test file (only compiled when `USE_TEMPLATES=ON`)
2. Test all precision types (double, mpreal, quad)
3. Test that all use same implementation
4. Verify no code duplication
5. Test precision flexibility with mpreal

**Estimated time**: 2 hours

**Total for Phase 5**: ~16 hours

---

## Phase 6: Integration and Testing

### Step 6.1: Create Unified Header

**File**: `include/polynomial_solver.h`

**Purpose**: Single include file that provides correct API based on configuration

**Tasks**:
1. Create master header that includes appropriate headers based on config
2. Provide type aliases for easy use:
   ```cpp
   // Always available
   using PolynomialD = Polynomial;
   
   #ifdef ENABLE_HIGH_PRECISION
       #ifdef USE_TEMPLATES
           using PolynomialHP = PolynomialMP;
       #else
           using PolynomialHP = PolynomialHP;  // Fallback version
       #endif
   #endif
   ```

**Estimated time**: 1 hour

### Step 6.2: Update Existing Tests

**Purpose**: Ensure all existing tests still pass

**Tasks**:
1. Run all tests with Tier 1 (default)
2. Verify no regressions
3. Check performance hasn't degraded

**Estimated time**: 1 hour

### Step 6.3: Create Comprehensive Test Suite

**File**: `tests/test_all_configurations.cpp`

**Purpose**: Test all three tiers

**Tasks**:
1. Test Tier 1 (double only)
2. Test Tier 2 (fallback, if enabled)
3. Test Tier 3 (templates, if enabled)
4. Compare results across tiers
5. Verify precision improvements

**Estimated time**: 3 hours

### Step 6.4: Create Example Programs

**Files**:
- `examples/example_double_precision.cpp`
- `examples/example_high_precision_fallback.cpp`
- `examples/example_high_precision_templates.cpp`

**Purpose**: Demonstrate usage of each tier

**Tasks**:
1. Create example for each tier
2. Show how to use each API
3. Demonstrate precision improvements
4. Add to documentation

**Estimated time**: 2 hours

**Total for Phase 6**: ~7 hours

---

## Phase 7: Documentation

### Step 7.1: Update Build Documentation

**File**: `docs/BUILD.md`

**Purpose**: Document build options and configurations

**Tasks**:
1. Document all CMake options
2. Explain three-tier system
3. Provide build examples for each tier
4. Add troubleshooting section

**Estimated time**: 2 hours

### Step 7.2: Create API Documentation

**File**: `docs/API.md`

**Purpose**: Document public APIs for each tier

**Tasks**:
1. Document Tier 1 API (double precision)
2. Document Tier 2 API (fallback high precision)
3. Document Tier 3 API (template high precision)
4. Show migration examples

**Estimated time**: 2 hours

### Step 7.3: Update README

**File**: `README.md`

**Purpose**: Update main documentation

**Tasks**:
1. Add high-precision features section
2. Update build instructions
3. Add quick start examples
4. Link to detailed documentation

**Estimated time**: 1 hour

**Total for Phase 7**: ~5 hours

---

## Summary Timeline

| Phase | Description | Estimated Time |
|-------|-------------|----------------|
| **Phase 1** | CMake configuration and dependency checking | 2 hours |
| **Phase 2** | Create type definition headers | 1 hour |
| **Phase 3** | Verify Tier 1 (minimum version) | 1.5 hours |
| **Phase 4** | Implement Tier 2 (fallback version) | 12 hours |
| **Phase 5** | Implement Tier 3 (template version) | 16 hours |
| **Phase 6** | Integration and testing | 7 hours |
| **Phase 7** | Documentation | 5 hours |
| **Total** | | **~44.5 hours** (~1 week) |

## Recommended Order

1. **Start with Phase 1-3** (infrastructure, ~4.5 hours)
2. **Then Phase 4** (fallback version, ~12 hours) - Gets you working high precision quickly
3. **Then Phase 5** (template version, ~16 hours) - Adds flexibility
4. **Finally Phase 6-7** (testing and docs, ~12 hours)

## Checkpoints

After each phase:
1. ✅ Commit changes
2. ✅ Push to gitee (then github)
3. ✅ Update progress in this document
4. ✅ Test that previous functionality still works

## Next Steps

Would you like me to:
1. Start with Phase 1 (CMake configuration)?
2. Create detailed file-by-file implementation plans?
3. Begin implementation of a specific phase?

