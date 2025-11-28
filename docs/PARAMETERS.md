# Solver Parameters Reference

## Overview

This document describes all parameters that influence the polynomial solver workflow, their default values, and how to tune them for different problems.

## Solver Parameters

### 1. Tolerance (`tolerance`)

**Description:** Box size threshold for convergence. When all dimensions of a box are smaller than this value, the box is considered converged.

**Type:** `double`  
**Default:** `1e-8`  
**Range:** `1e-15` to `1e-4`  
**Command-line:** `--tolerance` or `-t`

**Usage:**
```cpp
SubdivisionConfig config;
config.tolerance = 1e-10;  // Tighter tolerance for higher precision
```

**Guidelines:**
- **1e-6**: Fast, low precision (suitable for visualization)
- **1e-8**: Default, good balance (recommended for most problems)
- **1e-10**: High precision (for accurate initial guesses)
- **1e-12**: Very high precision (may be slow, use for critical applications)
- **< 1e-12**: Extreme precision (may hit numerical limits with double precision)

**Trade-offs:**
- Smaller tolerance → More accurate boxes, but slower (exponential growth in subdivisions)
- Larger tolerance → Faster, but less accurate initial boxes

**Note:** This is the **solver tolerance**, not the final root accuracy. Use `ResultRefiner` to achieve higher precision (1e-15) after solving.

---

### 2. Maximum Depth (`max_depth`)

**Description:** Maximum number of subdivision levels. Prevents infinite loops for degenerate or unsolvable problems.

**Type:** `unsigned int`  
**Default:** `100`  
**Range:** `10` to `1000`  
**Command-line:** `--max-depth` or `-d`

**Usage:**
```cpp
SubdivisionConfig config;
config.max_depth = 150;  // Allow deeper subdivisions
```

**Guidelines:**
- **50**: Fast, may miss roots in complex problems
- **100**: Default, sufficient for most problems
- **150**: For complex problems with many roots
- **200+**: For extremely complex problems (may be very slow)

**Depth vs. Box Size:**
```
depth = 0: box size = 1.0
depth = 10: box size ≈ 1e-3  (2^-10)
depth = 20: box size ≈ 1e-6  (2^-20)
depth = 30: box size ≈ 1e-9  (2^-30)
depth = 40: box size ≈ 1e-12 (2^-40)
```

**Trade-offs:**
- Larger max_depth → Can find roots in smaller regions, but slower
- Smaller max_depth → Faster, but may miss roots or fail to converge

---

### 3. Degeneracy Multiplier (`degeneracy_multiplier`)

**Description:** Multiplier for degeneracy detection. A box is considered degenerate if its size is less than `tolerance * degeneracy_multiplier` and it hasn't converged.

**Type:** `double`  
**Default:** `5.0`  
**Range:** `2.0` to `20.0`  
**Command-line:** `--degeneracy-multiplier` or `-m`

**Usage:**
```cpp
SubdivisionConfig config;
config.degeneracy_multiplier = 10.0;  // More aggressive degeneracy detection
```

**Guidelines:**
- **2.0**: Conservative, fewer false positives
- **5.0**: Default, good balance
- **10.0**: Aggressive, catches more degenerate cases
- **20.0**: Very aggressive, may flag non-degenerate cases

**What is Degeneracy?**

Degenerate cases include:
- **Multiple roots**: (x-0.5)² has a double root at x=0.5
- **Infinite solutions**: x² - x² = 0 (entire domain is a solution)
- **Curves/surfaces**: In 2D, a curve of solutions instead of isolated points

**Trade-offs:**
- Larger multiplier → Catches more degenerate cases, but may have false positives
- Smaller multiplier → Fewer false positives, but may miss degenerate cases

---

## Refinement Parameters

### 4. Target Tolerance (`target_tolerance`)

**Description:** Target error tolerance for condition-aware convergence and exclusion radius computation.

**Type:** `double`
**Default:** `1e-15`
**Range:** `1e-16` to `1e-10`
**Command-line:** `--target-tolerance`

**What it controls:**

1. **Condition-Aware Convergence** (PRIMARY PURPOSE):
   - Root is accepted only if `estimated_error < target_tolerance`
   - Estimated error = κ × |f(x)| / |f'(x)| where κ is the condition number
   - This prevents accepting inaccurate roots for ill-conditioned problems
   - Even if residual is small, root is rejected if estimated error is large

2. **Exclusion Radius Computation** (SECONDARY PURPOSE):
   - For simple roots: `radius ≈ target_tolerance / |f'(x)|`
   - For multiple roots: `radius ≈ (target_tolerance / |f^(m)(x)|)^(1/m)`
   - Roots within this radius are considered duplicates and merged

**Usage:**
```cpp
RefinementConfig config;
config.target_tolerance = 1e-15;  // For exclusion radius computation
```

**Guidelines:**
- **1e-15**: Maximum precision with double (recommended)
- **1e-12**: Slightly relaxed, larger exclusion radius
- **1e-10**: Relaxed, may merge distinct nearby roots

**Note:** This parameter affects how close two roots can be before being considered duplicates. Smaller values allow detecting roots that are closer together.

**Trade-offs:**
- Smaller tolerance → Smaller exclusion radius, can detect closer roots, may report duplicates
- Larger tolerance → Larger exclusion radius, merges nearby roots, fewer duplicates

---

### 5. Residual Tolerance (`residual_tolerance`)

**Description:** Residual threshold for triggering convergence check (NOT the acceptance criterion).

**Type:** `double`
**Default:** `1e-15`
**Range:** `1e-16` to `1e-8`
**Command-line:** `--residual-tolerance`

**What it controls:**
- Newton's method checks convergence when `|f(x)| < residual_tolerance`
- However, the root is NOT automatically accepted!
- The refiner uses **condition-aware convergence**:
  1. Check if `|f(x)| < residual_tolerance` (residual is small)
  2. Estimate condition number: κ ≈ |f''| / |f'|²
  3. Estimate actual error: error ≈ κ × |f(x)| / |f'(x)|
  4. Accept root only if `estimated_error < target_tolerance`

**Why Condition-Aware Convergence?**

Traditional residual-only convergence is INSUFFICIENT for ill-conditioned problems:
- For well-conditioned problems: small residual → small error ✓
- For ill-conditioned problems: small residual ≠ small error ✗

Example: (x - 0.5)³
- Residual: 6.2e-16 (machine epsilon) ✓
- Actual error: 8.5e-06 (10,000× larger!) ✗
- Condition number: 5.8e+29 (astronomical)

The condition-aware criterion rejects such roots, forcing the use of higher precision.

See [CONVERGENCE_CRITERIA_ANALYSIS.md](CONVERGENCE_CRITERIA_ANALYSIS.md) for details.

**Usage:**
```cpp
RefinementConfig config;
config.residual_tolerance = 1e-15;  // Converge when |f(x)| < 1e-15
```

**Guidelines:**
- **1e-15**: Strictest, maximum precision with double (recommended)
- **1e-12**: Slightly relaxed, faster convergence
- **1e-10**: Relaxed, for quick refinement

**Trade-offs:**
- Smaller tolerance → More iterations, stricter convergence, higher quality roots
- Larger tolerance → Fewer iterations, faster convergence, may be less accurate

**Example:**
```bash
# Use strictest residual tolerance (default)
./example_simple_cubic --residual-tolerance 1e-15

# Use more lenient residual tolerance for faster convergence
./example_simple_cubic --residual-tolerance 1e-12
```

---

### 6. Maximum Iterations (`max_iterations`)

**Description:** Maximum number of Newton iterations for refinement.

**Type:** `unsigned int`  
**Default:** `100`  
**Range:** `10` to `1000`

**Usage:**
```cpp
RefinementConfig config;
config.max_iterations = 100;  // Prevent infinite loops
```

**Guidelines:**
- **50**: Fast, may not converge for difficult problems
- **100**: Default, sufficient for most problems
- **200**: For difficult problems (high multiplicity, ill-conditioned)

---

### 7. Derivative Threshold (`derivative_threshold`)

**Description:** Threshold for detecting zero derivatives (for multiplicity detection).

**Type:** `double`  
**Default:** `1e-10`  
**Range:** `1e-14` to `1e-6`

**Usage:**
```cpp
// Internal parameter, not user-configurable yet
```

**Guidelines:**
- **1e-12**: Strict, fewer false positives for high multiplicity
- **1e-10**: Default, good balance
- **1e-8**: Relaxed, may detect multiplicity in ill-conditioned problems

**What is Multiplicity?**

A root at x=r has multiplicity m if:
- f(r) = 0
- f'(r) = 0
- ...
- f^(m-1)(r) = 0
- f^(m)(r) ≠ 0

**Example:** (x-0.5)³ has a triple root (m=3) at x=0.5

---

## Geometry Dump Parameters

### 8. Dump Geometry (`dump_geometry`)

**Description:** Enable geometry dump for visualization.

**Type:** `bool`
**Default:** `false`
**Command-line:** `--dump-geometry`

**Usage:**
```cpp
SubdivisionConfig config;
#ifdef ENABLE_GEOMETRY_DUMP
config.dump_geometry = true;
config.dump_prefix = "dumps/my_problem";
#endif
```

**Build Mode Behavior:**
- **Debug mode** (`CMAKE_BUILD_TYPE=Debug` or unspecified): The `ENABLE_GEOMETRY_DUMP` macro is defined. Geometry dumping is available and controlled by the `dump_geometry` runtime flag.
- **Release mode** (`CMAKE_BUILD_TYPE=Release`): The `ENABLE_GEOMETRY_DUMP` macro is not defined. All geometry dumping code is compiled out for maximum performance.

**Note:** Always wrap `dump_geometry` flag usage with `#ifdef ENABLE_GEOMETRY_DUMP` to ensure code compiles in both debug and release modes.

---

## Parameter Tuning Guide

### For Fast Solving (Visualization, Prototyping)

```cpp
SubdivisionConfig config;
config.tolerance = 1e-6;           // Relaxed tolerance
config.max_depth = 50;             // Shallow depth
config.degeneracy_multiplier = 5.0;
```

**Expected:** Fast solve, low precision boxes

---

### For Accurate Solving (Production, Scientific Computing)

```cpp
SubdivisionConfig config;
config.tolerance = 1e-10;          // Tight tolerance
config.max_depth = 150;            // Deep depth
config.degeneracy_multiplier = 5.0;

RefinementConfig refine_config;
refine_config.target_tolerance = 1e-15;   // Maximum precision
refine_config.residual_tolerance = 1e-15; // Strictest acceptance
```

**Expected:** Slower solve, high precision roots

**Command-line:**
```bash
./example_simple_cubic -t 1e-10 -d 150 --target-tolerance 1e-15 --residual-tolerance 1e-15
```

---

### For Ill-Conditioned Problems (Wilkinson, Multiple Roots)

```cpp
SubdivisionConfig config;
config.tolerance = 1e-8;           // Standard tolerance
config.max_depth = 100;
config.degeneracy_multiplier = 10.0;  // Aggressive degeneracy detection

RefinementConfig refine_config;
refine_config.max_iterations = 200;   // More iterations for convergence
```

**Expected:** Detects degeneracy, may need higher precision arithmetic

---

## Command-Line Usage

All examples support command-line parameters:

```bash
# Default parameters
./build/bin/example_simple_cubic

# Custom solver tolerance
./build/bin/example_simple_cubic -t 1e-10

# Custom solver tolerance and depth
./build/bin/example_simple_cubic -t 1e-12 -d 150

# Custom refinement tolerances
./build/bin/example_simple_cubic --target-tolerance 1e-12 --residual-tolerance 1e-12

# All parameters
./build/bin/example_simple_cubic -t 1e-10 -d 150 -m 10.0 \
    --target-tolerance 1e-12 --residual-tolerance 1e-12 --dump-geometry

# Show help
./build/bin/example_simple_cubic --help
```

---

## Summary Table

| Parameter | Default | Range | What It Controls | Command-Line |
|-----------|---------|-------|------------------|--------------|
| **Solver Parameters** |
| `tolerance` | 1e-8 | 1e-15 to 1e-4 | Box size threshold | `-t`, `--tolerance` |
| `max_depth` | 100 | 10 to 1000 | Max subdivisions | `-d`, `--max-depth` |
| `degeneracy_multiplier` | 5.0 | 2.0 to 20.0 | Degeneracy detection | `-m`, `--degeneracy-multiplier` |
| `dump_geometry` | false | true/false | Visualization | `--dump-geometry` |
| **Refinement Parameters** |
| `target_tolerance` | 1e-15 | 1e-16 to 1e-10 | Exclusion radius | `--target-tolerance` |
| `residual_tolerance` | 1e-15 | 1e-16 to 1e-8 | Convergence: \|f(x)\| < tol | `--residual-tolerance` |
| `max_iterations` | 100 | 10 to 1000 | Newton iterations | N/A |

