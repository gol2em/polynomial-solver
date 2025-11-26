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

**Description:** Target precision for refined roots. Newton's method iterates until the box size is smaller than this value.

**Type:** `double`  
**Default:** `1e-15`  
**Range:** `1e-16` to `1e-10`

**Usage:**
```cpp
RefinementConfig config;
config.target_tolerance = 1e-15;  // Maximum precision with double
```

**Guidelines:**
- **1e-15**: Maximum precision with double (recommended)
- **1e-12**: Slightly relaxed, faster convergence
- **1e-10**: Relaxed, for quick refinement

**Note:** This is limited by double precision (~16 decimal digits). For higher precision, use multiprecision arithmetic (future feature).

---

### 5. Residual Tolerance (`residual_tolerance`)

**Description:** Maximum acceptable residual |f(x)| for a refined root to be considered valid.

**Type:** `double`  
**Default:** `1e-12`  
**Range:** `1e-16` to `1e-8`

**Usage:**
```cpp
RefinementConfig config;
config.residual_tolerance = 1e-12;  // Verify root is accurate
```

**Guidelines:**
- **1e-14**: Very strict, may reject valid roots in ill-conditioned problems
- **1e-12**: Default, good balance
- **1e-10**: Relaxed, accepts more roots

**Warning:** Small residual doesn't guarantee accurate root for ill-conditioned problems! Check the condition number.

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
config.dump_geometry = true;
config.dump_prefix = "dumps/my_problem";
```

**Note:** Only available when compiled with `ENABLE_GEOMETRY_DUMP=ON` (default).

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
refine_config.target_tolerance = 1e-15;  // Maximum precision
refine_config.residual_tolerance = 1e-12;
```

**Expected:** Slower solve, high precision roots

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
./build/bin/example_cubic_1d

# Custom tolerance
./build/bin/example_cubic_1d -t 1e-10

# Custom tolerance and depth
./build/bin/example_cubic_1d -t 1e-12 -d 150

# All parameters
./build/bin/example_cubic_1d -t 1e-10 -d 150 -m 10.0 --dump-geometry

# Show help
./build/bin/example_cubic_1d --help
```

---

## Summary Table

| Parameter | Default | Range | Impact | Command-Line |
|-----------|---------|-------|--------|--------------|
| `tolerance` | 1e-8 | 1e-15 to 1e-4 | Solver precision | `-t` |
| `max_depth` | 100 | 10 to 1000 | Max subdivisions | `-d` |
| `degeneracy_multiplier` | 5.0 | 2.0 to 20.0 | Degeneracy detection | `-m` |
| `target_tolerance` | 1e-15 | 1e-16 to 1e-10 | Refinement precision | N/A |
| `residual_tolerance` | 1e-12 | 1e-16 to 1e-8 | Residual check | N/A |
| `max_iterations` | 100 | 10 to 1000 | Newton iterations | N/A |
| `dump_geometry` | false | true/false | Visualization | `--dump-geometry` |

