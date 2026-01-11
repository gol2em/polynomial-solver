# Templated Polynomial Solver Architecture

This document describes the templated polynomial solver architecture that provides
a unified interface across different scalar types (double, mpreal, etc.).

## Overview

The templated solver enables **precision-agnostic** polynomial root finding. The same
algorithms work with double precision for speed or arbitrary precision (via MPFR) for
accuracy on ill-conditioned problems like high-multiplicity roots.

## Components

### 1. PolynomialBase<Scalar> (`core/polynomial_base.h`)

Templated polynomial supporting both Power and Bernstein basis representations.

```cpp
#include "core/polynomial_base.h"
using Poly = PolynomialBase<double>;

// Create from power coefficients: 1 - 3x + 2x²
Poly p = Poly::fromPower({1.0, -3.0, 2.0});

// Operations
double val = p.evaluate(0.5);
Poly dp = p.differentiate(0);  // derivative w.r.t. variable 0
auto [left, right] = p.subdivide(0, 0.5);  // split at midpoint
```

### 2. SolverBase<Scalar> (`solver/solver_base.h`)

Bernstein subdivision solver that isolates roots into small boxes.

```cpp
#include "solver/solver_base.h"
using Solver = SolverBase<double>;
using Config = SubdivisionConfigBase<double>;

Config config;
config.tolerance = 1e-8;
config.max_depth = 100;

PolynomialSystemBase<double> system;
system.polynomials.push_back(poly.convertedToBernstein());

auto result = Solver::solve(system, config);
// result.boxes contains isolated root regions
```

### 3. ResultRefinerBase<Scalar> (`refinement/result_refiner_base.h`)

Newton-based root refinement with multiplicity detection and condition monitoring.

```cpp
#include "refinement/result_refiner_base.h"
using Refiner = ResultRefinerBase<double>;
using RefConfig = RefinementConfigBase<double>;

RefConfig config;
config.target_tolerance = 1e-15;
config.iteration_method = IterationMethodBase::ModifiedNewton;
config.multiplicity_method = MultiplicityMethodBase::Ostrowski;

auto result = Refiner::refineRoot1D(x0, poly, config);
// result.location, result.multiplicity, result.converged, etc.
```

## Iteration Methods

| Method | Description | Best For |
|--------|-------------|----------|
| `Newton` | Standard f/f' | Simple roots |
| `ModifiedNewton` | m·f/f' with known multiplicity | Multiple roots |
| `Halley` | Cubic convergence | Well-conditioned roots |
| `Schroder` | Quadratic, multiplicity-aware | Unknown multiplicity |

## Multiplicity Detection Methods

| Method | Description | Trade-off |
|--------|-------------|-----------|
| `Taylor` | Uses f, f', f'' ratios | Fast but less accurate |
| `Ostrowski` | Ratio of consecutive Newton steps | Most reliable |
| `SimpleThreshold` | |f'| < threshold | Simple heuristic |

## High-Multiplicity Root Workflow

For roots with multiplicity > 1, standard Newton converges slowly (linear instead of
quadratic). The recommended workflow:

### Step 1: Isolate with SolverBase<double>
```cpp
auto result = Solver::solve(system, config);
// Degenerate regions may indicate multiple roots
```

### Step 2: Attempt double-precision refinement
```cpp
auto refined = Refiner::refineBatch(result.boxes, poly, refine_config);
// Check refined.needs_higher_precision for problematic roots
```

### Step 3: HP refinement for flagged roots
```cpp
#ifdef ENABLE_HIGH_PRECISION
using HP = mpreal;
using RefinerHP = ResultRefinerBase<HP>;

mpfr::mpreal::set_default_prec(512);  // 512 bits ≈ 154 decimal digits

PolyHP poly_hp = PolyHP::fromPower(hp_coeffs);
auto hp_result = RefinerHP::refineRoot1D(HP(x0), poly_hp, hp_config);
// hp_result achieves full precision even for mult-6 roots
#endif
```

## Example: Mult-6 Root

```cpp
// p(x) = (x - 0.2)(x - 0.6)^6
// Expected: x = 0.2 (mult=1), x = 0.6 (mult=6)

// Double precision struggles with x = 0.6:
// - |f(x)| ~ 1e-80 at true root (underflows to 0)
// - |f'(x)| ~ 1e-70 (also underflows)
// - Newton step undefined or wrong

// HP with 512 bits succeeds:
// - Ostrowski detects mult = 6
// - Modified Newton: step = 6 * f/f'
// - Converges to x = 0.6 with error < 1e-15
```

## Type Aliases (Recommended)

```cpp
// Double precision
using Poly = PolynomialBase<double>;
using System = PolynomialSystemBase<double>;
using Solver = SolverBase<double>;
using Refiner = ResultRefinerBase<double>;
using SolverConfig = SubdivisionConfigBase<double>;
using RefineConfig = RefinementConfigBase<double>;

// High precision (requires ENABLE_HIGH_PRECISION)
#ifdef ENABLE_HIGH_PRECISION
using HP = mpreal;
using PolyHP = PolynomialBase<HP>;
using RefinerHP = ResultRefinerBase<HP>;
using RefineConfigHP = RefinementConfigBase<HP>;
#endif
```

## See Also

- `examples/multiplicity_1d_templated.cpp` - Complete HP escalation workflow
- `examples/simple_cubic_templated.cpp` - Basic 1D refinement
- `examples/wilkinson_1d_templated.cpp` - 20-root Wilkinson polynomial
- `examples/circle_ellipse_templated.cpp` - 2D Newton refinement

