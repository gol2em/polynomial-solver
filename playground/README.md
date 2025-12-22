# Playground

This directory contains the primary example and research materials for the polynomial solver library.

## Quick Start: Finding Zero Sets of Nonlinear Functions

The **`hessian_zero_set.cpp`** example demonstrates the complete workflow for finding zero sets of Hessian determinants (or any nonlinear function you define).

### Build and Run

```bash
cd playground
make hessian_zero_set
./hessian_zero_set
```

### Command Line Options

```
./hessian_zero_set [options]
  -r <half_width>   Region = [-r, r]², default: 1.5
  -s <subdivisions> Subdivisions per axis, default: 4
  -d <degree>       Polynomial degree for interpolation, default: 10
  -n <boxes>        Target boxes per subregion, default: 2000
  -t <tolerance>    Solver box tolerance, default: 1e-6
  -hp               Use high-precision refinement
  -b <bits>         Precision bits for HP mode, default: 128
  -o <file>         Output refined points to file (x y per line)
  -q                Quiet mode (output: total_boxes max_box_err max_refined_err)
```

### Example Output

```bash
# Normal mode with fewer boxes
$ ./hessian_zero_set -n 500 -s 2
Hessian Zero Set Finder
=======================
Region: [-1.5, 1.5]^2
Subdivisions: 2x2
Polynomial degree: 10 (Hessian det degree: 16)
Target boxes/subregion: 500
...
Total boxes: 9485
Refinement: double precision (h=1e-5, tol=1e-5)

=== Results ===
Refined: 9485/9485
Max box error:     2.457554e-03
Max refined error: 3.028793e-06

# Quiet mode for scripting (outputs: boxes box_err refined_err)
$ ./hessian_zero_set -n 500 -s 2 -q
9485 0.00245755 3.02879e-06

# High-precision refinement (default 128 bits → ~19 digits)
$ ./hessian_zero_set -n 500 -s 2 -hp -q
9485 0.00245755 2.82e-19

# Higher precision (256 bits → ~39 digits)
$ ./hessian_zero_set -n 500 -s 2 -hp -b 256 -q
9485 0.00245755 3.03e-39

# Output points to file for visualization
$ ./hessian_zero_set -n 500 -s 2 -hp -o dumps/hessian_points.txt -q
9485 0.00245755 2.82e-19
```

### Visualizing Results

Use the Python tool to visualize refined points:

```bash
# From the project root directory (uses root .venv)
cd playground
../.venv/bin/python visualize_refined_points.py dumps/hessian_points.txt --show-expected --analyze

# Save to file
../.venv/bin/python visualize_refined_points.py dumps/hessian_points.txt -o output.png
```

**Output files location:** `playground/dumps/`

Example comparison (f(x,y) = exp(-x²-y²), expected circle at r = 1/√2):

| Mode | Max Error | Point Uniformity (CV) |
|------|-----------|----------------------|
| Double precision | ~3e-6 | 4.79 |
| HP 128-bit | ~1e-16 (machine ε) | 1.69 |

The HP mode achieves machine epsilon accuracy for the radius.

### Customizing for Your Function

Modify the **USER-DEFINED FUNCTION** section at the top of the source file:

```cpp
// === USER-DEFINED FUNCTION ===
double f_user(double x, double y) {
    return std::exp(-(x*x + y*y));  // Replace with your function
}

double expected_radius() {
    return 1.0 / std::sqrt(2.0);    // Replace with your expected value (for validation)
}
```

### Workflow Summary

1. **Divide** domain into subregions for better local polynomial approximation
2. **Interpolate** f(x,y) as polynomial using Chebyshev nodes
3. **Compute** symbolic Hessian matrix via `Differentiation::hessian()`
4. **Compute** det(H) = H₁₁·H₂₂ - H₁₂² using polynomial arithmetic
5. **Solve** for zero set using subdivision solver (PP method)
6. **Refine** box centers onto curve via Newton's method with numerical gradient

### Key API Functions

```cpp
// Domain2D handles coordinate transformations between user space and [0,1]²
Domain2D domain = Domain2D::symmetric(1.5);  // [-1.5, 1.5]²
// or: Domain2D domain(x_min, x_max, y_min, y_max);

// Wrap user function for interpolation on [0,1]²
auto f_unit = domain.wrapFunction(f_user);
Polynomial f = Interpolation::interpolate2D(f_unit, degree, degree,
    0.0, 1.0, 0.0, 1.0, AbscissaeType::CHEBYSHEV);

// Compute symbolic Hessian matrix and determinant
auto H = Differentiation::hessian(f);
Polynomial det_H = H[0][0] * H[1][1] - H[0][1] * H[0][1];

// Configure and solve
SubdivisionConfig config = defaultSolverConfig();
config.tolerance = 1e-6;
config.degeneracy_multiplier = target_boxes / (degree * degree);

Solver solver;
auto result = solver.subdivisionSolve(PolynomialSystem({det_H}), config);

// Transform results back to user coordinates
double x, y;
domain.fromUnit(u_result, v_result, x, y);

// Refine (double precision)
refineCurveNumerical(g_func, x0, y0, CurveRefinementConfig{1e-5, 1e-5, 50});

// High-precision refinement with automatic parameter tuning
PrecisionContext ctx(128);  // Set working precision (128 bits → ~38 digits)
auto hp_config = CurveRefinementConfigHP::fromPrecisionBits(128);
auto f_hp = [](const mpreal& x, const mpreal& y) { return exp(-(x*x + y*y)); };
auto hess_det_hp = makeHessianDetFunctionHP(f_hp, hp_config.step_size_str);
refineCurveNumericalHP(hess_det_hp, x0, y0, hp_config);
```

## Directory Contents

- **`hessian_zero_set.cpp`** - Main example (see above)
- **`visualize_refined_points.py`** - Python visualization tool
- **`dumps/`** - Output directory for point files and visualizations
- **`docs/`** - Technical documentation on algorithms and methods
- **`research_tests/`** - Research and debugging test files

## Setup Guide

### Prerequisites

```bash
# High-precision arithmetic requires Boost + MPFR + GMP:
#   - Boost: Provides boost::multiprecision wrapper (header-only)
#   - MPFR: Actual arbitrary-precision floating-point library
#   - GMP:  Multi-precision integers (MPFR dependency)
sudo apt install libboost-dev libmpfr-dev libgmp-dev  # Debian/Ubuntu

# Python venv in project root (for visualization)
cd /path/to/polynomial-solver
python -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib
```

### Build with High-Precision Support

```bash
# From project root
mkdir -p build && cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
cd ../playground
make hessian_zero_set
```

### Verify HP is working

```bash
# Double precision (max error ~1e-6)
./hessian_zero_set -n 500 -s 2 -q
# Output: 9485 0.00245755 3.06e-06

# High precision (max error ~1e-16, machine epsilon for double output)
./hessian_zero_set -n 500 -s 2 -hp -q
# Output: 9485 0.00245755 2.82e-19

# The HP mode refines internally with arbitrary precision (MPFR),
# then converts to double for output. The ~1e-19 residual shows the
# curve is found with much higher accuracy than double can represent.
```

## Building Your Own Tests

The playground has an auto-generated Makefile that stays in sync with library dependencies:

```bash
# Compile any .cpp file
make your_test
./your_test
```

The Makefile is generated when you run `cmake -B build` in the project root. It automatically includes:
- All compiler flags and include paths
- Link libraries (polynomial_solver, MPFR, GMP, etc.)
- High-precision flags (if enabled)

## Research Tests

The `research_tests/` subdirectory contains tests for algorithm development:

```bash
cd playground/research_tests
make test_multiplicity_methods
./test_multiplicity_methods
```

See `research_tests/README.md` for details.
