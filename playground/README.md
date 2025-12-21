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
  -hp               Use high-precision refinement (~1e-30 vs ~1e-6 accuracy)
  -q                Quiet mode (output: total_boxes max_box_err max_refined_err)
```

### Example Output

```bash
# Normal mode
$ ./hessian_zero_set -r 1.0 -s 2 -n 1500
Hessian Zero Set Finder
=======================
Region: [-1, 1]^2
Subdivisions: 2x2
Polynomial degree: 10 (Hessian det degree: 16)
Target boxes/subregion: 1500
Solver tolerance: 1e-06
Expected radius: 0.707107

Region [0,0]: 6514 boxes
...
Total boxes: 26041
Refinement: double precision (h=1e-5, tol=1e-5)

=== Results ===
Refined: 20131/26041
Max box error:     1.407757e-03
Max refined error: 3.253795e-06

# Quiet mode for scripting (outputs: boxes box_err refined_err)
$ ./hessian_zero_set -q
48386 0.000666702 3.22552e-06

# High-precision refinement (~1e-30 accuracy)
$ ./hessian_zero_set -hp -q
48386 0.000666702 2.34e-31
```

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
// Interpolate nonlinear function as polynomial
Polynomial f = Interpolation::interpolate2D(func, degree, degree,
    u_min, u_max, v_min, v_max, AbscissaeType::CHEBYSHEV);

// Compute symbolic Hessian matrix
auto H = Differentiation::hessian(f);

// Compute determinant
Polynomial det_H = H[0][0] * H[1][1] - H[0][1] * H[0][1];

// Configure and solve
SubdivisionConfig config = defaultSolverConfig();
config.tolerance = 1e-6;
config.degeneracy_multiplier = target_boxes / (degree * degree);

Solver solver;
auto result = solver.subdivisionSolve(PolynomialSystem({det_H}), config);

// Refine (double precision)
refineCurveNumerical(g_func, x0, y0, CurveRefinementConfig{1e-5, 1e-5, 50});

// Refine (high precision, ~1e-30 accuracy)
refineCurveNumericalHP(g_hp_func, x0, y0, CurveRefinementConfigHP{"1e-20", "1e-30", 100});
```

## Directory Contents

- **`hessian_zero_set.cpp`** - Main example (see above)
- **`docs/`** - Technical documentation on algorithms and methods
- **`research_tests/`** - Research and debugging test files

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
