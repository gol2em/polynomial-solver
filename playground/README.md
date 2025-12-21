# Playground

This directory contains experimental tests for the polynomial solver library.

## Samples

The `samples/` folder contains simple, well-tested examples demonstrating successful usage:
- **Circle Test**: Finding the zero set of x² + y² - 1 = 0

See `samples/README.md` for details.

## Experimental Tests

### 1. Circle Test (`test_with_geometry_dump`)
Tests the solver on a simple circle equation: x² + y² - 1 = 0

```bash
make test_with_geometry_dump
./test_with_geometry_dump
python3 visualize_circle_boxes_2d.py
```

### 2. Hessian Zero Set Finder (`hessian_zero_set`)

**The recommended example for finding zero sets of Hessian determinants of nonlinear functions.**

This demonstrates the complete workflow:
1. Divide domain into subregions for better local approximation
2. Interpolate f(x,y) as polynomial in each subregion
3. Compute symbolic Hessian using `Differentiation::hessian()`
4. Compute det(H) = H[0][0]\*H[1][1] - H[0][1]² using polynomial arithmetic
5. Find zero set using subdivision solver
6. Refine box centers onto curve using Newton's method

**Usage:**
```bash
make hessian_zero_set
./hessian_zero_set [options]

Options:
  -r <half_width>   Region = [-r, r]², default: 1.5
  -s <subdivisions> Subdivisions per axis, default: 4
  -d <degree>       Polynomial degree, default: 10
  -t <tolerance>    Solver tolerance, default: 1e-6
  -m <max_depth>    Max subdivision depth, default: 15
  -q                Quiet mode (output: boxes box_err refined_err)
```

**Example function** (modify in source for your application):
- f(x,y) = exp(-(x² + y²))  (Gaussian surface)
- Zero set of det(H): circle of radius r = 1/√2 ≈ 0.7071

**Key API usage:**
```cpp
// 1. Interpolate nonlinear function as polynomial
Polynomial f = Interpolation::interpolate2D(func, degree, degree,
    0.0, 1.0, 0.0, 1.0, AbscissaeType::CHEBYSHEV);

// 2. Compute symbolic Hessian
auto H = Differentiation::hessian(f);

// 3. Compute determinant using polynomial arithmetic
Polynomial det_H = H[0][0] * H[1][1] - H[0][1] * H[0][1];

// 4. Solve
Solver solver;
auto result = solver.subdivisionSolve(PolynomialSystem({det_H}), config);

// 5. Refine using CurveRefiner (polynomial) or numerical gradient (function)
```

### 3. Legacy Hessian Test (`test_hessian_determinant`)
Older test for a complex analytical function. See source for details.

## Building

The playground has an **auto-generated Makefile** that stays in sync with the library's dependencies.

### Quick Start

```bash
cd playground
make test_program_name
./test_program_name
```

That's it! The Makefile is automatically generated when you run `cmake` in the main project.

### Workflow

1. **Configure the main project** (first time or after adding dependencies):
   ```bash
   cd /path/to/polynomial-solver
   cmake -B build
   ```
   This generates `playground/Makefile` with all the correct link flags.

2. **Create your test file**:
   ```bash
   cd playground
   # Create test.cpp with your experimental code
   ```

3. **Compile and run**:
   ```bash
   make test
   ./test
   ```

### Building All Programs

```bash
cd playground
make all
```

This compiles all `.cpp` files in the playground directory.

### Cleaning Up

```bash
make clean
```

### How It Works

- When you run `cmake` in the main project, it generates `playground/Makefile` from `Makefile.in`
- The generated Makefile includes all necessary:
  - Compiler flags
  - Include paths
  - Link libraries (polynomial_solver, MPFR, GMP, Boost, etc.)
  - High-precision flags (if enabled)
- You just use `make <name>` to compile `<name>.cpp` → `<name>` executable
- No need to manually specify `-I`, `-L`, `-l` flags!

## Research Tests

The `research_tests/` subdirectory also has an auto-generated Makefile:

```bash
cd playground/research_tests
make test_multiplicity_methods
./test_multiplicity_methods
```

**Note:** Research tests may require updates to compile with the current codebase.

## Why This Approach?

✅ **Simple**: Just `make test` to compile `test.cpp`
✅ **Automatic**: Dependencies tracked by CMake, no manual updates needed
✅ **Fast**: Direct compilation, no CMake overhead for quick experiments
✅ **Flexible**: Add any `.cpp` file and compile it immediately
