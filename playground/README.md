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

### 2. Hessian Determinant Zero Set (`test_hessian_determinant`)
Finds regions where a function is convex/concave/saddle by computing the zero set of det(Hessian).

**Function on [-1,1]²:**
- f1 = 10xy²x + √|x²y|
- f2 = atan2(0.001, sin(5y) - 2x)
- f3 = 10y³ + x³
- f4 = atan2(0.01, sin(5x) - 2y)
- f = f1 + f2 + f3 + f4

**Usage:**
```bash
make test_hessian_determinant
./test_hessian_determinant [degree]  # default degree=10
python3 visualize_hessian_det.py
```

**How it works:**
1. Transforms domain from [-1,1]² to [0,1]² via (x,y) = (2u-1, 2v-1)
2. Interpolates det(Hessian) with Bernstein polynomial of degree k×k
3. Finds zero set using the Projected Polyhedral method
4. Visualizes regions:
   - Red (det(H) > 0): Locally convex or concave
   - Blue (det(H) < 0): Saddle points
   - Black curve: Boundary (det(H) = 0)

**Parameters:**
- `degree`: Bernstein interpolation degree (higher = more accurate, slower)
- Recommended: 5-15 for good balance

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
