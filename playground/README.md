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

Use the Makefile to build individual test programs:

```bash
make test_program_name
```

The Makefile automatically:
- Compiles with C++11 standard
- Links against the library in `../build/lib`
- Includes headers from `../include`
- Enables geometry dump with `-DENABLE_GEOMETRY_DUMP`

## Cleaning Up

```bash
make clean
```
