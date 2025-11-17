# Polynomial Solver

A high-performance C++11 library for solving systems of multivariate polynomial equations using subdivision methods with Bernstein basis representation.

**GitHub**: https://github.com/gol2em/polynomial-solver.git

## Features

- **Multivariate Polynomial Support**: Handle polynomials in multiple variables with arbitrary degrees
- **Bernstein Basis Representation**: Numerically stable polynomial representation with power↔Bernstein conversion
- **Multiple Root Bounding Methods**:
  - **GraphHull**: Exact convex hull-based bounding in graph space
  - **ProjectedPolyhedral**: Direction-by-direction projection method
  - **None**: Uniform subdivision (baseline)
- **Intelligent Degeneracy Detection**: Automatically detects degenerate cases (multiple roots, infinite solutions)
- **Machine Epsilon Precision**: Achieves optimal error bounds (2.22×10⁻¹⁶) for linear systems
- **Robust Geometry**: Exact 2D/3D convex hull and hyperplane intersection algorithms
- **Comprehensive Test Suite**: 11 test suites covering all major functionality
- **Python Visualization**: Optional visualization tools for graphs and control points

## Quick Start

See [QUICKSTART.md](QUICKSTART.md) for detailed setup instructions.

### 1. Check Prerequisites

```bash
./check_prerequisites.sh
```

### 2. Build

```bash
./build.sh --test
```

### 3. Run Examples

```bash
./build/bin/test_method_comparison
```

## Project Structure

```
polynomial-solver/
├── include/              # Header files
│   ├── polynomial.h      # Multivariate polynomial class
│   ├── solver.h          # Main solver interface
│   ├── geometry.h        # Geometric algorithms
│   └── de_casteljau.h    # De Casteljau algorithm
├── src/                  # Implementation files
│   ├── polynomial.cpp
│   ├── solver.cpp
│   ├── geometry.cpp
│   ├── de_casteljau.cpp
│   └── main.cpp
├── tests/                # Test suite (11 tests)
├── python/               # Python visualization tools
├── docs/                 # Documentation
├── build.sh              # Automated build script
├── check_prerequisites.sh # Prerequisite checker
└── QUICKSTART.md         # Quick start guide
```

## Algorithm Overview

### Root Bounding Methods

#### GraphHull Method
- Computes convex hull of graph control points in R^{n+1}
- Intersects with hyperplane x_{n+1} = 0
- Projects to R^n parameter space
- Exact for linear functions (machine epsilon error)

#### ProjectedPolyhedral Method
- Processes each direction independently
- Projects to 2D (direction + function value)
- Computes convex hull and intersects with axis
- Simpler and more modular than GraphHull

### Subdivision Solver Workflow

1. Start with region [0,1]^n
2. Compute bounding box using selected method
3. If empty: discard (no roots)
4. If small enough: mark as converged
5. Else: subdivide and add to queue
6. Process boxes by depth (breadth-first)
7. Degeneracy detection for degenerate cases

## Test Suite

All 11 test suites pass:

| Test | Description |
|------|-------------|
| PolynomialConversionTest | Power ↔ Bernstein conversion |
| PolynomialGraphTest | Graph control points |
| PolynomialSystemExampleTest | System evaluation |
| SolverSubdivisionTest | Subdivision solver |
| GeometryConvexTest | Convex hull algorithms |
| 2DHyperplaneIntersectionTest | Hyperplane intersection |
| 2DConvexHullRobustTest | Robust 2D convex hull |
| 2DHullVertexOrderTest | Vertex ordering |
| DegenerateBoxesTest | Degenerate case handling |
| LinearGraphHullTest | Linear function verification |
| ProjectedPolyhedralTest | PP method verification |

## Documentation

- **[QUICKSTART.md](QUICKSTART.md)**: Quick start guide
- **[docs/GEOMETRY_ALGORITHMS.md](docs/GEOMETRY_ALGORITHMS.md)**: Geometric algorithms
- **[docs/GEOMETRY_ROBUSTNESS.md](docs/GEOMETRY_ROBUSTNESS.md)**: Robustness improvements
- **[docs/DEGENERATE_BOXES.md](docs/DEGENERATE_BOXES.md)**: Degeneracy handling

## Moving to Another Environment

This project is fully portable. To move to another environment:

1. **Copy the entire project directory**
2. **Run prerequisite checker**: `./check_prerequisites.sh`
3. **Install missing dependencies** (if any)
4. **Build**: `./build.sh --test`

All dependencies are standard C++11 libraries. No external libraries required.

## Author

- **User**: gol2em
- **Email**: 664862601@qq.com

## License

TBD
