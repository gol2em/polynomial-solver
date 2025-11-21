# Polynomial Solver

A high-performance C++11 library for solving systems of multivariate polynomial equations using subdivision methods with Bernstein basis representation.

**GitHub**: https://github.com/gol2em/polynomial-solver.git
**Gitee**: https://gitee.com/gol2em/polynomial-solver.git

## Features

- **Simple 2-Line Workflow**: Solve and refine roots with just two function calls
- **Multivariate Polynomial Support**: Handle polynomials in multiple variables with arbitrary degrees
- **Bernstein Basis Representation**: Numerically stable polynomial representation with power↔Bernstein conversion
- **Multiple Root Bounding Methods**:
  - **GraphHull**: Exact convex hull-based bounding in graph space
  - **ProjectedPolyhedral**: Direction-by-direction projection method
  - **None**: Uniform subdivision (baseline)
- **Intelligent Degeneracy Detection**: Automatically detects degenerate cases (multiple roots, infinite solutions)
- **High-Precision Result Refinement**: Newton's method with sign checking to achieve 1e-15 precision (1D)
- **Condition Number Estimation**: Automatically detects when higher precision arithmetic is needed
- **Multiplicity Detection**: Determines root multiplicity using derivative analysis
- **Machine Epsilon Precision**: Achieves optimal error bounds (2.22×10⁻¹⁶) for linear systems
- **Robust Geometry**: Exact 2D/3D convex hull and hyperplane intersection algorithms
- **Unified Library**: Single `libpolynomial_solver.a` library for easy linking
- **Comprehensive Test Suite**: 13 test suites covering all major functionality
- **Python Visualization**: Optional visualization tools for graphs and control points
- **Configurable Parameters**: All solver parameters accessible via command-line flags

## Quick Start

### 1. Build the Project

```bash
# Build with tests and examples
./build.sh --test

# Or build without tests
./build.sh
```

### 2. Run the Simple Example

The simplest way to get started is with the 2-line workflow example:

```bash
./build/bin/example_simple_cubic
```

This demonstrates solving a cubic polynomial `(x - 0.2)(x - 0.5)(x - 0.8)` with just two lines of code:

```cpp
#include <polynomial_solver.h>

// LINE 1: SOLVE (fast, double precision)
Solver solver;
auto result = solver.subdivisionSolve(system, config);

// LINE 2: REFINE (high precision, 1e-15)
ResultRefiner refiner;
auto refined = refiner.refine(result, system, refine_config);
```

### 3. Experiment with Parameters

All examples support command-line parameters for experimentation:

```bash
# Run with custom tolerance
./build/bin/example_simple_cubic -t 1e-10

# Run with custom parameters
./build/bin/example_simple_cubic -t 1e-12 -d 150 -m 3.0

# Show all options
./build/bin/example_simple_cubic --help
```

**Common Command-line Options**:
- `-t, --tolerance <value>`: Box size tolerance (default: 1e-8)
- `-d, --max-depth <value>`: Maximum subdivision depth (default: 100)
- `-m, --degeneracy-multiplier <value>`: Degeneracy detection multiplier (default: 5.0)
- `--dump-geometry`: Enable geometry dump for visualization
- `-h, --help`: Show help message

See [docs/PARAMETERS.md](docs/PARAMETERS.md) for detailed parameter documentation.

### 4. More Examples

**Ill-Conditioned Problems (Wilkinson Polynomial)**:
```bash
# Demonstrates condition number estimation
./build/bin/example_wilkinson_1d
```

**Multiple Roots**:
```bash
# Demonstrates multiplicity detection
./build/bin/example_multiplicity_1d
```

**2D Systems (Circle-Ellipse Intersection)**:
```bash
# Demonstrates 2D polynomial system solving
./build/bin/example_circle_ellipse
```

### 5. Using the Library in Your Code

**Include the unified header**:
```cpp
#include <polynomial_solver.h>
using namespace polynomial_solver;
```

**Link against the library**:
```cmake
target_link_libraries(your_target polynomial_solver)
```

**Basic usage**:
```cpp
// 1. Define polynomial system
std::vector<unsigned int> degrees = {3};
std::vector<double> power_coeffs = {-0.08, 0.66, -1.5, 1.0};
Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
PolynomialSystem system(std::vector<Polynomial>{poly});

// 2. Solve (fast, double precision)
Solver solver;
SubdivisionConfig config;
config.tolerance = 1e-8;
auto result = solver.subdivisionSolve(system, config, RootBoundingMethod::ProjectedPolyhedral);

// 3. Refine (high precision, 1e-15)
ResultRefiner refiner;
RefinementConfig refine_config;
refine_config.target_tolerance = 1e-15;
auto refined = refiner.refine(result, system, refine_config);

// 4. Check results
for (const auto& root : refined.roots) {
    if (root.needs_higher_precision) {
        // Use higher precision arithmetic (future feature)
    }
}
```

## Project Structure

```
polynomial-solver/
├── include/                      # Header files
│   ├── polynomial_solver.h       # Unified header (include this!)
│   ├── polynomial.h              # Multivariate polynomial class
│   ├── solver.h                  # Main solver interface
│   ├── result_refiner.h          # High-precision refinement
│   ├── differentiation.h         # Polynomial differentiation
│   ├── geometry.h                # Geometric algorithms
│   └── de_casteljau.h            # De Casteljau algorithm
├── src/                          # Implementation files
│   ├── polynomial.cpp
│   ├── solver.cpp
│   ├── result_refiner.cpp
│   ├── differentiation.cpp
│   ├── geometry.cpp
│   └── de_casteljau.cpp
├── build/lib/                    # Build output
│   └── libpolynomial_solver.a    # Unified library (link this!)
├── tests/                        # Test suite (13 tests)
├── examples/                     # Example programs
│   ├── simple_cubic.cpp          # 2-line workflow demo
│   ├── cubic_1d_roots.cpp        # Detailed 1D example
│   ├── multiplicity_1d_roots.cpp # Multiple roots
│   ├── wilkinson_1d_roots.cpp    # Ill-conditioned
│   └── circle_ellipse_intersection.cpp  # 2D system
├── tools/                        # Tools and utilities
│   └── refine_from_dumps.cpp     # Root refinement tool
├── docs/                         # Documentation
│   ├── PARAMETERS.md             # Parameter reference
│   ├── CONDITIONING_AND_PRECISION.md  # Condition numbers
│   ├── GEOMETRY_ALGORITHMS.md    # Geometric algorithms
│   ├── DEGENERATE_BOXES.md       # Degeneracy handling
│   └── result_refinement_design.md    # Refinement design
├── build.sh                      # Automated build script
└── README.md                     # This file
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

### Direct Contraction Implementation

The solver uses **direct contraction** to minimize error accumulation:

- **Traditional approach**: Restrict polynomials incrementally from current [0,1] to [a,b] repeatedly
  - Error grows linearly: ε ≈ k·n·ε·||b|| (k = number of contractions)

- **Direct contraction**: Restrict from original [0,1] to global [A,B] each time
  - Error stays constant: ε ≈ n·ε·||b|| (independent of k)
  - **2-8× error reduction** in practice
  - **Same computational cost** (2 de Casteljau subdivisions per contraction)
  - Better accuracy at extreme precision (tolerance < 10⁻¹²)

**Implementation**: Each subdivision node stores both current and original polynomials. During contraction, polynomials are recomputed from original using global box coordinates. During subdivision, the incremental approach is kept for efficiency.

## Examples

See [examples/README.md](examples/README.md) for detailed documentation.

### Circle-Ellipse Intersection

The circle-ellipse intersection example demonstrates solving the system:
- f₁(x,y) = x² + y² - 1 = 0 (unit circle)
- f₂(x,y) = x²/4 + 4y² - 1 = 0 (ellipse)

Expected root: (2/√5, 1/√5) ≈ (0.894427, 0.447214)

**Run the example:**
```bash
./build/bin/example_circle_ellipse
```

This generates geometry dumps in `dumps/` for all three strategies:
- ContractFirst: Pure contraction approach
- SubdivideFirst: Subdivision-first approach
- Simultaneous: Balanced approach

**Visualize the results:**
```bash
# Activate Python environment
source .venv/bin/activate

# Visualize all strategies
python examples/visualize_circle_ellipse.py

# Or visualize specific strategy with limited steps
python examples/visualize_circle_ellipse.py --strategy ContractFirst --max-steps 10
```

## Visualization API

The project provides a Python API for visualizing solver output. See [tools/README.md](tools/README.md)
for detailed documentation.

### Using the API

```python
from tools.solver_viz_api import visualize_solver

# Visualize all iterations
visualize_solver('dumps/strategy_ContractFirst_geometry.txt', 'output/')

# Visualize first 20 iterations
visualize_solver('dumps/example.txt', 'output/', max_steps=20)
```

### Command-line Usage

```bash
# Visualize using the API directly
python tools/solver_viz_api.py dumps/example.txt --output-dir output/ --max-steps 20
```

### Generating Dump Files

**Command-line control (examples):**

All examples support the `--dump-geometry` flag to enable geometry dumping:

```bash
# Run without geometry dumps (default, faster)
./build/bin/example_cubic_1d

# Run with geometry dumps for visualization
./build/bin/example_cubic_1d --dump-geometry
```

**Programmatic control:**

Enable geometry dumping in your code (only when `ENABLE_GEOMETRY_DUMP` macro is defined):

```cpp
SubdivisionConfig config;
#ifdef ENABLE_GEOMETRY_DUMP
config.dump_geometry = true;
config.dump_prefix = "dumps/my_dump";  // Creates dumps/my_dump_geometry.txt
#endif
```

**Build-time control:**

The `ENABLE_GEOMETRY_DUMP` macro is controlled by CMake:

```bash
# Enable geometry dump (default, for development)
cmake -DENABLE_GEOMETRY_DUMP=ON ..

# Disable geometry dump (for release builds, removes all dump code)
cmake -DENABLE_GEOMETRY_DUMP=OFF ..
```

When `ENABLE_GEOMETRY_DUMP=OFF`, all geometry dumping code is compiled out for maximum performance.

The solver automatically creates the `dumps/` directory if it doesn't exist.

**Standard tests** (like `test_strategies`) automatically generate dumps in `dumps/` directory.

### Visualization Output

Each iteration generates a PNG file with 3 subplots:

1. **Subplot 1 (Left)**: First equation
   - Polynomial surface (blue/red colormap)
   - Zero contour (blue line, the solution curve)
   - Control points (orange/green dots)
   - Convex hulls (yellow/cyan polygons)
   - Intersections with z=0 plane (red lines)
   - Bounding intervals (red thick lines with square markers)

2. **Subplot 2 (Middle)**: Second equation
   - Same visualization elements as subplot 1

3. **Subplot 3 (Right)**: Combined view
   - Both polynomial surfaces
   - Both zero contours (the intersection is the solution)
   - Current bounding box (red wireframe)
   - Final contracted box (green wireframe)

### Understanding the Visualization

- **Control Points**: Bernstein coefficients in graph space (x, y, f(x,y))
- **Convex Hull**: Convex hull of control points (used for root bounding)
- **Intersection**: Where convex hull intersects z=0 plane (potential root region)
- **Bounding Interval**: Projected interval on each axis (red thick line)
- **Zero Contour**: Where polynomial equals zero (the actual solution curve)
- **Solution**: Intersection of both zero contours in subplot 3

### Visualization Features

- **Automatic Z-axis Scaling**: Z-axis limits are automatically adjusted to show geometry clearly
- **Progressive Zoom**: Each iteration shows the previous iteration's box for context
- **Pruned Cases**: Boxes with empty bounding boxes are marked with red X
- **Decision Labels**: Each iteration shows the decision (CONTRACT/SUBDIVIDE/PRUNE)

### Requirements

```bash
# Install Python dependencies
pip install numpy matplotlib

# Or using uv (recommended)
uv pip install numpy matplotlib
```

## Test Suite

All 13 test suites pass:

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
| DifferentiationTest | Polynomial differentiation |
| ResultRefinerTest | High-precision Newton refinement |

## Documentation

### User Documentation
- **[docs/PARAMETERS.md](docs/PARAMETERS.md)**: Complete parameter reference with tuning guide
- **[docs/CONDITIONING_AND_PRECISION.md](docs/CONDITIONING_AND_PRECISION.md)**: Understanding condition numbers and when higher precision is needed

### Technical Documentation
- **[docs/GEOMETRY_ALGORITHMS.md](docs/GEOMETRY_ALGORITHMS.md)**: Geometric algorithms (convex hull, hyperplane intersection)
- **[docs/DEGENERATE_BOXES.md](docs/DEGENERATE_BOXES.md)**: Degeneracy detection and handling
- **[docs/result_refinement_design.md](docs/result_refinement_design.md)**: High-precision result refinement design

### API Documentation
- **[include/polynomial_solver.h](include/polynomial_solver.h)**: Unified header with comprehensive usage examples

## Moving to Another Environment

This project is fully portable. To move to another environment:

1. **Copy the entire project directory**
2. **Run prerequisite checker**: `./check_prerequisites.sh`
3. **Install missing dependencies** (if any)
4. **Build**: `./build.sh --test`

All dependencies are standard C++11 libraries. No external libraries required.

## Author

- **User**: gol2em
- **Email**: wenyd@lsec.cc.ac.cn

## License

TBD
