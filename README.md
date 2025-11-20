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
- **Comprehensive Test Suite**: 12 test suites covering all major functionality
- **Python Visualization**: Optional visualization tools for graphs and control points

## Quick Start

See [QUICKSTART.md](QUICKSTART.md) for detailed setup instructions.

### 1. Check Prerequisites

```bash
./check_prerequisites.sh
```

### 2. Build the Project

```bash
# Build with tests and examples
./build.sh --test

# Or build without tests
./build.sh
```

### 3. Run Examples

**2D Example: Circle-Ellipse Intersection**

```bash
# Run the example (generates geometry dumps)
./build/bin/example_circle_ellipse

# Or run the test version
./build/bin/test_strategies
```

This generates geometry dump files in `dumps/` directory:
- `dumps/strategy_ContractFirst_geometry.txt`
- `dumps/strategy_SubdivideFirst_geometry.txt`
- `dumps/strategy_Simultaneous_geometry.txt`

**1D Example: Cubic Polynomial with 3 Roots**

```bash
# Run the 1D example
./build/bin/example_cubic_1d
```

This generates 1D geometry dumps:
- `dumps/cubic_1d_ContractFirst_geometry.txt`
- `dumps/cubic_1d_SubdivideFirst_geometry.txt`
- `dumps/cubic_1d_Simultaneous_geometry.txt`

### 4. Visualize the Results

First, set up Python environment (one-time setup):

```bash
# Create virtual environment using uv
uv venv .venv

# Activate virtual environment
source .venv/bin/activate

# Install dependencies
uv pip install numpy matplotlib
```

Then visualize:

**2D Visualization (Circle-Ellipse):**

```bash
# Visualize all strategies
python examples/visualize_circle_ellipse.py

# Visualize first 20 iterations only
python examples/visualize_circle_ellipse.py --max-steps 20

# Visualize specific strategy
python examples/visualize_circle_ellipse.py --strategy ContractFirst
```

**1D Visualization (Cubic Polynomial):**

```bash
# Visualize all strategies
python examples/visualize_cubic_1d.py

# Visualize first 10 iterations only
python examples/visualize_cubic_1d.py --max-steps 10

# Visualize specific strategy
python examples/visualize_cubic_1d.py --strategy SubdivideFirst
```

Visualizations are saved to `visualizations/viz_*/` directories as PNG files.

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
├── tests/                # Test suite (12 tests)
├── examples/             # Example programs and scripts
│   ├── circle_ellipse_intersection.cpp  # C++ example
│   └── visualize_circle_ellipse.py      # Python visualization script
├── tools/                # Visualization tools and API
│   ├── visualize_solver.py              # Core visualization engine
│   ├── solver_viz_api.py                # Python API for visualization
│   └── README.md                        # API documentation
├── python/               # Python visualization tools (legacy)
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

Enable geometry dumping in your code:

```cpp
SubdivisionConfig config;
config.dump_geometry = true;
config.dump_prefix = "dumps/my_dump";  // Creates dumps/my_dump_geometry.txt
```

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

All 12 test suites pass:

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
- **Email**: wenyd@lsec.cc.ac.cn

## License

TBD
