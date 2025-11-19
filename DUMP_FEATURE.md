# Geometry Dump Feature

## Overview

The DUMP feature allows you to export detailed geometric information during the solving process for debugging and visualization purposes. This is particularly useful for understanding how the ProjectedPolyhedral (PP) method works step-by-step.

## Usage

### Command Line

```bash
# Enable dump with default prefix "dump"
./build/bin/polynomial_solver_app --dump

# Enable dump with custom prefix
./build/bin/polynomial_solver_app --dump my_debug

# Combine with other options
./build/bin/polynomial_solver_app --dump --tolerance 1e-6 --max-depth 50
```

### Programmatic API

```cpp
#include "solver.h"

SubdivisionConfig config;
config.dump_geometry = true;           // Enable dumping
config.dump_prefix = "my_output";      // Set output prefix (default: "dump")
config.tolerance = 1e-6;

Solver solver;
SubdivisionSolverResult result = solver.subdivisionSolve(
    system, config, RootBoundingMethod::ProjectedPolyhedral);
```

## Output Format

When enabled, the solver creates a text file `{prefix}_geometry.txt` containing:

### Header
- Method used (ProjectedPolyhedral, GraphHull, None)
- Dimension of the problem
- Number of equations

### For Each Iteration
- **Iteration number** and **depth**
- **Global box** bounds (current subdivision box in original coordinates)
- For each **direction** (0 to n-1):
  - For each **equation**:
    - **Projected_Points**: All control points projected to 2D (coordinate + function value)
    - **ConvexHull**: Vertices of the 2D convex hull
    - **Intersection**: Intersection of convex hull with axis (where function = 0)
    - **Interval**: 1D interval bounds for this direction
  - **Final_Interval**: Combined interval after intersecting all equations
- **Bounding Box (local)**: Final bounds in local [0,1]^n space
- **Bounding Box (global)**: Final bounds in original coordinate space

### Example Output

```
# Polynomial Solver Geometry Dump
# Method: ProjectedPolyhedral
# Dimension: 2
# Number of equations: 2

# Iteration 0, Depth 0
# Global box: [0, 1, 0, 1]

## Direction 0

### Equation 0
Projected_Points 9
0 -1
0 -1
0 0
0.5 -1
0.5 -1
0.5 0
1 0
1 0
1 1
ConvexHull 5
0 -1
0.5 -1
1 0
1 1
0 0
Intersection 2
0 0
1 0
Interval [0, 1]

### Equation 1
Projected_Points 9
...
```

## Test Case: 2D Ellipse Intersection

A complete example is provided in `tests/test_ellipse_dump.cpp`:

**System:**
- Equation 1: x² + y² - 1 = 0 (circle)
- Equation 2: x²/4 + 4y² - 1 = 0 (ellipse)
- Domain: [0,1]²

**Expected solution:** (√0.8, √0.2) ≈ (0.894, 0.447)

**Run the test:**
```bash
cd build
./bin/test_ellipse_dump
```

This generates `ellipse_test_geometry.txt` with detailed geometric data for each iteration.

## Visualization

A Python script `visualize_ellipse_dump.py` is provided to create step-by-step visualizations:

```bash
# Generate visualizations from dump file
python3 visualize_ellipse_dump.py build/ellipse_test_geometry.txt
```

This creates PNG images in `visualization_output/` directory showing:
1. **Left subplot**: Equation 1 in 3D with projected control points and convex hulls
2. **Middle subplot**: Equation 2 in 3D with projected control points and convex hulls
3. **Right subplot**: Both equations in 2D with current box and PP bounds

### Visualization Features

Each 3D subplot shows:
- Polynomial surface (colored by value)
- z=0 plane (gray)
- Zero contour (blue curve, intersection with z=0)
- Projected control points on background planes (orange dots)
- Convex hulls of projections (orange lines)
- Intersection points with axis (red circles)
- Current box boundary (black rectangle)

The 2D subplot shows:
- Zero contours of both equations (blue and red curves)
- Current subdivision box (black rectangle)
- PP bounding box (purple shaded region)
- Expected analytical solution (yellow star)

## Files Modified

- `include/solver.h`: Added `dump_geometry` and `dump_prefix` to `SubdivisionConfig`
- `src/solver.cpp`: Implemented `compute_projected_polyhedral_bounds_with_dump()`
- `src/main.cpp`: Added `--dump` command-line option
- `tests/test_ellipse_dump.cpp`: Test case for 2D ellipse intersection
- `tests/CMakeLists.txt`: Added test to build system
- `visualize_ellipse_dump.py`: Visualization script

## Performance Impact

- **When disabled** (default): No performance impact
- **When enabled**: Minimal overhead (~5-10% slower due to file I/O)
- Dump files can be large for deep subdivisions (25KB for 12 iterations in the ellipse test)

## Use Cases

1. **Debugging**: Understand why a particular box is discarded or subdivided
2. **Algorithm visualization**: Create animations of the solving process
3. **Research**: Analyze convergence behavior and geometric properties
4. **Teaching**: Demonstrate how subdivision methods work
5. **Optimization**: Identify bottlenecks or inefficiencies in bounding methods

