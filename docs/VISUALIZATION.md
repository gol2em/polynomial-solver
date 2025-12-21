# Visualization and Geometry Dump

This document describes the visualization tools and geometry dump feature for debugging and understanding the solver's behavior.

## Geometry Dump

### Enabling Geometry Dump

**Command-line (examples):**
```bash
# Run with geometry dumps for visualization
./build/bin/example_cubic_1d --dump-geometry
```

**Programmatic:**
```cpp
SubdivisionConfig config;
#ifdef ENABLE_GEOMETRY_DUMP
config.dump_geometry = true;
config.dump_prefix = "dumps/my_dump";  // Creates dumps/my_dump_geometry.txt
#endif
```

**Build-time control:**
```bash
# Debug mode: macro defined, feature controllable at runtime
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release mode: macro not defined, all dump code compiled out
cmake -DCMAKE_BUILD_TYPE=Release ..
```

The solver automatically creates the `dumps/` directory if needed.

## Python Visualization API

### Installation

```bash
pip install numpy matplotlib
# Or using uv
uv pip install numpy matplotlib
```

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
python tools/solver_viz_api.py dumps/example.txt --output-dir output/ --max-steps 20
```

## Visualization Output

Each iteration generates a PNG file with 3 subplots:

1. **Left**: First equation
   - Polynomial surface (blue/red colormap)
   - Zero contour (blue line)
   - Control points (orange/green dots)
   - Convex hulls (yellow/cyan polygons)
   - Intersections with z=0 plane (red lines)
   - Bounding intervals (red thick lines)

2. **Middle**: Second equation (same elements)

3. **Right**: Combined view
   - Both polynomial surfaces
   - Both zero contours
   - Current bounding box (red wireframe)
   - Final contracted box (green wireframe)

## Understanding the Visualization

- **Control Points**: Bernstein coefficients in graph space (x, y, f(x,y))
- **Convex Hull**: Convex hull of control points (used for root bounding)
- **Intersection**: Where convex hull intersects z=0 plane
- **Bounding Interval**: Projected interval on each axis
- **Zero Contour**: Where polynomial equals zero
- **Solution**: Intersection of both zero contours

## Visualization Features

- **Automatic Z-axis Scaling**: Adjusted to show geometry clearly
- **Progressive Zoom**: Shows previous iteration's box for context
- **Pruned Cases**: Boxes with empty bounds marked with red X
- **Decision Labels**: Shows CONTRACT/SUBDIVIDE/PRUNE decisions

## Circle-Ellipse Example

```bash
# Run example with geometry dumps
./build/bin/example_circle_ellipse

# Visualize results
source .venv/bin/activate
python examples/visualize_circle_ellipse.py

# Specific strategy with limited steps
python examples/visualize_circle_ellipse.py --strategy ContractFirst --max-steps 10
```

See [tools/README.md](../tools/README.md) for more details.

