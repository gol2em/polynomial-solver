# Sample Usage

This folder contains simple, well-tested examples demonstrating successful usage of the polynomial solver library.

## Circle Test

**File**: `circle_test.cpp`

Finds the zero set of the unit circle equation: x² + y² - 1 = 0

### Usage

```bash
# Compile
g++ -std=c++11 -Wall -I../../include -DENABLE_GEOMETRY_DUMP circle_test.cpp -o circle_test -L../../build/lib -lpolynomial_solver -lm

# Run with default degeneracy multiplier (5.0)
./circle_test

# Run with custom degeneracy multiplier
./circle_test 10.0

# Visualize
python3 visualize_circle_boxes_2d.py
```

### Expected Results

- **Converged boxes**: Points that converged to the circle (shown as green dots)
- **Unresolved boxes**: Boxes containing parts of the circle (shown as pink rectangles)
- **Degeneracy detected**: Yes (expected for curves)

With degeneracy_multiplier=5.0:
- 0 converged points
- 84 unresolved boxes
- Mean radius ≈ 0.999 (very close to unit circle)

With degeneracy_multiplier=20.0:
- 4 converged points
- 285 unresolved boxes
- Mean radius ≈ 0.9997 (even closer to unit circle)

### Key Features Demonstrated

1. **Single polynomial in 2D**: Finding zero set of a curve
2. **Degeneracy detection**: Properly handles curves (infinite roots)
3. **Converged vs unresolved**: Shows both types of results
4. **Visualization**: Overlays results on the actual curve

### Notes

- The polynomial is defined on [0,1]², so only the quarter circle in the first quadrant is found
- Converged boxes are points (lower = upper) where the solver determined the root location precisely
- Unresolved boxes are regions that need further subdivision but were stopped by degeneracy detection
- Both converged and unresolved boxes lie on or very close to the circle

