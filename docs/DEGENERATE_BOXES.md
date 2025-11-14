# Degenerate Bounding Box Handling

## Overview

The solver can encounter degenerate bounding boxes during the subdivision process. This document describes how these cases are detected and handled.

## Types of Degenerate Boxes

### 1. Empty Box
- **Definition**: No intersection between equation constraints
- **Handling**: Box is discarded (no roots in this region)
- **Detection**: `compute_graph_hull_bounds()` returns `false`

### 2. Single Point (0D)
- **Definition**: All dimensions converged (width < tolerance in all directions)
- **Handling**: 
  - Check if the center point is actually a root using `isApproximateRoot()`
  - If yes: Add to results as a root
  - If no: Discard (false positive from conservative bounding)
- **Detection**: All dimensions have `upper[i] - lower[i] < tolerance`
- **Note**: Only applied when using `RootBoundingMethod::GraphHull`

### 3. Line Segment (1D)
- **Definition**: Exactly one dimension active (width >= tolerance), others converged
- **Handling**: 
  - Mark as converged for now
  - **TODO**: Implement 1D subdivision along the active axis
- **Detection**: Exactly one dimension has `upper[i] - lower[i] >= tolerance`
- **Note**: Only applied when using `RootBoundingMethod::GraphHull`

### 4. Higher Dimensional Box
- **Definition**: Multiple dimensions active (width >= tolerance)
- **Handling**: Normal subdivision behavior
- **Detection**: Two or more dimensions have `upper[i] - lower[i] >= tolerance`

## Implementation Details

### Root Verification

When a box converges to a single point, we verify it's actually a root:

```cpp
bool PolynomialSystem::isApproximateRoot(const std::vector<double>& point, 
                                         double tolerance) const
```

- Evaluates all equations at the point
- Returns `true` if all equation values are within `tolerance`
- Default tolerance: `1e-6`
- For converged boxes, tolerance is set to `max(1e-6, max_box_width)`

### Dimension Analysis

Helper function to classify box dimension:

```cpp
int analyze_box_dimension(const std::vector<double>& lower,
                          const std::vector<double>& upper,
                          double tolerance,
                          std::size_t& active_axis)
```

Returns:
- `1`: Single point (all dimensions < tolerance)
- `2`: Line segment (exactly one dimension >= tolerance)
- `3`: Higher dimensional box (multiple dimensions >= tolerance)

For line segments, `active_axis` is set to the index of the active dimension.

## When Degenerate Handling Applies

Degenerate box handling is **only applied** when using `RootBoundingMethod::GraphHull`.

When using `RootBoundingMethod::None` (uniform subdivision):
- Converged boxes are added to results without verification
- No dimension analysis is performed
- This is because uniform subdivision doesn't claim to find roots, just subdivide the domain

## Example: Single Point Root

For the 1D system `p(x) = x - 0.5` on `[0,1]`:

1. GraphHull method contracts the box around the root
2. Box converges to a small region containing `x = 0.5`
3. Dimension analysis detects single point (0D)
4. Center point `x = 0.5` is evaluated: `p(0.5) = 0`
5. Point is verified as a root and added to results

## Example: False Positive Filtering

Due to conservative bounding, the GraphHull method might produce boxes that don't actually contain roots:

1. Convex hull bounding is conservative (over-approximates root region)
2. Box converges to a single point
3. Center point is evaluated
4. If not a root (value > tolerance), box is discarded
5. This prevents false positives in the results

## Future Work

### 1D Subdivision for Line Segments

When a box degenerates to a line segment (1D problem):

1. Identify the active axis
2. Subdivide only along that axis
3. Apply 1D root-finding techniques
4. Continue until convergence or max depth

This is currently marked as **TODO** in the code.

### Multiplicity Detection

For roots with multiplicity > 1:

1. Detect when multiple boxes converge to the same point
2. Analyze the Jacobian to determine multiplicity
3. Report multiplicity information to the user

### Infinite Root Sets

For 2D systems with infinite roots (e.g., curves):

1. Detect when boxes don't converge even at max depth
2. Check if the intrinsic dimension of the root set is 1 (curve)
3. Switch to curve tracing algorithms
4. Parameterize and report the curve

## Testing

See `tests/test_degenerate_boxes.cpp` for comprehensive tests:

1. **Single point root**: 1D system with isolated root
2. **isApproximateRoot**: Verification method testing
3. **evaluate**: System evaluation method testing

All tests verify that:
- Degenerate boxes are detected correctly
- Root verification works properly
- False positives are filtered out
- Method only applies with GraphHull bounding

## API Reference

### PolynomialSystem Methods

```cpp
// Evaluate all equations at a point
void evaluate(const std::vector<double>& point, 
              std::vector<double>& values) const;

// Check if a point is approximately a root
bool isApproximateRoot(const std::vector<double>& point, 
                       double tolerance = 1e-6) const;
```

### Solver Configuration

```cpp
struct SubdivisionConfig {
    double tolerance;      // Convergence tolerance (default: 1e-3)
    unsigned int max_depth;  // Maximum subdivision depth (default: 20)
    // ...
};

enum class RootBoundingMethod {
    None,      // Uniform subdivision (no degenerate handling)
    GraphHull  // Graph convex hull (with degenerate handling)
};
```

