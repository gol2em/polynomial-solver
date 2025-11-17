# Geometry Data Structures and Algorithms

This document provides detailed information about the geometry module's data structures and algorithms.

## Table of Contents

1. [Data Structures](#data-structures)
2. [Convex Hull Algorithms](#convex-hull-algorithms)
3. [Intersection Algorithms](#intersection-algorithms)
4. [Complexity Analysis](#complexity-analysis)
5. [Usage Examples](#usage-examples)

---

## Data Structures

### `ConvexPolyhedron`

Represents a convex polyhedron in n-dimensional space.

```cpp
struct ConvexPolyhedron {
    std::vector<std::vector<double>> vertices;  // List of vertices
    std::size_t intrinsic_dim;                  // Intrinsic dimension
    
    std::size_t ambient_dimension() const;      // Dimension of ambient space
};
```

**Fields:**
- `vertices`: List of vertex coordinates. Each vertex is a `std::vector<double>` of size n (ambient dimension).
- `intrinsic_dim`: The intrinsic (geometric) dimension of the object:
  - `0`: Point
  - `1`: Line segment or line
  - `2`: Polygon or plane
  - `3`: 3D polyhedron
  - etc.

**Key Properties:**
- For **2D convex polygons** (ambient_dim=2, intrinsic_dim=2):
  - Vertices are stored in **counter-clockwise order**
  - Adjacent vertices in the list form edges of the polygon
  - The last vertex connects back to the first vertex (implicit edge)
  
- For **line segments** (intrinsic_dim=1):
  - Exactly 2 vertices (endpoints)
  
- For **points** (intrinsic_dim=0):
  - Exactly 1 vertex

**Example:**
```cpp
// Square with vertices in CCW order
ConvexPolyhedron square;
square.vertices = {{0,0}, {1,0}, {1,1}, {0,1}};
square.intrinsic_dim = 2;
// Edges: (0,0)-(1,0), (1,0)-(1,1), (1,1)-(0,1), (0,1)-(0,0)
```

### `ConvexPolyhedronBox`

Represents an axis-aligned bounding box.

```cpp
struct ConvexPolyhedronBox {
    std::vector<double> min_coords;  // Minimum coordinates
    std::vector<double> max_coords;  // Maximum coordinates
    
    std::size_t dimension() const;   // Dimension of the box
};
```

**Usage:** Efficient representation for axis-aligned boxes, used for bounding and intersection tests.

---

## Convex Hull Algorithms

### 2D Convex Hull: Monotone Chain Algorithm

**Function:** `ConvexPolyhedron convex_hull(const std::vector<std::vector<double>>& points)`

**Algorithm:** Andrew's monotone chain algorithm (variant of Graham's scan)

**Steps:**
1. **Sort points** lexicographically by (x, y) coordinates
2. **Remove duplicates** (points with identical coordinates)
3. **Build lower hull** (left to right):
   - Add points while maintaining right turns (CCW orientation)
   - Remove points that create left turns or collinear segments
4. **Build upper hull** (right to left):
   - Same process in reverse direction
5. **Concatenate** lower and upper hulls (omitting duplicate endpoints)

**Output Guarantees:**
- Vertices are in **counter-clockwise order**
- Adjacent vertices form edges of the convex hull
- Last vertex implicitly connects to first vertex
- Handles degenerate cases:
  - 0 points → empty
  - 1 point → single point (intrinsic_dim=0)
  - 2 points → line segment (intrinsic_dim=1)
  - Collinear points → line segment with 2 endpoints
  - 3+ non-collinear points → polygon (intrinsic_dim=2)

**Complexity:**
- Time: **O(n log n)** (dominated by sorting)
- Space: **O(n)** for storing hull vertices

**Code Location:** `src/geometry.cpp`, lines 827-902

### 3D Convex Hull: Incremental Algorithm

**Algorithm:** Incremental convex hull construction

**Steps:**
1. Start with initial tetrahedron from first 4 non-coplanar points
2. For each remaining point:
   - Find visible faces (faces where point is above the plane)
   - Remove visible faces
   - Add new faces connecting point to horizon edges
3. Return final set of vertices

**Complexity:**
- Time: **O(n²)** worst case, **O(n log n)** expected for random points
- Space: **O(n)** for storing hull

**Code Location:** `src/geometry.cpp`, lines 600-750

---

## Intersection Algorithms

### Hyperplane Intersection with Last Coordinate Zero

**Function:** `bool intersect_convex_polyhedron_with_last_coordinate_zero(...)`

**Purpose:** Intersect a convex polyhedron with the hyperplane where the last coordinate equals zero (e.g., z=0 in 3D, y=0 in 2D).

**Algorithm:**

**For 2D (ambient_dim=2):**
1. Compute convex hull of input points
2. Collect vertices on the hyperplane (|y| ≤ ε)
3. **Check only adjacent pairs** (edges) for intersections:
   - For each edge (v[i], v[(i+1)%n]):
     - If edge straddles hyperplane (one above, one below)
     - Compute intersection point using linear interpolation
4. Compute convex hull of all collected points
5. Return result

**Optimization:** Uses the fact that 2D hull vertices are in order, so only O(n) edges need to be checked instead of O(n²) pairs.

**For 3D (ambient_dim=3):**
1. Compute convex hull of input points
2. Collect vertices on the hyperplane (|z| ≤ ε)
3. **Check all pairs** of vertices for intersections:
   - For each pair (v[i], v[j]) where i < j:
     - If segment straddles hyperplane
     - Compute intersection point
4. Compute convex hull of all collected points
5. Return result

**Note:** For 3D, we check all pairs because we don't have explicit edge information from the convex hull algorithm.

**Complexity:**
- **2D:** O(n) for edge checking + O(m log m) for final hull (m ≤ n)
- **3D:** O(n²) for pair checking + O(m log m) for final hull

**Epsilon:** Uses ε = 1e-12 (machine precision) for determining if a point lies on the hyperplane.

**Code Location:** `src/geometry.cpp`, lines 1316-1428

### Polygon Intersection (2D)

**Function:** `bool intersect_convex_polygons_2d(...)`

**Algorithm:** Sutherland-Hodgman polygon clipping

**Steps:**
1. Start with first polygon
2. For each edge of second polygon:
   - Clip first polygon against the half-plane defined by the edge
   - Keep only the portion on the "inside" of the edge
3. Return clipped polygon

**Handles:**
- Polygon-polygon intersection
- Segment-polygon intersection (via clipping)
- Point-polygon intersection (via containment test)

**Complexity:**
- Time: **O(nm)** where n, m are the number of vertices
- Space: **O(n+m)** for storing result

**Code Location:** `src/geometry.cpp`, lines 1027-1095

### General Polyhedra Intersection

**Function:** `bool intersect_convex_polyhedra(...)`

**Handles mixed-dimensional cases:**
- Point-point: Check if identical
- Point-segment: Check if point on segment
- Point-polygon: Check if point in polygon
- Segment-segment: Compute line segment intersection
- Segment-polygon: Clip segment against polygon
- Polygon-polygon: Use Sutherland-Hodgman algorithm

**For higher dimensions:** Falls back to axis-aligned bounding box intersection.

**Code Location:** `src/geometry.cpp`, lines 1097-1314

---

## Complexity Analysis

### Summary Table

| Operation | 2D | 3D | Notes |
|-----------|----|----|-------|
| Convex Hull | O(n log n) | O(n²) worst, O(n log n) expected | 2D uses monotone chain |
| Hyperplane Intersection | O(n) | O(n²) | 2D optimized for ordered vertices |
| Polygon Intersection | O(nm) | N/A | Sutherland-Hodgman clipping |
| Bounding Box | O(n) | O(n) | Simple min/max computation |

### Space Complexity

All algorithms use **O(n)** space for storing vertices and intermediate results.

---

## Usage Examples

### Example 1: Compute 2D Convex Hull

```cpp
#include "geometry.h"

std::vector<std::vector<double>> points = {
    {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}
};

ConvexPolyhedron hull = convex_hull(points);

// hull.vertices contains 4 points in CCW order: (0,0), (1,0), (1,1), (0,1)
// hull.intrinsic_dim = 2
```

### Example 2: Intersect Polygon with y=0

```cpp
// Square from (-1, -1) to (1, 1)
std::vector<std::vector<double>> square_pts = {
    {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}
};

ConvexPolyhedron square = convex_hull(square_pts);
ConvexPolyhedron intersection;

bool has_intersection = intersect_convex_polyhedron_with_last_coordinate_zero(
    square, intersection);

// has_intersection = true
// intersection.vertices = {(-1, 0), (1, 0)}
// intersection.intrinsic_dim = 1 (line segment)
```

### Example 3: Intersect Two Polygons

```cpp
ConvexPolyhedron poly1 = /* ... */;
ConvexPolyhedron poly2 = /* ... */;
ConvexPolyhedron result;

bool ok = intersect_convex_polyhedra({poly1, poly2}, result);

if (ok) {
    // result contains the intersection
    // result.intrinsic_dim tells you if it's a point, segment, or polygon
}
```

---

## Implementation Notes

1. **Numerical Stability:** All geometric predicates use epsilon tolerance (1e-12) for robustness.

2. **Degenerate Cases:** All algorithms handle degenerate cases (points, line segments, collinear points).

3. **Intrinsic Dimension:** Automatically computed using `compute_intrinsic_dimension()` to distinguish between points, lines, and polygons.

4. **2D Optimization:** The 2D hyperplane intersection is optimized to O(n) by exploiting the ordered vertex property of the monotone chain algorithm.

5. **Debug Mode:** Set `getGeometryConfig().debug = true` to enable debug output for geometric operations.

---

## References

- **Monotone Chain Algorithm:** Andrew, A. M. (1979). "Another efficient algorithm for convex hulls in two dimensions"
- **Sutherland-Hodgman Algorithm:** Sutherland, I. E.; Hodgman, G. W. (1974). "Reentrant polygon clipping"
- **3D Convex Hull:** Preparata, F. P.; Shamos, M. I. (1985). "Computational Geometry: An Introduction"

