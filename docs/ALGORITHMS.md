# Algorithm Overview

This document describes the algorithms used in the polynomial solver.

## Root Bounding Methods

### GraphHull Method

- Computes convex hull of graph control points in R^{n+1}
- Intersects with hyperplane x_{n+1} = 0
- Projects to R^n parameter space
- Exact for linear functions (machine epsilon error)

### ProjectedPolyhedral Method

- Processes each direction independently
- Projects to 2D (direction + function value)
- Computes convex hull and intersects with axis
- Simpler and more modular than GraphHull

## Subdivision Solver Workflow

1. Start with region [0,1]^n
2. Compute bounding box using selected method
3. If empty: discard (no roots)
4. If small enough: mark as converged
5. Else: subdivide and add to queue
6. Process boxes by depth (breadth-first)
7. Degeneracy detection for degenerate cases

## Direct Contraction Implementation

The solver uses **direct contraction** to minimize error accumulation:

### Traditional Approach
Restrict polynomials incrementally from current [0,1] to [a,b] repeatedly:
- Error grows linearly: ε ≈ k·n·ε·||b|| (k = number of contractions)

### Direct Contraction
Restrict from original [0,1] to global [A,B] each time:
- Error stays constant: ε ≈ n·ε·||b|| (independent of k)
- **2-8× error reduction** in practice
- **Same computational cost** (2 de Casteljau subdivisions per contraction)
- Better accuracy at extreme precision (tolerance < 10⁻¹²)

### Implementation Details

Each subdivision node stores both current and original polynomials:
- During contraction: polynomials recomputed from original using global coordinates
- During subdivision: incremental approach kept for efficiency

## Subdivision Strategies

Three strategies available:

1. **ContractFirst** (default): Pure contraction approach
2. **SubdivideFirst**: Subdivision-first approach  
3. **Simultaneous**: Balanced approach

## Degeneracy Handling

The solver detects degenerate cases:
- Multiple roots
- Infinite solutions
- Ill-conditioned systems

See [DEGENERATE_BOXES.md](DEGENERATE_BOXES.md) for details.

## Result Refinement

After subdivision, Newton refinement improves accuracy:

1. Use box centers as initial guesses
2. Apply Newton iteration with sign checking
3. Detect multiplicity for multiple roots
4. Estimate condition numbers

See [result_refinement_design.md](result_refinement_design.md) for details.

## References

- [GEOMETRY_ALGORITHMS.md](GEOMETRY_ALGORITHMS.md): Convex hull, hyperplane intersection
- [DEGENERATE_BOXES.md](DEGENERATE_BOXES.md): Degeneracy detection
- [CONDITIONING_AND_PRECISION.md](CONDITIONING_AND_PRECISION.md): Condition numbers

