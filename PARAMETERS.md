# Command-Line Parameters

This document describes the command-line parameters available for the polynomial solver.

## Main Application (`polynomial_solver_app`)

### Usage

```bash
./bin/polynomial_solver_app [options]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--debug` | flag | disabled | Enable debug output for geometric operations |
| `--tolerance <value>` | double | 1e-3 | Box width tolerance (convergence criterion) |
| `--max-depth <value>` | unsigned int | 20 | Maximum subdivision depth |
| `--max-boxes <value>` | unsigned int | 1000 | Max boxes per depth for degeneracy detection |
| `--help` | flag | - | Show help message and exit |

### Examples

**Show help:**
```bash
./bin/polynomial_solver_app --help
```

**Run with default parameters:**
```bash
./bin/polynomial_solver_app
```

**Run with custom parameters:**
```bash
./bin/polynomial_solver_app --debug --tolerance 1e-4 --max-depth 25 --max-boxes 500
```

## Test Suite (`test_solver_subdivision`)

### Usage

```bash
./bin/test_solver_subdivision [options]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--debug` | flag | disabled | Enable debug output for geometric operations |

### Examples

**Run tests without debug output:**
```bash
./bin/test_solver_subdivision
```

**Run tests with debug output:**
```bash
./bin/test_solver_subdivision --debug
```

**Run tests with debug output and filter:**
```bash
./bin/test_solver_subdivision --debug 2>&1 | grep "DEBUG:"
```

## Parameter Details

### `--debug`

Enables detailed debug output for geometric operations, particularly for:
- Segment-polygon clipping operations
- Convex hull computations
- Intersection calculations

This is useful for:
- Debugging geometric algorithms
- Understanding why certain intersections succeed or fail
- Visualizing the clipping process

**Note:** Debug output can be verbose. Consider piping to a file or using grep to filter.

### `--tolerance <value>`

The box width tolerance is the absolute minimum width a box can have in any dimension before it's considered converged.

- **Range:** > 0.0
- **Default:** 1e-3
- **Meaning:** The solver iteratively contracts boxes using root bounding methods. When all dimensions of a box are smaller than this tolerance, the box is considered converged and returned as a root candidate.

**Smaller values** (e.g., 1e-6):
- Higher precision
- More subdivisions
- Slower computation

**Larger values** (e.g., 1e-2):
- Lower precision
- Fewer subdivisions
- Faster computation

**Workflow:**
1. Start with a box (initially [0,1]^n)
2. Compute bounding box of roots using the selected method
3. If empty, discard the box
4. If all dimensions < tolerance, converge
5. Else, make the bounding box the new region and iterate from step 2
6. If not converged after iteration, subdivide and continue

**Smaller values** (e.g., 1e-6):
- Higher precision
- More iterations/subdivisions
- Slower computation

**Larger values** (e.g., 1e-2):
- Lower precision
- Fewer iterations/subdivisions
- Faster computation

**Note:** This is different from the geometric tolerance (1e-12), which is used for geometric comparisons and remains at machine precision.

### `--max-depth <value>`

The maximum subdivision depth limits how many times a box can be subdivided.

- **Range:** > 0
- **Default:** 20
- **Meaning:** Stop subdividing after this many levels, even if convergence criteria are not met

**Higher values** (e.g., 30):
- Allow deeper subdivision
- Can find smaller root regions
- May be slower

**Lower values** (e.g., 10):
- Limit subdivision depth
- Faster but may miss small root regions
- Useful for quick exploratory runs

### `--max-boxes <value>`

The maximum number of boxes allowed at any single depth level before declaring a degenerate case.

- **Range:** > 0
- **Default:** 1000
- **Meaning:** If the number of boxes at a given depth exceeds this threshold, the solver terminates and returns a degeneracy marker

**Purpose:** Detect degenerate cases such as:
- Multiple roots with high multiplicity
- Infinite solution sets (e.g., a curve of roots in 2D systems)
- Numerical instabilities causing excessive subdivision

**Higher values** (e.g., 5000):
- Allow more boxes before declaring degeneracy
- May take longer to detect degenerate cases
- Useful for systems with many isolated roots

**Lower values** (e.g., 100):
- Detect degeneracy earlier
- Faster termination for degenerate cases
- May prematurely terminate for systems with many roots

**Degeneracy marker:** When degeneracy is detected, the solver returns a special marker box with:
- `depth = max_depth + 1`
- `lower = [-1, -1, ...]` (invalid coordinates)
- `upper = [-1, -1, ...]` (invalid coordinates)

## Demo Script

A demonstration script is provided to show all parameter combinations:

```bash
cd build
../demo_parameters.sh
```

This script demonstrates:
1. Help message
2. Default parameters
3. Custom parameters
4. Tests without debug
5. Tests with debug output

