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
| `--crit <value>` | double | 0.8 | CRIT shrink factor for subdivision convergence |
| `--tolerance <value>` | double | 1e-3 | Minimum box width tolerance (convergence criterion) |
| `--max-depth <value>` | unsigned int | 20 | Maximum subdivision depth |
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
./bin/polynomial_solver_app --debug --crit 0.9 --tolerance 1e-4 --max-depth 25
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

### `--crit <value>`

The CRIT shrink factor determines when a contracted box is considered "small enough" to stop subdividing.

- **Range:** 0.0 to 1.0
- **Default:** 0.8
- **Meaning:** If the contracted box width is less than `crit * original_box_width`, the box is considered converged

**Higher values** (e.g., 0.9):
- More aggressive convergence criterion
- Fewer subdivisions
- Faster but potentially less accurate

**Lower values** (e.g., 0.5):
- More conservative convergence criterion
- More subdivisions
- Slower but potentially more accurate

### `--tolerance <value>`

The minimum box width tolerance is the absolute minimum width a box can have in any dimension before it's considered converged.

- **Range:** > 0.0
- **Default:** 1e-3
- **Meaning:** If any dimension of the box is smaller than this value, stop subdividing

**Smaller values** (e.g., 1e-6):
- Higher precision
- More subdivisions
- Slower computation

**Larger values** (e.g., 1e-2):
- Lower precision
- Fewer subdivisions
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

