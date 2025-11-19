# Build Configuration Options

This document describes the build-time configuration options for the Polynomial Solver project.

## Geometry Dump Feature

The geometry dump feature allows detailed debugging and visualization of the solver's internal geometric computations during the subdivision process.

### Overview

- **Purpose**: Export geometric data (control points, convex hulls, intersections) for debugging and visualization
- **Implementation**: Conditional compilation using preprocessor macros
- **Zero Overhead**: When disabled, dump code is completely removed at compile time

### Build Options

#### Development Mode (Default)
```bash
cmake -DENABLE_GEOMETRY_DUMP=ON ..
make
```

**Features:**
- Geometry dump functionality available
- `--dump [prefix]` command-line option works
- Slightly larger binary size (~5-10%)
- Minimal runtime overhead when dump is not used

#### Release Mode
```bash
cmake -DENABLE_GEOMETRY_DUMP=OFF ..
make
```

**Features:**
- Geometry dump code completely removed
- Smaller binary size
- Zero runtime overhead (no checks, no branches)
- `--dump` option will be ignored (if present in code)

### Usage

When built with `ENABLE_GEOMETRY_DUMP=ON`:

```bash
# Run solver with geometry dump
./bin/polynomial_solver_app --dump output_prefix

# Or in test code
SubdivisionConfig config;
config.dump_geometry = true;
config.dump_prefix = "debug_output";
```

This generates a text file with detailed geometric information that can be visualized using the provided Python scripts.

### Implementation Details

The feature uses C++11 preprocessor macros to conditionally compile dump I/O operations:

```cpp
bool compute_projected_polyhedral_bounds_with_dump(...) {
    // Core algorithm (always compiled)

#ifdef ENABLE_GEOMETRY_DUMP
    if (do_dump) {
        // Dump I/O operations (only compiled when enabled)
        dump << "data...";
    }
#endif

    // More core algorithm (always compiled)
}
```

**Key Design Decisions:**

1. **Single Implementation**: The algorithm exists in one place, always compiled the same way
2. **Conditional I/O**: Only the dump I/O operations are wrapped in `#ifdef` blocks
3. **Minimal Fragmentation**: Macros are used sparingly, only around I/O code
4. **Zero Overhead**: When disabled, dump I/O code is completely removed at compile time
5. **Synchronized Features**: Algorithm changes automatically apply to both dump and non-dump modes

### Performance Impact

| Configuration | Binary Size | Runtime Overhead | Use Case |
|---------------|-------------|------------------|----------|
| Dump Enabled  | +5-10%      | <1% (when not dumping) | Development, debugging |
| Dump Disabled | Baseline    | 0%               | Production, benchmarks |

### Recommendations

- **Development**: Use `ENABLE_GEOMETRY_DUMP=ON` (default)
- **Testing**: Use `ENABLE_GEOMETRY_DUMP=ON` for debugging failed tests
- **Benchmarking**: Use `ENABLE_GEOMETRY_DUMP=OFF` for accurate performance measurements
- **Production**: Use `ENABLE_GEOMETRY_DUMP=OFF` for deployed binaries

### Related Files

- `CMakeLists.txt`: Build configuration
- `src/solver.cpp`: Implementation with conditional compilation
- `visualize_ellipse_dump.py`: Visualization script for dump files
- `tests/test_ellipse_dump.cpp`: Example usage

### Future Considerations

If upgrading to C++17 or later, consider using `if constexpr` for cleaner conditional compilation:

```cpp
if constexpr (EnableDump) {
    // Dump code
}
```

This would provide better type checking and IDE support while maintaining zero runtime overhead.

