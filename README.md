# Polynomial Solver

A C++11 implementation of a polynomial solver using De Casteljau subdivision algorithm.

GitHub repository: https://github.com/gol2em/polynomial-solver.git

## Project Structure

```
polynomial-solver/
├── include/              # Header files
│   ├── polynomial.h      # Polynomial class interface
│   ├── geometry.h        # Geometry tools interface
│   ├── de_casteljau.h    # De Casteljau algorithm interface
│   └── solver.h          # Main solver interface
├── src/                  # Source files
│   ├── polynomial.cpp    # Polynomial implementation
│   ├── geometry.cpp      # Geometry tools implementation
│   ├── de_casteljau.cpp  # De Casteljau algorithm implementation
│   ├── solver.cpp        # Main solver implementation
│   └── main.cpp          # Application entry point
├── tests/                # Unit tests (conversion, graph, etc.)
├── python/               # Python visualization tools
│   ├── visualize.py                      # General visualization entry point (stub)
│   └── visualize_graph_from_test.py      # Visualize graph/control-point test data
├── build/                # Build directory (generated)
├── .venv/                # Python virtual environment managed by uv
└── CMakeLists.txt        # CMake configuration
```

## Modules

### 1. Polynomial Module
- Polynomial representation and manipulation
- Arithmetic operations
- Evaluation and differentiation

### 2. Geometry Module
- 2D point representation
- Interval operations
- Geometric utilities

### 3. De Casteljau Module
- De Casteljau subdivision algorithm
- Curve subdivision and evaluation
- Control point manipulation

### 4. Solver Module
- Main solver interface
- Root finding algorithms
- Root isolation and refinement

### 5. Python Visualization
- Polynomial plotting
- Root visualization
- Subdivision process visualization

## Build Instructions

### Prerequisites
- CMake 3.15 or higher
- C++11 compatible compiler (GCC 13.2.0 or higher)
- Python 3.11 or higher (for visualization)
- uv (for managing the `.venv` Python virtual environment)

### Building the Project

```bash
# Create build directory
mkdir -p build
cd build

# Configure with CMake
cmake ..

# Build
make

# Run the application
./bin/polynomial_solver_app
```

### Running Tests

```bash
cd build
ctest
```

## Python Visualization

### Environment

```bash
# From the project root
uv sync  # or equivalent to create/manage .venv
source .venv/bin/activate
```

### Visualize graph/control-point test cases

After building and running the C++ tests (from `build/`):

```bash
cd build
ctest --output-on-failure
```

This writes CSV files into `build/tests/` for the graph/control-point tests.
From the project root you can then run:

```bash
python python/visualize_graph_from_test.py
```

This produces PNGs in the `python/` directory showing:
- the 1D graph and its Bernstein control points, and
- the 2D surface and its 3D Bernstein control net.

( `python/visualize.py` is kept as a stub for future, more general
visualization tools.)

## Development Status

Core components are implemented and tested:
- Multivariate Bernstein polynomial representation and stable power→Bernstein conversion
- De Casteljau evaluation in 1D and 2D
- Graph control-net construction for scalar polynomials and systems
- Unit tests for coefficient conversion and graph/control-net correctness

The solver’s higher-level root-finding and subdivision strategies are still evolving.

## Author

- **User**: gol2em
- **Email**: 664862601@qq.com

## License

TBD

