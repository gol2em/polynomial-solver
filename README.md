# Polynomial Solver

A C++11 implementation of a polynomial solver using De Casteljau subdivision algorithm.

## Project Structure

```
polynomial-solver/
├── include/              # Header files
│   ├── polynomial/       # Polynomial class headers
│   ├── geometry/         # Geometry tools headers
│   ├── de_casteljau/     # De Casteljau algorithm headers
│   └── solver/           # Main solver interface headers
├── src/                  # Source files
│   ├── polynomial/       # Polynomial implementation
│   ├── geometry/         # Geometry tools implementation
│   ├── de_casteljau/     # De Casteljau algorithm implementation
│   └── solver/           # Main solver implementation
├── tests/                # Unit tests
├── python/               # Python visualization tools
├── build/                # Build directory (generated)
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

```bash
cd python
python3 visualize.py
```

## Development Status

This project is currently in the initial setup phase. Module interfaces have been defined,
but implementations are pending.

## Author

- **User**: gol2em
- **Email**: 664862601@qq.com

## License

TBD

