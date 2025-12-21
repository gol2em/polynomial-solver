# Polynomial Solver

A C++11 library for solving systems of multivariate polynomial equations using subdivision methods with Bernstein basis representation.

## Prerequisites

### Required
- C++11 compiler (GCC 4.8+, Clang 3.4+)
- CMake 3.15+
- Make

### Optional (High-Precision)
- Boost (header-only multiprecision wrapper)
- MPFR (arbitrary-precision floating-point)
- GMP (multi-precision integers, MPFR dependency)

```bash
# Check prerequisites
./check_prerequisites.sh

# Install high-precision libraries (Ubuntu/Debian)
sudo apt-get install libboost-dev libmpfr-dev libgmp-dev
```

## Quick Start

### Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Library Targets

After building, link against:
- **Library**: `build/lib/libpolynomial_solver.a`
- **Include**: `include/`

```cmake
target_include_directories(your_target PRIVATE /path/to/polynomial-solver/include)
target_link_libraries(your_target /path/to/polynomial-solver/build/lib/libpolynomial_solver.a)
```

### Basic Usage

```cpp
#include <polynomial_solver.h>
using namespace polynomial_solver;

// 1. Define polynomial: (x-0.2)(x-0.5)(x-0.8)
std::vector<unsigned int> degrees = {3};
std::vector<double> power_coeffs = {-0.08, 0.66, -1.5, 1.0};
Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
PolynomialSystem system({poly});

// 2. Solve
Solver solver;
SubdivisionConfig config;
config.tolerance = 1e-8;
auto result = solver.subdivisionSolve(system, config, RootBoundingMethod::ProjectedPolyhedral);

// 3. Refine
ResultRefiner refiner;
RefinementConfig refine_config;
refine_config.residual_tolerance = 1e-15;
auto refined = refiner.refine(result, system, refine_config);
```

### Run Examples

```bash
./build/bin/example_simple_cubic          # 1D cubic
./build/bin/example_wilkinson_1d          # Ill-conditioned
./build/bin/example_multiplicity_1d       # Multiple roots
./build/bin/example_circle_ellipse        # 2D system
```

Use `--help` for command-line options. See [docs/PARAMETERS.md](docs/PARAMETERS.md) for parameter details.

## High-Precision Build (Optional)

```bash
cmake .. -DENABLE_HIGH_PRECISION=ON
make -j$(nproc)
```

See [docs/HIGH_PRECISION.md](docs/HIGH_PRECISION.md) for details.

## Development

```bash
# Playground - rapid prototyping (has Makefile)
cd playground && make mytest && ./mytest

# Examples and tests are built via CMake
cd build
make examples           # Build all examples
make test               # Run all tests (via CTest)
ctest -V                # Run tests with verbose output
```

## Project Structure

```
polynomial-solver/
├── include/                      # Header files
│   ├── polynomial_solver.h       # Unified header (include this!)
│   ├── polynomial.h              # Multivariate polynomial class
│   ├── solver.h                  # Main solver interface
│   ├── result_refiner.h          # High-precision refinement
│   ├── differentiation.h         # Polynomial differentiation
│   ├── geometry.h                # Geometric algorithms
│   └── de_casteljau.h            # De Casteljau algorithm
├── src/                          # Implementation files
│   ├── polynomial.cpp
│   ├── solver.cpp
│   ├── result_refiner.cpp
│   ├── differentiation.cpp
│   ├── geometry.cpp
│   └── de_casteljau.cpp
├── build/lib/                    # Build output
│   └── libpolynomial_solver.a    # Unified library (link this!)
├── playground/                   # Quick testing (has Makefile)
│   ├── Makefile                  # Compile any .cpp file easily
│   └── README.md                 # Playground documentation
├── examples/                     # Example programs (CMake)
│   ├── simple_cubic.cpp          # 2-line workflow demo
│   ├── multiplicity_1d_roots.cpp # Multiple roots
│   ├── wilkinson_1d_roots.cpp    # Ill-conditioned
│   └── circle_ellipse_intersection.cpp  # 2D system
├── tests/                        # Test suite (CMake + CTest)
│   └── test_*.cpp                # Test files
├── tools/                        # Tools and utilities
│   └── refine_from_dumps.cpp     # Root refinement tool
├── docs/                         # Documentation
│   ├── ALGORITHMS.md             # Algorithm overview
│   ├── VISUALIZATION.md          # Visualization guide
│   ├── PARAMETERS.md             # Parameter reference
│   └── ...                       # More technical docs
├── check_prerequisites.sh        # Prerequisite checker
├── build.sh                      # Automated build script
└── README.md                     # This file
```

## Documentation

- [docs/ALGORITHMS.md](docs/ALGORITHMS.md) - Algorithm overview
- [docs/PARAMETERS.md](docs/PARAMETERS.md) - Parameter reference
- [docs/VISUALIZATION.md](docs/VISUALIZATION.md) - Visualization and geometry dump
- [docs/HIGH_PRECISION.md](docs/HIGH_PRECISION.md) - High-precision arithmetic
- [docs/CONDITIONING_AND_PRECISION.md](docs/CONDITIONING_AND_PRECISION.md) - Condition numbers
- [docs/GEOMETRY_ALGORITHMS.md](docs/GEOMETRY_ALGORITHMS.md) - Geometric algorithms
- [docs/DEGENERATE_BOXES.md](docs/DEGENERATE_BOXES.md) - Degeneracy handling
- [examples/README.md](examples/README.md) - Example programs

## Tests

```bash
cd build && ctest
```

## License

TBD
