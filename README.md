# Polynomial Solver

A C++11 research library for solving systems of multivariate polynomial equations using subdivision methods with Bernstein basis representation.

**Status**: Research/experimental. Not verified for production use.

## Overview

This library implements:
- Bernstein subdivision for root isolation
- Newton-based refinement with several iteration methods
- Templated architecture supporting double and arbitrary precision (mpreal)
- Multiplicity estimation using Ostrowski's method

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

### Basic Usage (Original API)

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

### Templated API (Recommended for new code)

```cpp
#include "core/polynomial_base.h"
#include "solver/solver_base.h"
#include "refinement/result_refiner_base.h"
using namespace polynomial_solver;

// Type aliases
using Poly = PolynomialBase<double>;
using Solver = SolverBase<double>;
using Refiner = ResultRefinerBase<double>;

// 1. Define polynomial
Poly poly = Poly::fromPower({-0.08, 0.66, -1.5, 1.0});

// 2. Solve
PolynomialSystemBase<double> system;
system.polynomials.push_back(poly.convertedToBernstein());
SubdivisionConfigBase<double> config;
config.tolerance = 1e-8;
auto result = Solver::solve(system, config);

// 3. Refine
RefinementConfigBase<double> refine_config;
refine_config.target_tolerance = 1e-15;
auto refined = Refiner::refineBatch(result.boxes, poly, refine_config);
```

See [docs/TEMPLATED_SOLVER.md](docs/TEMPLATED_SOLVER.md) for the full templated API guide.

### Run Examples

```bash
# Original API examples
./build/bin/example_simple_cubic          # 1D cubic
./build/bin/example_wilkinson_1d          # Ill-conditioned
./build/bin/example_multiplicity_1d       # Multiple roots
./build/bin/example_circle_ellipse        # 2D system

# Templated API examples
./build/bin/example_simple_cubic_templated      # Basic workflow
./build/bin/example_wilkinson_1d_templated      # 20-root Wilkinson
./build/bin/example_multiplicity_1d_templated   # HP escalation for mult-6
./build/bin/example_circle_ellipse_templated    # 2D Newton refinement
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
├── include/
│   ├── polynomial_solver.h       # Main umbrella header
│   │
│   ├── core/                     # Fundamental types
│   │   ├── polynomial.h          # Original Polynomial class
│   │   ├── polynomial_base.h     # Templated PolynomialBase<Scalar>
│   │   ├── geometry.h            # Convex hull, intersection
│   │   ├── geometry_base.h       # Templated geometry
│   │   ├── de_casteljau.h        # De Casteljau subdivision
│   │   ├── differentiation.h     # Polynomial differentiation
│   │   └── interpolation.h       # Polynomial interpolation
│   │
│   ├── solver/                   # Solving algorithms
│   │   ├── solver.h              # Original Solver class
│   │   ├── solver_base.h         # Templated SolverBase<Scalar>
│   │   ├── bounding_strategy.h   # Root bounding methods
│   │   └── newton_multidim.h     # Multi-dimensional Newton
│   │
│   ├── refinement/               # Root refinement
│   │   ├── result_refiner.h      # Original ResultRefiner
│   │   └── result_refiner_base.h # Templated ResultRefinerBase<Scalar>
│   │
│   └── hp/                       # High-precision support
│       ├── high_precision_types.h
│       ├── polynomial_hp.h
│       ├── result_refiner_hp.h
│       └── precision_conversion.h
│
├── src/                          # Implementation files
├── examples/                     # Example programs
│   ├── simple_cubic.cpp          # Original API
│   ├── simple_cubic_templated.cpp    # Templated API
│   ├── multiplicity_1d_templated.cpp # HP escalation workflow
│   └── ...
├── tests/                        # Test suite (33 tests)
├── tools/                        # Utilities
├── playground/                   # Rapid prototyping
├── docs/                         # Documentation
└── build/lib/libpolynomial_solver.a  # Link this!
```

## Documentation

- [docs/TEMPLATED_SOLVER.md](docs/TEMPLATED_SOLVER.md) - **Templated API guide** (recommended)
- [docs/ALGORITHMS.md](docs/ALGORITHMS.md) - Algorithm overview
- [docs/PARAMETERS.md](docs/PARAMETERS.md) - Parameter reference
- [docs/HIGH_PRECISION.md](docs/HIGH_PRECISION.md) - High-precision arithmetic
- [docs/CONDITIONING_AND_PRECISION.md](docs/CONDITIONING_AND_PRECISION.md) - Condition numbers
- [docs/VISUALIZATION.md](docs/VISUALIZATION.md) - Visualization and geometry dump
- [docs/GEOMETRY_ALGORITHMS.md](docs/GEOMETRY_ALGORITHMS.md) - Geometric algorithms
- [docs/DEGENERATE_BOXES.md](docs/DEGENERATE_BOXES.md) - Degeneracy handling
- [examples/README.md](examples/README.md) - Example programs

## Tests

```bash
cd build && ctest        # Run all 33 tests
ctest -V                 # Verbose output
ctest -R Refiner         # Run tests matching "Refiner"
```

## License

TBD
