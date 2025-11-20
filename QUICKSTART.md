# Quick Start Guide

Get up and running with Polynomial Solver in 5 minutes.

## Prerequisites

You need:
- C++ compiler with C++11 support (GCC 4.8+, Clang 3.4+)
- CMake 3.15+
- Make
- Git

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/gol2em/polynomial-solver.git
cd polynomial-solver
```

### Step 2: Check Prerequisites

```bash
./check_prerequisites.sh
```

If any required tools are missing, install them:

**Ubuntu/Debian:**
```bash
sudo apt-get install build-essential cmake git
```

**Fedora/RHEL:**
```bash
sudo dnf install gcc-c++ cmake git
```

**Arch Linux:**
```bash
sudo pacman -S base-devel cmake git
```

### Step 3: Build

```bash
./build.sh --test
```

This will:
1. Create a `build/` directory
2. Configure with CMake
3. Compile the project
4. Run all tests

Expected output:
```
100% tests passed, 0 tests failed out of 12
```

## Running Examples

### Example 1: Circle-Ellipse Intersection (Standard Verification)

Run the standard verification workflow:

```bash
# Run the test with all 3 strategies
./build/bin/test_strategies
```

This generates geometry dumps in `dumps/` directory for visualization.

**Visualize the results:**

```bash
# Set up Python environment (one-time)
uv venv .venv
source .venv/bin/activate
uv pip install numpy matplotlib

# Visualize all strategies
python examples/visualize_circle_ellipse.py

# Or visualize specific strategy with limited steps
python examples/visualize_circle_ellipse.py --strategy ContractFirst --max-steps 10
```

Output: PNG visualizations in `visualizations/viz_*/` showing step-by-step solving process.

### Example 2: Method Comparison

Compare GraphHull vs ProjectedPolyhedral methods:

```bash
./build/bin/test_method_comparison
```

Output shows both methods produce identical results for linear and quadratic systems.

### Example 3: Run All Tests

```bash
cd build
ctest --output-on-failure
```

## Basic Usage

### Solving a 1D Polynomial

Create a file `example.cpp`:

```cpp
#include "polynomial.h"
#include "solver.h"
#include <iostream>

using namespace polynomial_solver;

int main() {
    // Solve: p(x) = x - 0.5
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{-0.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    
    PolynomialSystem system({p});
    
    SubdivisionConfig config;
    config.tolerance = 1e-8;
    
    Solver solver;
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "Found " << result.num_resolved << " root(s):" << std::endl;
    for (size_t i = 0; i < result.num_resolved; ++i) {
        std::cout << "  x = " << result.boxes[i].center[0] << std::endl;
    }
    
    return 0;
}
```

Compile and run:

```bash
cd build
g++ -std=c++11 -I../include ../example.cpp -L./lib -lsolver -lpolynomial -lgeometry -lde_casteljau -o example
./example
```

Output:
```
Found 1 root(s):
  x = 0.5
```

### Solving a 2D System

```cpp
// Solve: p1(x,y) = x - 0.5, p2(x,y) = y - 0.3
std::vector<unsigned int> degrees{1u, 1u};

// p1(x,y) = x - 0.5
std::vector<double> power1(4, 0.0);
power1[0] = -0.5;  // constant
power1[2] = 1.0;   // x coefficient
Polynomial p1 = Polynomial::fromPower(degrees, power1);

// p2(x,y) = y - 0.3
std::vector<double> power2(4, 0.0);
power2[0] = -0.3;  // constant
power2[1] = 1.0;   // y coefficient
Polynomial p2 = Polynomial::fromPower(degrees, power2);

PolynomialSystem system({p1, p2});

// Solve...
SubdivisionSolverResult result = solver.subdivisionSolve(
    system, config, RootBoundingMethod::ProjectedPolyhedral);

std::cout << "Root: (" << result.boxes[0].center[0] 
          << ", " << result.boxes[0].center[1] << ")" << std::endl;
```

## Configuration Options

### Solver Configuration

```cpp
SubdivisionConfig config;
config.tolerance = 1e-8;              // Convergence tolerance
config.max_depth = 100;               // Maximum subdivision depth
config.degeneracy_multiplier = 5.0;   // Degeneracy detection threshold
```

### Root Bounding Methods

```cpp
// Method 1: ProjectedPolyhedral (recommended)
RootBoundingMethod::ProjectedPolyhedral

// Method 2: GraphHull (exact for 1D/2D)
RootBoundingMethod::GraphHull

// Method 3: None (uniform subdivision)
RootBoundingMethod::None
```

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Explore [docs/](docs/) for algorithm details
- Check [tests/](tests/) for more examples
- Try [examples/](examples/) for complete workflows
- Use visualization tools in [tools/](tools/)

## Troubleshooting

### Build fails with "C++11 required"

Make sure your compiler supports C++11:
```bash
g++ --version  # Should be 4.8 or higher
```

### CMake version too old

Update CMake:
```bash
# Ubuntu/Debian
sudo apt-get install cmake

# Or download from https://cmake.org/download/
```

### Tests fail

Run with verbose output:
```bash
cd build
ctest --output-on-failure --verbose
```

## Support

For issues or questions:
- Email: wenyd@lsec.cc.ac.cn
- GitHub: https://github.com/gol2em/polynomial-solver

