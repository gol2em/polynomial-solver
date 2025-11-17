# Examples

This directory contains example programs demonstrating the polynomial solver library.

## Available Examples

### circle_ellipse_intersection.cpp

Demonstrates solving a system of polynomial equations representing the intersection
of a circle and an ellipse.

**Problem:**
- f1(x,y) = x² + y² - 1 = 0 (unit circle)
- f2(x,y) = x²/4 + 4y² - 1 = 0 (ellipse)
- Domain: [0, 1] × [0, 1]

**Features demonstrated:**
- Creating polynomials from power basis coefficients
- Building polynomial systems
- Testing all three subdivision strategies
- Enabling geometry dumps for visualization
- Comparing strategy performance

**Expected output:**
- Root at approximately (0.894, 0.447)
- Geometry dumps in `dumps/example_*.txt`
- Convergence with machine epsilon precision

**Build and run:**
```bash
./build.sh
./build/bin/example_circle_ellipse
```

**Visualize results:**
```bash
source .venv/bin/activate
python examples/visualize_circle_ellipse.py
```

### visualize_circle_ellipse.py

Python script for visualizing circle-ellipse intersection results using the visualization API.

**Usage:**
```bash
# Visualize all strategies
python examples/visualize_circle_ellipse.py

# Visualize first 20 iterations
python examples/visualize_circle_ellipse.py --max-steps 20

# Visualize specific strategy
python examples/visualize_circle_ellipse.py --strategy ContractFirst
```

## Building Examples

Examples are built automatically when you run `./build.sh`. To build only examples:

```bash
cd build
cmake ..
make examples
```

To disable building examples:

```bash
cmake -DBUILD_EXAMPLES=OFF ..
make
```

## Adding New Examples

To add a new example:

1. Create a new `.cpp` file in this directory
2. Add it to `examples/CMakeLists.txt`:
   ```cmake
   add_executable(example_myexample myexample.cpp)
   target_link_libraries(example_myexample polynomial solver geometry de_casteljau)
   target_include_directories(example_myexample PRIVATE ${CMAKE_SOURCE_DIR}/include)
   add_dependencies(examples example_myexample)
   ```
3. Rebuild: `./build.sh`

## Example Template

```cpp
#include "polynomial.h"
#include "solver.h"
#include <iostream>

using namespace polynomial_solver;

int main() {
    // 1. Define polynomials
    std::vector<unsigned int> degrees{2u, 2u};
    std::vector<double> coeffs{/* ... */};
    Polynomial p = Polynomial::fromPower(degrees, coeffs);
    
    // 2. Create system
    PolynomialSystem system({p});
    
    // 3. Configure solver
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = 1e-6;
    config.max_depth = 100;
    config.dump_geometry = true;
    config.dump_prefix = "dumps/myexample";
    
    // 4. Solve
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    // 5. Print results
    std::cout << "Found " << result.num_resolved << " roots\n";
    
    return 0;
}
```

## Example Ideas

- **Linear system**: Test machine epsilon precision
- **Cubic curves**: Intersection of two cubic curves
- **Degenerate cases**: Multiple roots, infinite solutions
- **High-degree polynomials**: Performance testing
- **3D problems**: Intersection of surfaces
- **Comparison**: Different root bounding methods

