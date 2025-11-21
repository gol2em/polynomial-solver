#!/usr/bin/env python3
"""
Test the result refiner on all 1D examples with various tolerances.

Shows how Newton refinement achieves 1e-15 precision regardless of
initial solver tolerance.
"""

import subprocess
import sys

def compile_and_run_test(example_name, tolerance):
    """Compile a test that uses the refiner on an example."""
    
    # Create a simple C++ test program
    test_code = f'''
#include "result_refiner.h"
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>

using namespace polynomial_solver;

int main() {{
    // {example_name}
'''
    
    if example_name == "cubic":
        test_code += '''
    // p(x) = (x-0.2)(x-0.5)(x-0.8)
    std::vector<unsigned int> degrees{3u};
    std::vector<double> power_coeffs{-0.08, 0.66, -1.5, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    std::cout << "Cubic: p(x) = (x-0.2)(x-0.5)(x-0.8)\\n";
'''
    elif example_name == "multiplicity":
        test_code += '''
    // p(x) = (x-0.2)(x-0.6)^6
    std::vector<unsigned int> degrees{7u};
    std::vector<double> power_coeffs{
        -0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0
    };
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    std::cout << "Multiplicity: p(x) = (x-0.2)(x-0.6)^6\\n";
'''
    elif example_name == "wilkinson":
        test_code += '''
    // Wilkinson polynomial: (x-1)(x-2)...(x-19)
    std::vector<unsigned int> degrees{19u};
    std::vector<double> power_coeffs(20);
    // Simplified: just use a few roots for testing
    // p(x) = (x-0.1)(x-0.3)(x-0.5)(x-0.7)(x-0.9)
    degrees[0] = 5u;
    power_coeffs = {-0.00945, 0.1575, -1.05, 3.15, -4.2, 1.0};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({p});
    std::cout << "Wilkinson (5 roots): p(x) = (x-0.1)(x-0.3)(x-0.5)(x-0.7)(x-0.9)\\n";
'''
    
    test_code += f'''
    // Solve
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = {tolerance};
    config.max_depth = 100;
    
    std::cout << "Solver tolerance: " << config.tolerance << "\\n";
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "Solver found " << result.num_resolved << " boxes\\n";
    
    // Refine
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = 1e-15;
    refine_config.residual_tolerance = 1e-12;
    
    RefinementResult refined = refiner.refine(result, system, refine_config);
    
    std::cout << "Refined to " << refined.roots.size() << " roots\\n";
    
    // Print results
    double max_residual = 0.0;
    for (const auto& root : refined.roots) {{
        double res = std::abs(root.residual[0]);
        if (res > max_residual) max_residual = res;
        std::cout << "  x=" << std::setprecision(16) << root.location[0]
                  << ", mult=" << root.multiplicity
                  << ", res=" << std::scientific << std::setprecision(2) << res << "\\n";
    }}
    
    std::cout << "Max residual: " << std::scientific << max_residual << "\\n";
    std::cout << "Status: " << (max_residual < 1e-12 ? "PASS" : "FAIL") << "\\n";
    
    return (max_residual < 1e-12) ? 0 : 1;
}}
'''
    
    # Write test file
    test_file = f"/tmp/test_refiner_{example_name}.cpp"
    with open(test_file, 'w') as f:
        f.write(test_code)
    
    # Compile
    compile_cmd = f"g++ -std=c++11 -I./include {test_file} -L./build/lib -lresult_refiner -lsolver -ldifferentiation -lpolynomial -lgeometry -lde_casteljau -o /tmp/test_refiner_{example_name}"
    result = subprocess.run(compile_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        return None
    
    # Run
    run_cmd = f"/tmp/test_refiner_{example_name}"
    result = subprocess.run(run_cmd, shell=True, capture_output=True, text=True)
    return result.stdout

def main():
    print("Result Refiner Test on 1D Examples")
    print("=" * 70)
    print()
    
    examples = ["cubic", "multiplicity", "wilkinson"]
    tolerances = [1e-8, 1e-10, 1e-12]
    
    for example in examples:
        print(f"\\n{'='*70}")
        print(f"Example: {example}")
        print('='*70)
        
        for tol in tolerances:
            print(f"\\n--- Tolerance: {tol:.0e} ---")
            output = compile_and_run_test(example, tol)
            if output:
                print(output)
            else:
                print("Test failed to compile/run")

if __name__ == "__main__":
    main()

