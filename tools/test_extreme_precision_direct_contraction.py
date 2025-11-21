#!/usr/bin/env python3
"""
Test direct contraction implementation with extreme precisions on 1D examples.

This script tests various tolerance levels to find the optimal precision
that maximizes root resolution while maintaining numerical stability.
"""

import subprocess
import re
import sys
import json
from pathlib import Path

def build_test_executable(tolerance):
    """Build a custom test executable with specific tolerance."""
    # Create a temporary test file
    test_code = f"""
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace polynomial_solver;

int main() {{
    // Test 1: Cubic (x-0.2)(x-0.5)(x-0.8)
    std::vector<unsigned int> degrees1{{3u}};
    std::vector<double> power_coeffs1{{-0.08, 0.66, -1.5, 1.0}};
    Polynomial p1 = Polynomial::fromPower(degrees1, power_coeffs1);
    PolynomialSystem system1({{p1}});
    
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = {tolerance};
    config.max_depth = 100;
    config.contraction_threshold = 0.8;
    config.strategy = SubdivisionStrategy::SubdivideFirst;
    
    SubdivisionSolverResult result1 = solver.subdivisionSolve(
        system1, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "CUBIC:" << result1.num_resolved << ":" 
              << (result1.boxes.size() - result1.num_resolved) << ":"
              << (result1.degeneracy_detected ? "yes" : "no");
    
    if (result1.num_resolved > 0) {{
        double max_width = 0.0;
        for (size_t i = 0; i < result1.num_resolved; ++i) {{
            double width = result1.boxes[i].upper[0] - result1.boxes[i].lower[0];
            if (width > max_width) max_width = width;
        }}
        std::cout << ":" << std::scientific << std::setprecision(3) << max_width;
    }} else {{
        std::cout << ":N/A";
    }}
    std::cout << std::endl;
    
    // Test 2: Multiplicity (x-0.2)(x-0.6)^6
    std::vector<unsigned int> degrees2{{7u}};
    std::vector<double> power_coeffs2{{
        -0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0
    }};
    Polynomial p2 = Polynomial::fromPower(degrees2, power_coeffs2);
    PolynomialSystem system2({{p2}});
    
    SubdivisionSolverResult result2 = solver.subdivisionSolve(
        system2, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "MULTIPLICITY:" << result2.num_resolved << ":" 
              << (result2.boxes.size() - result2.num_resolved) << ":"
              << (result2.degeneracy_detected ? "yes" : "no") << std::endl;
    
    return 0;
}}
"""
    
    # Write test file
    test_file = Path("tools/temp_precision_test.cpp")
    test_file.write_text(test_code)
    
    # Compile
    cmd = "g++ -std=c++11 -O2 -I./include tools/temp_precision_test.cpp " \
          "build/lib/libsolver.a build/lib/libpolynomial.a " \
          "build/lib/libgeometry.a build/lib/libde_casteljau.a " \
          "-o tools/temp_precision_test"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        return False
    return True

def run_precision_test(tolerance):
    """Run test with specific tolerance."""
    if not build_test_executable(tolerance):
        return None
    
    result = subprocess.run("./tools/temp_precision_test", 
                          shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        return None
    
    lines = result.stdout.strip().split('\n')
    results = {}
    
    for line in lines:
        if line.startswith("CUBIC:"):
            parts = line.split(':')
            results['cubic'] = {
                'resolved': int(parts[1]),
                'unresolved': int(parts[2]),
                'degeneracy': parts[3],
                'max_width': parts[4] if len(parts) > 4 else 'N/A'
            }
        elif line.startswith("MULTIPLICITY:"):
            parts = line.split(':')
            results['multiplicity'] = {
                'resolved': int(parts[1]),
                'unresolved': int(parts[2]),
                'degeneracy': parts[3]
            }
    
    return results

def main():
    print("="*80)
    print("Extreme Precision Test: Direct Contraction Implementation")
    print("="*80)
    print()
    print("Testing 1D examples with various tolerance levels to find optimal precision")
    print()
    
    # Test tolerance levels from 10^-8 to 10^-16
    tolerances = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16]
    
    print(f"{'Tolerance':<12} {'Cubic':<25} {'Multiplicity':<25}")
    print(f"{'':12} {'Resolved/Unres/Width':<25} {'Resolved/Unres':<25}")
    print("-"*80)
    
    results_log = []
    
    for tol in tolerances:
        print(f"{tol:<12.0e} ", end='', flush=True)
        
        results = run_precision_test(tol)
        
        if results is None:
            print("FAILED")
            continue
        
        # Cubic results
        cubic = results.get('cubic', {})
        cubic_str = f"{cubic.get('resolved', 0)}/{cubic.get('unresolved', 0)}"
        if cubic.get('max_width') != 'N/A':
            cubic_str += f"/{cubic.get('max_width')}"
        print(f"{cubic_str:<25} ", end='')
        
        # Multiplicity results
        mult = results.get('multiplicity', {})
        mult_str = f"{mult.get('resolved', 0)}/{mult.get('unresolved', 0)}"
        print(f"{mult_str:<25}")
        
        results_log.append({
            'tolerance': tol,
            'cubic': cubic,
            'multiplicity': mult
        })
    
    # Cleanup
    Path("tools/temp_precision_test.cpp").unlink(missing_ok=True)
    Path("tools/temp_precision_test").unlink(missing_ok=True)
    
    print()
    print("="*80)
    print("Analysis")
    print("="*80)
    
    # Find best tolerance for cubic (all roots resolved, smallest width)
    best_cubic = None
    for r in results_log:
        cubic = r['cubic']
        if cubic.get('resolved', 0) >= 3:  # All 3 roots
            if best_cubic is None:
                best_cubic = r
            elif cubic.get('max_width', 'N/A') != 'N/A':
                if best_cubic['cubic'].get('max_width', 'N/A') == 'N/A':
                    best_cubic = r
    
    if best_cubic:
        print(f"\n✅ Best for Cubic (simple roots):")
        print(f"   Tolerance: {best_cubic['tolerance']:.0e}")
        print(f"   Resolved: {best_cubic['cubic']['resolved']}/3")
        print(f"   Max width: {best_cubic['cubic'].get('max_width', 'N/A')}")
    
    # Find best tolerance for multiplicity (simple root resolved, fewest unresolved)
    best_mult = None
    for r in results_log:
        mult = r['multiplicity']
        if mult.get('resolved', 0) >= 1:  # Simple root resolved
            if best_mult is None or mult.get('unresolved', 999) < best_mult['multiplicity'].get('unresolved', 999):
                best_mult = r
    
    if best_mult:
        print(f"\n✅ Best for Multiplicity (degenerate case):")
        print(f"   Tolerance: {best_mult['tolerance']:.0e}")
        print(f"   Resolved: {best_mult['multiplicity']['resolved']}/2")
        print(f"   Unresolved: {best_mult['multiplicity']['unresolved']}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

