#!/usr/bin/env python3
"""
Detailed accuracy analysis for direct contraction at extreme precisions.

Tests 1D examples and computes actual root errors to find optimal precision.
"""

import subprocess
import re
import sys
import numpy as np
from pathlib import Path

def create_test_code(tolerance):
    """Create C++ test code with specific tolerance."""
    return f"""
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

double evaluate_cubic(double x) {{
    return (x - 0.2) * (x - 0.5) * (x - 0.8);
}}

int main() {{
    // Test: Cubic (x-0.2)(x-0.5)(x-0.8)
    std::vector<unsigned int> degrees{{3u}};
    std::vector<double> power_coeffs{{-0.08, 0.66, -1.5, 1.0}};
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({{p}});
    
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = {tolerance};
    config.max_depth = 100;
    config.contraction_threshold = 0.8;
    config.strategy = SubdivisionStrategy::SubdivideFirst;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "STATS:" << result.num_resolved << ":" 
              << (result.boxes.size() - result.num_resolved) << ":"
              << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    
    // Known roots
    double true_roots[] = {{0.2, 0.5, 0.8}};
    
    // Compute errors for resolved roots
    for (size_t i = 0; i < result.num_resolved && i < 10; ++i) {{
        double center = result.boxes[i].center[0];
        double width = result.boxes[i].upper[0] - result.boxes[i].lower[0];
        double residual = std::abs(evaluate_cubic(center));
        
        // Find closest true root
        double min_error = 1e10;
        for (int j = 0; j < 3; ++j) {{
            double error = std::abs(center - true_roots[j]);
            if (error < min_error) min_error = error;
        }}
        
        std::cout << "ROOT:" << i << ":" 
                  << std::scientific << std::setprecision(16) << center << ":"
                  << width << ":" << min_error << ":" << residual << ":"
                  << result.boxes[i].depth << std::endl;
    }}
    
    return 0;
}}
"""

def run_test(tolerance):
    """Run test with specific tolerance and parse results."""
    # Write test file
    test_file = Path("tools/temp_accuracy_test.cpp")
    test_file.write_text(create_test_code(tolerance))
    
    # Compile
    cmd = "g++ -std=c++11 -O2 -I./include tools/temp_accuracy_test.cpp " \
          "build/lib/libsolver.a build/lib/libpolynomial.a " \
          "build/lib/libgeometry.a build/lib/libde_casteljau.a " \
          "-o tools/temp_accuracy_test 2>&1"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return None
    
    # Run test
    result = subprocess.run("./tools/temp_accuracy_test", 
                          shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        return None
    
    # Parse results
    lines = result.stdout.strip().split('\n')
    stats = None
    roots = []
    
    for line in lines:
        if line.startswith("STATS:"):
            parts = line.split(':')
            stats = {
                'resolved': int(parts[1]),
                'unresolved': int(parts[2]),
                'degeneracy': parts[3]
            }
        elif line.startswith("ROOT:"):
            parts = line.split(':')
            roots.append({
                'index': int(parts[1]),
                'center': float(parts[2]),
                'width': float(parts[3]),
                'error': float(parts[4]),
                'residual': float(parts[5]),
                'depth': int(parts[6])
            })
    
    return {'stats': stats, 'roots': roots}

def main():
    print("="*90)
    print("Detailed Accuracy Analysis: Direct Contraction at Extreme Precisions")
    print("="*90)
    print()
    print("Testing cubic polynomial (x-0.2)(x-0.5)(x-0.8) with true roots at 0.2, 0.5, 0.8")
    print()
    
    # Test tolerance levels
    tolerances = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17]
    
    print(f"{'Tolerance':<12} {'Resolved':<10} {'Max Error':<15} {'Max Width':<15} {'Max Residual':<15} {'Avg Depth':<10}")
    print("-"*90)
    
    all_results = []
    
    for tol in tolerances:
        result = run_test(tol)
        
        if result is None or result['stats'] is None:
            print(f"{tol:<12.0e} FAILED")
            continue
        
        stats = result['stats']
        roots = result['roots']
        
        if len(roots) == 0:
            print(f"{tol:<12.0e} {stats['resolved']:<10} {'N/A':<15} {'N/A':<15} {'N/A':<15} {'N/A':<10}")
            continue
        
        max_error = max(r['error'] for r in roots)
        max_width = max(r['width'] for r in roots)
        max_residual = max(r['residual'] for r in roots)
        avg_depth = sum(r['depth'] for r in roots) / len(roots)
        
        print(f"{tol:<12.0e} {stats['resolved']:<10} {max_error:<15.3e} {max_width:<15.3e} {max_residual:<15.3e} {avg_depth:<10.1f}")
        
        all_results.append({
            'tolerance': tol,
            'stats': stats,
            'max_error': max_error,
            'max_width': max_width,
            'max_residual': max_residual,
            'avg_depth': avg_depth,
            'roots': roots
        })
    
    # Cleanup
    Path("tools/temp_accuracy_test.cpp").unlink(missing_ok=True)
    Path("tools/temp_accuracy_test").unlink(missing_ok=True)
    
    print()
    print("="*90)
    print("Analysis & Recommendations")
    print("="*90)
    
    # Find best results
    best_by_error = min((r for r in all_results if r['stats']['resolved'] >= 3), 
                        key=lambda x: x['max_error'], default=None)
    
    if best_by_error:
        print(f"\n✅ BEST ACCURACY (minimum error):")
        print(f"   Tolerance:     {best_by_error['tolerance']:.0e}")
        print(f"   Max error:     {best_by_error['max_error']:.3e}")
        print(f"   Max width:     {best_by_error['max_width']:.3e}")
        print(f"   Max residual:  {best_by_error['max_residual']:.3e}")
        print(f"   Avg depth:     {best_by_error['avg_depth']:.1f}")
        print(f"   Resolved:      {best_by_error['stats']['resolved']}/3")
    
    # Find practical optimum (good accuracy, reasonable depth)
    practical = None
    for r in all_results:
        if r['stats']['resolved'] >= 3 and r['max_error'] < 1e-14:
            if practical is None or r['avg_depth'] < practical['avg_depth']:
                practical = r
    
    if practical:
        print(f"\n✅ PRACTICAL OPTIMUM (accuracy vs. depth):")
        print(f"   Tolerance:     {practical['tolerance']:.0e}")
        print(f"   Max error:     {practical['max_error']:.3e}")
        print(f"   Max width:     {practical['max_width']:.3e}")
        print(f"   Avg depth:     {practical['avg_depth']:.1f}")
    
    # Machine epsilon comparison
    print(f"\n📊 MACHINE EPSILON COMPARISON:")
    print(f"   Machine ε:     2.22e-16")
    if best_by_error:
        ratio = best_by_error['max_error'] / 2.22e-16
        print(f"   Best error:    {best_by_error['max_error']:.3e} ({ratio:.1f}× machine ε)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

