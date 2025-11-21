#!/usr/bin/env python3
"""
Test Wilkinson polynomial with extreme precisions using direct contraction.

The Wilkinson polynomial is notoriously ill-conditioned, so this tests
the limits of direct contraction accuracy.
"""

import subprocess
import sys
from pathlib import Path

def create_wilkinson_test(tolerance):
    """Create C++ test for Wilkinson polynomial."""
    return f"""
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace polynomial_solver;

int main() {{
    // Wilkinson polynomial: (x-1)(x-2)...(x-19)
    // Coefficients computed using Python
    std::vector<unsigned int> degrees{{19u}};
    std::vector<double> power_coeffs{{
        1.216451004088320e+17,
        -6.040978405331200e+17,
        1.337731838289920e+18,
        -1.732776261120000e+18,
        1.428329123020800e+18,
        -7.890371113950720e+17,
        3.056447448985600e+17,
        -8.535575456000000e+16,
        1.723949560000000e+16,
        -2.573775360000000e+15,
        2.849834880000000e+14,
        -2.345852800000000e+13,
        1.414014720000000e+12,
        -6.187239200e+10,
        1.931559200e+09,
        -4.151347200e+07,
        5.749760000e+05,
        -4.845000000e+03,
        1.900000000e+01,
        -1.0
    }};
    
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    // Restrict to domain [1, 19] (where the roots are)
    std::vector<double> lower{{1.0}};
    std::vector<double> upper{{19.0}};
    p = p.restrict(lower, upper);

    PolynomialSystem system({{p}});

    Solver solver;
    SubdivisionConfig config;
    config.tolerance = {tolerance};
    config.max_depth = 100;
    config.contraction_threshold = 0.8;
    config.strategy = SubdivisionStrategy::SubdivideFirst;

    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);

    // Map results back to [1, 19] domain
    for (size_t i = 0; i < result.boxes.size(); ++i) {{
        result.boxes[i].lower[0] = 1.0 + result.boxes[i].lower[0] * 18.0;
        result.boxes[i].upper[0] = 1.0 + result.boxes[i].upper[0] * 18.0;
        result.boxes[i].center[0] = 1.0 + result.boxes[i].center[0] * 18.0;
        result.boxes[i].max_error[0] = result.boxes[i].max_error[0] * 18.0;
    }}
    
    std::cout << "STATS:" << result.num_resolved << ":" 
              << (result.boxes.size() - result.num_resolved) << ":"
              << (result.degeneracy_detected ? "yes" : "no") << std::endl;
    
    // Analyze resolved roots
    int roots_in_range[19] = {{0}};  // Count roots near each integer 1-19
    
    for (size_t i = 0; i < result.num_resolved; ++i) {{
        double center = result.boxes[i].center[0];
        double width = result.boxes[i].upper[0] - result.boxes[i].lower[0];
        
        // Find closest integer root
        int closest = (int)(center + 0.5);
        if (closest >= 1 && closest <= 19) {{
            double error = std::abs(center - closest);
            roots_in_range[closest-1]++;
            
            std::cout << "ROOT:" << i << ":" << closest << ":"
                      << std::scientific << std::setprecision(16) 
                      << center << ":" << error << ":" << width << ":"
                      << result.boxes[i].depth << std::endl;
        }}
    }}
    
    // Count unique roots found
    int unique_roots = 0;
    for (int i = 0; i < 19; ++i) {{
        if (roots_in_range[i] > 0) unique_roots++;
    }}
    std::cout << "UNIQUE:" << unique_roots << std::endl;
    
    return 0;
}}
"""

def run_wilkinson_test(tolerance):
    """Run Wilkinson test with specific tolerance."""
    # Write test file
    test_file = Path("tools/temp_wilkinson_test.cpp")
    test_file.write_text(create_wilkinson_test(tolerance))
    
    # Compile
    cmd = "g++ -std=c++11 -O2 -I./include tools/temp_wilkinson_test.cpp " \
          "build/lib/libsolver.a build/lib/libpolynomial.a " \
          "build/lib/libgeometry.a build/lib/libde_casteljau.a " \
          "-o tools/temp_wilkinson_test 2>&1"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return None
    
    # Run test
    result = subprocess.run("./tools/temp_wilkinson_test", 
                          shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        return None
    
    # Parse results
    lines = result.stdout.strip().split('\n')
    stats = None
    roots = []
    unique = 0
    
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
                'true_root': int(parts[2]),
                'center': float(parts[3]),
                'error': float(parts[4]),
                'width': float(parts[5]),
                'depth': int(parts[6])
            })
        elif line.startswith("UNIQUE:"):
            unique = int(line.split(':')[1])
    
    return {'stats': stats, 'roots': roots, 'unique': unique}

def main():
    print("="*95)
    print("Wilkinson Polynomial: Extreme Precision Test with Direct Contraction")
    print("="*95)
    print()
    print("Testing: (x-1)(x-2)(x-3)...(x-19) with 19 integer roots")
    print("This is a notoriously ill-conditioned problem.")
    print()
    
    # Test tolerance levels
    tolerances = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15]
    
    print(f"{'Tolerance':<12} {'Resolved':<10} {'Unique':<8} {'Unresolved':<12} {'Max Error':<15} {'Avg Depth':<10}")
    print("-"*95)
    
    all_results = []
    
    for tol in tolerances:
        result = run_wilkinson_test(tol)
        
        if result is None or result['stats'] is None:
            print(f"{tol:<12.0e} FAILED")
            continue
        
        stats = result['stats']
        roots = result['roots']
        unique = result['unique']
        
        if len(roots) == 0:
            print(f"{tol:<12.0e} {stats['resolved']:<10} {unique:<8} {stats['unresolved']:<12} {'N/A':<15} {'N/A':<10}")
        else:
            max_error = max(r['error'] for r in roots)
            avg_depth = sum(r['depth'] for r in roots) / len(roots)
            print(f"{tol:<12.0e} {stats['resolved']:<10} {unique:<8} {stats['unresolved']:<12} {max_error:<15.3e} {avg_depth:<10.1f}")
        
        all_results.append({
            'tolerance': tol,
            'stats': stats,
            'unique': unique,
            'roots': roots
        })
    
    # Cleanup
    Path("tools/temp_wilkinson_test.cpp").unlink(missing_ok=True)
    Path("tools/temp_wilkinson_test").unlink(missing_ok=True)
    
    print()
    print("="*95)
    print("Analysis")
    print("="*95)
    
    # Find best result (most unique roots resolved)
    best = max(all_results, key=lambda x: x['unique'])
    
    print(f"\n✅ BEST RESULT:")
    print(f"   Tolerance:     {best['tolerance']:.0e}")
    print(f"   Unique roots:  {best['unique']}/19 ({100*best['unique']/19:.1f}%)")
    print(f"   Total boxes:   {best['stats']['resolved']}")
    print(f"   Unresolved:    {best['stats']['unresolved']}")
    
    if len(best['roots']) > 0:
        max_error = max(r['error'] for r in best['roots'])
        avg_error = sum(r['error'] for r in best['roots']) / len(best['roots'])
        print(f"   Max error:     {max_error:.3e}")
        print(f"   Avg error:     {avg_error:.3e}")
        
        # Show which roots were found
        found_roots = sorted(set(r['true_root'] for r in best['roots']))
        print(f"   Roots found:   {found_roots}")
    
    # Compare different tolerances
    print(f"\n📊 TOLERANCE COMPARISON:")
    for r in all_results:
        if r['unique'] > 0:
            print(f"   {r['tolerance']:.0e}: {r['unique']}/19 roots, {r['stats']['unresolved']} unresolved boxes")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

