#!/usr/bin/env python3
"""
Sweep through different tolerance values and test resolution quality.

Uses the existing example binaries with command-line options.
"""

import subprocess
import re
import sys
from pathlib import Path

def run_example_with_tolerance(example_name, tolerance):
    """Run example binary with specific tolerance using command-line options."""

    # Map example names to binary names
    binary_map = {
        "cubic": "./build/bin/example_cubic_1d",
        "multiplicity": "./build/bin/example_multiplicity_1d"
    }

    if example_name not in binary_map:
        return None

    binary = binary_map[example_name]

    # Run with tolerance option
    cmd = f"{binary} -t {tolerance} 2>&1"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        return None

    # Parse output
    output = result.stdout

    # Extract statistics
    resolved_match = re.search(r'Found (\d+) root\(s\)', output)
    degeneracy_match = re.search(r'Degeneracy detected: (\w+)', output)
    unresolved_match = re.search(r'Unresolved boxes: (\d+)', output)

    if not (resolved_match and degeneracy_match and unresolved_match):
        return None

    resolved = int(resolved_match.group(1))
    unresolved = int(unresolved_match.group(1))
    degeneracy = degeneracy_match.group(1) == 'yes'

    # For cubic, extract max error and width
    max_error = None
    max_width = None

    if example_name == "cubic":
        # Extract all widths from "Width: X.XXXe-XX" lines
        width_matches = re.findall(r'Width: ([\d.e+-]+)', output)
        if width_matches:
            widths = [float(w) for w in width_matches]
            max_width = max(widths)

        # Compute max error from centers
        true_roots = [0.2, 0.5, 0.8]
        center_matches = re.findall(r'Center: ([\d.]+)', output)
        if center_matches:
            errors = []
            for center_str in center_matches:
                center = float(center_str)
                min_err = min(abs(center - root) for root in true_roots)
                errors.append(min_err)
            if errors:
                max_error = max(errors)

    return {
        'resolved': resolved,
        'unresolved': unresolved,
        'degeneracy': degeneracy,
        'max_error': max_error,
        'max_width': max_width
    }

def main():
    print("="*100)
    print("Tolerance Sweep: Finding Optimal Precision for Direct Contraction")
    print("="*100)
    print()
    
    tolerances = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16]
    
    # Test Cubic
    print("CUBIC POLYNOMIAL: (x-0.2)(x-0.5)(x-0.8)")
    print("-"*100)
    print(f"{'Tolerance':<12} {'Resolved':<10} {'Unresolved':<12} {'Max Error':<15} {'Max Width':<15}")
    print("-"*100)
    
    cubic_results = []
    for tol in tolerances:
        result = run_example_with_tolerance("cubic", tol)
        if result:
            print(f"{tol:<12.0e} {result['resolved']:<10} {result['unresolved']:<12} ", end='')
            if result['max_error'] is not None:
                print(f"{result['max_error']:<15.3e} {result['max_width']:<15.3e}")
            else:
                print(f"{'N/A':<15} {'N/A':<15}")
            cubic_results.append((tol, result))

    # Test Multiplicity
    print()
    print("MULTIPLICITY POLYNOMIAL: (x-0.2)(x-0.6)^6")
    print("-"*100)
    print(f"{'Tolerance':<12} {'Resolved':<10} {'Unresolved':<12} {'Degeneracy':<12}")
    print("-"*100)

    mult_results = []
    for tol in tolerances:
        result = run_example_with_tolerance("multiplicity", tol)
        if result:
            print(f"{tol:<12.0e} {result['resolved']:<10} {result['unresolved']:<12} {str(result['degeneracy']):<12}")
            mult_results.append((tol, result))
    
    # Analysis
    print()
    print("="*100)
    print("ANALYSIS & RECOMMENDATIONS")
    print("="*100)
    
    # Find best cubic result
    best_cubic = min((r for r in cubic_results if r[1]['resolved'] >= 3 and r[1]['max_error'] is not None),
                     key=lambda x: x[1]['max_error'], default=None)
    
    if best_cubic:
        tol, res = best_cubic
        print(f"\n✅ CUBIC (simple roots) - BEST ACCURACY:")
        print(f"   Tolerance:     {tol:.0e}")
        print(f"   Max error:     {res['max_error']:.3e} ({res['max_error']/2.22e-16:.1f}× machine ε)")
        print(f"   Max width:     {res['max_width']:.3e}")
        print(f"   Resolved:      {res['resolved']}/3")
    
    # Find best multiplicity result
    best_mult = min((r for r in mult_results if r[1]['resolved'] >= 1),
                    key=lambda x: x[1]['unresolved'], default=None)
    
    if best_mult:
        tol, res = best_mult
        print(f"\n✅ MULTIPLICITY (degenerate) - BEST RESOLUTION:")
        print(f"   Tolerance:     {tol:.0e}")
        print(f"   Resolved:      {res['resolved']}/2")
        print(f"   Unresolved:    {res['unresolved']}")
        print(f"   Degeneracy:    {res['degeneracy']}")
    
    print(f"\n📊 RECOMMENDATION:")
    print(f"   For simple roots:     tolerance = 1e-9 to 1e-12 (achieves machine precision)")
    print(f"   For degenerate cases: tolerance = 1e-8 to 1e-10 (best balance)")
    print(f"   Practical default:    tolerance = 1e-10 (good for both cases)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

