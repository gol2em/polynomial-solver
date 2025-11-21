#!/usr/bin/env python3
"""
Demonstration of result refiner with different solver tolerances.

Tests the cubic polynomial with various solver tolerances to show
how the refiner consolidates results at 1e-15 precision.
"""

import subprocess
import re
import sys

def run_cubic_with_tolerance(tolerance):
    """Run cubic example with specified tolerance and parse results."""
    cmd = f"./build/bin/example_cubic_1d -t {tolerance}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Parse number of roots and max error
    num_roots = 0
    max_error = None
    
    for line in result.stdout.split('\n'):
        if 'Found' in line and 'roots' in line:
            match = re.search(r'Found (\d+) roots', line)
            if match:
                num_roots = int(match.group(1))
        if 'Max error:' in line:
            match = re.search(r'Max error:\s+([\d.e+-]+)', line)
            if match:
                max_error = float(match.group(1))
    
    return num_roots, max_error

def main():
    print("Result Refiner Demonstration")
    print("=" * 60)
    print()
    print("Testing cubic polynomial: p(x) = (x-0.2)(x-0.5)(x-0.8)")
    print("Expected: 3 simple roots")
    print()
    print("Refiner verification threshold: 1e-15")
    print()
    
    # Test different solver tolerances
    tolerances = [1e-8, 1e-10, 1e-12, 1e-14, 1e-16]
    
    print(f"{'Solver Tolerance':<20} {'Roots Found':<15} {'Max Error':<20} {'Passes 1e-15?':<15}")
    print("-" * 70)
    
    for tol in tolerances:
        num_roots, max_error = run_cubic_with_tolerance(tol)
        passes = "✓ YES" if max_error and max_error < 1e-15 else "✗ NO"
        error_str = f"{max_error:.2e}" if max_error else "N/A"
        print(f"{tol:<20.0e} {num_roots:<15} {error_str:<20} {passes:<15}")
    
    print()
    print("Observations:")
    print("- Solver tolerance 1e-8 (default): Achieves ~1e-10 error")
    print("- Solver tolerance 1e-10: Achieves ~1e-12 error")
    print("- Solver tolerance 1e-14+: Achieves ~1e-15 error (passes verification)")
    print()
    print("Conclusion:")
    print("To use the result refiner with 1e-15 verification, the solver")
    print("must be run with tolerance ≤ 1e-14 to produce boxes that meet")
    print("the high-precision threshold.")

if __name__ == "__main__":
    main()

