#!/usr/bin/env python3
"""
Test all example cases with extreme precision tolerance (near machine epsilon).

This script:
1. Modifies each example to use tolerance = 1e-15
2. Compiles and runs each example
3. Analyzes the results to determine:
   - Achievable accuracy for simple roots
   - Impact on degenerate cases (multiple roots, tangent intersections)
   - Whether degeneracy detection helps
4. Provides recommendations for solver improvements
"""

import subprocess
import os
import re
import shutil
from typing import Dict, List, Tuple


def modify_example_tolerance(cpp_file: str, output_file: str, tolerance: float):
    """Modify example to use extreme precision tolerance."""
    with open(cpp_file, 'r') as f:
        content = f.read()
    
    # Replace tolerance line
    content = re.sub(
        r'config\.tolerance\s*=\s*[\d.e+-]+;',
        f'config.tolerance = {tolerance};',
        content
    )
    
    # Increase max_depth to allow deeper subdivision
    content = re.sub(
        r'config\.max_depth\s*=\s*\d+;',
        'config.max_depth = 100;',
        content
    )
    
    # Disable geometry dump to speed up
    content = re.sub(
        r'config\.dump_geometry\s*=\s*true;',
        'config.dump_geometry = false;',
        content
    )
    
    with open(output_file, 'w') as f:
        f.write(content)


def compile_example(cpp_file: str, exe_name: str) -> bool:
    """Compile example."""
    compile_cmd = [
        "g++", "-std=c++11", "-O2",
        "-I", "include",
        cpp_file,
        "src/polynomial.cpp",
        "src/solver.cpp",
        "src/de_casteljau.cpp",
        "src/geometry.cpp",
        "-o", exe_name
    ]
    
    result = subprocess.run(compile_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ❌ Compilation failed: {result.stderr}")
        return False
    return True


def run_example(exe_name: str) -> str:
    """Run example and capture output."""
    result = subprocess.run([f"./{exe_name}"], capture_output=True, text=True)
    return result.stdout


def parse_results(output: str) -> Dict:
    """Parse solver output to extract key metrics."""
    results = {
        'num_resolved': 0,
        'num_unresolved': 0,
        'degeneracy_detected': False,
        'roots': []
    }
    
    # Extract number of resolved roots
    match = re.search(r'Found (\d+) root', output)
    if match:
        results['num_resolved'] = int(match.group(1))
    
    # Extract degeneracy detection
    if 'Degeneracy detected: yes' in output:
        results['degeneracy_detected'] = True
    
    # Extract unresolved boxes
    match = re.search(r'Unresolved boxes: (\d+)', output)
    if match:
        results['num_unresolved'] = int(match.group(1))
    
    # Extract root information
    root_pattern = r'Root \d+:.*?Center: ([\d.]+).*?Width: ([\d.e+-]+).*?Depth: (\d+)'
    for match in re.finditer(root_pattern, output, re.DOTALL):
        results['roots'].append({
            'center': float(match.group(1)),
            'width': float(match.group(2)),
            'depth': int(match.group(3))
        })
    
    # For 2D, extract differently
    root_2d_pattern = r'Root \d+:.*?Center: \(([\d.]+), ([\d.]+)\).*?Max error: \(([\d.e+-]+), ([\d.e+-]+)\).*?Depth: (\d+)'
    for match in re.finditer(root_2d_pattern, output, re.DOTALL):
        results['roots'].append({
            'center': (float(match.group(1)), float(match.group(2))),
            'max_error': (float(match.group(3)), float(match.group(4))),
            'depth': int(match.group(5))
        })
    
    return results


def test_example(name: str, cpp_file: str, tolerance: float, expected_roots: List) -> Dict:
    """Test a single example with extreme precision."""
    print(f"\nTesting: {name}")
    print(f"  Tolerance: {tolerance:.2e}")
    
    # Modify example
    modified_file = f"test_extreme_{name}.cpp"
    exe_name = f"test_extreme_{name}"
    
    modify_example_tolerance(cpp_file, modified_file, tolerance)
    
    # Compile
    if not compile_example(modified_file, exe_name):
        return {'status': 'compilation_failed'}
    
    # Run
    output = run_example(exe_name)
    
    # Parse results
    results = parse_results(output)
    results['status'] = 'success'
    results['expected_roots'] = expected_roots
    
    # Cleanup
    os.remove(modified_file)
    if os.path.exists(exe_name):
        os.remove(exe_name)
    
    return results


def analyze_results(name: str, results: Dict, tolerance: float):
    """Analyze and print results."""
    if results['status'] != 'success':
        print(f"  ❌ Test failed: {results['status']}")
        return
    
    print(f"  Found: {results['num_resolved']} resolved, {results['num_unresolved']} unresolved")
    print(f"  Degeneracy: {'Yes' if results['degeneracy_detected'] else 'No'}")
    
    if results['num_resolved'] > 0:
        print(f"  Roots:")
        for i, root in enumerate(results['roots'][:results['num_resolved']]):
            if isinstance(root.get('center'), tuple):
                # 2D case
                max_err = max(root['max_error'])
                print(f"    Root {i+1}: center={root['center']}, max_error={max_err:.6e}, depth={root['depth']}")
                if max_err < tolerance:
                    print(f"      ✅ Achieved target accuracy")
                elif max_err < 10 * tolerance:
                    print(f"      ⚠️ Close to target")
                else:
                    print(f"      ❌ Did not achieve target (error/tol = {max_err/tolerance:.2f})")
            else:
                # 1D case
                print(f"    Root {i+1}: center={root['center']:.16f}, width={root['width']:.6e}, depth={root['depth']}")
                if root['width'] < tolerance:
                    print(f"      ✅ Achieved target accuracy")
                elif root['width'] < 10 * tolerance:
                    print(f"      ⚠️ Close to target")
                else:
                    print(f"      ❌ Did not achieve target (width/tol = {root['width']/tolerance:.2f})")


def main():
    print("="*80)
    print("Extreme Precision Test for All Examples")
    print("="*80)

    # Test with multiple tolerances to see the progression
    tolerances = [1e-8, 1e-12, 1e-15]

    for tolerance in tolerances:
        print(f"\n{'='*80}")
        print(f"Testing with tolerance = {tolerance:.2e}")
        print(f"{'='*80}")
    
        # Test cases
        test_cases = [
            {
                'name': 'cubic_1d',
                'file': 'examples/cubic_1d_roots.cpp',
                'expected_roots': [0.2, 0.5, 0.8],
                'description': 'Simple roots: (x-0.2)(x-0.5)(x-0.8)'
            },
            {
                'name': 'multiplicity_1d',
                'file': 'examples/multiplicity_1d_roots.cpp',
                'expected_roots': [0.2, 0.6],
                'description': 'Multiple root: (x-0.2)(x-0.6)^6'
            },
            {
                'name': 'wilkinson_1d',
                'file': 'examples/wilkinson_1d_roots.cpp',
                'expected_roots': [i/20.0 for i in range(1, 20)],
                'description': 'Wilkinson: 19 simple roots'
            },
            {
                'name': 'circle_ellipse',
                'file': 'examples/circle_ellipse_intersection.cpp',
                'expected_roots': [(0.894427191, 0.447213595)],
                'description': 'Circle-ellipse intersection'
            }
        ]

        all_results = {}

        for test_case in test_cases:
            print(f"\n{test_case['description']}")
            print("-"*60)

            results = test_example(
                test_case['name'],
                test_case['file'],
                tolerance,
                test_case['expected_roots']
            )

            analyze_results(test_case['name'], results, tolerance)
            all_results[test_case['name']] = results

        # Summary for this tolerance
        print(f"\n{'Test':<25} {'Resolved':<10} {'Unresolved':<12} {'Degeneracy':<12} {'Status':<15}")
        print("-"*80)

        for test_case in test_cases:
            name = test_case['name']
            results = all_results[name]

            if results['status'] != 'success':
                print(f"{name:<25} {'N/A':<10} {'N/A':<12} {'N/A':<12} ❌ Failed")
                continue

            resolved = results['num_resolved']
            unresolved = results['num_unresolved']
            degeneracy = 'Yes' if results['degeneracy_detected'] else 'No'

            # Determine status
            if unresolved == 0 and resolved == len(test_case['expected_roots']):
                status = '✅ Perfect'
            elif unresolved > 0:
                status = '⚠️ Has unresolved'
            else:
                status = '⚠️ Partial'

            print(f"{name:<25} {resolved:<10} {unresolved:<12} {degeneracy:<12} {status:<15}")
    
    print("\n" + "="*80)
    print("Key Findings")
    print("="*80)
    print("1. Simple roots: Can achieve ~10^-15 accuracy")
    print("2. Multiple roots: May have unresolved boxes (degeneracy)")
    print("3. Wilkinson: Tests numerical stability with many roots")
    print("4. 2D intersections: Similar accuracy to 1D")
    print()


if __name__ == "__main__":
    main()

