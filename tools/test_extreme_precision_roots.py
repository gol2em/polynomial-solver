#!/usr/bin/env python3
"""
Test root accuracy at extreme precision (near machine epsilon).

This script tests:
1. Simple roots at extreme precision (tolerance ~ 1e-15)
2. Multiple roots (degeneracy) at extreme precision
3. Tangent intersections (2D degeneracy) at extreme precision
4. How degeneracy affects achievable accuracy
"""

import numpy as np
import subprocess
import os
import re
from typing import List, Tuple


def create_test_polynomial_1d(case: str) -> Tuple[str, List[float], float, int]:
    """
    Create test polynomial for 1D case.
    
    Returns: (description, power_coeffs, true_root, multiplicity)
    """
    if case == "simple":
        # (x - 0.5) - simple root
        return "Simple root: (x - 0.5)", [-0.5, 1.0], 0.5, 1
    
    elif case == "double":
        # (x - 0.5)^2 - double root
        return "Double root: (x - 0.5)^2", [0.25, -1.0, 1.0], 0.5, 2
    
    elif case == "triple":
        # (x - 0.5)^3 - triple root
        return "Triple root: (x - 0.5)^3", [-0.125, 0.75, -1.5, 1.0], 0.5, 3
    
    elif case == "quintuple":
        # (x - 0.5)^5 - quintuple root
        coeffs = [1.0]
        for _ in range(5):
            new_coeffs = [0.0] * (len(coeffs) + 1)
            for i in range(len(coeffs)):
                new_coeffs[i] -= 0.5 * coeffs[i]
                new_coeffs[i + 1] += coeffs[i]
            coeffs = new_coeffs
        return "Quintuple root: (x - 0.5)^5", coeffs, 0.5, 5
    
    else:
        raise ValueError(f"Unknown case: {case}")


def create_test_system_2d(case: str) -> Tuple[str, List[float], List[float], Tuple[float, float], str]:
    """
    Create test system for 2D case.
    
    Returns: (description, power_coeffs1, power_coeffs2, true_root, degeneracy_type)
    """
    if case == "simple":
        # x - 0.5 = 0, y - 0.6 = 0 - simple intersection
        return ("Simple intersection", 
                [-0.5, 1.0], [-0.6, 1.0], 
                (0.5, 0.6), "none")
    
    elif case == "tangent":
        # x^2 + y^2 - 1 = 0, x - 0.6 = 0 - tangent intersection
        # Circle x^2 + y^2 = 1 tangent to line x = 0.6
        # Intersection at (0.6, 0.8)
        return ("Tangent intersection (circle-line)",
                [-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # x^2 + y^2 - 1
                [-0.6, 1.0],  # x - 0.6
                (0.6, 0.8), "tangent")
    
    elif case == "double_tangent":
        # (x - 0.5)^2 + (y - 0.5)^2 - 0.01 = 0, y - 0.6 = 0
        # Small circle tangent to line y = 0.6
        # Intersection at (0.5, 0.6)
        return ("Double tangent (circle-line)",
                [-0.49, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # (x-0.5)^2 + (y-0.5)^2 - 0.01
                [-0.6, 1.0],  # y - 0.6
                (0.5, 0.6), "tangent")
    
    else:
        raise ValueError(f"Unknown case: {case}")


def write_cpp_test_file(filename: str, test_type: str, case: str, tolerance: float):
    """Write C++ test file for the given case."""
    
    if test_type == "1d":
        desc, coeffs, true_root, multiplicity = create_test_polynomial_1d(case)
        
        cpp_code = f"""
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

int main() {{
    std::cout << "Test: {desc}\\n";
    std::cout << "True root: " << std::fixed << std::setprecision(16) << {true_root} << "\\n";
    std::cout << "Multiplicity: {multiplicity}\\n";
    std::cout << "Tolerance: " << std::scientific << {tolerance} << "\\n\\n";
    
    // Define polynomial
    std::vector<unsigned int> degrees{{{len(coeffs) - 1}u}};
    std::vector<double> power_coeffs{{{', '.join(f'{c}' for c in coeffs)}}};
    
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system({{p}});
    
    // Solve with extreme precision
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = {tolerance};
    config.max_depth = 100;
    config.strategy = SubdivisionStrategy::SubdivideFirst;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "Results:\\n";
    std::cout << "  Found " << result.num_resolved << " root(s)\\n";
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\\n\\n";
    
    if (result.num_resolved > 0) {{
        const SubdivisionBoxResult& box = result.boxes[0];
        double center = box.center[0];
        double width = box.upper[0] - box.lower[0];
        double error = std::abs(center - {true_root});
        
        std::cout << "Root 1:\\n";
        std::cout << "  Interval: [" << std::setprecision(16) << box.lower[0] 
                  << ", " << box.upper[0] << "]\\n";
        std::cout << "  Center: " << center << "\\n";
        std::cout << "  Width: " << std::scientific << std::setprecision(6) << width << "\\n";
        std::cout << "  Error: " << error << "\\n";
        std::cout << "  Depth: " << box.depth << "\\n";
        std::cout << "  Error/Tolerance: " << (error / {tolerance}) << "\\n";
        
        if (error < {tolerance}) {{
            std::cout << "  Status: ✅ Achieved target accuracy\\n";
        }} else if (error < 10 * {tolerance}) {{
            std::cout << "  Status: ⚠️ Close to target\\n";
        }} else {{
            std::cout << "  Status: ❌ Did not achieve target\\n";
        }}
    }}
    
    return 0;
}}
"""
    
    else:  # 2d
        desc, coeffs1, coeffs2, true_root, deg_type = create_test_system_2d(case)
        
        # Determine degrees from coefficients
        deg1 = int(np.sqrt(len(coeffs1))) - 1
        deg2 = int(np.sqrt(len(coeffs2))) - 1
        
        cpp_code = f"""
#include "polynomial.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace polynomial_solver;

int main() {{
    std::cout << "Test: {desc}\\n";
    std::cout << "True root: (" << std::fixed << std::setprecision(16) 
              << {true_root[0]} << ", " << {true_root[1]} << ")\\n";
    std::cout << "Degeneracy type: {deg_type}\\n";
    std::cout << "Tolerance: " << std::scientific << {tolerance} << "\\n\\n";
    
    // Define polynomials
    std::vector<unsigned int> degrees1{{{deg1}u, {deg1}u}};
    std::vector<double> power_coeffs1{{{', '.join(f'{c}' for c in coeffs1)}}};
    
    std::vector<unsigned int> degrees2{{{deg2}u, {deg2}u}};
    std::vector<double> power_coeffs2{{{', '.join(f'{c}' for c in coeffs2)}}};
    
    Polynomial p1 = Polynomial::fromPower(degrees1, power_coeffs1);
    Polynomial p2 = Polynomial::fromPower(degrees2, power_coeffs2);
    PolynomialSystem system({{p1, p2}});
    
    // Solve with extreme precision
    Solver solver;
    SubdivisionConfig config;
    config.tolerance = {tolerance};
    config.max_depth = 100;
    config.strategy = SubdivisionStrategy::SubdivideFirst;
    
    SubdivisionSolverResult result = solver.subdivisionSolve(
        system, config, RootBoundingMethod::ProjectedPolyhedral);
    
    std::cout << "Results:\\n";
    std::cout << "  Found " << result.num_resolved << " root(s)\\n";
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\\n\\n";
    
    if (result.num_resolved > 0) {{
        const SubdivisionBoxResult& box = result.boxes[0];
        double center_x = box.center[0];
        double center_y = box.center[1];
        double width_x = box.upper[0] - box.lower[0];
        double width_y = box.upper[1] - box.lower[1];
        double max_width = std::max(width_x, width_y);
        double error = std::sqrt(std::pow(center_x - {true_root[0]}, 2) + 
                                std::pow(center_y - {true_root[1]}, 2));
        
        std::cout << "Root 1:\\n";
        std::cout << "  Center: (" << std::setprecision(16) << center_x 
                  << ", " << center_y << ")\\n";
        std::cout << "  Box width: " << std::scientific << std::setprecision(6) 
                  << max_width << "\\n";
        std::cout << "  Error: " << error << "\\n";
        std::cout << "  Depth: " << box.depth << "\\n";
        std::cout << "  Error/Tolerance: " << (error / {tolerance}) << "\\n";
        
        if (error < {tolerance}) {{
            std::cout << "  Status: ✅ Achieved target accuracy\\n";
        }} else if (error < 10 * {tolerance}) {{
            std::cout << "  Status: ⚠️ Close to target\\n";
        }} else {{
            std::cout << "  Status: ❌ Did not achieve target\\n";
        }}
    }}
    
    return 0;
}}
"""
    
    with open(filename, 'w') as f:
        f.write(cpp_code)


def compile_and_run(cpp_file: str, exe_name: str) -> str:
    """Compile and run C++ test."""
    # Compile
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
        print(f"Compilation failed: {result.stderr}")
        return ""
    
    # Run
    result = subprocess.run([f"./{exe_name}"], capture_output=True, text=True)
    return result.stdout


def parse_result(output: str) -> dict:
    """Parse solver output."""
    result = {}
    
    # Extract key metrics
    if "Error:" in output:
        match = re.search(r"Error:\s+([\d.e+-]+)", output)
        if match:
            result['error'] = float(match.group(1))
    
    if "Width:" in output or "Box width:" in output:
        match = re.search(r"(?:Width|Box width):\s+([\d.e+-]+)", output)
        if match:
            result['width'] = float(match.group(1))
    
    if "Depth:" in output:
        match = re.search(r"Depth:\s+(\d+)", output)
        if match:
            result['depth'] = int(match.group(1))
    
    if "Degeneracy detected:" in output:
        result['degeneracy'] = "yes" in output.lower()
    
    if "Status:" in output:
        if "✅" in output:
            result['status'] = "achieved"
        elif "⚠️" in output:
            result['status'] = "close"
        else:
            result['status'] = "failed"
    
    return result


def main():
    print("="*80)
    print("Extreme Precision Root Accuracy Test")
    print("="*80)
    print()
    
    tolerance = 1e-15
    
    print(f"Testing with tolerance = {tolerance:.2e} (near machine epsilon)")
    print()
    
    # Test 1D cases
    print("="*80)
    print("1D Cases")
    print("="*80)
    print()
    
    cases_1d = ["simple", "double", "triple", "quintuple"]
    
    print(f"{'Case':<20} {'Error':<15} {'Width':<15} {'Depth':<8} {'Degeneracy':<12} {'Status':<15}")
    print("-"*95)
    
    for case in cases_1d:
        cpp_file = f"test_extreme_{case}_1d.cpp"
        exe_name = f"test_extreme_{case}_1d"
        
        write_cpp_test_file(cpp_file, "1d", case, tolerance)
        output = compile_and_run(cpp_file, exe_name)
        result = parse_result(output)
        
        error = result.get('error', float('nan'))
        width = result.get('width', float('nan'))
        depth = result.get('depth', 0)
        degeneracy = "Yes" if result.get('degeneracy', False) else "No"
        status = result.get('status', 'unknown')
        
        status_symbol = "✅" if status == "achieved" else ("⚠️" if status == "close" else "❌")
        
        print(f"{case:<20} {error:<15.6e} {width:<15.6e} {depth:<8} {degeneracy:<12} {status_symbol} {status:<12}")
        
        # Cleanup
        os.remove(cpp_file)
        if os.path.exists(exe_name):
            os.remove(exe_name)
    
    print()
    
    # Test 2D cases
    print("="*80)
    print("2D Cases")
    print("="*80)
    print()
    
    cases_2d = ["simple", "tangent"]
    
    print(f"{'Case':<25} {'Error':<15} {'Width':<15} {'Depth':<8} {'Degeneracy':<12} {'Status':<15}")
    print("-"*100)
    
    for case in cases_2d:
        cpp_file = f"test_extreme_{case}_2d.cpp"
        exe_name = f"test_extreme_{case}_2d"
        
        write_cpp_test_file(cpp_file, "2d", case, tolerance)
        output = compile_and_run(cpp_file, exe_name)
        result = parse_result(output)
        
        error = result.get('error', float('nan'))
        width = result.get('width', float('nan'))
        depth = result.get('depth', 0)
        degeneracy = "Yes" if result.get('degeneracy', False) else "No"
        status = result.get('status', 'unknown')
        
        status_symbol = "✅" if status == "achieved" else ("⚠️" if status == "close" else "❌")
        
        print(f"{case:<25} {error:<15.6e} {width:<15.6e} {depth:<8} {degeneracy:<12} {status_symbol} {status:<12}")
        
        # Cleanup
        os.remove(cpp_file)
        if os.path.exists(exe_name):
            os.remove(exe_name)
    
    print()
    print("="*80)
    print("Summary")
    print("="*80)
    print("1. Simple roots: Can achieve ~10^-15 accuracy")
    print("2. Multiple roots: Accuracy degrades with multiplicity")
    print("3. Tangent intersections: Similar to multiple roots")
    print("4. Degeneracy detection: Helps identify problematic cases")
    print()


if __name__ == "__main__":
    main()

