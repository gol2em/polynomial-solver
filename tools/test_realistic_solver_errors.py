#!/usr/bin/env python3
"""
Analyze accumulated errors for realistic solver scenarios.

Test cases:
1. 1D polynomials (degree 3, 5, 10)
2. 2D polynomials (degree (3,3), (5,5), (10,10))
3. Realistic solver parameters (max_depth = 30-50, tolerance = 1e-8)

Measure:
- Error after typical solver run (depth 20-30)
- Error vs. tolerance
- Whether error handling is necessary
"""

import numpy as np
import mpmath
from typing import List, Tuple

mpmath.mp.dps = 50


def de_casteljau_subdivide_double(coeffs: np.ndarray, t: float) -> Tuple[np.ndarray, np.ndarray]:
    """Subdivide at t using double precision."""
    n = len(coeffs) - 1
    b = np.zeros((n + 1, n + 1))
    b[0, :] = coeffs
    
    for r in range(1, n + 1):
        for i in range(n - r + 1):
            b[r, i] = (1 - t) * b[r-1, i] + t * b[r-1, i+1]
    
    left = np.array([b[j, 0] for j in range(n + 1)])
    right = np.array([b[n - j, j] for j in range(n + 1)])
    
    return left, right


def de_casteljau_subdivide_hp(coeffs: List, t) -> Tuple[List, List]:
    """Subdivide at t using high precision."""
    n = len(coeffs) - 1
    b = [[mpmath.mpf(0) for _ in range(n + 1)] for _ in range(n + 1)]
    for i in range(n + 1):
        b[0][i] = mpmath.mpf(coeffs[i])
    
    t_mp = mpmath.mpf(t)
    for r in range(1, n + 1):
        for i in range(n - r + 1):
            b[r][i] = (1 - t_mp) * b[r-1][i] + t_mp * b[r-1][i+1]
    
    left = [b[j][0] for j in range(n + 1)]
    right = [b[n - j][j] for j in range(n + 1)]
    
    return left, right


def simulate_solver_run(coeffs: np.ndarray, max_depth: int, num_contracts: int = 5) -> Tuple[np.ndarray, List]:
    """
    Simulate a realistic solver run with contract + subdivide operations.
    
    Returns: (final_coeffs_double, final_coeffs_hp)
    """
    current_double = coeffs.copy()
    current_hp = [mpmath.mpf(c) for c in coeffs]
    
    for depth in range(max_depth):
        # Contract: restrict to [0.1, 0.9] (simulate bound tightening)
        for _ in range(num_contracts):
            # Restrict to [0.1, 1.0]
            left_d, right_d = de_casteljau_subdivide_double(current_double, 0.1)
            left_hp, right_hp = de_casteljau_subdivide_hp(current_hp, mpmath.mpf(0.1))
            current_double = right_d
            current_hp = right_hp
            
            # Restrict to [0, 0.9] (relative to [0.1, 1.0])
            t_rel = (0.9 - 0.1) / (1.0 - 0.1)
            left_d, right_d = de_casteljau_subdivide_double(current_double, t_rel)
            left_hp, right_hp = de_casteljau_subdivide_hp(current_hp, mpmath.mpf(t_rel))
            current_double = left_d
            current_hp = left_hp
        
        # Subdivide: split at midpoint
        left_d, right_d = de_casteljau_subdivide_double(current_double, 0.5)
        left_hp, right_hp = de_casteljau_subdivide_hp(current_hp, mpmath.mpf(0.5))
        current_double = left_d  # Take left child
        current_hp = left_hp
    
    return current_double, current_hp


def compute_error(result: np.ndarray, reference: List) -> Tuple[float, float, float]:
    """Compute absolute, relative, and normalized errors."""
    abs_errors = [abs(float(reference[i]) - result[i]) for i in range(len(result))]
    max_abs_error = max(abs_errors)
    
    norm_ref = float(mpmath.sqrt(sum(c**2 for c in reference)))
    rel_error = max_abs_error / norm_ref if norm_ref > 0 else 0
    
    # Normalized error: relative to initial coefficient magnitude
    initial_norm = 1.0  # Assume initial ||b|| ≈ 1
    normalized_error = max_abs_error / initial_norm
    
    return max_abs_error, rel_error, normalized_error


def test_1d_polynomials():
    """Test 1D polynomials of various degrees."""
    print("="*80)
    print("1D Polynomial Error Analysis")
    print("="*80)
    print()
    
    degrees = [3, 5, 10]
    max_depth = 30
    tolerance = 1e-8
    
    print(f"Solver parameters:")
    print(f"  Max depth: {max_depth}")
    print(f"  Tolerance: {tolerance:.2e}")
    print(f"  Contract iterations per depth: 5")
    print()
    
    print(f"{'Degree':<8} {'Max Abs Error':<20} {'Rel Error':<20} {'Error/Tolerance':<20} {'Status':<15}")
    print("-"*90)
    
    for degree in degrees:
        # Create test polynomial: x^degree
        coeffs = np.array([(i / degree) ** degree for i in range(degree + 1)])
        
        # Simulate solver run
        result_double, result_hp = simulate_solver_run(coeffs, max_depth)
        
        # Compute errors
        max_abs_error, rel_error, _ = compute_error(result_double, result_hp)
        error_ratio = max_abs_error / tolerance
        
        # Determine status
        if error_ratio < 0.01:
            status = "✅ Excellent"
        elif error_ratio < 0.1:
            status = "✅ Good"
        elif error_ratio < 1.0:
            status = "⚠️ Acceptable"
        else:
            status = "❌ Problematic"
        
        print(f"{degree:<8} {max_abs_error:<20.6e} {rel_error:<20.6e} {error_ratio:<20.6e} {status:<15}")
    
    print()


def test_2d_polynomials():
    """Test 2D polynomials (tensor product)."""
    print("="*80)
    print("2D Polynomial Error Analysis")
    print("="*80)
    print()
    
    degrees = [(3, 3), (5, 5), (10, 10)]
    max_depth = 30
    tolerance = 1e-8
    
    print(f"Solver parameters:")
    print(f"  Max depth: {max_depth}")
    print(f"  Tolerance: {tolerance:.2e}")
    print(f"  Contract iterations per depth: 5")
    print(f"  Note: 2D requires 2× operations per dimension")
    print()
    
    print(f"{'Degree':<12} {'Max Abs Error':<20} {'Rel Error':<20} {'Error/Tolerance':<20} {'Status':<15}")
    print("-"*95)
    
    for deg_x, deg_y in degrees:
        # For 2D, we need to apply operations along both dimensions
        # Approximate by doubling the number of operations
        
        # Create 1D slice (simplified model)
        max_degree = max(deg_x, deg_y)
        coeffs = np.array([(i / max_degree) ** max_degree for i in range(max_degree + 1)])
        
        # Simulate with 2× operations (both dimensions)
        result_double, result_hp = simulate_solver_run(coeffs, max_depth, num_contracts=10)
        
        # Compute errors
        max_abs_error, rel_error, _ = compute_error(result_double, result_hp)
        error_ratio = max_abs_error / tolerance
        
        # Determine status
        if error_ratio < 0.01:
            status = "✅ Excellent"
        elif error_ratio < 0.1:
            status = "✅ Good"
        elif error_ratio < 1.0:
            status = "⚠️ Acceptable"
        else:
            status = "❌ Problematic"
        
        print(f"({deg_x},{deg_y}){'':<6} {max_abs_error:<20.6e} {rel_error:<20.6e} {error_ratio:<20.6e} {status:<15}")
    
    print()


def test_depth_sensitivity():
    """Test how error grows with depth."""
    print("="*80)
    print("Error vs. Solver Depth (Degree 10, 1D)")
    print("="*80)
    print()
    
    degree = 10
    depths = [10, 20, 30, 40, 50]
    tolerance = 1e-8
    
    coeffs = np.array([(i / degree) ** degree for i in range(degree + 1)])
    
    print(f"{'Depth':<8} {'Max Abs Error':<20} {'Rel Error':<20} {'Error/Tolerance':<20} {'Status':<15}")
    print("-"*90)
    
    for depth in depths:
        result_double, result_hp = simulate_solver_run(coeffs, depth)
        max_abs_error, rel_error, _ = compute_error(result_double, result_hp)
        error_ratio = max_abs_error / tolerance
        
        if error_ratio < 0.01:
            status = "✅ Excellent"
        elif error_ratio < 0.1:
            status = "✅ Good"
        elif error_ratio < 1.0:
            status = "⚠️ Acceptable"
        else:
            status = "❌ Problematic"
        
        print(f"{depth:<8} {max_abs_error:<20.6e} {rel_error:<20.6e} {error_ratio:<20.6e} {status:<15}")
    
    print()


if __name__ == "__main__":
    test_1d_polynomials()
    test_2d_polynomials()
    test_depth_sensitivity()
    
    print("="*80)
    print("Summary and Recommendations")
    print("="*80)
    print()
    print("Key Findings:")
    print("1. Check if errors are << tolerance (error/tolerance < 0.01)")
    print("2. Determine if error handling is necessary")
    print("3. Identify problematic cases (if any)")
    print()

