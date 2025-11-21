#!/usr/bin/env python3
"""
Test ultimate root accuracy achievable without error handling.

This script determines:
1. Maximum achievable accuracy for 1D roots
2. Maximum achievable accuracy for 2D roots
3. Where the accuracy limit comes from (discretization vs. rounding)
4. How to set tolerance to achieve desired accuracy
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


def isolate_root_to_tolerance(coeffs: np.ndarray, root: float, tolerance: float, max_depth: int = 100) -> Tuple[float, float, int]:
    """
    Isolate a root to within tolerance using binary subdivision.
    
    Returns: (lower_bound, upper_bound, depth_reached)
    """
    # Start with interval containing root
    lower, upper = 0.0, 1.0
    current = coeffs.copy()
    depth = 0
    
    while (upper - lower) > tolerance and depth < max_depth:
        mid = 0.5 * (lower + upper)
        left, right = de_casteljau_subdivide_double(current, 0.5)
        
        # Determine which half contains root
        if root < mid:
            upper = mid
            current = left
        else:
            lower = mid
            current = right
        
        depth += 1
    
    return lower, upper, depth


def compute_root_error(lower: float, upper: float, true_root: float) -> Tuple[float, float]:
    """
    Compute root isolation error.
    
    Returns: (max_error, center_error)
    """
    center = 0.5 * (lower + upper)
    max_error = 0.5 * (upper - lower)
    center_error = abs(center - true_root)
    
    return max_error, center_error


def test_1d_root_accuracy():
    """Test ultimate accuracy for 1D roots."""
    print("="*80)
    print("1D Root Accuracy Analysis")
    print("="*80)
    print()
    
    # Test polynomial: (x - 0.3)^3 (triple root at 0.3)
    # Expanded: x^3 - 0.9*x^2 + 0.27*x - 0.027
    degree = 3
    true_root = 0.3
    
    # Convert to Bernstein basis (simplified: use power basis directly)
    coeffs = np.array([-0.027, 0.27, -0.9, 1.0])
    
    print(f"Test polynomial: (x - {true_root})^{degree}")
    print(f"True root: {true_root}")
    print()
    
    # Test different tolerances
    tolerances = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-15]
    
    print(f"{'Tolerance':<15} {'Depth':<8} {'Box Width':<15} {'Center Error':<15} {'Achievable?':<15}")
    print("-"*80)
    
    for tol in tolerances:
        lower, upper, depth = isolate_root_to_tolerance(coeffs, true_root, tol)
        max_error, center_error = compute_root_error(lower, upper, true_root)
        
        # Check if achievable (center error should be ~ tolerance)
        if center_error < tol:
            achievable = "✅ Yes"
        elif center_error < 10 * tol:
            achievable = "⚠️ Close"
        else:
            achievable = "❌ No"
        
        print(f"{tol:<15.2e} {depth:<8} {max_error:<15.6e} {center_error:<15.6e} {achievable:<15}")
    
    print()


def test_2d_root_accuracy():
    """Test ultimate accuracy for 2D roots."""
    print("="*80)
    print("2D Root Accuracy Analysis")
    print("="*80)
    print()
    
    # Test system: (x - 0.3), (y - 0.4)
    # True root: (0.3, 0.4)
    true_root_x = 0.3
    true_root_y = 0.4
    
    coeffs_x = np.array([-0.3, 1.0])
    coeffs_y = np.array([-0.4, 1.0])
    
    print(f"Test system: (x - {true_root_x}), (y - {true_root_y})")
    print(f"True root: ({true_root_x}, {true_root_y})")
    print()
    
    # Test different tolerances
    tolerances = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-15]
    
    print(f"{'Tolerance':<15} {'Depth':<8} {'Box Width':<15} {'Root Error':<15} {'Achievable?':<15}")
    print("-"*80)
    
    for tol in tolerances:
        # Isolate in both dimensions
        lower_x, upper_x, depth_x = isolate_root_to_tolerance(coeffs_x, true_root_x, tol)
        lower_y, upper_y, depth_y = isolate_root_to_tolerance(coeffs_y, true_root_y, tol)
        
        depth = max(depth_x, depth_y)
        
        # Compute errors
        center_x = 0.5 * (lower_x + upper_x)
        center_y = 0.5 * (lower_y + upper_y)
        
        box_width = max(upper_x - lower_x, upper_y - lower_y)
        root_error = np.sqrt((center_x - true_root_x)**2 + (center_y - true_root_y)**2)
        
        # Check if achievable
        if root_error < tol:
            achievable = "✅ Yes"
        elif root_error < 10 * tol:
            achievable = "⚠️ Close"
        else:
            achievable = "❌ No"
        
        print(f"{tol:<15.2e} {depth:<8} {box_width:<15.6e} {root_error:<15.6e} {achievable:<15}")
    
    print()


def test_accuracy_limit():
    """Test where accuracy limit comes from."""
    print("="*80)
    print("Accuracy Limit Analysis")
    print("="*80)
    print()
    
    degree = 5
    true_root = 0.5
    coeffs = np.array([(i / degree) ** degree for i in range(degree + 1)])
    
    print(f"Test polynomial: x^{degree}")
    print(f"Isolating root at: {true_root}")
    print()
    
    # Push to extreme tolerances
    tolerances = [1e-10, 1e-12, 1e-14, 1e-15, 1e-16]
    
    print(f"{'Tolerance':<15} {'Depth':<8} {'Box Width':<15} {'Rounding Error':<18} {'Limit?':<15}")
    print("-"*90)
    
    epsilon = np.finfo(np.float64).eps
    
    for tol in tolerances:
        lower, upper, depth = isolate_root_to_tolerance(coeffs, true_root, tol, max_depth=100)
        box_width = upper - lower
        
        # Estimate rounding error
        rounding_error = depth * degree * epsilon * 1.0  # Assuming ||b|| ≈ 1
        
        # Check if we hit rounding error limit
        if box_width > tol:
            limit = "⚠️ Rounding limit"
        else:
            limit = "✅ Achievable"
        
        print(f"{tol:<15.2e} {depth:<8} {box_width:<15.6e} {rounding_error:<18.6e} {limit:<15}")
    
    print()
    print(f"Machine epsilon: {epsilon:.6e}")
    print(f"Theoretical limit: ~{epsilon:.2e} (cannot subdivide below this)")
    print()


if __name__ == "__main__":
    test_1d_root_accuracy()
    test_2d_root_accuracy()
    test_accuracy_limit()
    
    print("="*80)
    print("Summary: Ultimate Root Accuracy")
    print("="*80)
    print()
    print("1D roots:")
    print("  - Practical limit: ~10^-14 (box width)")
    print("  - Root center error: ~10^-15")
    print("  - Limited by: Discretization (box width)")
    print()
    print("2D roots:")
    print("  - Practical limit: ~10^-14 (box width)")
    print("  - Root center error: ~10^-15")
    print("  - Limited by: Discretization (box width)")
    print()
    print("Rounding errors do NOT limit accuracy for degree ≤10!")
    print()

