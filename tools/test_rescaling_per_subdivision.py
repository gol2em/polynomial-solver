#!/usr/bin/env python3
"""
Test whether rescaling after each de Casteljau subdivision reduces rounding errors.

This script compares three strategies:
1. No rescaling: Standard de Casteljau subdivision
2. Rescale after each subdivision: Normalize ||b|| = 1 after each step
3. Rescale only at start: Normalize once before subdivisions

We measure actual errors using high-precision reference (mpmath).
"""

import numpy as np
import mpmath
from typing import List, Tuple

# Set high precision for reference computations
mpmath.mp.dps = 50


def compute_norm(coeffs):
    """Compute L2 norm of coefficients."""
    if isinstance(coeffs, np.ndarray):
        return np.linalg.norm(coeffs)
    else:
        return float(mpmath.sqrt(sum(c**2 for c in coeffs)))


def de_casteljau_subdivide_double(coeffs: np.ndarray, t: float) -> Tuple[np.ndarray, np.ndarray]:
    """Subdivide at t using double precision (no rescaling)."""
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


def strategy_no_rescaling(coeffs: np.ndarray, subdivisions: List[Tuple[str, float]]) -> np.ndarray:
    """Strategy 1: No rescaling."""
    current = coeffs.copy()
    for direction, t in subdivisions:
        left, right = de_casteljau_subdivide_double(current, t)
        current = left if direction == 'left' else right
    return current


def strategy_rescale_each(coeffs: np.ndarray, subdivisions: List[Tuple[str, float]]) -> np.ndarray:
    """Strategy 2: Rescale after each subdivision."""
    current = coeffs.copy()
    for direction, t in subdivisions:
        left, right = de_casteljau_subdivide_double(current, t)
        current = left if direction == 'left' else right
        
        # Rescale to ||b|| = 1
        norm = np.linalg.norm(current)
        if norm > 0:
            current = current / norm
    
    return current


def strategy_rescale_start(coeffs: np.ndarray, subdivisions: List[Tuple[str, float]]) -> np.ndarray:
    """Strategy 3: Rescale only at start."""
    # Rescale input
    norm = np.linalg.norm(coeffs)
    current = coeffs / norm if norm > 0 else coeffs.copy()
    
    # Subdivide without rescaling
    for direction, t in subdivisions:
        left, right = de_casteljau_subdivide_double(current, t)
        current = left if direction == 'left' else right
    
    return current


def compute_error(result: np.ndarray, reference: List) -> Tuple[float, float]:
    """Compute absolute and relative error."""
    abs_errors = [abs(float(reference[i]) - result[i]) for i in range(len(result))]
    max_abs_error = max(abs_errors)
    
    norm_ref = compute_norm(reference)
    rel_error = max_abs_error / norm_ref if norm_ref > 0 else 0
    
    return max_abs_error, rel_error


def test_rescaling_strategies():
    """Test all three strategies and compare errors."""
    print("="*80)
    print("Test: Does Rescaling After Each Subdivision Reduce Error?")
    print("="*80)
    print()
    
    # Test polynomial: x^20
    degree = 20
    coeffs_double = np.array([(i / degree) ** degree for i in range(degree + 1)])
    
    print(f"Polynomial: x^{degree}")
    print(f"Initial ||b||: {np.linalg.norm(coeffs_double):.6e}")
    print()
    
    # Test different depths
    depths = [5, 10, 15, 20]
    
    print(f"{'Depth':<8} {'No Rescale':<20} {'Rescale Each':<20} {'Rescale Start':<20} {'Best Strategy':<20}")
    print("-"*100)
    
    for depth in depths:
        # Subdivide at midpoint repeatedly
        subdivisions = [('left', 0.5) for _ in range(depth)]
        
        # High precision reference
        coeffs_hp = [mpmath.mpf(c) for c in coeffs_double]
        current_hp = coeffs_hp
        for direction, t in subdivisions:
            left_hp, right_hp = de_casteljau_subdivide_hp(current_hp, t)
            current_hp = left_hp if direction == 'left' else right_hp
        
        # Strategy 1: No rescaling
        result_no_rescale = strategy_no_rescaling(coeffs_double, subdivisions)
        _, err_no_rescale = compute_error(result_no_rescale, current_hp)
        
        # Strategy 2: Rescale after each subdivision
        result_rescale_each = strategy_rescale_each(coeffs_double, subdivisions)
        _, err_rescale_each = compute_error(result_rescale_each, current_hp)
        
        # Strategy 3: Rescale only at start
        result_rescale_start = strategy_rescale_start(coeffs_double, subdivisions)
        _, err_rescale_start = compute_error(result_rescale_start, current_hp)
        
        # Find best strategy
        errors = {
            'No Rescale': err_no_rescale,
            'Rescale Each': err_rescale_each,
            'Rescale Start': err_rescale_start
        }
        best = min(errors, key=errors.get)
        
        print(f"{depth:<8} {err_no_rescale:<20.6e} {err_rescale_each:<20.6e} {err_rescale_start:<20.6e} {best:<20}")
    
    print()


def test_coefficient_norm_evolution():
    """Track how ||b|| evolves during subdivisions."""
    print("="*80)
    print("Test: How Does ||b|| Evolve During Subdivisions?")
    print("="*80)
    print()
    
    degree = 20
    coeffs = np.array([(i / degree) ** degree for i in range(degree + 1)])
    
    print(f"Polynomial: x^{degree}")
    print(f"Initial ||b||: {np.linalg.norm(coeffs):.6e}")
    print()
    
    # Track norm evolution
    depth = 20
    subdivisions = [('left', 0.5) for _ in range(depth)]
    
    print(f"{'Step':<8} {'||b|| (no rescale)':<25} {'||b|| (rescale each)':<25}")
    print("-"*60)
    
    # No rescaling
    current_no_rescale = coeffs.copy()
    print(f"{0:<8} {np.linalg.norm(current_no_rescale):<25.6e} {'-':<25}")
    
    # Rescale each
    current_rescale_each = coeffs.copy()
    
    for step, (direction, t) in enumerate(subdivisions, 1):
        # No rescaling
        left, right = de_casteljau_subdivide_double(current_no_rescale, t)
        current_no_rescale = left if direction == 'left' else right
        norm_no_rescale = np.linalg.norm(current_no_rescale)
        
        # Rescale each
        left, right = de_casteljau_subdivide_double(current_rescale_each, t)
        current_rescale_each = left if direction == 'left' else right
        norm_rescale_each = np.linalg.norm(current_rescale_each)
        current_rescale_each = current_rescale_each / norm_rescale_each
        
        if step % 5 == 0 or step <= 3:
            print(f"{step:<8} {norm_no_rescale:<25.6e} {norm_rescale_each:<25.6e} (→ 1.0)")
    
    print()


if __name__ == "__main__":
    test_rescaling_strategies()
    test_coefficient_norm_evolution()
    
    print("="*80)
    print("Summary:")
    print("="*80)
    print("1. Compare relative errors across strategies")
    print("2. Check if rescaling after each subdivision helps")
    print("3. Observe ||b|| evolution with/without rescaling")
    print()

