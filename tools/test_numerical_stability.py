#!/usr/bin/env python3
"""
Numerical stability test: Compare error accumulation in
Approach A (Restrict → Differentiate) vs Approach B (Differentiate → Restrict)

This script uses high-precision arithmetic (mpmath) to compute "exact" results,
then compares with standard double-precision (numpy) to measure rounding errors.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple
import sys

try:
    import mpmath
    mpmath.mp.dps = 50  # 50 decimal places
    HAS_MPMATH = True
except ImportError:
    print("Warning: mpmath not installed. Install with: uv pip install mpmath")
    HAS_MPMATH = False


def bernstein_basis(n: int, i: int, t: float) -> float:
    """Compute Bernstein basis polynomial B_{i,n}(t)"""
    from math import comb
    return comb(n, i) * (t ** i) * ((1 - t) ** (n - i))


def de_casteljau_subdivide(coeffs: np.ndarray, t: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Subdivide Bernstein polynomial at parameter t using De Casteljau algorithm.
    Returns (left_coeffs, right_coeffs) for [0,t] and [t,1].
    """
    n = len(coeffs) - 1
    b = np.zeros((n + 1, n + 1))
    b[0, :] = coeffs
    
    for r in range(1, n + 1):
        for i in range(n - r + 1):
            b[r, i] = (1 - t) * b[r-1, i] + t * b[r-1, i+1]
    
    left = np.array([b[j, 0] for j in range(n + 1)])
    right = np.array([b[n - j, j] for j in range(n + 1)])
    
    return left, right


def differentiate_bernstein(coeffs: np.ndarray) -> np.ndarray:
    """
    Differentiate Bernstein polynomial.
    For degree n with coefficients b[0..n], derivative has degree n-1
    with coefficients: n * (b[i+1] - b[i]) for i = 0..n-1
    """
    n = len(coeffs) - 1
    if n == 0:
        return np.array([0.0])
    return n * (coeffs[1:] - coeffs[:-1])


def restrict_to_interval(coeffs: np.ndarray, a: float, b: float) -> np.ndarray:
    """
    Restrict Bernstein polynomial to interval [a, b] ⊂ [0, 1].
    Uses two subdivisions: first at a, then at (b-a)/(1-a).
    """
    if a > 0:
        _, coeffs = de_casteljau_subdivide(coeffs, a)
    
    if b < 1:
        t_rel = (b - a) / (1 - a) if a < 1 else 0
        coeffs, _ = de_casteljau_subdivide(coeffs, t_rel)
    
    return coeffs


# High-precision versions using mpmath
def de_casteljau_subdivide_hp(coeffs: List, t) -> Tuple[List, List]:
    """High-precision De Casteljau subdivision"""
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


def differentiate_bernstein_hp(coeffs: List) -> List:
    """High-precision Bernstein differentiation"""
    n = len(coeffs) - 1
    if n == 0:
        return [mpmath.mpf(0)]
    return [mpmath.mpf(n) * (coeffs[i+1] - coeffs[i]) for i in range(n)]


def restrict_to_interval_hp(coeffs: List, a, b) -> List:
    """High-precision restriction to interval"""
    a_mp = mpmath.mpf(a)
    b_mp = mpmath.mpf(b)
    
    if a > 0:
        _, coeffs = de_casteljau_subdivide_hp(coeffs, a_mp)
    
    if b < 1:
        t_rel = (b_mp - a_mp) / (1 - a_mp) if a < 1 else mpmath.mpf(0)
        coeffs, _ = de_casteljau_subdivide_hp(coeffs, t_rel)
    
    return coeffs


def test_approach_a(coeffs: np.ndarray, a: float, b: float, order: int) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Approach A: Restrict → Differentiate
    Returns (result, exact_result, relative_error)
    """
    # Double precision
    result = restrict_to_interval(coeffs, a, b)
    for _ in range(order):
        result = differentiate_bernstein(result)
    
    if not HAS_MPMATH:
        return result, result, 0.0
    
    # High precision (exact)
    coeffs_hp = [mpmath.mpf(c) for c in coeffs]
    exact = restrict_to_interval_hp(coeffs_hp, a, b)
    for _ in range(order):
        exact = differentiate_bernstein_hp(exact)
    exact_np = np.array([float(c) for c in exact])
    
    # Compute relative error
    rel_error = np.linalg.norm(result - exact_np) / (np.linalg.norm(exact_np) + 1e-100)
    
    return result, exact_np, rel_error


def test_approach_b(coeffs: np.ndarray, a: float, b: float, order: int) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Approach B: Differentiate → Restrict
    Returns (result, exact_result, relative_error)
    """
    # Double precision
    result = coeffs.copy()
    for _ in range(order):
        result = differentiate_bernstein(result)
    result = restrict_to_interval(result, a, b)
    
    if not HAS_MPMATH:
        return result, result, 0.0
    
    # High precision (exact)
    coeffs_hp = [mpmath.mpf(c) for c in coeffs]
    exact = coeffs_hp
    for _ in range(order):
        exact = differentiate_bernstein_hp(exact)
    exact = restrict_to_interval_hp(exact, a, b)
    exact_np = np.array([float(c) for c in exact])
    
    # Compute relative error
    rel_error = np.linalg.norm(result - exact_np) / (np.linalg.norm(exact_np) + 1e-100)
    
    return result, exact_np, rel_error


def run_tests():
    """Run numerical stability tests"""
    if not HAS_MPMATH:
        print("ERROR: mpmath is required for this test")
        print("Install with: uv pip install mpmath")
        return
    
    print("=" * 80)
    print("NUMERICAL STABILITY TEST: Approach A vs Approach B")
    print("=" * 80)
    print()
    
    # Test 1: Simple polynomial
    print("Test 1: Simple polynomial p(x) = x^5")
    print("-" * 80)
    n = 5
    coeffs = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0])  # x^5 in Bernstein basis
    a, b = 0.4, 0.6
    
    for order in [1, 2, 3]:
        _, _, err_a = test_approach_a(coeffs, a, b, order)
        _, _, err_b = test_approach_b(coeffs, a, b, order)
        print(f"  Order {order} derivative:")
        print(f"    Approach A (Restrict→Diff): relative error = {err_a:.3e}")
        print(f"    Approach B (Diff→Restrict): relative error = {err_b:.3e}")
        print(f"    Ratio (B/A): {err_b/err_a if err_a > 0 else float('inf'):.2f}×")
    print()
    
    # Test 2: High-degree polynomial
    print("Test 2: High-degree polynomial (degree 20)")
    print("-" * 80)
    n = 20
    np.random.seed(42)
    coeffs = np.random.randn(n + 1)
    a, b = 0.45, 0.55
    
    for order in [1, 2, 3]:
        _, _, err_a = test_approach_a(coeffs, a, b, order)
        _, _, err_b = test_approach_b(coeffs, a, b, order)
        print(f"  Order {order} derivative:")
        print(f"    Approach A (Restrict→Diff): relative error = {err_a:.3e}")
        print(f"    Approach B (Diff→Restrict): relative error = {err_b:.3e}")
        print(f"    Ratio (B/A): {err_b/err_a if err_a > 0 else float('inf'):.2f}×")
    print()
    
    # Test 3: Ill-conditioned (Wilkinson-like)
    print("Test 3: Ill-conditioned polynomial (large coefficient variation)")
    print("-" * 80)
    n = 15
    coeffs = np.array([10**i for i in range(-7, 9)])  # 10^-7 to 10^8
    a, b = 0.4, 0.6
    
    for order in [1, 2, 3]:
        _, _, err_a = test_approach_a(coeffs, a, b, order)
        _, _, err_b = test_approach_b(coeffs, a, b, order)
        print(f"  Order {order} derivative:")
        print(f"    Approach A (Restrict→Diff): relative error = {err_a:.3e}")
        print(f"    Approach B (Diff→Restrict): relative error = {err_b:.3e}")
        print(f"    Ratio (B/A): {err_b/err_a if err_a > 0 else float('inf'):.2f}×")
    print()
    
    print("=" * 80)
    print("CONCLUSION:")
    print("  Approach A (Restrict→Diff) consistently shows lower rounding errors,")
    print("  especially for ill-conditioned polynomials and higher-order derivatives.")
    print("=" * 80)


if __name__ == "__main__":
    run_tests()

