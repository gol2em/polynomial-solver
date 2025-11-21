#!/usr/bin/env python3
"""
Analyze coefficient norms for example polynomials.

This script computes:
1. Power basis coefficients
2. Bernstein basis coefficients (via conversion)
3. L2 norm, max norm, condition number
4. Predicted rounding error
"""

import numpy as np
import mpmath
from typing import List, Tuple
import matplotlib.pyplot as plt

# Set high precision for accurate conversion
mpmath.mp.dps = 50


def power_to_bernstein_1d(power_coeffs: np.ndarray) -> np.ndarray:
    """
    Convert 1D polynomial from power basis to Bernstein basis.
    
    Power basis: p(x) = sum_i c_i * x^i
    Bernstein basis: p(x) = sum_i b_i * B_i^n(x)
    
    Conversion formula: b_i = sum_{j=i}^n c_j * C(j,i) / C(n,i)
    where C(n,k) = n! / (k! * (n-k)!)
    """
    n = len(power_coeffs) - 1
    bernstein_coeffs = np.zeros(n + 1)
    
    # Compute binomial coefficients
    def binomial(n, k):
        if k < 0 or k > n:
            return 0
        result = 1
        for i in range(min(k, n - k)):
            result = result * (n - i) // (i + 1)
        return result
    
    for i in range(n + 1):
        for j in range(i, n + 1):
            bernstein_coeffs[i] += power_coeffs[j] * binomial(j, i) / binomial(n, i)
    
    return bernstein_coeffs


def analyze_polynomial(name: str, power_coeffs: np.ndarray, depth: int = 10):
    """
    Analyze a polynomial's coefficient properties and predicted error.
    """
    print(f"\n{'='*80}")
    print(f"Polynomial: {name}")
    print(f"{'='*80}")
    
    degree = len(power_coeffs) - 1
    print(f"Degree: {degree}")
    
    # Power basis analysis
    print(f"\nPower Basis Coefficients:")
    print(f"  Count: {len(power_coeffs)}")
    print(f"  Max: {np.max(np.abs(power_coeffs)):.6e}")
    print(f"  Min (non-zero): {np.min(np.abs(power_coeffs[power_coeffs != 0])):.6e}")
    print(f"  L2 norm: {np.linalg.norm(power_coeffs):.6e}")
    print(f"  Max/Min ratio: {np.max(np.abs(power_coeffs)) / np.min(np.abs(power_coeffs[power_coeffs != 0])):.6e}")
    
    # Convert to Bernstein basis
    bernstein_coeffs = power_to_bernstein_1d(power_coeffs)
    
    print(f"\nBernstein Basis Coefficients:")
    print(f"  Max: {np.max(np.abs(bernstein_coeffs)):.6e}")
    print(f"  Min (non-zero): {np.min(np.abs(bernstein_coeffs[bernstein_coeffs != 0])):.6e}")
    print(f"  L2 norm: {np.linalg.norm(bernstein_coeffs):.6e}")
    print(f"  Max/Min ratio: {np.max(np.abs(bernstein_coeffs)) / np.min(np.abs(bernstein_coeffs[bernstein_coeffs != 0])):.6e}")
    
    # Error prediction
    norm_b = np.linalg.norm(bernstein_coeffs)
    epsilon = np.finfo(np.float64).eps
    predicted_error = depth * degree * epsilon * norm_b
    
    print(f"\nError Prediction (depth={depth}):")
    print(f"  ||b|| = {norm_b:.6e}")
    print(f"  Predicted error: {predicted_error:.6e}")
    print(f"  Relative error: {predicted_error / norm_b:.6e}")
    
    # Rescaling benefit
    if norm_b > 1.0:
        rescaled_error = depth * degree * epsilon * 1.0
        improvement = predicted_error / rescaled_error
        print(f"\nRescaling Benefit:")
        print(f"  Error after rescaling: {rescaled_error:.6e}")
        print(f"  Improvement factor: {improvement:.2f}×")
    else:
        print(f"\nRescaling Benefit:")
        print(f"  Already well-scaled (||b|| ≈ 1)")
    
    return {
        'name': name,
        'degree': degree,
        'power_norm': np.linalg.norm(power_coeffs),
        'bernstein_norm': norm_b,
        'predicted_error': predicted_error,
        'power_coeffs': power_coeffs,
        'bernstein_coeffs': bernstein_coeffs
    }


def wilkinson_coefficients(n: int) -> np.ndarray:
    """
    Compute coefficients of (x - 1/n)(x - 2/n)...(x - (n-1)/n).
    """
    scale = 1.0 / n
    coeffs = np.array([-scale, 1.0])
    
    for k in range(2, n):
        new_coeffs = np.zeros(len(coeffs) + 1)
        root = k * scale
        
        for i in range(len(coeffs)):
            new_coeffs[i] -= root * coeffs[i]
            new_coeffs[i + 1] += coeffs[i]
        
        coeffs = new_coeffs
    
    return coeffs


def main():
    print("Coefficient Norm Analysis for Example Polynomials")
    print("="*80)
    
    results = []
    
    # 1. Cubic polynomial: (x - 0.2)(x - 0.5)(x - 0.8)
    cubic_coeffs = np.array([-0.08, 0.66, -1.5, 1.0])
    results.append(analyze_polynomial("Cubic: (x-0.2)(x-0.5)(x-0.8)", cubic_coeffs))
    
    # 2. Multiplicity polynomial: (x - 0.2)(x - 0.6)^6
    multiplicity_coeffs = np.array([
        -0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0
    ])
    results.append(analyze_polynomial("Multiplicity: (x-0.2)(x-0.6)^6", multiplicity_coeffs))
    
    # 3. Wilkinson polynomial (scaled to [0,1])
    wilkinson_coeffs = wilkinson_coefficients(20)
    results.append(analyze_polynomial("Wilkinson: (x-1/20)...(x-19/20)", wilkinson_coeffs))
    
    # 4. Synthetic ill-conditioned: scale cubic by 10^6
    ill_conditioned_coeffs = cubic_coeffs * 1e6
    results.append(analyze_polynomial("Ill-conditioned: 10^6 × cubic", ill_conditioned_coeffs))
    
    # Summary comparison
    print(f"\n{'='*80}")
    print("Summary Comparison")
    print(f"{'='*80}\n")
    
    print(f"{'Polynomial':<35} {'Degree':<8} {'||b|| (Bernstein)':<20} {'Predicted Error':<20}")
    print("-"*80)
    for r in results:
        print(f"{r['name']:<35} {r['degree']:<8} {r['bernstein_norm']:<20.6e} {r['predicted_error']:<20.6e}")
    
    print(f"\n{'='*80}")
    print("Key Findings:")
    print(f"{'='*80}")
    
    # Find worst case
    worst = max(results, key=lambda r: r['bernstein_norm'])
    print(f"\n1. Worst-case coefficient norm: {worst['name']}")
    print(f"   ||b|| = {worst['bernstein_norm']:.6e}")
    print(f"   Predicted error: {worst['predicted_error']:.6e}")
    
    # Find best case
    best = min(results, key=lambda r: r['bernstein_norm'])
    print(f"\n2. Best-case coefficient norm: {best['name']}")
    print(f"   ||b|| = {best['bernstein_norm']:.6e}")
    print(f"   Predicted error: {best['predicted_error']:.6e}")
    
    # Ratio
    ratio = worst['bernstein_norm'] / best['bernstein_norm']
    print(f"\n3. Worst/Best ratio: {ratio:.2f}×")
    print(f"   Error ratio: {worst['predicted_error'] / best['predicted_error']:.2f}×")
    
    print(f"\n4. Rescaling recommendations:")
    for r in results:
        if r['bernstein_norm'] > 1000:
            print(f"   ⚠️  {r['name']}: RESCALE (||b|| = {r['bernstein_norm']:.2e})")
        elif r['bernstein_norm'] > 10:
            print(f"   ⚡ {r['name']}: Consider rescaling (||b|| = {r['bernstein_norm']:.2e})")
        else:
            print(f"   ✅ {r['name']}: Well-scaled (||b|| = {r['bernstein_norm']:.2e})")
    
    print()


if __name__ == "__main__":
    main()

