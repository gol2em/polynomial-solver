#!/usr/bin/env python3
"""
High-precision visualization of Hessian determinant zero contour.
Uses mpmath for multiprecision floating-point arithmetic.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp, sqrt, sin, atan2, mpf
import sys

# Set precision (number of decimal digits)
mp.dps = 50  # 50 decimal places

print("=" * 70)
print("High-Precision Hessian Determinant Visualization")
print("=" * 70)
print(f"\nPrecision: {mp.dps} decimal places")
print(f"Working precision: ~{int(mp.dps * 3.32)} bits\n")

def f_original_mp(x, y):
    """Original function with multiprecision arithmetic."""
    x = mpf(x)
    y = mpf(y)
    
    # f1 = 10*x^2*y^2 + sqrt(|x^2*y|)
    f1 = 10 * x * y * x * y + sqrt(abs(x * x * y))
    
    # f2 = atan2(0.001, sin(5*y) - 2*x)
    f2 = atan2(mpf('0.001'), sin(5 * y) - 2 * x)
    
    # f3 = 10*y^3 + x^3
    f3 = 10 * y * y * y + x * x * x
    
    # f4 = atan2(0.01, sin(5*x) - 2*y)
    f4 = atan2(mpf('0.01'), sin(5 * x) - 2 * y)
    
    return f1 + f2 + f3 + f4

def compute_hessian_determinant_mp(x, y, h=None):
    """Compute Hessian determinant with multiprecision using finite differences."""
    if h is None:
        h = mpf(10) ** (-mp.dps // 2)  # Adaptive step size based on precision
    
    x = mpf(x)
    y = mpf(y)
    h = mpf(h)
    
    # Second derivatives using central differences
    # f_xx = (f(x+h, y) - 2*f(x, y) + f(x-h, y)) / h^2
    f_xx = (f_original_mp(x + h, y) - 2 * f_original_mp(x, y) + f_original_mp(x - h, y)) / (h * h)
    
    # f_yy = (f(x, y+h) - 2*f(x, y) + f(x, y-h)) / h^2
    f_yy = (f_original_mp(x, y + h) - 2 * f_original_mp(x, y) + f_original_mp(x, y - h)) / (h * h)
    
    # f_xy = (f(x+h, y+h) - f(x+h, y-h) - f(x-h, y+h) + f(x-h, y-h)) / (4*h^2)
    f_xy = (f_original_mp(x + h, y + h) - f_original_mp(x + h, y - h) - 
            f_original_mp(x - h, y + h) + f_original_mp(x - h, y - h)) / (4 * h * h)
    
    # det(H) = f_xx * f_yy - f_xy^2
    det_H = f_xx * f_yy - f_xy * f_xy
    
    return float(det_H)

def compute_hessian_grid_high_precision(resolution=1000):
    """Compute Hessian determinant on a grid with high precision."""
    print(f"Computing high-precision Hessian determinant on {resolution}×{resolution} grid...")
    print(f"Total evaluations: {resolution * resolution:,}")
    print(f"This may take several minutes...\n")
    
    x = np.linspace(-1, 1, resolution)
    y = np.linspace(-1, 1, resolution)
    X, Y = np.meshgrid(x, y)
    
    det_H = np.zeros_like(X)
    
    total = resolution * resolution
    for i in range(resolution):
        for j in range(resolution):
            det_H[i, j] = compute_hessian_determinant_mp(X[i, j], Y[i, j])
        
        # Progress indicator
        if (i + 1) % (resolution // 20) == 0:
            progress = (i + 1) * resolution / total * 100
            print(f"  Progress: {progress:.1f}%")
    
    print("  Done!\n")
    return X, Y, det_H

def visualize_high_precision_contour(X, Y, det_H):
    """Visualize the zero contour with high precision."""
    print("Creating high-precision visualization...")
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 11))
    
    # Plot filled contours for context
    contourf = ax.contourf(X, Y, det_H, levels=100, cmap='RdBu_r', alpha=0.6)
    
    # Plot zero contour with high visibility
    contour_zero = ax.contour(X, Y, det_H, levels=[0], colors='black', linewidths=4)
    
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xlabel('x', fontsize=18)
    ax.set_ylabel('y', fontsize=18)
    ax.set_title(f'High-Precision Hessian Determinant Zero Contour\n' + 
                 f'Precision: {mp.dps} decimal places (~{int(mp.dps * 3.32)} bits)',
                 fontsize=20, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_aspect('equal')
    
    # Colorbar
    cbar = plt.colorbar(contourf, ax=ax, label='det(H)', pad=0.02)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('det(H)', fontsize=16)
    
    # Add text annotation
    ax.text(0.02, 0.98, 
            f'Zero contour: det(H) = 0\n' +
            f'Grid: {X.shape[0]}×{X.shape[1]}\n' +
            f'Precision: {mp.dps} digits',
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    filename = 'dumps/hessian_det_high_precision.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {filename}")
    plt.close()

def main():
    # Compute high-precision grid
    resolution = 1000  # 1000x1000 grid for smooth contour
    X, Y, det_H = compute_hessian_grid_high_precision(resolution=resolution)
    
    # Visualize
    visualize_high_precision_contour(X, Y, det_H)
    
    print("\n" + "=" * 70)
    print("Summary:")
    print("=" * 70)
    print(f"Precision: {mp.dps} decimal places")
    print(f"Grid resolution: {resolution}×{resolution}")
    print(f"Total evaluations: {resolution * resolution:,}")
    print(f"\nHessian determinant range:")
    print(f"  Min: {np.min(det_H):.6e}")
    print(f"  Max: {np.max(det_H):.6e}")
    print(f"\nThe zero contour shows the boundary between:")
    print(f"  • det(H) > 0: Locally convex or concave regions")
    print(f"  • det(H) < 0: Saddle point regions")

if __name__ == '__main__':
    main()

