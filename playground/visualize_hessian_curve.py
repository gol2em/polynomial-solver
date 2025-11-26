#!/usr/bin/env python3
"""
Construct the Hessian determinant zero curve using contour methods.
Compare with the polynomial solver results.
"""

import numpy as np
import matplotlib.pyplot as plt

def f_original(x, y):
    """Original function on [-1,1]^2."""
    f1 = 10*x**2*y**2 + np.sqrt(np.abs(x**2*y))
    f2 = np.arctan2(0.001, np.sin(5*y) - 2*x)
    f3 = 10*y**3 + x**3
    f4 = np.arctan2(0.01, np.sin(5*x) - 2*y)
    return f1 + f2 + f3 + f4

def compute_hessian_determinant_numerical(x, y, h=1e-6):
    """
    Compute Hessian determinant using numerical derivatives.
    
    Hessian matrix H = [[f_xx, f_xy],
                        [f_xy, f_yy]]
    
    det(H) = f_xx * f_yy - f_xy^2
    
    Second derivatives using finite differences:
    - f_xx = (f(x+h,y) - 2*f(x,y) + f(x-h,y)) / h^2
    - f_yy = (f(x,y+h) - 2*f(x,y) + f(x,y-h)) / h^2
    - f_xy = (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)) / (4*h^2)
    """
    # Second derivative with respect to x
    f_xx = (f_original(x+h, y) - 2*f_original(x, y) + f_original(x-h, y)) / (h*h)
    
    # Second derivative with respect to y
    f_yy = (f_original(x, y+h) - 2*f_original(x, y) + f_original(x, y-h)) / (h*h)
    
    # Mixed derivative (cross derivative)
    f_xy = (f_original(x+h, y+h) - f_original(x+h, y-h) - 
            f_original(x-h, y+h) + f_original(x-h, y-h)) / (4*h*h)
    
    # Determinant of Hessian
    det_H = f_xx * f_yy - f_xy * f_xy
    
    return det_H, f_xx, f_yy, f_xy

def compute_hessian_determinant_grid(resolution=500):
    """Compute Hessian determinant on a grid."""
    x = np.linspace(-1, 1, resolution)
    y = np.linspace(-1, 1, resolution)
    X, Y = np.meshgrid(x, y)

    print(f"Computing Hessian determinant on {resolution}x{resolution} grid...")
    print(f"  Total points: {resolution * resolution:,}")

    # Vectorized computation
    det_H = np.zeros_like(X)
    total = resolution
    progress_interval = max(1, resolution // 20)  # Show 20 progress updates

    for i in range(resolution):
        if i % progress_interval == 0:
            percent = 100.0 * i / resolution
            print(f"  Progress: {i}/{resolution} ({percent:.1f}%)")
        for j in range(resolution):
            det_H[i, j], _, _, _ = compute_hessian_determinant_numerical(X[i, j], Y[i, j])

    print("  Progress: 100.0% - Done!")
    return X, Y, det_H

def load_solver_results():
    """Load results from polynomial solver."""
    converged = []
    unresolved = []
    
    try:
        with open('dumps/hessian_det_converged.txt', 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    u, v = float(parts[0]), float(parts[1])
                    # Transform to [-1,1]^2
                    x, y = 2*u - 1, 2*v - 1
                    converged.append([x, y])
    except FileNotFoundError:
        print("Warning: converged points file not found")
    
    try:
        with open('dumps/hessian_det_boxes.txt', 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    u_min, u_max, v_min, v_max = [float(x) for x in parts[:4]]
                    u_center = (u_min + u_max) / 2
                    v_center = (v_min + v_max) / 2
                    # Transform to [-1,1]^2
                    x, y = 2*u_center - 1, 2*v_center - 1
                    unresolved.append([x, y])
    except FileNotFoundError:
        print("Warning: unresolved boxes file not found")
    
    return np.array(converged) if converged else np.array([]).reshape(0, 2), \
           np.array(unresolved) if unresolved else np.array([]).reshape(0, 2)

def main():
    print("=" * 60)
    print("Hessian Determinant Curve Construction")
    print("=" * 60)
    
    # Show how Hessian is computed at a test point
    print("\n1. Hessian Computation Example at (x=0, y=0):")
    print("-" * 60)
    x_test, y_test = 0.0, 0.0
    det_H, f_xx, f_yy, f_xy = compute_hessian_determinant_numerical(x_test, y_test)
    print(f"Point: (x={x_test}, y={y_test})")
    print(f"f(x,y) = {f_original(x_test, y_test):.6f}")
    print(f"\nSecond derivatives:")
    print(f"  f_xx = ∂²f/∂x² = {f_xx:.6f}")
    print(f"  f_yy = ∂²f/∂y² = {f_yy:.6f}")
    print(f"  f_xy = ∂²f/∂x∂y = {f_xy:.6f}")
    print(f"\nHessian matrix H = [[{f_xx:.4f}, {f_xy:.4f}],")
    print(f"                     [{f_xy:.4f}, {f_yy:.4f}]]")
    print(f"\ndet(H) = f_xx * f_yy - f_xy² = {f_xx:.4f} * {f_yy:.4f} - ({f_xy:.4f})²")
    print(f"       = {f_xx * f_yy:.4f} - {f_xy * f_xy:.4f}")
    print(f"       = {det_H:.4f}")
    
    # Compute on grid with higher resolution for smooth curve
    print(f"\n2. Computing Hessian determinant on grid...")
    print("-" * 60)
    X, Y, det_H_grid = compute_hessian_determinant_grid(resolution=800)
    
    print(f"\nHessian determinant statistics:")
    print(f"  Min: {np.min(det_H_grid):.2e}")
    print(f"  Max: {np.max(det_H_grid):.2e}")
    print(f"  Mean: {np.mean(det_H_grid):.2e}")

    # Load solver results
    print(f"\n3. Loading polynomial solver results...")
    print("-" * 60)
    converged, unresolved = load_solver_results()
    print(f"  Converged points: {len(converged)}")
    print(f"  Unresolved points: {len(unresolved)}")

    # Create visualization
    print(f"\n4. Creating visualization...")
    print("-" * 60)

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))

    # Left plot: Contour plot with zero level curve
    ax1 = axes[0]

    # Plot filled contours for context
    levels = np.linspace(np.min(det_H_grid), np.max(det_H_grid), 100)
    contourf = ax1.contourf(X, Y, det_H_grid, levels=levels, cmap='RdBu_r', alpha=0.6)

    # Plot zero level curve in black (the actual zero set) with higher smoothness
    contour_zero = ax1.contour(X, Y, det_H_grid, levels=[0], colors='black', linewidths=4, antialiased=True)

    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax1.set_xlabel('x', fontsize=14)
    ax1.set_ylabel('y', fontsize=14)
    ax1.set_title('Hessian Determinant Zero Curve (Contour Method)', fontsize=16, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')

    # Add colorbar
    cbar = plt.colorbar(contourf, ax=ax1)
    cbar.set_label('det(H)', fontsize=12)

    # Right plot: Comparison with solver results
    ax2 = axes[1]

    # Plot contour zero curve with smooth antialiasing
    ax2.contour(X, Y, det_H_grid, levels=[0], colors='black', linewidths=3, antialiased=True, label='Contour (ground truth)')

    # Plot solver results
    if len(converged) > 0:
        ax2.scatter(converged[:, 0], converged[:, 1], c='lime', s=3, alpha=0.6,
                   label=f'Solver converged ({len(converged)})', edgecolors='none')

    if len(unresolved) > 0:
        ax2.scatter(unresolved[:, 0], unresolved[:, 1], c='orange', s=3, alpha=0.6,
                   label=f'Solver unresolved ({len(unresolved)})', edgecolors='none')

    ax2.set_xlim(-1, 1)
    ax2.set_ylim(-1, 1)
    ax2.set_xlabel('x', fontsize=14)
    ax2.set_ylabel('y', fontsize=14)
    ax2.set_title('Comparison: Contour vs Polynomial Solver', fontsize=16, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right', fontsize=11)
    ax2.set_aspect('equal')

    plt.tight_layout()
    plt.savefig('dumps/hessian_det_curve_comparison.png', dpi=250, bbox_inches='tight')
    print(f"  ✅ Saved: dumps/hessian_det_curve_comparison.png (high resolution)")

    # Create a detailed view of just the contour curve
    fig2, ax = plt.subplots(1, 1, figsize=(14, 14))

    # Plot with better color scheme and more levels for smoothness
    contourf2 = ax.contourf(X, Y, det_H_grid, levels=100, cmap='RdBu_r', alpha=0.7)
    contour_zero2 = ax.contour(X, Y, det_H_grid, levels=[0], colors='black', linewidths=5, antialiased=True)

    # Add labels to show regions
    ax.text(0.5, 0.5, 'det(H) > 0\n(convex/concave)', fontsize=12, ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(-0.5, -0.5, 'det(H) < 0\n(saddle)', fontsize=12, ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xlabel('x', fontsize=16)
    ax.set_ylabel('y', fontsize=16)
    ax.set_title('Hessian Determinant Zero Set Curve\n(Separates Convex/Concave from Saddle Regions)',
                 fontsize=18, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    cbar2 = plt.colorbar(contourf2, ax=ax)
    cbar2.set_label('det(H)', fontsize=14)

    plt.savefig('dumps/hessian_det_curve_detailed.png', dpi=250, bbox_inches='tight')
    print(f"  ✅ Saved: dumps/hessian_det_curve_detailed.png (high resolution)")

    # Also create an ultra-high resolution version for publication quality
    plt.savefig('dumps/hessian_det_curve_detailed_hires.png', dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: dumps/hessian_det_curve_detailed_hires.png (ultra-high resolution)")

    print("\n" + "=" * 60)
    print("Summary:")
    print("=" * 60)
    print("The Hessian determinant zero curve has been constructed using:")
    print("  • Numerical finite differences for second derivatives")
    print("  • Matplotlib contour method to find zero level set")
    print("  • Comparison with polynomial solver results")
    print("\nThe curve separates regions where:")
    print("  • det(H) > 0: Function is locally convex or concave")
    print("  • det(H) < 0: Function has saddle point behavior")
    print("  • det(H) = 0: Boundary between these regions")

if __name__ == '__main__':
    main()

