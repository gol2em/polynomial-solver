#!/usr/bin/env python3
"""
Compare contour method with polynomial solver method.
Focus on unresolved boxes with reasonable threshold.
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

def compute_hessian_determinant(x, y, h=1e-6):
    """Compute Hessian determinant using finite differences."""
    f_xx = (f_original(x+h, y) - 2*f_original(x, y) + f_original(x-h, y)) / (h*h)
    f_yy = (f_original(x, y+h) - 2*f_original(x, y) + f_original(x, y-h)) / (h*h)
    f_xy = (f_original(x+h, y+h) - f_original(x+h, y-h) - 
            f_original(x-h, y+h) + f_original(x-h, y-h)) / (4*h*h)
    return f_xx * f_yy - f_xy * f_xy

def load_all_boxes():
    """Load all unresolved boxes without filtering."""
    boxes = []

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
                    residual = abs(compute_hessian_determinant(x, y))
                    boxes.append([x, y, residual])
    except FileNotFoundError:
        print("Warning: boxes file not found")

    return np.array(boxes) if boxes else np.array([]).reshape(0, 3)

def compute_hessian_grid(resolution=800):
    """Compute Hessian determinant on grid."""
    x = np.linspace(-1, 1, resolution)
    y = np.linspace(-1, 1, resolution)
    X, Y = np.meshgrid(x, y)
    
    print(f"Computing Hessian determinant on {resolution}x{resolution} grid...")
    det_H = np.zeros_like(X)
    progress_interval = max(1, resolution // 20)
    
    for i in range(resolution):
        if i % progress_interval == 0:
            percent = 100.0 * i / resolution
            print(f"  Progress: {percent:.1f}%")
        for j in range(resolution):
            det_H[i, j] = compute_hessian_determinant(X[i, j], Y[i, j])
    
    print("  Done!")
    return X, Y, det_H

def main():
    print("=" * 70)
    print("Comparing Contour Method vs Polynomial Solver")
    print("=" * 70)

    # Compute contour once
    print("\n1. Computing ground truth contour...")
    print("-" * 70)
    X, Y, det_H = compute_hessian_grid(resolution=800)

    # Load all boxes without filtering
    print(f"\n2. Loading all unresolved boxes...")
    print("-" * 70)

    boxes = load_all_boxes()
    print(f"  Total unresolved boxes: {len(boxes)}")

    if len(boxes) > 0:
        print(f"  Residual range: [{np.min(boxes[:, 2]):.2e}, {np.max(boxes[:, 2]):.2e}]")
        print(f"  Mean residual: {np.mean(boxes[:, 2]):.2e}")
        print(f"  Median residual: {np.median(boxes[:, 2]):.2e}")

        # Show distribution
        percentiles = [50, 75, 90, 95, 99]
        print(f"\n  Residual percentiles:")
        for p in percentiles:
            val = np.percentile(boxes[:, 2], p)
            print(f"    {p}th percentile: {val:.2e}")

    # Create side-by-side comparison (no overlay)
    print(f"\n3. Creating visualization...")
    print("-" * 70)

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))

    # Left: Contour method only
    ax1 = axes[0]
    contourf1 = ax1.contourf(X, Y, det_H, levels=100, cmap='RdBu_r', alpha=0.6)
    ax1.contour(X, Y, det_H, levels=[0], colors='black', linewidths=3)
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax1.set_xlabel('x', fontsize=16)
    ax1.set_ylabel('y', fontsize=16)
    ax1.set_title('Ground Truth (Contour Method)\n800×800 grid evaluation',
                 fontsize=18, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    cbar1 = plt.colorbar(contourf1, ax=ax1, label='det(H)')
    cbar1.ax.tick_params(labelsize=12)

    # Right: Polynomial solver only
    ax2 = axes[1]
    if len(boxes) > 0:
        # Use log scale for better visualization of wide residual range
        scatter = ax2.scatter(boxes[:, 0], boxes[:, 1], c=np.log10(boxes[:, 2] + 1e-10),
                             cmap='viridis', s=2, alpha=0.7, edgecolors='none')
        cbar2 = plt.colorbar(scatter, ax=ax2, label='log10(|det(H)|)')
        cbar2.ax.tick_params(labelsize=12)
    ax2.set_xlim(-1, 1)
    ax2.set_ylim(-1, 1)
    ax2.set_xlabel('x', fontsize=16)
    ax2.set_ylabel('y', fontsize=16)
    ax2.set_title(f'Polynomial Solver Method\n{len(boxes)} unresolved boxes',
                 fontsize=18, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')

    plt.tight_layout()
    filename = 'dumps/method_comparison_side_by_side.png'
    plt.savefig(filename, dpi=250, bbox_inches='tight')
    print(f"  ✅ Saved: {filename}")
    plt.close()

    print("\n" + "=" * 70)
    print("Summary:")
    print("=" * 70)
    print("Comparison complete!")
    print(f"\nTotal unresolved boxes: {len(boxes)}")
    if len(boxes) > 0:
        print(f"Residual statistics:")
        print(f"  Min: {np.min(boxes[:, 2]):.2e}")
        print(f"  Max: {np.max(boxes[:, 2]):.2e}")
        print(f"  Mean: {np.mean(boxes[:, 2]):.2e}")
        print(f"  Median: {np.median(boxes[:, 2]):.2e}")
    print("\nThe polynomial solver method:")
    print("  • Uses piecewise polynomial interpolation")
    print("  • Subdivides domain for better local accuracy")
    print("  • Finds boxes containing the zero set")
    print("  • Shows ALL unresolved boxes")
    print("\nThe contour method:")
    print("  • Direct numerical evaluation on 800×800 grid")
    print("  • Finds exact zero level curve")
    print("  • Serves as ground truth for comparison")
    print("\nNote: Run test_hessian_determinant with different subdivision")
    print("      parameters to generate results with varying accuracy.")

if __name__ == '__main__':
    main()

