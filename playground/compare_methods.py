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

def load_all_boxes(threshold=None):
    """Load unresolved boxes with optional threshold filtering."""
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
                    if threshold is None or residual <= threshold:
                        boxes.append([x, y, residual])
    except FileNotFoundError:
        print("Warning: boxes file not found")

    return np.array(boxes) if boxes else np.array([]).reshape(0, 3)

def load_converged_points(threshold=1.0):
    """Load converged points with strict threshold filtering."""
    points = []

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
                    residual = abs(compute_hessian_determinant(x, y))
                    if residual <= threshold:
                        points.append([x, y, residual])
    except FileNotFoundError:
        print("Warning: converged points file not found")

    return np.array(points) if points else np.array([]).reshape(0, 3)

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

    # Thresholds
    converged_threshold = 1.0  # Strict threshold for converged points
    unresolved_threshold = None  # No threshold for unresolved boxes (show all)

    print(f"\nThresholds:")
    print(f"  Converged points: residual <= {converged_threshold}")
    print(f"  Unresolved boxes: no threshold (show all)")

    # Compute contour once
    print("\n1. Computing ground truth contour...")
    print("-" * 70)
    X, Y, det_H = compute_hessian_grid(resolution=800)

    # Load filtered boxes and converged points
    print(f"\n2. Loading solver results...")
    print("-" * 70)

    boxes = load_all_boxes(threshold=unresolved_threshold)
    converged = load_converged_points(threshold=converged_threshold)

    print(f"  Converged points: {len(converged)}")
    print(f"  Unresolved boxes: {len(boxes)}")
    print(f"  Total: {len(converged) + len(boxes)}")

    if len(converged) > 0:
        print(f"\n  Converged residuals:")
        print(f"    Min: {np.min(converged[:, 2]):.2e}")
        print(f"    Max: {np.max(converged[:, 2]):.2e}")
        print(f"    Median: {np.median(converged[:, 2]):.2e}")

    if len(boxes) > 0:
        print(f"\n  Unresolved residuals:")
        print(f"    Min: {np.min(boxes[:, 2]):.2e}")
        print(f"    Max: {np.max(boxes[:, 2]):.2e}")
        print(f"    Median: {np.median(boxes[:, 2]):.2e}")

        # Show distribution
        percentiles = [50, 75, 90, 95, 99]
        print(f"\n  Unresolved percentiles:")
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

    # Right: Polynomial solver (converged + unresolved)
    ax2 = axes[1]

    # Plot converged points first (blue, larger, more visible)
    if len(converged) > 0:
        ax2.scatter(converged[:, 0], converged[:, 1], c='royalblue', s=8, alpha=0.7,
                   edgecolors='none', label=f'Converged (≤{converged_threshold}): {len(converged)}')

    # Plot unresolved boxes (orange, larger, more visible)
    if len(boxes) > 0:
        ax2.scatter(boxes[:, 0], boxes[:, 1], c='darkorange', s=8, alpha=0.7,
                   edgecolors='none', label=f'Unresolved (all): {len(boxes)}')

    ax2.set_xlim(-1, 1)
    ax2.set_ylim(-1, 1)
    ax2.set_xlabel('x', fontsize=16)
    ax2.set_ylabel('y', fontsize=16)
    total_points = len(converged) + len(boxes)
    ax2.set_title(f'Polynomial Solver Method\n{total_points} points',
                 fontsize=18, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right', fontsize=11)
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
    print(f"\nTotal points: {len(converged) + len(boxes)}")
    print(f"  Converged: {len(converged)}")
    print(f"  Unresolved: {len(boxes)}")

    if len(converged) > 0:
        print(f"\nConverged residual statistics:")
        print(f"  Min: {np.min(converged[:, 2]):.2e}")
        print(f"  Max: {np.max(converged[:, 2]):.2e}")
        print(f"  Median: {np.median(converged[:, 2]):.2e}")

    if len(boxes) > 0:
        print(f"\nUnresolved residual statistics:")
        print(f"  Min: {np.min(boxes[:, 2]):.2e}")
        print(f"  Max: {np.max(boxes[:, 2]):.2e}")
        print(f"  Median: {np.median(boxes[:, 2]):.2e}")
    print("\nThe polynomial solver method:")
    print("  • Uses piecewise polynomial interpolation")
    print("  • Subdivides domain for better local accuracy")
    print("  • Finds boxes containing the zero set")
    print(f"  • Converged points filtered: residual ≤ {converged_threshold}")
    print(f"  • Unresolved boxes: no filtering (all boxes shown)")
    print("\nThe contour method:")
    print("  • Direct numerical evaluation on 800×800 grid")
    print("  • Finds exact zero level curve")
    print("  • Serves as ground truth for comparison")
    print("\nNote: Run test_hessian_determinant with different subdivision")
    print("      parameters to generate results with varying accuracy.")

if __name__ == '__main__':
    main()

