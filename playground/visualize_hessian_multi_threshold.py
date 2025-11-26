#!/usr/bin/env python3
"""
Visualize Hessian determinant zero set with multiple thresholds.
Compares different residual thresholds side-by-side.
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

def hessian_determinant(u, v):
    """Compute Hessian determinant at (u,v) in [0,1]^2."""
    # Transform to [-1,1]^2
    x = 2*u - 1
    y = 2*v - 1
    
    # Original function components
    f1 = 10*x**2*y**2 + np.sqrt(np.abs(x**2*y))
    f2 = np.arctan2(0.001, np.sin(5*y) - 2*x)
    f3 = 10*y**3 + x**3
    f4 = np.arctan2(0.01, np.sin(5*x) - 2*y)
    
    # Compute numerical derivatives with finite differences
    h = 1e-6
    
    def f(u_val, v_val):
        x_val = 2*u_val - 1
        y_val = 2*v_val - 1
        f1_val = 10*x_val**2*y_val**2 + np.sqrt(np.abs(x_val**2*y_val))
        f2_val = np.arctan2(0.001, np.sin(5*y_val) - 2*x_val)
        f3_val = 10*y_val**3 + x_val**3
        f4_val = np.arctan2(0.01, np.sin(5*x_val) - 2*y_val)
        return f1_val + f2_val + f3_val + f4_val
    
    # Second derivatives
    fuu = (f(u+h, v) - 2*f(u, v) + f(u-h, v)) / (h*h)
    fvv = (f(u, v+h) - 2*f(u, v) + f(u, v-h)) / (h*h)
    fuv = (f(u+h, v+h) - f(u+h, v-h) - f(u-h, v+h) + f(u-h, v-h)) / (4*h*h)
    
    return fuu * fvv - fuv * fuv

def filter_points_by_threshold(converged_file, unresolved_file, threshold):
    """Filter points based on residual threshold."""
    # Read converged points
    converged = []
    if os.path.exists(converged_file):
        with open(converged_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    u, v = float(parts[0]), float(parts[1])
                    residual = abs(hessian_determinant(u, v))
                    if residual <= threshold:
                        converged.append([u, v])
    
    # Read unresolved boxes and use their centers
    unresolved = []
    if os.path.exists(unresolved_file):
        with open(unresolved_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    u_min, u_max, v_min, v_max = [float(x) for x in parts[:4]]
                    u_center = (u_min + u_max) / 2
                    v_center = (v_min + v_max) / 2
                    residual = abs(hessian_determinant(u_center, v_center))
                    if residual <= threshold:
                        unresolved.append([u_center, v_center])
    
    return np.array(converged) if converged else np.array([]).reshape(0, 2), \
           np.array(unresolved) if unresolved else np.array([]).reshape(0, 2)

def main():
    # Test multiple thresholds - can be adjusted
    # Very strict: [0.1, 0.5, 1, 2]
    # Strict: [0.5, 1, 5, 10]
    # Moderate: [1, 5, 10, 50]
    # Loose: [10, 50, 100, 500]
    thresholds = [0.5, 1, 5, 10]

    # Point size - can be adjusted (1=tiny, 3=small, 5=medium, 10=large)
    point_size = 5

    converged_file = 'dumps/hessian_det_converged.txt'
    unresolved_file = 'dumps/hessian_det_boxes.txt'

    # Create individual plots for each threshold
    for threshold in thresholds:
        print(f"\nProcessing threshold = {threshold}...")
        converged, unresolved = filter_points_by_threshold(
            converged_file, unresolved_file, threshold)

        total_points = len(converged) + len(unresolved)
        print(f"  Converged points: {len(converged)}")
        print(f"  Unresolved points: {len(unresolved)}")
        print(f"  Total points: {total_points}")

        # Create individual figure
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))

        # Plot both as scatter points with adjustable size
        if len(converged) > 0:
            # Transform to [-1,1]^2 for display
            x_conv = 2*converged[:, 0] - 1
            y_conv = 2*converged[:, 1] - 1
            ax.scatter(x_conv, y_conv, c='lime', s=point_size, alpha=0.7, label='Converged', edgecolors='none')

        if len(unresolved) > 0:
            x_unres = 2*unresolved[:, 0] - 1
            y_unres = 2*unresolved[:, 1] - 1
            ax.scatter(x_unres, y_unres, c='orange', s=point_size, alpha=0.7, label='Unresolved', edgecolors='none')

        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xlabel('x', fontsize=14)
        ax.set_ylabel('y', fontsize=14)
        ax.set_title(f'Hessian Determinant Zero Set (Threshold = {threshold}, {total_points} points)',
                     fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=12)
        ax.set_aspect('equal')

        # Save individual plot
        filename = f'dumps/hessian_det_threshold_{threshold}.png'
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        print(f"  ✅ Saved: {filename}")
        plt.close()

    print(f"\n✅ All individual plots saved!")

if __name__ == '__main__':
    main()

