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
    # Strict: [0.1, 0.5, 1, 5]
    # Moderate: [1, 5, 10, 50]
    # Loose: [10, 50, 100, 500]
    thresholds = [0.5, 1, 5, 10]
    
    converged_file = 'dumps/hessian_det_converged.txt'
    unresolved_file = 'dumps/hessian_det_boxes.txt'
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    axes = axes.flatten()
    
    for idx, threshold in enumerate(thresholds):
        ax = axes[idx]
        
        print(f"\nProcessing threshold = {threshold}...")
        converged, unresolved = filter_points_by_threshold(
            converged_file, unresolved_file, threshold)
        
        total_points = len(converged) + len(unresolved)
        print(f"  Converged points: {len(converged)}")
        print(f"  Unresolved points: {len(unresolved)}")
        print(f"  Total points: {total_points}")
        
        # Plot both as scatter points
        if len(converged) > 0:
            # Transform to [-1,1]^2 for display
            x_conv = 2*converged[:, 0] - 1
            y_conv = 2*converged[:, 1] - 1
            ax.scatter(x_conv, y_conv, c='lime', s=1, alpha=0.6, label='Converged')
        
        if len(unresolved) > 0:
            x_unres = 2*unresolved[:, 0] - 1
            y_unres = 2*unresolved[:, 1] - 1
            ax.scatter(x_unres, y_unres, c='orange', s=1, alpha=0.6, label='Unresolved')
        
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('y', fontsize=12)
        ax.set_title(f'Threshold = {threshold} ({total_points} points)', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
        ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig('dumps/hessian_det_multi_threshold.png', dpi=150, bbox_inches='tight')
    print(f"\n✅ Saved: dumps/hessian_det_multi_threshold.png")

if __name__ == '__main__':
    main()

