#!/usr/bin/env python3
"""
Visualize the Hessian determinant zero set.
Shows the regions where det(Hessian) > 0, < 0, and = 0.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys

def parse_boxes(filename):
    """Parse boxes from result file."""
    boxes = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            coords = [float(x) for x in line.split()]
            if len(coords) == 4:
                boxes.append(coords)
    return boxes

def parse_samples(filename):
    """Parse Hessian determinant samples."""
    u_vals = []
    v_vals = []
    det_vals = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 3:
                u_vals.append(float(parts[0]))
                v_vals.append(float(parts[1]))
                det_vals.append(float(parts[2]))

    return np.array(u_vals), np.array(v_vals), np.array(det_vals)

def parse_converged_points(filename):
    """Parse converged points from file."""
    points = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    u, v = float(parts[0]), float(parts[1])
                    points.append((u, v))
    except FileNotFoundError:
        pass
    return points

def main():
    # Parse data
    print("Loading data...")
    boxes = parse_boxes('dumps/hessian_det_boxes.txt')
    converged_points = parse_converged_points('dumps/hessian_det_converged.txt')
    u_vals, v_vals, det_vals = parse_samples('dumps/hessian_det_samples.txt')

    print(f"Loaded {len(boxes)} unresolved boxes")
    print(f"Loaded {len(converged_points)} converged points")
    print(f"Loaded {len(det_vals)} samples")
    
    # Determine grid size
    n_samples = int(np.sqrt(len(det_vals))) - 1
    u_grid = u_vals.reshape((n_samples + 1, n_samples + 1))
    v_grid = v_vals.reshape((n_samples + 1, n_samples + 1))
    det_grid = det_vals.reshape((n_samples + 1, n_samples + 1))
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Left plot: Hessian determinant contour with zero set boxes
    ax1 = axes[0]
    
    # Plot contour of det(Hessian)
    levels = np.linspace(det_grid.min(), det_grid.max(), 50)
    contour = ax1.contourf(u_grid, v_grid, det_grid, levels=levels, cmap='RdBu_r', alpha=0.7)
    
    # Highlight zero contour
    zero_contour = ax1.contour(u_grid, v_grid, det_grid, levels=[0], colors='black', linewidths=2)
    ax1.clabel(zero_contour, inline=True, fontsize=10, fmt='det(H)=0')
    
    # Draw unresolved boxes
    for box in boxes:
        u_min, u_max, v_min, v_max = box
        width = u_max - u_min
        height = v_max - v_min
        rect = Rectangle((u_min, v_min), width, height,
                        linewidth=0.5, edgecolor='yellow', facecolor='none', alpha=0.8)
        ax1.add_patch(rect)

    # Draw converged points
    if converged_points:
        conv_u = [p[0] for p in converged_points]
        conv_v = [p[1] for p in converged_points]
        ax1.scatter(conv_u, conv_v, c='lime', s=20, marker='o',
                   edgecolors='black', linewidths=0.5, alpha=0.9, zorder=5, label='Converged')

    ax1.set_xlabel('u (transformed from x)', fontsize=12)
    ax1.set_ylabel('v (transformed from y)', fontsize=12)
    title = f'Hessian Determinant Zero Set\n{len(boxes)} unresolved boxes'
    if converged_points:
        title += f', {len(converged_points)} converged points'
    ax1.set_title(title + ' on [0,1]²', fontsize=14)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar1 = plt.colorbar(contour, ax=ax1)
    cbar1.set_label('det(Hessian)', fontsize=11)
    
    # Right plot: Same but in original coordinates [-1,1]²
    ax2 = axes[1]
    
    # Transform to original coordinates
    x_grid = 2 * u_grid - 1
    y_grid = 2 * v_grid - 1
    
    # Plot contour
    contour2 = ax2.contourf(x_grid, y_grid, det_grid, levels=levels, cmap='RdBu_r', alpha=0.7)
    
    # Highlight zero contour
    zero_contour2 = ax2.contour(x_grid, y_grid, det_grid, levels=[0], colors='black', linewidths=2)
    ax2.clabel(zero_contour2, inline=True, fontsize=10, fmt='det(H)=0')
    
    # Draw unresolved boxes in original coordinates
    for box in boxes:
        u_min, u_max, v_min, v_max = box
        x_min = 2 * u_min - 1
        x_max = 2 * u_max - 1
        y_min = 2 * v_min - 1
        y_max = 2 * v_max - 1
        width = x_max - x_min
        height = y_max - y_min
        rect = Rectangle((x_min, y_min), width, height,
                        linewidth=0.5, edgecolor='yellow', facecolor='none', alpha=0.8)
        ax2.add_patch(rect)

    # Draw converged points in original coordinates
    if converged_points:
        conv_x = [2 * p[0] - 1 for p in converged_points]
        conv_y = [2 * p[1] - 1 for p in converged_points]
        ax2.scatter(conv_x, conv_y, c='lime', s=20, marker='o',
                   edgecolors='black', linewidths=0.5, alpha=0.9, zorder=5, label='Converged')

    ax2.set_xlabel('x (original)', fontsize=12)
    ax2.set_ylabel('y (original)', fontsize=12)
    title2 = f'Hessian Determinant Zero Set\n{len(boxes)} unresolved boxes'
    if converged_points:
        title2 += f', {len(converged_points)} converged points'
    ax2.set_title(title2 + ' on [-1,1]²', fontsize=14)
    ax2.set_xlim(-1, 1)
    ax2.set_ylim(-1, 1)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar2 = plt.colorbar(contour2, ax=ax2)
    cbar2.set_label('det(Hessian)', fontsize=11)
    
    # Add interpretation text
    fig.text(0.5, 0.02, 
             'Red regions (det(H) > 0): Locally convex or concave | Blue regions (det(H) < 0): Saddle points | Black curve: Boundary',
             ha='center', fontsize=11, style='italic')
    
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    
    # Save
    output_file = 'dumps/hessian_det_visualization.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n✅ Saved: {output_file}")
    
    plt.show()

if __name__ == '__main__':
    main()

