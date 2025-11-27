#!/usr/bin/env python3
"""
Single Polynomial 2D Visualizer

Visualizes 2D single polynomial solving process from geometry dump files.
Shows only the polynomial surface with control points and bounding box.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
import os

def parse_dump_file(filename):
    """Parse the geometry dump file and extract iterations."""
    iterations = []
    current_iter = None
    current_dir = None
    current_eq = None
    state = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if line.startswith('# Iteration'):
                    parts = line.split(',')
                    iter_num = int(parts[0].split()[2])
                    depth = int(parts[1].split()[1])
                    current_iter = {'iter': iter_num, 'depth': depth, 'directions': []}
                    iterations.append(current_iter)
                elif line.startswith('# Decision:'):
                    if current_iter:
                        current_iter['decision'] = line.split(':', 1)[1].strip()
                elif line.startswith('# Global box:'):
                    if current_iter:
                        box_str = line.split(':', 1)[1].strip()
                        box_vals = [float(x) for x in box_str.strip('[]').split(',')]
                        current_iter['global_box'] = box_vals
                elif line.startswith('## Direction'):
                    dir_num = int(line.split()[2])
                    current_dir = {'direction': dir_num, 'equations': []}
                    if current_iter:
                        current_iter['directions'].append(current_dir)
                elif line.startswith('### Equation'):
                    eq_num = int(line.split()[2])
                    current_eq = {'equation': eq_num}
                    if current_dir:
                        current_dir['equations'].append(current_eq)
                continue

            if line.startswith('Control_Points_3D'):
                count = int(line.split()[1])
                state = ('control_points_3d', count)
                current_eq['control_points_3d'] = []
            elif line.startswith('Projected_Points'):
                count = int(line.split()[1])
                state = ('projected_points', count)
                current_eq['projected_points'] = []
            elif state:
                state_type, count = state
                if state_type == 'control_points_3d' and current_eq:
                    vals = [float(x) for x in line.split()]
                    current_eq['control_points_3d'].append(vals)
                    if len(current_eq['control_points_3d']) >= count:
                        state = None
                elif state_type == 'projected_points' and current_eq:
                    vals = [float(x) for x in line.split()]
                    current_eq['projected_points'].append(vals)
                    if len(current_eq['projected_points']) >= count:
                        state = None

    return iterations

def plot_single_poly_2d(ax, control_points_3d, global_box, title):
    """Plot a single 2D polynomial surface with control points."""
    # Extract x, y, z from control points
    pts = np.array(control_points_3d)
    x = pts[:, 0]
    y = pts[:, 1]
    z = pts[:, 2]
    
    # Reshape to grid (assuming 3x3 for degree 2)
    n = int(np.sqrt(len(pts)))
    X = x.reshape(n, n)
    Y = y.reshape(n, n)
    Z = z.reshape(n, n)
    
    # Plot surface
    ax.plot_surface(X, Y, Z, alpha=0.6, cmap='viridis', edgecolor='none')
    
    # Plot control points
    ax.scatter(x, y, z, c='red', s=50, marker='o', label='Control Points')
    
    # Plot f=0 plane
    x_range = [global_box[0], global_box[1]]
    y_range = [global_box[2], global_box[3]]
    xx, yy = np.meshgrid(x_range, y_range)
    zz = np.zeros_like(xx)
    ax.plot_surface(xx, yy, zz, alpha=0.2, color='gray')
    
    # Draw bounding box
    x0, x1, y0, y1 = global_box
    # Bottom rectangle (z=min(z))
    z_min = np.min(z)
    z_max = np.max(z)
    
    # Draw vertical edges of bounding box
    for xi in [x0, x1]:
        for yi in [y0, y1]:
            ax.plot([xi, xi], [yi, yi], [z_min, z_max], 'k-', linewidth=1, alpha=0.5)
    
    # Draw horizontal edges at bottom and top
    for zi in [z_min, z_max]:
        ax.plot([x0, x1], [y0, y0], [zi, zi], 'k-', linewidth=1, alpha=0.5)
        ax.plot([x0, x1], [y1, y1], [zi, zi], 'k-', linewidth=1, alpha=0.5)
        ax.plot([x0, x0], [y0, y1], [zi, zi], 'k-', linewidth=1, alpha=0.5)
        ax.plot([x1, x1], [y0, y1], [zi, zi], 'k-', linewidth=1, alpha=0.5)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('f(X,Y)')
    ax.set_title(title)
    ax.legend()

def visualize_single_poly_2d(dump_file, output_dir='viz_output', max_steps=None):
    """Visualize single polynomial 2D solving process."""
    iterations = parse_dump_file(dump_file)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"Visualizing {len(iterations)} iterations...")
    
    for i, iteration in enumerate(iterations):
        if max_steps and i >= max_steps:
            break
        
        iter_num = iteration['iter']
        depth = iteration['depth']
        decision = iteration.get('decision', 'UNKNOWN')
        global_box = iteration['global_box']
        
        # Get control points from direction 0, equation 0
        if not iteration['directions'] or not iteration['directions'][0]['equations']:
            print(f"  Skipping iteration {iter_num}: no data")
            continue
        
        control_points_3d = iteration['directions'][0]['equations'][0].get('control_points_3d')
        if not control_points_3d:
            print(f"  Skipping iteration {iter_num}: no control points")
            continue
        
        # Create figure
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        title = f'Iteration {iter_num} (Depth {depth}): {decision}\nBox: [{global_box[0]:.2f}, {global_box[1]:.2f}] × [{global_box[2]:.2f}, {global_box[3]:.2f}]'
        plot_single_poly_2d(ax, control_points_3d, global_box, title)
        
        filename = f'{output_dir}/step_{iter_num:03d}_depth_{depth}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [{i+1}/{len(iterations)}] Saved: {filename}")
    
    print(f"\nVisualization complete! Output: {output_dir}/")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python visualize_single_poly_2d.py <dump_file> [output_dir] [max_steps]")
        sys.exit(1)
    
    dump_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'viz_output'
    max_steps = int(sys.argv[3]) if len(sys.argv) > 3 else None
    
    visualize_single_poly_2d(dump_file, output_dir, max_steps)

