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

def evaluate_polynomial_circle(X, Y):
    """Evaluate x^2 + y^2 - 1"""
    return X**2 + Y**2 - 1.0

def plot_single_poly_2d(ax, poly_func, box, title, control_points_dir0, control_points_dir1,
                        hull_dir0, hull_dir1, intersection_dir0, intersection_dir1,
                        interval_dir0=None, interval_dir1=None, view_box=None, z_limit=None,
                        control_points_3d=None):
    """Plot 3D graph of polynomial with control points and projections.

    This matches the format of plot_3d_graph from visualize_solver.py exactly.

    Args:
        poly_func: Function to evaluate polynomial (e.g., lambda X, Y: X**2 + Y**2 - 1)
        box: Current bounding box [x_min, x_max, y_min, y_max]
        title: Plot title
        control_points_dir0: Projected points in direction 0 (x-direction)
        control_points_dir1: Projected points in direction 1 (y-direction)
        hull_dir0: Convex hull in direction 0
        hull_dir1: Convex hull in direction 1
        intersection_dir0: Intersection points in direction 0
        intersection_dir1: Intersection points in direction 1
        interval_dir0: Bounding interval in direction 0
        interval_dir1: Bounding interval in direction 1
        view_box: Viewing scope (if None, uses box)
        z_limit: Z-axis limit (if None, computed from control points)
        control_points_3d: 3D control points (Bernstein coefficients) in local [0,1]^2 space
    """
    # Box format from solver: [x_min, x_max, y_min, y_max]
    x_min, x_max, y_min, y_max = box[0], box[1], box[2], box[3]

    # Use view_box for visualization scope
    if view_box:
        vis_x_min, vis_x_max = view_box[0], view_box[1]
        vis_y_min, vis_y_max = view_box[2], view_box[3]
    else:
        vis_x_min, vis_x_max = x_min, x_max
        vis_y_min, vis_y_max = y_min, y_max

    # Create mesh for surface over visible domain
    x = np.linspace(vis_x_min, vis_x_max, 50)
    y = np.linspace(vis_y_min, vis_y_max, 50)
    X, Y = np.meshgrid(x, y)
    Z = poly_func(X, Y)

    # Plot surface
    surf = ax.plot_surface(X, Y, Z, alpha=0.3, cmap='coolwarm', edgecolor='none')

    # Plot z=0 plane
    Z_plane = np.zeros_like(X)
    ax.plot_surface(X, Y, Z_plane, alpha=0.2, color='gray')

    # Plot zero contour (intersection with z=0)
    contour = ax.contour(X, Y, Z, levels=[0], colors='blue', linewidths=3)

    # Plot 3D control points (Bernstein coefficients)
    # These are in local [0,1]^2 space, need to transform to global [x_min, x_max] x [y_min, y_max]
    if control_points_3d:
        pts_3d = np.array(control_points_3d)
        # Transform from local [0,1]^2 to global box coordinates
        pts_global = pts_3d.copy()
        pts_global[:, 0] = x_min + pts_3d[:, 0] * (x_max - x_min)
        pts_global[:, 1] = y_min + pts_3d[:, 1] * (y_max - y_min)
        # pts_global[:, 2] is the function value (z), already in correct space
        ax.scatter(pts_global[:, 0], pts_global[:, 1], pts_global[:, 2],
                  c='green', s=50, marker='o', edgecolors='darkgreen', linewidths=1.5,
                  alpha=0.8, label='Control points', zorder=5)

    # Plot projections on background planes
    # Direction 0 projects along x-axis: (x, f(x,y)) - plot on x-z plane at y=vis_y_max (back wall)
    if control_points_dir0:
        pts = np.array(control_points_dir0)
        # Transform from [0,1] to [x_min, x_max]
        pts_transformed = pts.copy()
        pts_transformed[:, 0] = x_min + pts[:, 0] * (x_max - x_min)
        ax.scatter(pts_transformed[:, 0], np.full(len(pts_transformed), vis_y_max), pts_transformed[:, 1],
                  c='orange', s=20, alpha=0.6, label='X-proj points')
        if hull_dir0:
            hull_pts = np.array(hull_dir0)
            hull_pts_transformed = hull_pts.copy()
            hull_pts_transformed[:, 0] = x_min + hull_pts[:, 0] * (x_max - x_min)
            hull_pts_closed = np.vstack([hull_pts_transformed, hull_pts_transformed[0]])
            ax.plot(hull_pts_closed[:, 0], np.full(len(hull_pts_closed), vis_y_max), hull_pts_closed[:, 1],
                   '-', color='orange', linewidth=2, label='X-proj hull')
        if intersection_dir0:
            int_pts = np.array(intersection_dir0)
            int_pts_transformed = int_pts.copy()
            int_pts_transformed[:, 0] = x_min + int_pts[:, 0] * (x_max - x_min)
            ax.scatter(int_pts_transformed[:, 0], np.full(len(int_pts_transformed), vis_y_max), int_pts_transformed[:, 1],
                      s=150, c='red', marker='o', edgecolors='darkred', linewidths=3,
                      alpha=0.8, label='X-axis intersect', zorder=10)
            for pt in int_pts_transformed:
                ax.plot([pt[0], pt[0]], [vis_y_max, vis_y_max], [pt[1], 0],
                       'r--', linewidth=2, alpha=0.7)

    # Direction 1 projects along y-axis: (y, f(x,y)) - plot on y-z plane at x=vis_x_min (left wall)
    if control_points_dir1:
        pts = np.array(control_points_dir1)
        # Transform from [0,1] to [y_min, y_max]
        pts_transformed = pts.copy()
        pts_transformed[:, 0] = y_min + pts[:, 0] * (y_max - y_min)
        ax.scatter(np.full(len(pts_transformed), vis_x_min), pts_transformed[:, 0], pts_transformed[:, 1],
                  c='cyan', s=20, alpha=0.6, label='Y-proj points')
        if hull_dir1:
            hull_pts = np.array(hull_dir1)
            hull_pts_transformed = hull_pts.copy()
            hull_pts_transformed[:, 0] = y_min + hull_pts[:, 0] * (y_max - y_min)
            hull_pts_closed = np.vstack([hull_pts_transformed, hull_pts_transformed[0]])
            ax.plot(np.full(len(hull_pts_closed), vis_x_min), hull_pts_closed[:, 0], hull_pts_closed[:, 1],
                   '-', color='cyan', linewidth=2, label='Y-proj hull')
        if intersection_dir1:
            int_pts = np.array(intersection_dir1)
            int_pts_transformed = int_pts.copy()
            int_pts_transformed[:, 0] = y_min + int_pts[:, 0] * (y_max - y_min)
            ax.scatter(np.full(len(int_pts_transformed), vis_x_min), int_pts_transformed[:, 0], int_pts_transformed[:, 1],
                      s=150, c='red', marker='o', edgecolors='darkred', linewidths=3,
                      alpha=0.8, label='Y-axis intersect', zorder=10)
            for pt in int_pts_transformed:
                ax.plot([vis_x_min, vis_x_min], [pt[0], pt[0]], [pt[1], 0],
                       'r--', linewidth=2, alpha=0.7)

    # Plot current box on z=0 plane
    corners = [
        [x_min, y_min, 0], [x_max, y_min, 0],
        [x_max, y_max, 0], [x_min, y_max, 0],
        [x_min, y_min, 0]
    ]
    corners = np.array(corners)
    ax.plot(corners[:, 0], corners[:, 1], corners[:, 2], 'k-', linewidth=2, label='Current box')

    # Highlight bounding intervals on axes (z=0 plane)
    if interval_dir0:
        x_int_min = x_min + interval_dir0[0] * (x_max - x_min)
        x_int_max = x_min + interval_dir0[1] * (x_max - x_min)
        ax.plot([x_int_min, x_int_max], [vis_y_max, vis_y_max], [0, 0],
               'r-', linewidth=5, alpha=0.8, label=f'X-interval [{x_int_min:.3f}, {x_int_max:.3f}]')
        ax.scatter([x_int_min, x_int_max], [vis_y_max, vis_y_max], [0, 0],
                  s=100, c='red', marker='s', edgecolors='darkred', linewidths=2, zorder=10)

    if interval_dir1:
        y_int_min = y_min + interval_dir1[0] * (y_max - y_min)
        y_int_max = y_min + interval_dir1[1] * (y_max - y_min)
        ax.plot([vis_x_min, vis_x_min], [y_int_min, y_int_max], [0, 0],
               'r-', linewidth=5, alpha=0.8, label=f'Y-interval [{y_int_min:.3f}, {y_int_max:.3f}]')
        ax.scatter([vis_x_min, vis_x_min], [y_int_min, y_int_max], [0, 0],
                  s=100, c='red', marker='s', edgecolors='darkred', linewidths=2, zorder=10)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)
    ax.set_xlim(vis_x_min, vis_x_max)
    ax.set_ylim(vis_y_min, vis_y_max)

    # Set z limits
    if z_limit is not None:
        ax.set_zlim(-z_limit, z_limit)
    else:
        # Compute from control points
        z_values = []
        if control_points_dir0:
            z_values.extend([pt[1] for pt in control_points_dir0])
        if control_points_dir1:
            z_values.extend([pt[1] for pt in control_points_dir1])

        if z_values:
            max_abs_z = max(abs(z) for z in z_values)
            computed_z_limit = max(max_abs_z * 1.2, 0.1)
            ax.set_zlim(-computed_z_limit, computed_z_limit)
        else:
            z_min = min(Z.min(), -1.5)
            z_max = max(Z.max(), 1.5)
            ax.set_zlim(z_min, z_max)

    # Add legend
    ax.legend(loc='upper left', fontsize=7, framealpha=0.8)

def visualize_single_poly_2d(dump_file, output_dir='viz_output', max_steps=None):
    """Visualize single polynomial 2D solving process."""
    iterations = parse_dump_file(dump_file)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Visualizing {len(iterations)} iterations...")

    # Compute global z-limit from all control points
    z_limit = 1.0
    for iteration in iterations:
        if iteration['directions']:
            for direction in iteration['directions']:
                if direction['equations']:
                    for equation in direction['equations']:
                        if 'projected_points' in equation:
                            pts = equation['projected_points']
                            for pt in pts:
                                z_limit = max(z_limit, abs(pt[1]))
    z_limit *= 1.2  # Add 20% margin

    # Track previous box for view_box
    prev_box = None

    for i, iteration in enumerate(iterations):
        if max_steps and i >= max_steps:
            break

        iter_num = iteration['iter']
        depth = iteration['depth']
        decision = iteration.get('decision', 'UNKNOWN')
        global_box = iteration['global_box']

        # Get data from both directions
        if not iteration['directions'] or len(iteration['directions']) < 2:
            print(f"  Skipping iteration {iter_num}: insufficient directions")
            continue

        dir0 = iteration['directions'][0]
        dir1 = iteration['directions'][1]

        if not dir0['equations'] or not dir1['equations']:
            print(f"  Skipping iteration {iter_num}: no equations")
            continue

        # Extract data for equation 0 from both directions
        dir0_eq0 = dir0['equations'][0]
        dir1_eq0 = dir1['equations'][0]

        control_points_3d = dir0_eq0.get('control_points_3d')
        if not control_points_3d:
            print(f"  Skipping iteration {iter_num}: no control points")
            continue

        # Extract all the data needed for plotting
        control_points_dir0 = dir0_eq0.get('projected_points')
        control_points_dir1 = dir1_eq0.get('projected_points')
        hull_dir0 = dir0_eq0.get('convex_hull')
        hull_dir1 = dir1_eq0.get('convex_hull')
        intersection_dir0 = dir0_eq0.get('intersection')
        intersection_dir1 = dir1_eq0.get('intersection')
        interval_dir0 = dir0_eq0.get('interval')
        interval_dir1 = dir1_eq0.get('interval')

        # Use previous box as view_box (to show context)
        view_box = prev_box if prev_box else global_box

        # Create figure
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        title = f'Iteration {iter_num} (Depth {depth}): {decision}\nBox: [{global_box[0]:.3f}, {global_box[1]:.3f}] × [{global_box[2]:.3f}, {global_box[3]:.3f}]'

        # Call the plotting function with all parameters
        plot_single_poly_2d(ax, evaluate_polynomial_circle, global_box, title,
                           control_points_dir0, control_points_dir1,
                           hull_dir0, hull_dir1,
                           intersection_dir0, intersection_dir1,
                           interval_dir0, interval_dir1,
                           view_box=view_box, z_limit=z_limit,
                           control_points_3d=control_points_3d)

        filename = f'{output_dir}/step_{iter_num:03d}_depth_{depth}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  [{i+1}/{len(iterations)}] Saved: {filename}")

        # Update previous box
        prev_box = global_box

    print(f"\nVisualization complete! Output: {output_dir}/")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python visualize_single_poly_2d.py <dump_file> [output_dir] [max_steps]")
        sys.exit(1)
    
    dump_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'viz_output'
    max_steps = int(sys.argv[3]) if len(sys.argv) > 3 else None
    
    visualize_single_poly_2d(dump_file, output_dir, max_steps)

