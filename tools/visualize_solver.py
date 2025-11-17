#!/usr/bin/env python3
"""
Polynomial Solver Geometry Visualizer

Visualizes the polynomial solving process step by step from geometry dump files.
Creates 3-subplot figures for each iteration showing:
- Subplot 1: First equation with control points and projections
- Subplot 2: Second equation with control points and projections
- Subplot 3: Combined view with both equations and bounding boxes

Usage:
    python visualize_ellipse_dump.py <dump_file> [options]

Arguments:
    dump_file           Path to the geometry dump file (required)

Options:
    --max-steps N       Maximum number of iterations to visualize (default: all)
    --output-dir DIR    Output directory for PNG files (default: visualization_output)
    --help, -h          Show this help message

Examples:
    # Visualize all iterations
    python visualize_ellipse_dump.py strategy_ContractFirst_geometry.txt

    # Visualize first 20 iterations
    python visualize_ellipse_dump.py strategy_ContractFirst_geometry.txt --max-steps 20

    # Custom output directory
    python visualize_ellipse_dump.py dump.txt --max-steps 10 --output-dir my_viz
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
import os
import argparse

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
            if not line:
                continue

            # Check for section markers first (before single #)
            if line.startswith('### Equation'):
                eq_num = int(line.split()[2])
                current_eq = {
                    'equation': eq_num,
                    'projected_points': [],
                    'convex_hull': [],
                    'intersection': [],
                    'interval': None
                }
                if current_dir:
                    current_dir['equations'].append(current_eq)
                state = None
            elif line.startswith('## Direction'):
                dir_num = int(line.split()[2])
                current_dir = {
                    'direction': dir_num,
                    'equations': []
                }
                if current_iter:
                    current_iter['directions'].append(current_dir)
                else:
                    print(f"Warning: Direction {dir_num} found but no current iteration")
                current_eq = None
                state = None
            elif line.startswith('#'):
                if line.startswith('# Iteration'):
                    # Before starting new iteration, finalize previous iteration's decision
                    if current_iter and current_iter['decision'] is None:
                        # No explicit decision means it contracted and continued iterating
                        current_iter['decision'] = 'CONTRACTED (continuing iteration)'

                    # Start new iteration
                    parts = line.split(',')
                    iter_num = int(parts[0].split()[2])
                    depth = int(parts[1].split()[1])
                    current_iter = {
                        'iteration': iter_num,
                        'depth': depth,
                        'decision': None,
                        'global_box': None,
                        'directions': [],
                        'accumulated_z_norm': {}  # Will store accumulated factors per equation
                    }
                    iterations.append(current_iter)
                elif line.startswith('# Decision:') and current_iter:
                    # Parse decision (bounding strategy info)
                    decision_str = line.split(':', 1)[1].strip()
                    current_iter['bounding_info'] = decision_str
                elif line.startswith('# FINAL_DECISION:'):
                    # Parse final decision (what actually happened)
                    # This belongs to the current iteration if it has no decision yet,
                    # otherwise it belongs to the previous iteration (late SUBDIVIDE decision)
                    decision_str = line.split(':', 1)[1].strip()
                    if current_iter and current_iter['decision'] is None:
                        current_iter['decision'] = decision_str
                    elif len(iterations) > 1 and iterations[-2]['decision'] is None:
                        iterations[-2]['decision'] = decision_str
                elif line.startswith('# Global box:') and current_iter:
                    # Parse global box
                    box_str = line.split('[')[1].split(']')[0]
                    values = [float(x.strip()) for x in box_str.split(',')]
                    current_iter['global_box'] = values
                elif line.startswith('# Bounding Box (global):') and current_iter:
                    # Parse final bounding box
                    box_str = line.split('[')[1].split(']')[0]
                    values = [float(x.strip()) for x in box_str.split(',')]
                    current_iter['final_box'] = values
                # Skip other comment lines
            elif line.startswith('Projected_Points'):
                state = 'projected'
                count = int(line.split()[1])
            elif line.startswith('ConvexHull'):
                parts = line.split()
                if len(parts) > 1 and parts[1] == 'EMPTY':
                    state = None
                    count = 0
                else:
                    state = 'hull'
                    count = int(parts[1])
            elif line.startswith('Intersection'):
                parts = line.split()
                if len(parts) > 1 and parts[1] == 'EMPTY':
                    state = None
                    count = 0
                else:
                    state = 'intersection'
                    count = int(parts[1])
            elif line.startswith('Interval') or line.startswith('Direction_Interval'):
                if 'EMPTY' in line:
                    # Empty interval, skip
                    state = None
                elif current_eq and '[' in line:
                    interval_str = line.split('[')[1].split(']')[0]
                    values = [float(x.strip()) for x in interval_str.split(',')]
                    current_eq['interval'] = values
                    state = None
                else:
                    state = None
            elif line.startswith('Final_Interval'):
                state = None
            else:
                # Parse coordinate data
                parts = line.split()
                if len(parts) == 2 and state and current_eq:
                    try:
                        x, y = float(parts[0]), float(parts[1])
                        if state == 'projected':
                            current_eq['projected_points'].append([x, y])
                        elif state == 'hull':
                            current_eq['convex_hull'].append([x, y])
                        elif state == 'intersection':
                            current_eq['intersection'].append([x, y])
                    except ValueError:
                        pass

    # Finalize last iteration's decision if needed
    if current_iter and current_iter['decision'] is None:
        current_iter['decision'] = 'CONTRACTED (continuing iteration)'

    return iterations

def evaluate_polynomial_1(x, y):
    """Evaluate f1(x,y) = x^2 + y^2 - 1 (original polynomial)"""
    return x**2 + y**2 - 1

def evaluate_polynomial_2(x, y):
    """Evaluate f2(x,y) = x^2/4 + 4*y^2 - 1 (original polynomial)"""
    return x**2 / 4 + 4 * y**2 - 1

def create_mesh_grid(box):
    """Create mesh grid for the given box [x_min, x_max, y_min, y_max]"""
    x_min, x_max, y_min, y_max = box
    x = np.linspace(x_min, x_max, 50)
    y = np.linspace(y_min, y_max, 50)
    X, Y = np.meshgrid(x, y)
    return X, Y

def plot_3d_graph(ax, poly_func, box, title, control_points_dir0, control_points_dir1,
                  hull_dir0, hull_dir1, intersection_dir0, intersection_dir1,
                  interval_dir0=None, interval_dir1=None, view_box=None):
    """Plot 3D graph of polynomial with control points and projections.

    Args:
        box: Current bounding box to draw
        view_box: Viewing scope (if None, uses box)

    Note: Control points are already normalized in the dump file (z-values in [-1,1]).
          The polynomial surface is evaluated from the original polynomial (not normalized).
    """
    # Box format from solver: [x_min, x_max, y_min, y_max]
    x_min, x_max, y_min, y_max = box[0], box[1], box[2], box[3]

    # Use view_box for visualization scope (shows previous iteration's box)
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
    # Note: We evaluate the original polynomial, not the subdivided one
    # So we should NOT apply z_norm_factor here (it's for Bernstein coefficients on subdivided box)
    Z = poly_func(X, Y)

    # Plot surface
    surf = ax.plot_surface(X, Y, Z, alpha=0.3, cmap='coolwarm', edgecolor='none')

    # Plot z=0 plane
    Z_plane = np.zeros_like(X)
    ax.plot_surface(X, Y, Z_plane, alpha=0.2, color='gray')

    # Plot zero contour (intersection with z=0)
    contour = ax.contour(X, Y, Z, levels=[0], colors='blue', linewidths=3)

    # Plot projections on background planes
    # Direction 0 projects along x-axis: (x, f(x,y)) - plot on x-z plane at y=vis_y_max (back wall)
    # NOTE: Projected points are in normalized [0,1] coordinates, need affine transformation
    if control_points_dir0:
        pts = np.array(control_points_dir0)
        # Transform from [0,1] to [x_min, x_max]
        pts_transformed = pts.copy()
        pts_transformed[:, 0] = x_min + pts[:, 0] * (x_max - x_min)
        # pts[:, 1] is function value (already in correct space)
        ax.scatter(pts_transformed[:, 0], np.full(len(pts_transformed), vis_y_max), pts_transformed[:, 1],
                  c='orange', s=20, alpha=0.6, label='X-proj points')
        if hull_dir0:
            hull_pts = np.array(hull_dir0)
            # Transform from [0,1] to [x_min, x_max]
            hull_pts_transformed = hull_pts.copy()
            hull_pts_transformed[:, 0] = x_min + hull_pts[:, 0] * (x_max - x_min)
            # Close the hull for visualization
            hull_pts_closed = np.vstack([hull_pts_transformed, hull_pts_transformed[0]])
            ax.plot(hull_pts_closed[:, 0], np.full(len(hull_pts_closed), vis_y_max), hull_pts_closed[:, 1],
                   '-', color='orange', linewidth=2, label='X-proj hull')
        if intersection_dir0:
            int_pts = np.array(intersection_dir0)
            # Transform from [0,1] to [x_min, x_max]
            int_pts_transformed = int_pts.copy()
            int_pts_transformed[:, 0] = x_min + int_pts[:, 0] * (x_max - x_min)
            # Highlight intersection points with axis (z=0)
            ax.scatter(int_pts_transformed[:, 0], np.full(len(int_pts_transformed), vis_y_max), int_pts_transformed[:, 1],
                      s=150, c='red', marker='o', edgecolors='darkred', linewidths=3,
                      alpha=0.8, label='X-axis intersect', zorder=10)
            # Draw lines from intersection points to x-axis on z=0 plane
            for pt in int_pts_transformed:
                ax.plot([pt[0], pt[0]], [vis_y_max, vis_y_max], [pt[1], 0],
                       'r--', linewidth=2, alpha=0.7)

    # Direction 1 projects along y-axis: (y, f(x,y)) - plot on y-z plane at x=vis_x_min (left wall)
    # NOTE: Projected points are in normalized [0,1] coordinates, need affine transformation
    if control_points_dir1:
        pts = np.array(control_points_dir1)
        # Transform from [0,1] to [y_min, y_max]
        pts_transformed = pts.copy()
        pts_transformed[:, 0] = y_min + pts[:, 0] * (y_max - y_min)
        # pts[:, 1] is function value (already in correct space)
        ax.scatter(np.full(len(pts_transformed), vis_x_min), pts_transformed[:, 0], pts_transformed[:, 1],
                  c='cyan', s=20, alpha=0.6, label='Y-proj points')
        if hull_dir1:
            hull_pts = np.array(hull_dir1)
            # Transform from [0,1] to [y_min, y_max]
            hull_pts_transformed = hull_pts.copy()
            hull_pts_transformed[:, 0] = y_min + hull_pts[:, 0] * (y_max - y_min)
            # Close the hull for visualization
            hull_pts_closed = np.vstack([hull_pts_transformed, hull_pts_transformed[0]])
            ax.plot(np.full(len(hull_pts_closed), vis_x_min), hull_pts_closed[:, 0], hull_pts_closed[:, 1],
                   '-', color='cyan', linewidth=2, label='Y-proj hull')
        if intersection_dir1:
            int_pts = np.array(intersection_dir1)
            # Transform from [0,1] to [y_min, y_max]
            int_pts_transformed = int_pts.copy()
            int_pts_transformed[:, 0] = y_min + int_pts[:, 0] * (y_max - y_min)
            # Highlight intersection points with axis (z=0)
            ax.scatter(np.full(len(int_pts_transformed), vis_x_min), int_pts_transformed[:, 0], int_pts_transformed[:, 1],
                      s=150, c='red', marker='o', edgecolors='darkred', linewidths=3,
                      alpha=0.8, label='Y-axis intersect', zorder=10)
            # Draw lines from intersection points to y-axis on z=0 plane
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
    # NOTE: Intervals are in normalized [0,1] coordinates, need affine transformation
    if interval_dir0:
        # X-direction interval on x-axis (y=vis_y_max, z=0)
        # Transform from [0,1] to [x_min, x_max]
        x_int_min = x_min + interval_dir0[0] * (x_max - x_min)
        x_int_max = x_min + interval_dir0[1] * (x_max - x_min)
        ax.plot([x_int_min, x_int_max], [vis_y_max, vis_y_max], [0, 0],
               'r-', linewidth=5, alpha=0.8, label=f'X-interval [{x_int_min:.3f}, {x_int_max:.3f}]')
        # Mark endpoints
        ax.scatter([x_int_min, x_int_max], [vis_y_max, vis_y_max], [0, 0],
                  s=100, c='red', marker='s', edgecolors='darkred', linewidths=2, zorder=10)

    if interval_dir1:
        # Y-direction interval on y-axis (x=vis_x_min, z=0)
        # Transform from [0,1] to [y_min, y_max]
        y_int_min = y_min + interval_dir1[0] * (y_max - y_min)
        y_int_max = y_min + interval_dir1[1] * (y_max - y_min)
        ax.plot([vis_x_min, vis_x_min], [y_int_min, y_int_max], [0, 0],
               'r-', linewidth=5, alpha=0.8, label=f'Y-interval [{y_int_min:.3f}, {y_int_max:.3f}]')
        # Mark endpoints
        ax.scatter([vis_x_min, vis_x_min], [y_int_min, y_int_max], [0, 0],
                  s=100, c='red', marker='s', edgecolors='darkred', linewidths=2, zorder=10)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)
    ax.set_xlim(vis_x_min, vis_x_max)
    ax.set_ylim(vis_y_min, vis_y_max)

    # Set z limits to show the surface well
    z_min = min(Z.min(), -1.5)
    z_max = max(Z.max(), 1.5)
    ax.set_zlim(z_min, z_max)

    # Add legend
    ax.legend(loc='upper left', fontsize=7, framealpha=0.8)

def plot_3d_combined(ax, box, final_box, decision, view_box=None):
    """Plot both equations in 3D with surfaces, contours, and bounding boxes.

    Args:
        box: Current bounding box to draw
        final_box: PP bounds to draw
        view_box: Viewing scope (if None, uses [0,1]^2)
    """
    # Box format from solver: [x_min, x_max, y_min, y_max]
    x_min, x_max, y_min, y_max = box[0], box[1], box[2], box[3]

    # Use view_box for visualization scope
    if view_box:
        vis_x_min, vis_x_max = view_box[0], view_box[1]
        vis_y_min, vis_y_max = view_box[2], view_box[3]
    else:
        # Default to [0,1]^2 domain
        vis_x_min, vis_x_max = 0.0, 1.0
        vis_y_min, vis_y_max = 0.0, 1.0

    # Create mesh over full domain
    x = np.linspace(vis_x_min, vis_x_max, 50)
    y = np.linspace(vis_y_min, vis_y_max, 50)
    X, Y = np.meshgrid(x, y)
    # Evaluate original polynomials (not subdivided, so no normalization needed)
    Z1 = evaluate_polynomial_1(X, Y)
    Z2 = evaluate_polynomial_2(X, Y)

    # Plot both surfaces with transparency
    ax.plot_surface(X, Y, Z1, alpha=0.2, color='blue', edgecolor='none')
    ax.plot_surface(X, Y, Z2, alpha=0.2, color='red', edgecolor='none')

    # Plot z=0 plane
    Z_plane = np.zeros_like(X)
    ax.plot_surface(X, Y, Z_plane, alpha=0.15, color='gray')

    # Plot zero contours on z=0 plane
    ax.contour(X, Y, Z1, levels=[0], colors='blue', linewidths=3, offset=0)
    ax.contour(X, Y, Z2, levels=[0], colors='red', linewidths=3, offset=0)

    # Plot current box on z=0 plane
    corners = [
        [x_min, y_min, 0], [x_max, y_min, 0],
        [x_max, y_max, 0], [x_min, y_max, 0],
        [x_min, y_min, 0]
    ]
    corners = np.array(corners)
    ax.plot(corners[:, 0], corners[:, 1], corners[:, 2], 'k-', linewidth=3, label='Current box')

    # Plot final bounding box (PP bounds) if available
    if final_box:
        # Box format from solver: [x_min, x_max, y_min, y_max]
        fx_min, fx_max, fy_min, fy_max = final_box[0], final_box[1], final_box[2], final_box[3]
        pp_corners = [
            [fx_min, fy_min, 0], [fx_max, fy_min, 0],
            [fx_max, fy_max, 0], [fx_min, fy_max, 0],
            [fx_min, fy_min, 0]
        ]
        pp_corners = np.array(pp_corners)
        ax.plot(pp_corners[:, 0], pp_corners[:, 1], pp_corners[:, 2],
               color='purple', linewidth=3, linestyle='--', label='PP bounds')

    # Plot expected solution (analytical)
    # Solve: x^2 + y^2 = 1 and x^2/4 + 4*y^2 = 1
    y_sol = np.sqrt(0.2)
    x_sol = np.sqrt(0.8)
    ax.plot([x_sol], [y_sol], [0], 'k*', markersize=15,
            markeredgewidth=2, markerfacecolor='yellow', label='Expected solution')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(f'Both Equations + Bounds\n{decision}')
    ax.set_xlim(vis_x_min, vis_x_max)
    ax.set_ylim(vis_y_min, vis_y_max)

    # Set z limits
    z_min = min(Z1.min(), Z2.min(), -1.5)
    z_max = max(Z1.max(), Z2.max(), 1.5)
    ax.set_zlim(z_min, z_max)

    ax.legend(loc='upper left', fontsize=8)

def visualize_iteration(iteration, prev_iteration=None, output_dir='visualization_output'):
    """Create visualization for one iteration."""
    os.makedirs(output_dir, exist_ok=True)

    iter_num = iteration['iteration']
    depth = iteration['depth']
    bounding_info = iteration.get('bounding_info', '')
    decision = iteration.get('decision', 'UNKNOWN')
    global_box = iteration['global_box']
    final_box = iteration.get('final_box', None)

    # Use previous iteration's bounding box as the viewing scope (if available)
    # This shows the contraction more clearly
    if prev_iteration and prev_iteration.get('final_box'):
        view_box = prev_iteration['final_box']
    else:
        view_box = global_box

    # Check if this is a pruned case
    is_pruned = 'PRUNED' in decision

    # For pruned cases, we can visualize even with incomplete data
    # For non-pruned cases, we need complete data
    if not is_pruned:
        if len(iteration['directions']) < 2:
            print(f"  Skipping iteration {iter_num}: not enough directions")
            return
        if len(iteration['directions'][0]['equations']) < 2:
            print(f"  Skipping iteration {iter_num}: not enough equations in direction 0")
            return
        if len(iteration['directions'][1]['equations']) < 2:
            print(f"  Skipping iteration {iter_num}: not enough equations in direction 1")
            return

    # Create figure with 3 subplots (all 3D)
    fig = plt.figure(figsize=(20, 6))
    title = f'Iteration {iter_num} (Depth {depth}): {decision}'
    fig.suptitle(title, fontsize=16, fontweight='bold')

    # Handle pruned cases specially
    if is_pruned:
        # For pruned cases, show the box and indicate it was empty
        ax = fig.add_subplot(111, projection='3d')

        # Draw the current box
        x_min, x_max, y_min, y_max = global_box[0], global_box[1], global_box[2], global_box[3]

        # Use view_box for plotting limits
        vis_x_min, vis_x_max = view_box[0], view_box[1]
        vis_y_min, vis_y_max = view_box[2], view_box[3]

        # Create mesh grid for the surfaces
        X, Y = create_mesh_grid(view_box)
        Z1 = evaluate_polynomial_1(X, Y)
        Z2 = evaluate_polynomial_2(X, Y)

        # Plot surfaces with transparency
        ax.plot_surface(X, Y, Z1, alpha=0.3, cmap='Blues', label='f1')
        ax.plot_surface(X, Y, Z2, alpha=0.3, cmap='Reds', label='f2')

        # Plot zero contours
        ax.contour(X, Y, Z1, levels=[0], colors='blue', linewidths=2, offset=0)
        ax.contour(X, Y, Z2, levels=[0], colors='red', linewidths=2, offset=0)

        # Draw the current box in RED with X marks to indicate pruned
        corners = [
            [x_min, y_min, 0], [x_max, y_min, 0],
            [x_max, y_max, 0], [x_min, y_max, 0],
            [x_min, y_min, 0]
        ]
        corners = np.array(corners)
        ax.plot(corners[:, 0], corners[:, 1], corners[:, 2], 'r-', linewidth=4, label='Pruned box')

        # Add X marks at corners
        for i in range(4):
            ax.scatter([corners[i, 0]], [corners[i, 1]], [corners[i, 2]],
                      c='red', marker='x', s=200, linewidths=4)

        # Add text annotation
        center_x = (x_min + x_max) / 2
        center_y = (y_min + y_max) / 2
        ax.text(center_x, center_y, 0, 'EMPTY\nBOUNDING BOX',
               fontsize=14, color='red', weight='bold',
               ha='center', va='center',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('Pruned Box (No Roots)')
        ax.set_xlim(vis_x_min, vis_x_max)
        ax.set_ylim(vis_y_min, vis_y_max)
        ax.legend(loc='upper left', fontsize=8)

        plt.tight_layout()

        # Save figure
        output_file = os.path.join(output_dir, f'step_{iter_num:03d}_depth_{depth}.png')
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        print(f"Saved: {output_file}")
        return

    # Extract data for each equation and direction (non-pruned case)
    dir0_eq0 = iteration['directions'][0]['equations'][0]
    dir0_eq1 = iteration['directions'][0]['equations'][1]
    dir1_eq0 = iteration['directions'][1]['equations'][0]
    dir1_eq1 = iteration['directions'][1]['equations'][1]

    # Subplot 1: Equation 1
    ax1 = fig.add_subplot(131, projection='3d')
    plot_3d_graph(ax1, evaluate_polynomial_1, global_box,
                 f'Equation 1: x² + y² - 1 = 0',
                 dir0_eq0['projected_points'], dir1_eq0['projected_points'],
                 dir0_eq0['convex_hull'], dir1_eq0['convex_hull'],
                 dir0_eq0['intersection'], dir1_eq0['intersection'],
                 dir0_eq0.get('interval'), dir1_eq0.get('interval'),
                 view_box=view_box)

    # Subplot 2: Equation 2
    ax2 = fig.add_subplot(132, projection='3d')
    plot_3d_graph(ax2, evaluate_polynomial_2, global_box,
                 f'Equation 2: x²/4 + 4y² - 1 = 0',
                 dir0_eq1['projected_points'], dir1_eq1['projected_points'],
                 dir0_eq1['convex_hull'], dir1_eq1['convex_hull'],
                 dir0_eq1['intersection'], dir1_eq1['intersection'],
                 dir0_eq1.get('interval'), dir1_eq1.get('interval'),
                 view_box=view_box)

    # Subplot 3: Combined 3D view
    ax3 = fig.add_subplot(133, projection='3d')
    plot_3d_combined(ax3, global_box, final_box, decision, view_box=view_box)

    plt.tight_layout()

    # Save figure
    output_file = os.path.join(output_dir, f'step_{iter_num:03d}_depth_{depth}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved: {output_file}')
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Visualize polynomial solver geometry dump files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Visualize all iterations
  python visualize_ellipse_dump.py strategy_ContractFirst_geometry.txt

  # Visualize first 20 iterations
  python visualize_ellipse_dump.py strategy_ContractFirst_geometry.txt --max-steps 20

  # Custom output directory
  python visualize_ellipse_dump.py dump.txt --max-steps 10 --output-dir my_viz
        '''
    )

    parser.add_argument('dump_file', help='Path to the geometry dump file')
    parser.add_argument('--max-steps', type=int, default=None, metavar='N',
                        help='Maximum number of iterations to visualize (default: all)')
    parser.add_argument('--output-dir', default='visualization_output', metavar='DIR',
                        help='Output directory for PNG files (default: visualization_output)')

    args = parser.parse_args()

    if not os.path.exists(args.dump_file):
        print(f"Error: File not found: {args.dump_file}")
        sys.exit(1)

    print(f"Parsing dump file: {args.dump_file}")
    iterations = parse_dump_file(args.dump_file)
    print(f"Found {len(iterations)} iterations")

    # Limit iterations if max_steps specified
    if args.max_steps is not None and args.max_steps < len(iterations):
        print(f"Limiting to first {args.max_steps} iterations")
        iterations = iterations[:args.max_steps]

    # Debug: print structure
    for i, iteration in enumerate(iterations):
        print(f"  Iteration {i}: depth={iteration['depth']}, "
              f"directions={len(iteration['directions'])}, "
              f"global_box={iteration.get('global_box', 'None')}")
        for j, direction in enumerate(iteration['directions']):
            print(f"    Direction {j}: equations={len(direction['equations'])}")

    print(f"\nGenerating visualizations to: {args.output_dir}/")
    for i, iteration in enumerate(iterations):
        # Pass previous iteration for viewing scope
        prev_iteration = iterations[i-1] if i > 0 else None
        visualize_iteration(iteration, prev_iteration=prev_iteration, output_dir=args.output_dir)

    print(f"\nDone! Generated {len(iterations)} visualization(s)")
    print(f"Output directory: {args.output_dir}/")

if __name__ == '__main__':
    main()


