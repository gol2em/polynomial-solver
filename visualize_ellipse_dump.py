#!/usr/bin/env python3
"""
Visualize the ellipse intersection solving process step by step.
Reads the geometry dump file and creates 3-subplot figures for each iteration.
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
                    # Start new iteration
                    parts = line.split(',')
                    iter_num = int(parts[0].split()[2])
                    depth = int(parts[1].split()[1])
                    current_iter = {
                        'iteration': iter_num,
                        'depth': depth,
                        'decision': 'CONTRACT',
                        'global_box': None,
                        'directions': []
                    }
                    iterations.append(current_iter)
                elif line.startswith('# Decision:') and current_iter:
                    # Parse decision
                    decision_str = line.split(':', 1)[1].strip()
                    current_iter['decision'] = decision_str
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
                state = 'hull'
                count = int(line.split()[1])
            elif line.startswith('Intersection'):
                state = 'intersection'
                count = int(line.split()[1])
            elif line.startswith('Interval'):
                if current_eq:
                    interval_str = line.split('[')[1].split(']')[0]
                    values = [float(x.strip()) for x in interval_str.split(',')]
                    current_eq['interval'] = values
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

    return iterations

def evaluate_polynomial_1(x, y):
    """Evaluate f1(x,y) = x^2 + y^2 - 1"""
    return x**2 + y**2 - 1

def evaluate_polynomial_2(x, y):
    """Evaluate f2(x,y) = x^2/4 + 4*y^2 - 1"""
    return x**2 / 4 + 4 * y**2 - 1

def create_mesh_grid(box):
    """Create mesh grid for the given box [x_min, x_max, y_min, y_max]"""
    x_min, x_max, y_min, y_max = box
    x = np.linspace(x_min, x_max, 50)
    y = np.linspace(y_min, y_max, 50)
    X, Y = np.meshgrid(x, y)
    return X, Y

def plot_3d_graph(ax, poly_func, box, title, control_points_dir0, control_points_dir1,
                  hull_dir0, hull_dir1, intersection_dir0, intersection_dir1, show_full_domain=False):
    """Plot 3D graph of polynomial with control points and projections."""
    x_min, x_max, y_min, y_max = box

    # For visualization, always show [0,1]^2 domain if requested
    if show_full_domain:
        vis_x_min, vis_x_max = 0.0, 1.0
        vis_y_min, vis_y_max = 0.0, 1.0
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

    # Plot projections on background planes
    # Direction 0 projects along x-axis: (x, f(x,y)) - plot on x-z plane at y=vis_y_min
    if control_points_dir0:
        pts = np.array(control_points_dir0)
        # pts[:, 0] is x coordinate, pts[:, 1] is function value
        ax.scatter(pts[:, 0], np.full(len(pts), vis_y_min), pts[:, 1],
                  c='orange', s=20, alpha=0.6)
        if hull_dir0:
            hull_pts = np.array(hull_dir0)
            # Close the hull for visualization
            hull_pts_closed = np.vstack([hull_pts, hull_pts[0]])
            ax.plot(hull_pts_closed[:, 0], np.full(len(hull_pts_closed), vis_y_min), hull_pts_closed[:, 1],
                   '-', color='orange', linewidth=2)
        if intersection_dir0:
            int_pts = np.array(intersection_dir0)
            ax.plot(int_pts[:, 0], np.full(len(int_pts), vis_y_min), int_pts[:, 1],
                   'o', color='red', markersize=8, markeredgewidth=2, markerfacecolor='none')

    # Direction 1 projects along y-axis: (y, f(x,y)) - plot on y-z plane at x=vis_x_min
    if control_points_dir1:
        pts = np.array(control_points_dir1)
        # pts[:, 0] is y coordinate, pts[:, 1] is function value
        ax.scatter(np.full(len(pts), vis_x_min), pts[:, 0], pts[:, 1],
                  c='cyan', s=20, alpha=0.6)
        if hull_dir1:
            hull_pts = np.array(hull_dir1)
            # Close the hull for visualization
            hull_pts_closed = np.vstack([hull_pts, hull_pts[0]])
            ax.plot(np.full(len(hull_pts_closed), vis_x_min), hull_pts_closed[:, 0], hull_pts_closed[:, 1],
                   '-', color='cyan', linewidth=2)
        if intersection_dir1:
            int_pts = np.array(intersection_dir1)
            ax.plot(np.full(len(int_pts), vis_x_min), int_pts[:, 0], int_pts[:, 1],
                   'o', color='red', markersize=8, markeredgewidth=2, markerfacecolor='none')

    # Plot current box on z=0 plane
    corners = [
        [x_min, y_min, 0], [x_max, y_min, 0],
        [x_max, y_max, 0], [x_min, y_max, 0],
        [x_min, y_min, 0]
    ]
    corners = np.array(corners)
    ax.plot(corners[:, 0], corners[:, 1], corners[:, 2], 'k-', linewidth=2)

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

def plot_3d_combined(ax, box, final_box, decision):
    """Plot both equations in 3D with surfaces, contours, and bounding boxes."""
    x_min, x_max, y_min, y_max = box

    # Always show [0,1]^2 domain
    vis_x_min, vis_x_max = 0.0, 1.0
    vis_y_min, vis_y_max = 0.0, 1.0

    # Create mesh over full domain
    x = np.linspace(vis_x_min, vis_x_max, 50)
    y = np.linspace(vis_y_min, vis_y_max, 50)
    X, Y = np.meshgrid(x, y)
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
        fx_min, fx_max, fy_min, fy_max = final_box
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

def visualize_iteration(iteration, output_dir='visualization_output'):
    """Create visualization for one iteration."""
    os.makedirs(output_dir, exist_ok=True)

    iter_num = iteration['iteration']
    depth = iteration['depth']
    decision = iteration.get('decision', 'CONTRACT')
    global_box = iteration['global_box']
    final_box = iteration.get('final_box', None)

    # Check if we have enough data
    if len(iteration['directions']) < 2:
        print(f"  Skipping iteration {iter_num}: not enough directions")
        return
    if len(iteration['directions'][0]['equations']) < 2:
        print(f"  Skipping iteration {iter_num}: not enough equations in direction 0")
        return
    if len(iteration['directions'][1]['equations']) < 2:
        print(f"  Skipping iteration {iter_num}: not enough equations in direction 1")
        return

    # Extract data for each equation and direction
    dir0_eq0 = iteration['directions'][0]['equations'][0]
    dir0_eq1 = iteration['directions'][0]['equations'][1]
    dir1_eq0 = iteration['directions'][1]['equations'][0]
    dir1_eq1 = iteration['directions'][1]['equations'][1]

    # Create figure with 3 subplots (all 3D)
    fig = plt.figure(figsize=(20, 6))
    fig.suptitle(f'Iteration {iter_num}: Depth {depth} - {decision}',
                fontsize=16, fontweight='bold')

    # Subplot 1: Equation 1
    ax1 = fig.add_subplot(131, projection='3d')
    plot_3d_graph(ax1, evaluate_polynomial_1, global_box,
                 f'Equation 1: x² + y² - 1 = 0',
                 dir0_eq0['projected_points'], dir1_eq0['projected_points'],
                 dir0_eq0['convex_hull'], dir1_eq0['convex_hull'],
                 dir0_eq0['intersection'], dir1_eq0['intersection'],
                 show_full_domain=True)

    # Subplot 2: Equation 2
    ax2 = fig.add_subplot(132, projection='3d')
    plot_3d_graph(ax2, evaluate_polynomial_2, global_box,
                 f'Equation 2: x²/4 + 4y² - 1 = 0',
                 dir0_eq1['projected_points'], dir1_eq1['projected_points'],
                 dir0_eq1['convex_hull'], dir1_eq1['convex_hull'],
                 dir0_eq1['intersection'], dir1_eq1['intersection'],
                 show_full_domain=True)

    # Subplot 3: Combined 3D view
    ax3 = fig.add_subplot(133, projection='3d')
    plot_3d_combined(ax3, global_box, final_box, decision)

    plt.tight_layout()

    # Save figure
    output_file = os.path.join(output_dir, f'step_{iter_num:03d}_depth_{depth}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved: {output_file}')
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize_ellipse_dump.py <dump_file>")
        sys.exit(1)

    dump_file = sys.argv[1]
    if not os.path.exists(dump_file):
        print(f"Error: File not found: {dump_file}")
        sys.exit(1)

    print(f"Parsing dump file: {dump_file}")
    iterations = parse_dump_file(dump_file)
    print(f"Found {len(iterations)} iterations")

    # Debug: print structure
    for i, iteration in enumerate(iterations):
        print(f"  Iteration {i}: depth={iteration['depth']}, "
              f"directions={len(iteration['directions'])}, "
              f"global_box={iteration.get('global_box', 'None')}")
        for j, direction in enumerate(iteration['directions']):
            print(f"    Direction {j}: equations={len(direction['equations'])}")

    print("\nGenerating visualizations...")
    for iteration in iterations:
        visualize_iteration(iteration)

    print(f"\nDone! Generated {len(iterations)} visualization(s)")
    print(f"Output directory: visualization_output/")

if __name__ == '__main__':
    main()


