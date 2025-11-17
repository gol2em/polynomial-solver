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
                        'global_box': None,
                        'directions': []
                    }
                    iterations.append(current_iter)
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
                  hull_dir0, hull_dir1, intersection_dir0, intersection_dir1):
    """Plot 3D graph of polynomial with control points and projections."""
    x_min, x_max, y_min, y_max = box
    
    # Create mesh for surface
    X, Y = create_mesh_grid(box)
    Z = poly_func(X, Y)
    
    # Plot surface
    surf = ax.plot_surface(X, Y, Z, alpha=0.3, cmap='coolwarm', edgecolor='none')
    
    # Plot z=0 plane
    Z_plane = np.zeros_like(X)
    ax.plot_surface(X, Y, Z_plane, alpha=0.2, color='gray')
    
    # Plot control points in 3D
    if control_points_dir0 and control_points_dir1:
        # Reconstruct 3D control points from 2D projections
        # This is approximate - we use the projected points
        pass
    
    # Plot projections on background planes
    # X-direction projection (on x-z plane at y=y_min)
    if control_points_dir0:
        pts = np.array(control_points_dir0)
        ax.scatter(pts[:, 0], np.full(len(pts), y_min), pts[:, 1],
                  c='orange', s=20, alpha=0.6, label='Proj X')
        if hull_dir0:
            hull_pts = np.array(hull_dir0)
            ax.plot(hull_pts[:, 0], np.full(len(hull_pts), y_min), hull_pts[:, 1],
                   'o-', color='orange', linewidth=2, markersize=4)
        if intersection_dir0:
            int_pts = np.array(intersection_dir0)
            ax.plot(int_pts[:, 0], np.full(len(int_pts), y_min), int_pts[:, 1],
                   'o', color='red', markersize=8, markeredgewidth=2, markerfacecolor='none')

    # Y-direction projection (on y-z plane at x=x_min)
    if control_points_dir1:
        pts = np.array(control_points_dir1)
        ax.scatter(np.full(len(pts), x_min), pts[:, 0], pts[:, 1],
                  c='orange', s=20, alpha=0.6, label='Proj Y')
        if hull_dir1:
            hull_pts = np.array(hull_dir1)
            ax.plot(np.full(len(hull_pts), x_min), hull_pts[:, 0], hull_pts[:, 1],
                   'o-', color='orange', linewidth=2, markersize=4)
        if intersection_dir1:
            int_pts = np.array(intersection_dir1)
            ax.plot(np.full(len(int_pts), x_min), int_pts[:, 0], int_pts[:, 1],
                   'o', color='red', markersize=8, markeredgewidth=2, markerfacecolor='none')

    # Plot zero contour (intersection with z=0)
    contour = ax.contour(X, Y, Z, levels=[0], colors='blue', linewidths=3)

    # Plot current box
    # Draw box edges
    corners = [
        [x_min, y_min, 0], [x_max, y_min, 0],
        [x_max, y_max, 0], [x_min, y_max, 0],
        [x_min, y_min, 0]
    ]
    corners = np.array(corners)
    ax.plot(corners[:, 0], corners[:, 1], corners[:, 2], 'k-', linewidth=2, label='Current box')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Set z limits to show the surface well
    z_min = min(Z.min(), -1)
    z_max = max(Z.max(), 1)
    ax.set_zlim(z_min, z_max)

def plot_2d_combined(ax, box, final_box, iteration, depth):
    """Plot both equations in 2D with contours and bounding boxes."""
    x_min, x_max, y_min, y_max = box

    # Create mesh
    X, Y = create_mesh_grid(box)
    Z1 = evaluate_polynomial_1(X, Y)
    Z2 = evaluate_polynomial_2(X, Y)

    # Plot contours
    ax.contour(X, Y, Z1, levels=[0], colors='blue', linewidths=2, label='Eq1')
    ax.contour(X, Y, Z2, levels=[0], colors='red', linewidths=2, label='Eq2')

    # Plot current box
    rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                         fill=False, edgecolor='black', linewidth=2, label='Current box')
    ax.add_patch(rect)

    # Plot final bounding box (PP bounds)
    if final_box:
        fx_min, fx_max, fy_min, fy_max = final_box
        rect_final = plt.Rectangle((fx_min, fy_min), fx_max - fx_min, fy_max - fy_min,
                                   fill=True, facecolor='purple', alpha=0.2,
                                   edgecolor='purple', linewidth=2, label='PP bounds')
        ax.add_patch(rect_final)

    # Plot expected solution (analytical)
    # Solve: x^2 + y^2 = 1 and x^2/4 + 4*y^2 = 1
    # From eq1: x^2 = 1 - y^2
    # Substitute into eq2: (1-y^2)/4 + 4*y^2 = 1
    # => 1/4 - y^2/4 + 4*y^2 = 1
    # => 15*y^2/4 = 3/4
    # => y^2 = 1/5 = 0.2
    # => y = sqrt(0.2) ≈ 0.447 (in [0,1])
    # => x^2 = 1 - 0.2 = 0.8
    # => x = sqrt(0.8) ≈ 0.894 (in [0,1])
    y_sol = np.sqrt(0.2)
    x_sol = np.sqrt(0.8)
    ax.plot(x_sol, y_sol, 'k*', markersize=15, label='Expected solution',
            markeredgewidth=2, markerfacecolor='yellow')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Both Equations + Bounds')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=8)

def visualize_iteration(iteration, output_dir='visualization_output'):
    """Create visualization for one iteration."""
    os.makedirs(output_dir, exist_ok=True)

    iter_num = iteration['iteration']
    depth = iteration['depth']
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

    # Create figure with 3 subplots
    fig = plt.figure(figsize=(18, 6))
    fig.suptitle(f'Step {iter_num + 1}: Depth {depth} - SUBDIVIDE', fontsize=16, fontweight='bold')

    # Subplot 1: Equation 1
    ax1 = fig.add_subplot(131, projection='3d')
    plot_3d_graph(ax1, evaluate_polynomial_1, global_box,
                 f'Equation 1: x² + y² - 1 = 0',
                 dir0_eq0['projected_points'], dir1_eq0['projected_points'],
                 dir0_eq0['convex_hull'], dir1_eq0['convex_hull'],
                 dir0_eq0['intersection'], dir1_eq0['intersection'])

    # Subplot 2: Equation 2
    ax2 = fig.add_subplot(132, projection='3d')
    plot_3d_graph(ax2, evaluate_polynomial_2, global_box,
                 f'Equation 2: x²/4 + 4y² - 1 = 0',
                 dir0_eq1['projected_points'], dir1_eq1['projected_points'],
                 dir0_eq1['convex_hull'], dir1_eq1['convex_hull'],
                 dir0_eq1['intersection'], dir1_eq1['intersection'])

    # Subplot 3: Combined 2D view
    ax3 = fig.add_subplot(133)
    plot_2d_combined(ax3, global_box, final_box, iter_num, depth)

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


