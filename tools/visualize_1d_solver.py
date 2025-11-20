#!/usr/bin/env python3
"""
1D Polynomial Solver Visualizer

Visualizes the 1D polynomial solving process step by step.
For 1D problems, shows:
- The polynomial curve (evaluated from Bernstein control points)
- Graph control points in 2D: (x, p(x))
- Convex hull of control points
- Intersection with y=0 axis (the root bounds)
- Current bounding box (interval)

The bounding box is computed as the intersection of the convex hull with y=0.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def parse_1d_dump_file(filename):
    """Parse 1D geometry dump file."""
    iterations = []
    current_iter = None
    state = None  # Track which section we're reading: 'projected', 'hull', 'intersection'
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('#'):
                if line.startswith('# Iteration'):
                    parts = line.split(',')
                    iter_num = int(parts[0].split()[2])
                    depth = int(parts[1].split()[1])
                    current_iter = {
                        'iteration': iter_num,
                        'depth': depth,
                        'global_box': None,
                        'projected_points': [],
                        'convex_hull': [],
                        'intersection': [],
                        'interval': None,
                        'decision': None
                    }
                    iterations.append(current_iter)
                    state = None
                elif line.startswith('# Global box:') and current_iter:
                    box_str = line.split(':', 1)[1].strip()
                    box_str = box_str.strip('[]')
                    coords = [float(x) for x in box_str.split(',')]
                    current_iter['global_box'] = coords
                elif line.startswith('# Decision:') and current_iter:
                    current_iter['decision'] = line.split(':', 1)[1].strip()
                elif line.startswith('# FINAL_DECISION:') and current_iter:
                    current_iter['decision'] = line.split(':', 1)[1].strip()
                continue
            
            if not current_iter:
                continue
                
            if line.startswith('Projected_Points'):
                state = 'projected'
            elif line.startswith('ConvexHull'):
                state = 'hull'
            elif line.startswith('Intersection'):
                state = 'intersection'
            elif line.startswith('Interval') or line.startswith('Final_Interval'):
                interval_str = line.split('[', 1)[1].split(']')[0]
                coords = [float(x) for x in interval_str.split(',')]
                current_iter['interval'] = coords
                state = None
            else:
                # Parse point data
                parts = line.split()
                if len(parts) == 2:
                    try:
                        x, y = float(parts[0]), float(parts[1])
                        if state == 'projected':
                            current_iter['projected_points'].append([x, y])
                        elif state == 'hull':
                            current_iter['convex_hull'].append([x, y])
                        elif state == 'intersection':
                            current_iter['intersection'].append([x, y])
                    except ValueError:
                        pass
    
    return iterations

def evaluate_bernstein_1d(control_points, num_samples=200):
    """
    Evaluate 1D polynomial from Bernstein control points.

    Args:
        control_points: Array of [t, p(t)] control points in [0,1] parameter space
        num_samples: Number of points to evaluate

    Returns:
        t_vals, p_vals: Arrays of parameter values and polynomial values
    """
    if control_points is None or len(control_points) == 0:
        return np.array([]), np.array([])
    
    control_points = np.array(control_points)
    n = len(control_points) - 1  # degree
    
    # Extract y-coordinates (polynomial values at Bernstein points)
    b_coeffs = control_points[:, 1]
    
    # Evaluate using De Casteljau algorithm
    t_vals = np.linspace(0, 1, num_samples)
    p_vals = np.zeros(num_samples)
    
    for i, t in enumerate(t_vals):
        # De Casteljau algorithm
        coeffs = b_coeffs.copy()
        for r in range(1, n + 1):
            for j in range(n - r + 1):
                coeffs[j] = (1 - t) * coeffs[j] + t * coeffs[j + 1]
        p_vals[i] = coeffs[0]
    
    return t_vals, p_vals

def visualize_1d_iteration(iteration, output_dir='visualizations', show_function=True, expected_roots=None):
    """
    Visualize a single iteration of 1D polynomial solving.

    Args:
        iteration: Dictionary containing iteration data
        output_dir: Output directory for PNG files
        show_function: Whether to plot the actual polynomial curve
        expected_roots: List of expected root locations to mark on the plot

    Returns:
        Path to saved PNG file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    # Get data
    projected_points = np.array(iteration['projected_points']) if iteration['projected_points'] else None
    convex_hull = np.array(iteration['convex_hull']) if iteration['convex_hull'] else None
    intersection = np.array(iteration['intersection']) if iteration['intersection'] else None
    interval = iteration['interval']
    global_box = iteration['global_box']
    decision = iteration.get('decision', 'BOUNDING')
    
    # Plot the actual polynomial curve
    if show_function and projected_points is not None and len(projected_points) > 0:
        t_vals, p_vals = evaluate_bernstein_1d(projected_points, num_samples=300)
        # Map from [0,1] parameter space to global box
        x_vals = global_box[0] + t_vals * (global_box[1] - global_box[0])
        ax.plot(x_vals, p_vals, 'b-', linewidth=2, label='Polynomial p(x)', zorder=4, alpha=0.8)

    # Plot control points (in global coordinates)
    if projected_points is not None and len(projected_points) > 0:
        # Map control points from [0,1] to global box
        cp_x = global_box[0] + projected_points[:, 0] * (global_box[1] - global_box[0])
        cp_y = projected_points[:, 1]
        ax.scatter(cp_x, cp_y, c='orange', s=80, marker='o',
                  label='Control Points', zorder=5, edgecolors='black', linewidths=1)

    # Plot convex hull (in global coordinates)
    if convex_hull is not None and len(convex_hull) > 1:
        # Map hull from [0,1] to global box
        hull_x = global_box[0] + convex_hull[:, 0] * (global_box[1] - global_box[0])
        hull_y = convex_hull[:, 1]
        hull_points = np.column_stack([hull_x, hull_y])
        hull_closed = np.vstack([hull_points, hull_points[0]])
        ax.plot(hull_closed[:, 0], hull_closed[:, 1],
               'g-', linewidth=2, label='Convex Hull', zorder=3)
        ax.fill(hull_closed[:, 0], hull_closed[:, 1],
               alpha=0.15, color='green', zorder=2)

    # Plot intersection with y=0 (these are already in [0,1] parameter space)
    if intersection is not None and len(intersection) > 0:
        # Map intersection from [0,1] to global box
        int_x = global_box[0] + intersection[:, 0] * (global_box[1] - global_box[0])
        int_y = intersection[:, 1]
        ax.scatter(int_x, int_y, c='red', s=200, marker='x', linewidths=3,
                  label='Hull ∩ y=0', zorder=6)

    # Plot bounding box (interval) - this should match the intersection
    # The interval IS the bounding box computed from hull intersection with y=0
    if interval is not None and len(interval) == 2:
        # Map interval from [0,1] to global box
        bbox_min = global_box[0] + interval[0] * (global_box[1] - global_box[0])
        bbox_max = global_box[0] + interval[1] * (global_box[1] - global_box[0])

        ax.axvline(bbox_min, color='purple', linestyle='--', linewidth=2.5, alpha=0.8, zorder=7)
        ax.axvline(bbox_max, color='purple', linestyle='--', linewidth=2.5, alpha=0.8, zorder=7)

        # Get y-limits for shading
        y_min, y_max = ax.get_ylim() if ax.get_ylim() != (0, 1) else (-0.2, 0.2)
        ax.axvspan(bbox_min, bbox_max, alpha=0.15, color='purple', zorder=1,
                  label=f'Bounding Box [{bbox_min:.6f}, {bbox_max:.6f}]')

    # Plot y=0 line
    ax.axhline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.6, zorder=3)

    # Plot expected roots if provided - mark as points on the x-axis
    if expected_roots is not None:
        for root in expected_roots:
            ax.scatter([root], [0], c='darkgreen', s=150, marker='o',
                      edgecolors='black', linewidths=2, zorder=9, alpha=0.8)
        # Add to legend (only once)
        ax.scatter([], [], c='darkgreen', s=150, marker='o',
                  edgecolors='black', linewidths=2, label='Expected Roots')

    # Set labels and title
    ax.set_xlabel('x', fontsize=14, fontweight='bold')
    ax.set_ylabel('p(x)', fontsize=14, fontweight='bold')

    title = f'Iteration {iteration["iteration"]} (Depth {iteration["depth"]})'
    ax.set_title(title, fontsize=16, fontweight='bold')

    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')

    # Set axis limits with some padding
    if projected_points is not None and len(projected_points) > 0:
        # Use global box for x limits
        x_range = global_box[1] - global_box[0]
        ax.set_xlim(global_box[0] - 0.05*x_range, global_box[1] + 0.05*x_range)

        # Use control points for y limits
        y_min, y_max = projected_points[:, 1].min(), projected_points[:, 1].max()
        y_range = y_max - y_min
        if y_range < 1e-10:  # Handle flat case
            y_range = 0.1
        ax.set_ylim(y_min - 0.2*y_range, y_max + 0.2*y_range)

    plt.tight_layout()

    # Save figure
    output_file = os.path.join(output_dir, f'iteration_{iteration["iteration"]:03d}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    return output_file

def visualize_1d_solver(dump_file, output_dir, max_steps=None, show_function=True, expected_roots=None):
    """
    Visualize 1D polynomial solver geometry dump.

    Args:
        dump_file: Path to 1D geometry dump file
        output_dir: Output directory for PNG files
        max_steps: Maximum number of iterations to visualize (default: all)
        show_function: Whether to plot the actual polynomial curve (default: True)
        expected_roots: List of expected root locations to mark on the plot (default: None)

    Returns:
        Number of iterations visualized
    """
    os.makedirs(output_dir, exist_ok=True)

    # Parse dump file
    iterations = parse_1d_dump_file(dump_file)

    # Limit iterations if max_steps specified
    if max_steps is not None:
        iterations = iterations[:max_steps]

    # Visualize each iteration
    for iteration in iterations:
        visualize_1d_iteration(iteration, output_dir, show_function=show_function, expected_roots=expected_roots)

    return len(iterations)

def get_1d_iteration_count(dump_file):
    """
    Get the number of iterations in a 1D dump file.

    Args:
        dump_file: Path to 1D geometry dump file

    Returns:
        Number of iterations
    """
    iterations = parse_1d_dump_file(dump_file)
    return len(iterations)

def visualize_1d_single_iteration(dump_file, iteration_num, output_dir, show_function=True, expected_roots=None):
    """
    Visualize a single iteration from a 1D dump file.

    Args:
        dump_file: Path to 1D geometry dump file
        iteration_num: Iteration number to visualize (0-indexed)
        output_dir: Output directory for PNG file
        show_function: Whether to plot the actual polynomial curve (default: True)
        expected_roots: List of expected root locations to mark on the plot (default: None)

    Returns:
        Path to saved PNG file, or None if iteration not found
    """
    iterations = parse_1d_dump_file(dump_file)

    if iteration_num < 0 or iteration_num >= len(iterations):
        print(f"Error: Iteration {iteration_num} not found (file has {len(iterations)} iterations)")
        return None

    return visualize_1d_iteration(iterations[iteration_num], output_dir, show_function=show_function, expected_roots=expected_roots)



