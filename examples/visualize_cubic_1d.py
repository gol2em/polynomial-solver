#!/usr/bin/env python3
"""
1D Polynomial Solver Visualizer

Visualizes the 1D polynomial solving process step by step.
For 1D problems, shows:
- The polynomial graph (control points in 2D: x vs p(x))
- Convex hull of control points
- Intersection with y=0 axis (the roots)
- Current bounding interval

Usage:
    python visualize_cubic_1d.py [--strategy STRATEGY] [--max-steps N]

Options:
    --strategy STRATEGY  Strategy to visualize: ContractFirst, SubdivideFirst, Simultaneous, or all (default: all)
    --max-steps N        Maximum number of iterations to visualize (default: all)
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Add tools directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'tools'))

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

def visualize_1d_iteration(iteration, output_dir='visualizations'):
    """Visualize a single iteration of 1D polynomial solving."""
    os.makedirs(output_dir, exist_ok=True)
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Get data
    projected_points = np.array(iteration['projected_points']) if iteration['projected_points'] else None
    convex_hull = np.array(iteration['convex_hull']) if iteration['convex_hull'] else None
    intersection = np.array(iteration['intersection']) if iteration['intersection'] else None
    interval = iteration['interval']
    global_box = iteration['global_box']
    
    # Plot control points
    if projected_points is not None and len(projected_points) > 0:
        ax.scatter(projected_points[:, 0], projected_points[:, 1], 
                  c='blue', s=100, marker='o', label='Control Points', zorder=5)
    
    # Plot convex hull
    if convex_hull is not None and len(convex_hull) > 1:
        hull_closed = np.vstack([convex_hull, convex_hull[0]])
        ax.plot(hull_closed[:, 0], hull_closed[:, 1], 
               'g-', linewidth=2, label='Convex Hull', zorder=3)
        ax.fill(hull_closed[:, 0], hull_closed[:, 1], 
               alpha=0.2, color='green', zorder=2)
    
    # Plot intersection with y=0
    if intersection is not None and len(intersection) > 0:
        ax.scatter(intersection[:, 0], intersection[:, 1], 
                  c='red', s=150, marker='x', linewidths=3,
                  label='Intersection (y=0)', zorder=6)
    
    # Plot bounding interval
    if interval is not None and len(interval) == 2:
        ax.axvline(interval[0], color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.axvline(interval[1], color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.axhspan(-0.1, 0.1, xmin=(interval[0]-global_box[0])/(global_box[1]-global_box[0]),
                  xmax=(interval[1]-global_box[0])/(global_box[1]-global_box[0]),
                  alpha=0.3, color='red', label=f'Interval [{interval[0]:.6f}, {interval[1]:.6f}]')
    
    # Plot y=0 line
    ax.axhline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    
    # Set labels and title
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('p(x)', fontsize=12)
    ax.set_title(f'Iteration {iteration["iteration"]}, Depth {iteration["depth"]}\n' +
                f'Domain: [{global_box[0]:.6f}, {global_box[1]:.6f}]',
                fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Set axis limits with some padding
    if projected_points is not None and len(projected_points) > 0:
        x_min, x_max = projected_points[:, 0].min(), projected_points[:, 0].max()
        y_min, y_max = projected_points[:, 1].min(), projected_points[:, 1].max()
        x_range = x_max - x_min
        y_range = y_max - y_min
        ax.set_xlim(x_min - 0.1*x_range, x_max + 0.1*x_range)
        ax.set_ylim(y_min - 0.2*y_range, y_max + 0.2*y_range)
    
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, f'iteration_{iteration["iteration"]:03d}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    return output_file

def main():
    parser = argparse.ArgumentParser(description='Visualize 1D polynomial solver geometry dumps')
    parser.add_argument('--strategy',
                       choices=['ContractFirst', 'SubdivideFirst', 'Simultaneous', 'all'],
                       default='all',
                       help='Strategy to visualize (default: all)')
    parser.add_argument('--max-steps', type=int, default=None,
                       help='Maximum number of iterations to visualize')

    args = parser.parse_args()

    # Determine which strategies to visualize
    if args.strategy == 'all':
        strategies = ['ContractFirst', 'SubdivideFirst', 'Simultaneous']
    else:
        strategies = [args.strategy]

    print("=" * 60)
    print("1D Cubic Polynomial Visualization")
    print("=" * 60)
    print()

    for strategy in strategies:
        dump_file = f'dumps/cubic_1d_{strategy}_geometry.txt'
        output_dir = f'visualizations/viz_cubic_1d_{strategy}'

        if not os.path.exists(dump_file):
            print(f"Warning: {dump_file} not found, skipping {strategy}")
            continue

        print(f"Processing {strategy} strategy...")
        print(f"  Input: {dump_file}")
        print(f"  Output: {output_dir}/")

        # Parse dump file
        iterations = parse_1d_dump_file(dump_file)

        if args.max_steps is not None:
            iterations = iterations[:args.max_steps]

        print(f"  Iterations: {len(iterations)}")

        # Visualize each iteration
        for i, iteration in enumerate(iterations):
            output_file = visualize_1d_iteration(iteration, output_dir)
            if (i + 1) % 5 == 0 or (i + 1) == len(iterations):
                print(f"    Generated {i + 1}/{len(iterations)} visualizations")

        print(f"  ✓ Completed {strategy}")
        print()

    print("=" * 60)
    print("Visualization complete!")
    print("=" * 60)

if __name__ == '__main__':
    main()


