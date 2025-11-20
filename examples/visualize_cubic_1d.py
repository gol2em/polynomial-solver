#!/usr/bin/env python3
"""
1D Polynomial Solver Visualizer

Visualizes the 1D polynomial solving process step by step.
For 1D problems, shows:
- The polynomial curve p(x) (evaluated from Bernstein control points)
- Graph control points in 2D: (x, p(x))
- Convex hull of control points
- Intersection with y=0 axis (the root bounds)
- Current bounding box (interval computed from hull ∩ y=0)

Usage:
    python visualize_cubic_1d.py [--strategy STRATEGY] [--max-steps N]

Options:
    --strategy STRATEGY  Strategy to visualize: ContractFirst, SubdivideFirst, Simultaneous, or all (default: all)
    --max-steps N        Maximum number of iterations to visualize (default: all)
"""

import sys
import os
import argparse

# Add tools directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'tools'))

from solver_viz_api import visualize_solver

def main():
    parser = argparse.ArgumentParser(description='Visualize 1D polynomial solver geometry dumps')
    parser.add_argument('--strategy',
                       choices=['ContractFirst', 'SubdivideFirst', 'Simultaneous', 'all'],
                       default='all',
                       help='Strategy to visualize (default: all)')
    parser.add_argument('--method',
                       choices=['ProjectedPolyhedral', 'GraphHull', 'all'],
                       default='ProjectedPolyhedral',
                       help='Bounding method to visualize (default: ProjectedPolyhedral)')
    parser.add_argument('--max-steps', type=int, default=None,
                       help='Maximum number of iterations to visualize')

    args = parser.parse_args()

    # Determine which strategies to visualize
    if args.strategy == 'all':
        strategies = ['ContractFirst', 'SubdivideFirst', 'Simultaneous']
    else:
        strategies = [args.strategy]

    # Determine which methods to visualize
    if args.method == 'all':
        methods = ['ProjectedPolyhedral', 'GraphHull']
    else:
        methods = [args.method]

    print("=" * 60)
    print("1D Cubic Polynomial Visualization")
    print("=" * 60)
    print()

    for method in methods:
        for strategy in strategies:
            dump_file = f'dumps/cubic_1d_{method}_{strategy}_geometry.txt'
            output_dir = f'visualizations/viz_cubic_1d_{method}_{strategy}'

            if not os.path.exists(dump_file):
                print(f"Warning: {dump_file} not found, skipping {method}/{strategy}")
                continue

            print(f"Processing {method}/{strategy}...")
            print(f"  Input: {dump_file}")
            print(f"  Output: {output_dir}/")

            # Use the visualization API (auto-detects 1D)
            num_visualized = visualize_solver(dump_file, output_dir, max_steps=args.max_steps)

            print(f"  ✓ Completed {method}/{strategy}: {num_visualized} iterations")
            print()

    print("=" * 60)
    print("Visualization complete!")
    print("=" * 60)

if __name__ == '__main__':
    main()


