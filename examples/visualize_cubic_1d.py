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
    python visualize_cubic_1d.py [--max-steps N]

Options:
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
    parser.add_argument('--max-steps', type=int, default=None,
                       help='Maximum number of iterations to visualize')

    args = parser.parse_args()

    dump_file = 'dumps/cubic_1d_geometry.txt'
    output_dir = 'visualizations/viz_cubic_1d'

    # Expected roots for p(x) = (x - 0.2)(x - 0.5)(x - 0.8)
    expected_roots = [0.2, 0.5, 0.8]

    if not os.path.exists(dump_file):
        print(f"Error: {dump_file} not found")
        print("Run './build/bin/example_cubic_1d' first to generate the dump file")
        return 1

    print("=" * 60)
    print("1D Cubic Polynomial Visualization")
    print("=" * 60)
    print()
    print(f"Problem: p(x) = (x - 0.2)(x - 0.5)(x - 0.8)")
    print(f"Expected roots: {expected_roots}")
    print()
    print(f"Input: {dump_file}")
    print(f"Output: {output_dir}/")
    print()

    # Use the visualization API (auto-detects 1D)
    num_visualized = visualize_solver(dump_file, output_dir, max_steps=args.max_steps, expected_roots=expected_roots)

    print()
    print("=" * 60)
    print(f"Visualization complete! Generated {num_visualized} images")
    print(f"View images in: {output_dir}/")
    print("=" * 60)

    return 0

if __name__ == '__main__':
    sys.exit(main())


