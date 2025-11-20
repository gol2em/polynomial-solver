#!/usr/bin/env python3
"""
Visualize 1D polynomial solver for Wilkinson polynomial.

Problem: p(x) = (x-1)(x-2)(x-3)...(x-20)
Expected roots: x = 1, 2, 3, ..., 20
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

    dump_file = 'dumps/wilkinson_1d_geometry.txt'
    output_dir = 'visualizations/viz_wilkinson_1d'

    # Expected roots for Wilkinson polynomial (scaled to [0,1])
    # Roots at 1/20, 2/20, ..., 19/20
    expected_roots = [i/20.0 for i in range(1, 20)]

    print("=" * 60)
    print("Wilkinson Polynomial Visualization")
    print("=" * 60)
    print()
    print("Problem: p(x) = (x-1/20)(x-2/20)(x-3/20)...(x-19/20)")
    print("Expected roots: [1/20, 2/20, 3/20, ..., 19/20]")
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

if __name__ == '__main__':
    main()

