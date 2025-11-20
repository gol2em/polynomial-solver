#!/usr/bin/env python3
"""
Visualize 1D polynomial solver for high multiplicity root case.

Problem: p(x) = (x - 0.2)(x - 0.6)^6
Expected roots: x = 0.2 (multiplicity 1), x = 0.6 (multiplicity 6)
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

    dump_file = 'dumps/multiplicity_1d_geometry.txt'
    output_dir = 'visualizations/viz_multiplicity_1d'
    
    # Expected roots for p(x) = (x - 0.2)(x - 0.6)^6
    expected_roots = [0.2, 0.6]

    print("=" * 60)
    print("1D Polynomial with High Multiplicity Root Visualization")
    print("=" * 60)
    print()
    print("Problem: p(x) = (x - 0.2)(x - 0.6)^6")
    print("Expected roots: [0.2, 0.6] (multiplicity 6 at x=0.6)")
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

