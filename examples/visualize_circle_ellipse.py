#!/usr/bin/env python3
"""
Circle-Ellipse Intersection Visualization Example

This script demonstrates how to visualize the polynomial solver output
for the circle-ellipse intersection problem.

Usage:
    python examples/visualize_circle_ellipse.py [--max-steps N]

Requirements:
    - Python 3.6+
    - numpy
    - matplotlib

Install dependencies:
    pip install numpy matplotlib
    # or using uv:
    uv pip install numpy matplotlib
"""

import sys
import os

# Add tools directory to Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'tools'))

from solver_viz_api import visualize_solver

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Visualize circle-ellipse intersection solver output'
    )
    parser.add_argument(
        '--max-steps', 
        type=int, 
        default=None,
        help='Maximum number of iterations to visualize (default: all)'
    )
    parser.add_argument(
        '--strategy',
        choices=['ContractFirst', 'SubdivideFirst', 'Simultaneous', 'all'],
        default='all',
        help='Which strategy to visualize (default: all)'
    )
    
    args = parser.parse_args()
    
    # Define strategies to visualize
    if args.strategy == 'all':
        strategies = ['ContractFirst', 'SubdivideFirst', 'Simultaneous']
    else:
        strategies = [args.strategy]
    
    print("=" * 60)
    print("Circle-Ellipse Intersection Visualization")
    print("=" * 60)
    print()
    
    # Visualize each strategy
    for strategy in strategies:
        dump_file = f'dumps/strategy_{strategy}_geometry.txt'
        output_dir = f'visualizations/viz_{strategy}'
        
        # Check if dump file exists
        if not os.path.exists(dump_file):
            print(f"⚠️  Dump file not found: {dump_file}")
            print(f"   Run './build/bin/test_strategies' first to generate dumps")
            print()
            continue
        
        print(f"Strategy: {strategy}")
        print(f"  Dump file: {dump_file}")
        print(f"  Output dir: {output_dir}")
        
        try:
            count = visualize_solver(dump_file, output_dir, max_steps=args.max_steps)
            print(f"  ✓ Visualized {count} iterations")
        except Exception as e:
            print(f"  ✗ Error: {e}")
        
        print()
    
    print("=" * 60)
    print("Visualization complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()

