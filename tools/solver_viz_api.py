"""
Polynomial Solver Visualization API

This module provides a unified API for visualizing polynomial solver geometry dumps.
Automatically detects 1D vs 2D problems and uses the appropriate visualizer.

Example usage:
    from tools.solver_viz_api import visualize_solver

    # Visualize all iterations (auto-detects 1D or 2D)
    visualize_solver('dumps/example.txt', 'output/')

    # Visualize first 20 iterations
    visualize_solver('dumps/example.txt', 'output/', max_steps=20)
"""

import sys
import os

# Add tools directory to path
sys.path.insert(0, os.path.dirname(__file__))

from visualize_solver import parse_dump_file, visualize_iteration
from visualize_1d_solver import visualize_1d_solver, get_1d_iteration_count, visualize_1d_single_iteration

def detect_dimension(dump_file):
    """
    Detect whether a dump file is 1D or 2D by reading the header.

    Args:
        dump_file: Path to geometry dump file

    Returns:
        Dimension (1 or 2), or None if not found
    """
    with open(dump_file, 'r') as f:
        for line in f:
            if line.startswith('# Dimension:'):
                dim_str = line.split(':')[1].strip()
                return int(dim_str)
    return None

def visualize_solver(dump_file, output_dir, max_steps=None, expected_roots=None):
    """
    Visualize polynomial solver geometry dump (auto-detects 1D or 2D).

    Args:
        dump_file: Path to geometry dump file (e.g., 'dumps/strategy_ContractFirst_geometry.txt')
        output_dir: Output directory for PNG files (e.g., 'visualizations/viz_ContractFirst/')
        max_steps: Maximum number of iterations to visualize (default: all)
        expected_roots: List of expected root locations to mark on 1D plots (default: None)

    Returns:
        Number of iterations visualized

    Example:
        >>> visualize_solver('dumps/example.txt', 'output/', max_steps=10, expected_roots=[0.2, 0.5, 0.8])
        Detected 1D problem
        Visualized 10 iterations
        10
    """
    # Detect dimension
    dimension = detect_dimension(dump_file)

    if dimension == 1:
        print(f"Detected 1D problem")
        return visualize_1d_solver(dump_file, output_dir, max_steps=max_steps, expected_roots=expected_roots)
    elif dimension == 2:
        print(f"Detected 2D problem")
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Parse dump file
        print(f"Parsing {dump_file}...")
        iterations = parse_dump_file(dump_file)

        # Limit iterations if max_steps specified
        if max_steps is not None:
            iterations = iterations[:max_steps]

        print(f"Visualizing {len(iterations)} iterations...")

        # Visualize each iteration
        for i, iteration in enumerate(iterations):
            prev_iteration = iterations[i-1] if i > 0 else None
            visualize_iteration(iteration, prev_iteration=prev_iteration, output_dir=output_dir)
            print(f"  [{i+1}/{len(iterations)}] Iteration {iteration['iteration']}")

        print(f"\nVisualized {len(iterations)} iterations")
        print(f"Output directory: {output_dir}")

        return len(iterations)
    else:
        raise ValueError(f"Unsupported dimension: {dimension} (only 1D and 2D are supported)")


def get_iteration_count(dump_file):
    """
    Get the number of iterations in a dump file without visualizing (auto-detects 1D or 2D).

    Args:
        dump_file: Path to geometry dump file

    Returns:
        Number of iterations in the dump file
    """
    dimension = detect_dimension(dump_file)

    if dimension == 1:
        return get_1d_iteration_count(dump_file)
    elif dimension == 2:
        iterations = parse_dump_file(dump_file)
        return len(iterations)
    else:
        raise ValueError(f"Unsupported dimension: {dimension}")


def visualize_single_iteration(dump_file, iteration_num, output_dir, expected_roots=None):
    """
    Visualize a single iteration from a dump file (auto-detects 1D or 2D).

    Args:
        dump_file: Path to geometry dump file
        iteration_num: Iteration number to visualize (0-based index)
        output_dir: Output directory for PNG file
        expected_roots: List of expected root locations to mark on 1D plots (default: None)

    Example:
        >>> visualize_single_iteration('dumps/example.txt', 5, 'output/', expected_roots=[0.2, 0.5, 0.8])
    """
    dimension = detect_dimension(dump_file)

    if dimension == 1:
        return visualize_1d_single_iteration(dump_file, iteration_num, output_dir, expected_roots=expected_roots)
    elif dimension == 2:
        iterations = parse_dump_file(dump_file)

        if iteration_num >= len(iterations):
            raise ValueError(f"Iteration {iteration_num} not found (only {len(iterations)} iterations)")

        prev_iteration = iterations[iteration_num-1] if iteration_num > 0 else None
        visualize_iteration(iterations[iteration_num], prev_iteration=prev_iteration, output_dir=output_dir)
        print(f"Visualized iteration {iteration_num}")
    else:
        raise ValueError(f"Unsupported dimension: {dimension}")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize polynomial solver geometry dumps')
    parser.add_argument('dump_file', help='Path to geometry dump file')
    parser.add_argument('--output-dir', default='visualization_output', 
                        help='Output directory for PNG files (default: visualization_output)')
    parser.add_argument('--max-steps', type=int, default=None,
                        help='Maximum number of iterations to visualize (default: all)')
    
    args = parser.parse_args()
    
    visualize_solver(args.dump_file, args.output_dir, args.max_steps)

