"""
Polynomial Solver Visualization API

This module provides a simple API for visualizing polynomial solver geometry dumps.

Example usage:
    from tools.solver_viz_api import visualize_solver
    
    # Visualize all iterations
    visualize_solver('dumps/example.txt', 'output/')
    
    # Visualize first 20 iterations
    visualize_solver('dumps/example.txt', 'output/', max_steps=20)
"""

import sys
import os

# Add tools directory to path to import visualize_solver
sys.path.insert(0, os.path.dirname(__file__))

from visualize_solver import parse_dump_file, visualize_iteration

def visualize_solver(dump_file, output_dir, max_steps=None):
    """
    Visualize polynomial solver geometry dump.
    
    Args:
        dump_file: Path to geometry dump file (e.g., 'dumps/strategy_ContractFirst_geometry.txt')
        output_dir: Output directory for PNG files (e.g., 'visualizations/viz_ContractFirst/')
        max_steps: Maximum number of iterations to visualize (default: all)
        
    Returns:
        Number of iterations visualized
        
    Example:
        >>> visualize_solver('dumps/example.txt', 'output/', max_steps=10)
        Visualized 10 iterations
        10
    """
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


def get_iteration_count(dump_file):
    """
    Get the number of iterations in a dump file without visualizing.
    
    Args:
        dump_file: Path to geometry dump file
        
    Returns:
        Number of iterations in the dump file
    """
    iterations = parse_dump_file(dump_file)
    return len(iterations)


def visualize_single_iteration(dump_file, iteration_num, output_dir):
    """
    Visualize a single iteration from a dump file.

    Args:
        dump_file: Path to geometry dump file
        iteration_num: Iteration number to visualize (0-based index)
        output_dir: Output directory for PNG file

    Example:
        >>> visualize_single_iteration('dumps/example.txt', 5, 'output/')
    """
    iterations = parse_dump_file(dump_file)

    if iteration_num >= len(iterations):
        raise ValueError(f"Iteration {iteration_num} not found (only {len(iterations)} iterations)")

    prev_iteration = iterations[iteration_num-1] if iteration_num > 0 else None
    visualize_iteration(iterations[iteration_num], prev_iteration=prev_iteration, output_dir=output_dir)
    print(f"Visualized iteration {iteration_num}")


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

