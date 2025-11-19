#!/usr/bin/env python3
"""
Visualize and compare the three subdivision strategies.
Creates side-by-side visualizations showing the first few iterations of each strategy.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
import os

# Import the parser from the existing visualization script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from visualize_ellipse_dump import parse_dump_file

def compare_strategies():
    """Compare the three strategies side by side."""
    strategies = ['ContractFirst', 'SubdivideFirst', 'Simultaneous']
    dump_files = [f'build/strategy_{s}_geometry.txt' for s in strategies]
    
    # Check if files exist
    for f in dump_files:
        if not os.path.exists(f):
            print(f"Error: {f} not found. Run ./bin/test_strategies first.")
            return
    
    # Parse all dump files
    all_iterations = {}
    for strategy, dump_file in zip(strategies, dump_files):
        iterations = parse_dump_file(dump_file)
        all_iterations[strategy] = iterations
        print(f"{strategy}: {len(iterations)} iterations")
    
    # Create output directory
    output_dir = 'build/strategy_comparison'
    os.makedirs(output_dir, exist_ok=True)
    
    # Find the maximum number of iterations to visualize (limit to first 5)
    max_iters = min(5, max(len(iters) for iters in all_iterations.values()))
    
    print(f"\nGenerating comparison visualizations for first {max_iters} iterations...")
    
    # For each iteration, create a comparison figure
    for iter_idx in range(max_iters):
        fig = plt.figure(figsize=(24, 8))
        fig.suptitle(f'Iteration {iter_idx} - Strategy Comparison', fontsize=16, fontweight='bold')
        
        for col_idx, strategy in enumerate(strategies):
            iterations = all_iterations[strategy]
            
            if iter_idx >= len(iterations):
                # This strategy finished earlier
                ax = fig.add_subplot(1, 3, col_idx + 1, projection='3d')
                ax.text(0.5, 0.5, 0.5, f'{strategy}\n(Finished)', 
                       ha='center', va='center', fontsize=14, transform=ax.transAxes)
                ax.set_title(f'{strategy} (Finished)')
                continue
            
            iteration = iterations[iter_idx]
            
            # Create 3D subplot for this strategy
            ax = fig.add_subplot(1, 3, col_idx + 1, projection='3d')
            
            # Plot the iteration (simplified version)
            plot_strategy_iteration(ax, iteration, strategy)
        
        plt.tight_layout()
        output_file = os.path.join(output_dir, f'comparison_iter_{iter_idx:03d}.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file}")
    
    print(f"\nComparison visualizations saved to: {output_dir}/")
    print(f"Total files: {max_iters}")

def plot_strategy_iteration(ax, iteration, strategy_name):
    """Plot a single iteration for one strategy (simplified)."""
    # Extract global box
    global_box = iteration.get('global_box', [[0, 1], [0, 1]])
    x_min, x_max = global_box[0]
    y_min, y_max = global_box[1]
    
    # Set up the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z (function value)')
    ax.set_title(f"{strategy_name}\nDepth {iteration['depth']}, Decision: {iteration.get('decision', 'N/A')}")
    
    # Plot the box
    plot_box_3d(ax, x_min, x_max, y_min, y_max)
    
    # Plot ellipse surfaces (simplified)
    x = np.linspace(x_min, x_max, 30)
    y = np.linspace(y_min, y_max, 30)
    X, Y = np.meshgrid(x, y)
    
    # f1 = x^2 + y^2 - 1
    Z1 = X**2 + Y**2 - 1
    ax.plot_surface(X, Y, Z1, alpha=0.2, color='blue', label='f1')
    
    # f2 = x^2/4 + 4*y^2 - 1
    Z2 = X**2/4 + 4*Y**2 - 1
    ax.plot_surface(X, Y, Z2, alpha=0.2, color='red', label='f2')
    
    # Plot projected control points and hulls for first direction
    if iteration['directions']:
        dir0 = iteration['directions'][0]
        for eq in dir0['equations']:
            # Plot projected points
            if eq['projected_points']:
                pts = np.array(eq['projected_points'])
                # Transform from normalized [0,1] to actual box coordinates
                pts_x = x_min + pts[:, 0] * (x_max - x_min)
                pts_z = pts[:, 1]
                pts_y = np.full_like(pts_x, y_max)  # Project on back wall
                ax.scatter(pts_x, pts_y, pts_z, c='green', s=20, alpha=0.6)
    
    # Set view limits
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    z_range = max(abs(Z1.max()), abs(Z1.min()), abs(Z2.max()), abs(Z2.min()))
    ax.set_zlim(-z_range*1.1, z_range*1.1)

def plot_box_3d(ax, x_min, x_max, y_min, y_max):
    """Plot a 3D box outline."""
    # Draw box edges at z=0
    ax.plot([x_min, x_max], [y_min, y_min], [0, 0], 'k-', linewidth=1, alpha=0.3)
    ax.plot([x_min, x_max], [y_max, y_max], [0, 0], 'k-', linewidth=1, alpha=0.3)
    ax.plot([x_min, x_min], [y_min, y_max], [0, 0], 'k-', linewidth=1, alpha=0.3)
    ax.plot([x_max, x_max], [y_min, y_max], [0, 0], 'k-', linewidth=1, alpha=0.3)

if __name__ == '__main__':
    compare_strategies()

