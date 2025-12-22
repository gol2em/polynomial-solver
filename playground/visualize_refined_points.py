#!/usr/bin/env python3
"""
Visualize refined points from hessian_zero_set output.

Usage:
    uv run visualize_refined_points.py dumps/hessian_points.txt
    uv run visualize_refined_points.py dumps/hessian_points.txt -o output.png
    uv run visualize_refined_points.py dumps/hessian_points.txt --show-expected

The input file should have one point per line: x y
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def load_points(filepath: str) -> np.ndarray:
    """Load points from file (x y per line)."""
    points = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                points.append([float(parts[0]), float(parts[1])])
    return np.array(points)


def compute_arc_length_spacing(points: np.ndarray) -> np.ndarray:
    """
    Compute spacing between consecutive points along the curve.
    Points should be sorted by angle for circular curves.
    """
    # Sort by angle
    angles = np.arctan2(points[:, 1], points[:, 0])
    sorted_indices = np.argsort(angles)
    sorted_points = points[sorted_indices]
    
    # Compute distances between consecutive points (wrap around)
    n = len(sorted_points)
    distances = np.zeros(n)
    for i in range(n):
        j = (i + 1) % n
        distances[i] = np.linalg.norm(sorted_points[j] - sorted_points[i])
    
    return distances


def main():
    parser = argparse.ArgumentParser(description='Visualize refined curve points')
    parser.add_argument('input', help='Input file with points (x y per line)')
    parser.add_argument('-o', '--output', help='Output image file (default: show interactively)')
    parser.add_argument('--show-expected', action='store_true', 
                        help='Show expected circle at r=1/sqrt(2)')
    parser.add_argument('--analyze', action='store_true',
                        help='Show point distribution analysis')
    args = parser.parse_args()

    # Load points
    points = load_points(args.input)
    print(f"Loaded {len(points)} points from {args.input}")

    # Compute statistics
    radii = np.sqrt(points[:, 0]**2 + points[:, 1]**2)
    expected_r = 1.0 / np.sqrt(2.0)
    
    print(f"Expected radius: {expected_r:.10f}")
    print(f"Mean radius:     {np.mean(radii):.10f}")
    print(f"Std radius:      {np.std(radii):.2e}")
    print(f"Max error:       {np.max(np.abs(radii - expected_r)):.2e}")

    # Create figure
    fig, axes = plt.subplots(1, 2 if args.analyze else 1, figsize=(12 if args.analyze else 6, 6))
    if not args.analyze:
        axes = [axes]

    # Plot 1: Scatter plot of points
    ax = axes[0]
    ax.scatter(points[:, 0], points[:, 1], s=1, alpha=0.5, label=f'Refined points ({len(points)})')
    
    if args.show_expected:
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(expected_r * np.cos(theta), expected_r * np.sin(theta), 
                'r-', linewidth=1.5, label=f'Expected (r={expected_r:.4f})')
    
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Hessian Zero Set: det(H) = 0')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Distribution analysis
    if args.analyze:
        ax2 = axes[1]
        spacings = compute_arc_length_spacing(points)
        ax2.hist(spacings, bins=50, edgecolor='black', alpha=0.7)
        ax2.axvline(np.mean(spacings), color='r', linestyle='--', 
                    label=f'Mean: {np.mean(spacings):.4f}')
        ax2.axvline(np.median(spacings), color='g', linestyle='--',
                    label=f'Median: {np.median(spacings):.4f}')
        ax2.set_xlabel('Arc-length spacing')
        ax2.set_ylabel('Count')
        ax2.set_title('Point Distribution (spacing between consecutive points)')
        ax2.legend()
        
        # Print uniformity analysis
        cv = np.std(spacings) / np.mean(spacings)  # Coefficient of variation
        print(f"\nPoint distribution:")
        print(f"  Mean spacing:   {np.mean(spacings):.6f}")
        print(f"  Std spacing:    {np.std(spacings):.6f}")
        print(f"  CV (std/mean):  {cv:.4f}  (0 = perfectly uniform)")
        print(f"  Min spacing:    {np.min(spacings):.6f}")
        print(f"  Max spacing:    {np.max(spacings):.6f}")

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150, bbox_inches='tight')
        print(f"Saved to {args.output}")
    else:
        plt.show()


if __name__ == '__main__':
    main()

