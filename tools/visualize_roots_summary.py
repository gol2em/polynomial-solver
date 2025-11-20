#!/usr/bin/env python3
"""
Visualize roots found and degeneracy locations from geometry dump files.

This tool provides a high-level overview of:
- All boxes (converged, pruned, subdivided)
- Where degeneracy was detected
- Root locations
"""

import sys
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pathlib import Path


def parse_geometry_dump(dump_file):
    """Parse geometry dump file to extract box information."""
    
    boxes = {
        'converged': [],
        'pruned': [],
        'subdivided': []
    }
    
    dimension = None
    current_box = None
    current_depth = None
    
    with open(dump_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Parse dimension
            if line.startswith('# Dimension:'):
                dimension = int(line.split(':')[1].strip())
            
            # Parse iteration header
            if line.startswith('# Iteration'):
                match = re.search(r'Depth (\d+)', line)
                if match:
                    current_depth = int(match.group(1))
            
            # Parse global box
            if line.startswith('# Global box:'):
                box_str = line.split(':', 1)[1].strip()
                box_str = box_str.strip('[]')
                coords = [float(x) for x in box_str.split(',')]
                current_box = {
                    'coords': coords,
                    'depth': current_depth
                }
            
            # Parse final decision
            if line.startswith('# FINAL_DECISION:'):
                if current_box is not None:
                    if 'CONVERGED' in line:
                        # Extract iteration count
                        iter_match = re.search(r'iter=(\d+)', line)
                        current_box['iters'] = int(iter_match.group(1)) if iter_match else 0
                        boxes['converged'].append(current_box)
                    elif 'PRUNED' in line:
                        boxes['pruned'].append(current_box)
                    elif 'SUBDIVIDE' in line:
                        boxes['subdivided'].append(current_box)
                    current_box = None
    
    return dimension, boxes


def visualize_1d_summary(boxes, expected_roots=None, output_file=None):
    """Visualize 1D root finding summary."""
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot all boxes at different y-levels based on depth
    max_depth = max(
        max([b['depth'] for b in boxes['converged']] or [0]),
        max([b['depth'] for b in boxes['pruned']] or [0]),
        max([b['depth'] for b in boxes['subdivided']] or [0])
    )
    
    # Plot subdivided boxes (gray)
    for box in boxes['subdivided']:
        x_min, x_max = box['coords'][0], box['coords'][1]
        y = box['depth']
        rect = patches.Rectangle((x_min, y - 0.3), x_max - x_min, 0.6,
                                 linewidth=1, edgecolor='gray', facecolor='lightgray', alpha=0.3)
        ax.add_patch(rect)
    
    # Plot pruned boxes (red)
    for box in boxes['pruned']:
        x_min, x_max = box['coords'][0], box['coords'][1]
        y = box['depth']
        rect = patches.Rectangle((x_min, y - 0.3), x_max - x_min, 0.6,
                                 linewidth=2, edgecolor='red', facecolor='pink', alpha=0.5)
        ax.add_patch(rect)
    
    # Plot converged boxes (green)
    for box in boxes['converged']:
        x_min, x_max = box['coords'][0], box['coords'][1]
        y = box['depth']
        x_center = (x_min + x_max) / 2
        rect = patches.Rectangle((x_min, y - 0.3), x_max - x_min, 0.6,
                                 linewidth=2, edgecolor='darkgreen', facecolor='lightgreen', alpha=0.7)
        ax.add_patch(rect)
        
        # Mark root center
        ax.plot(x_center, y, 'go', markersize=10, markeredgecolor='darkgreen', markeredgewidth=2)
        
        # Add iteration count label
        iters = box.get('iters', 0)
        ax.text(x_center, y + 0.5, f'iter={iters}', ha='center', va='bottom', fontsize=8)
    
    # Plot expected roots if provided
    if expected_roots:
        for root in expected_roots:
            ax.axvline(root, color='blue', linestyle='--', linewidth=2, alpha=0.6, zorder=1)
            ax.plot(root, -0.5, 'b*', markersize=15, markeredgecolor='darkblue', markeredgewidth=1.5)
    
    # Set axis properties
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-1, max_depth + 1)
    ax.set_xlabel('x', fontsize=14)
    ax.set_ylabel('Subdivision Depth', fontsize=14)
    ax.set_title('Root Finding Summary: Boxes and Degeneracy', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add legend
    legend_elements = [
        patches.Patch(facecolor='lightgreen', edgecolor='darkgreen', label=f'Converged ({len(boxes["converged"])})'),
        patches.Patch(facecolor='pink', edgecolor='red', label=f'Pruned ({len(boxes["pruned"])})'),
        patches.Patch(facecolor='lightgray', edgecolor='gray', label=f'Subdivided ({len(boxes["subdivided"])})'),
    ]
    if expected_roots:
        legend_elements.append(plt.Line2D([0], [0], color='blue', linestyle='--', linewidth=2, 
                                         label=f'Expected Roots ({len(expected_roots)})'))
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to: {output_file}")
    else:
        plt.show()
    
    plt.close()


def visualize_2d_summary(boxes, expected_roots=None, output_file=None):
    """Visualize 2D root finding summary."""

    fig, ax = plt.subplots(figsize=(12, 12))

    # Plot all boxes
    # Plot subdivided boxes (gray) - lowest layer
    for box in boxes['subdivided']:
        x_min, x_max, y_min, y_max = box['coords']
        rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                                 linewidth=0.5, edgecolor='gray', facecolor='lightgray', alpha=0.2)
        ax.add_patch(rect)

    # Plot pruned boxes (red)
    for box in boxes['pruned']:
        x_min, x_max, y_min, y_max = box['coords']
        rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                                 linewidth=1.5, edgecolor='red', facecolor='pink', alpha=0.4)
        ax.add_patch(rect)

    # Plot converged boxes (green) - top layer
    for box in boxes['converged']:
        x_min, x_max, y_min, y_max = box['coords']
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                                 linewidth=2, edgecolor='darkgreen', facecolor='lightgreen', alpha=0.6)
        ax.add_patch(rect)

        # Mark root center
        ax.plot(x_center, y_center, 'go', markersize=10, markeredgecolor='darkgreen', markeredgewidth=2)

        # Add iteration count label
        iters = box.get('iters', 0)
        ax.text(x_center, y_center, f'{iters}', ha='center', va='center',
               fontsize=8, fontweight='bold', color='darkgreen')

    # Plot expected roots if provided
    if expected_roots:
        for root in expected_roots:
            ax.plot(root[0], root[1], 'b*', markersize=15, markeredgecolor='darkblue',
                   markeredgewidth=1.5, zorder=100)

    # Set axis properties
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('x', fontsize=14)
    ax.set_ylabel('y', fontsize=14)
    ax.set_title('Root Finding Summary: Boxes and Degeneracy', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    # Add legend
    legend_elements = [
        patches.Patch(facecolor='lightgreen', edgecolor='darkgreen', label=f'Converged ({len(boxes["converged"])})'),
        patches.Patch(facecolor='pink', edgecolor='red', label=f'Pruned ({len(boxes["pruned"])})'),
        patches.Patch(facecolor='lightgray', edgecolor='gray', label=f'Subdivided ({len(boxes["subdivided"])})'),
    ]
    if expected_roots:
        legend_elements.append(plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='blue',
                                         markersize=12, label=f'Expected Roots ({len(expected_roots)})'))
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to: {output_file}")
    else:
        plt.show()

    plt.close()


def print_summary(dimension, boxes):
    """Print text summary of results."""
    print("\n" + "=" * 60)
    print("ROOT FINDING SUMMARY")
    print("=" * 60)
    print(f"Dimension: {dimension}D")
    print(f"\nTotal boxes processed:")
    print(f"  Converged:  {len(boxes['converged']):4d} (found roots)")
    print(f"  Pruned:     {len(boxes['pruned']):4d} (no roots)")
    print(f"  Subdivided: {len(boxes['subdivided']):4d} (intermediate)")
    print(f"\nTotal:        {sum(len(v) for v in boxes.values()):4d}")

    if boxes['converged']:
        print(f"\nConverged boxes details:")
        for i, box in enumerate(boxes['converged'], 1):
            if dimension == 1:
                x_min, x_max = box['coords']
                x_center = (x_min + x_max) / 2
                width = x_max - x_min
                print(f"  Root {i}: x = {x_center:.10f} (width: {width:.3e}, depth: {box['depth']}, iters: {box.get('iters', 0)})")
            else:
                x_min, x_max, y_min, y_max = box['coords']
                x_center = (x_min + x_max) / 2
                y_center = (y_min + y_max) / 2
                width_x = x_max - x_min
                width_y = y_max - y_min
                print(f"  Root {i}: ({x_center:.10f}, {y_center:.10f}) " +
                      f"(width: {width_x:.3e} x {width_y:.3e}, depth: {box['depth']}, iters: {box.get('iters', 0)})")

    # Check for degeneracy indicators
    if boxes['converged']:
        depths = [b['depth'] for b in boxes['converged']]
        max_depth = max(depths)
        if max_depth > 10:
            print(f"\n⚠️  WARNING: Maximum depth {max_depth} suggests possible degeneracy")

        # Check for multiple roots at similar locations
        if dimension == 1 and len(boxes['converged']) > 1:
            centers = [(b['coords'][0] + b['coords'][1]) / 2 for b in boxes['converged']]
            centers.sort()
            for i in range(len(centers) - 1):
                dist = centers[i+1] - centers[i]
                if dist < 1e-6:
                    print(f"⚠️  WARNING: Multiple roots very close together (distance: {dist:.3e})")

    print("=" * 60 + "\n")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Visualize root finding summary from geometry dump')
    parser.add_argument('dump_file', help='Path to geometry dump file')
    parser.add_argument('--output', '-o', help='Output image file (default: show plot)')
    parser.add_argument('--expected-roots', nargs='+', type=float,
                       help='Expected root locations (1D: x1 x2 ..., 2D: x1 y1 x2 y2 ...)')

    args = parser.parse_args()

    # Parse dump file
    dimension, boxes = parse_geometry_dump(args.dump_file)

    # Print summary
    print_summary(dimension, boxes)

    # Parse expected roots
    expected_roots = None
    if args.expected_roots:
        if dimension == 1:
            expected_roots = args.expected_roots
        else:  # 2D
            if len(args.expected_roots) % 2 != 0:
                print("Error: For 2D, expected roots must be pairs of (x, y)")
                return 1
            expected_roots = [(args.expected_roots[i], args.expected_roots[i+1])
                            for i in range(0, len(args.expected_roots), 2)]

    # Visualize
    if dimension == 1:
        visualize_1d_summary(boxes, expected_roots, args.output)
    else:
        visualize_2d_summary(boxes, expected_roots, args.output)

    return 0


if __name__ == '__main__':
    sys.exit(main())

