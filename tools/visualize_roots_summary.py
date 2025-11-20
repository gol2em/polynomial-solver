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


def parse_result_dump(dump_file):
    """Parse result dump file to extract box information."""

    boxes = {
        'resolved': [],
        'unresolved': []
    }

    dimension = None
    degeneracy_detected = False
    current_section = None
    current_box = {}

    with open(dump_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Parse dimension
            if line.startswith('# Dimension:'):
                dimension = int(line.split(':')[1].strip())

            # Parse degeneracy
            if line.startswith('# Degeneracy detected:'):
                degeneracy_detected = 'yes' in line

            # Parse section headers
            if line.startswith('## Resolved Boxes'):
                # Save previous box before switching sections
                if current_box and current_section:
                    boxes[current_section].append(current_box)
                    current_box = {}
                current_section = 'resolved'
                continue
            elif line.startswith('## Unresolved Boxes'):
                # Save previous box before switching sections
                if current_box and current_section:
                    boxes[current_section].append(current_box)
                    current_box = {}
                current_section = 'unresolved'
                continue

            # Parse box data
            if line.startswith('Box '):
                # Save previous box if exists
                if current_box and current_section:
                    boxes[current_section].append(current_box)
                current_box = {}
            elif line.startswith('Lower:'):
                values = line.split(':')[1].strip().split()
                current_box['lower'] = [float(v) for v in values]
            elif line.startswith('Upper:'):
                values = line.split(':')[1].strip().split()
                current_box['upper'] = [float(v) for v in values]
            elif line.startswith('Center:'):
                values = line.split(':')[1].strip().split()
                current_box['center'] = [float(v) for v in values]
            elif line.startswith('Depth:'):
                current_box['depth'] = int(line.split(':')[1].strip())
            elif line.startswith('Converged:'):
                current_box['converged'] = 'yes' in line

    # Save last box
    if current_box and current_section:
        boxes[current_section].append(current_box)

    return dimension, boxes, degeneracy_detected


def visualize_1d_summary(boxes, expected_roots=None, output_file=None):
    """Visualize 1D root finding summary."""

    fig, ax = plt.subplots(figsize=(14, 8))

    # Plot all boxes at different y-levels based on depth
    max_depth = max(
        max([b['depth'] for b in boxes['resolved']] or [0]),
        max([b['depth'] for b in boxes['unresolved']] or [0])
    )

    # Plot unresolved boxes (red)
    for box in boxes['unresolved']:
        x_min, x_max = box['lower'][0], box['upper'][0]
        y = box['depth']
        x_center = (x_min + x_max) / 2
        rect = patches.Rectangle((x_min, y - 0.3), x_max - x_min, 0.6,
                                 linewidth=2, edgecolor='red', facecolor='pink', alpha=0.5)
        ax.add_patch(rect)

        # Mark center
        ax.plot(x_center, y, 'ro', markersize=8, markeredgecolor='darkred', markeredgewidth=2)

    # Plot resolved boxes (green)
    for box in boxes['resolved']:
        x_min, x_max = box['lower'][0], box['upper'][0]
        y = box['depth']
        x_center = (x_min + x_max) / 2
        rect = patches.Rectangle((x_min, y - 0.3), x_max - x_min, 0.6,
                                 linewidth=2, edgecolor='darkgreen', facecolor='lightgreen', alpha=0.7)
        ax.add_patch(rect)

        # Mark root center
        ax.plot(x_center, y, 'go', markersize=10, markeredgecolor='darkgreen', markeredgewidth=2)

        # Add width label
        width = x_max - x_min
        ax.text(x_center, y + 0.5, f'w={width:.2e}', ha='center', va='bottom', fontsize=8)

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
    ax.set_title('Root Finding Summary: Resolved and Unresolved Boxes', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Add legend
    legend_elements = [
        patches.Patch(facecolor='lightgreen', edgecolor='darkgreen', label=f'Resolved ({len(boxes["resolved"])})'),
        patches.Patch(facecolor='pink', edgecolor='red', label=f'Unresolved ({len(boxes["unresolved"])})'),
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

    # Plot unresolved boxes (red)
    for box in boxes['unresolved']:
        x_min, x_max = box['lower'][0], box['upper'][0]
        y_min, y_max = box['lower'][1], box['upper'][1]
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                                 linewidth=1.5, edgecolor='red', facecolor='pink', alpha=0.4)
        ax.add_patch(rect)

        # Mark center
        ax.plot(x_center, y_center, 'ro', markersize=8, markeredgecolor='darkred', markeredgewidth=2)

    # Plot resolved boxes (green) - top layer
    for box in boxes['resolved']:
        x_min, x_max = box['lower'][0], box['upper'][0]
        y_min, y_max = box['lower'][1], box['upper'][1]
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                                 linewidth=2, edgecolor='darkgreen', facecolor='lightgreen', alpha=0.6)
        ax.add_patch(rect)

        # Mark root center
        ax.plot(x_center, y_center, 'go', markersize=10, markeredgecolor='darkgreen', markeredgewidth=2)

        # Add depth label
        depth = box['depth']
        ax.text(x_center, y_center, f'd={depth}', ha='center', va='center',
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
    ax.set_title('Root Finding Summary: Resolved and Unresolved Boxes', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    # Add legend
    legend_elements = [
        patches.Patch(facecolor='lightgreen', edgecolor='darkgreen', label=f'Resolved ({len(boxes["resolved"])})'),
        patches.Patch(facecolor='pink', edgecolor='red', label=f'Unresolved ({len(boxes["unresolved"])})'),
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


def print_summary(dimension, boxes, degeneracy_detected):
    """Print text summary of results."""
    print("\n" + "=" * 60)
    print("ROOT FINDING SUMMARY")
    print("=" * 60)
    print(f"Dimension: {dimension}D")
    print(f"Degeneracy detected: {'YES' if degeneracy_detected else 'NO'}")
    print(f"\nTotal boxes:")
    print(f"  Resolved:   {len(boxes['resolved']):4d} (found roots)")
    print(f"  Unresolved: {len(boxes['unresolved']):4d} (degeneracy/max depth)")
    print(f"\nTotal:        {sum(len(v) for v in boxes.values()):4d}")

    if boxes['resolved']:
        print(f"\nResolved boxes details:")
        for i, box in enumerate(boxes['resolved'], 1):
            if dimension == 1:
                x_min, x_max = box['lower'][0], box['upper'][0]
                x_center = box['center'][0]
                width = x_max - x_min
                print(f"  Root {i}: x = {x_center:.10f} (width: {width:.3e}, depth: {box['depth']})")
            else:
                x_min, x_max = box['lower'][0], box['upper'][0]
                y_min, y_max = box['lower'][1], box['upper'][1]
                x_center = box['center'][0]
                y_center = box['center'][1]
                width_x = x_max - x_min
                width_y = y_max - y_min
                print(f"  Root {i}: ({x_center:.10f}, {y_center:.10f}) " +
                      f"(width: {width_x:.3e} x {width_y:.3e}, depth: {box['depth']})")

    if boxes['unresolved']:
        print(f"\nUnresolved boxes details:")
        for i, box in enumerate(boxes['unresolved'], 1):
            if dimension == 1:
                x_min, x_max = box['lower'][0], box['upper'][0]
                x_center = box['center'][0]
                width = x_max - x_min
                print(f"  Box {i}: x = {x_center:.10f} (width: {width:.3e}, depth: {box['depth']})")
            else:
                x_min, x_max = box['lower'][0], box['upper'][0]
                y_min, y_max = box['lower'][1], box['upper'][1]
                x_center = box['center'][0]
                y_center = box['center'][1]
                width_x = x_max - x_min
                width_y = y_max - y_min
                print(f"  Box {i}: ({x_center:.10f}, {y_center:.10f}) " +
                      f"(width: {width_x:.3e} x {width_y:.3e}, depth: {box['depth']})")

    # Check for degeneracy indicators
    if boxes['resolved']:
        depths = [b['depth'] for b in boxes['resolved']]
        max_depth = max(depths)
        if max_depth > 10:
            print(f"\n⚠️  WARNING: Maximum depth {max_depth} suggests possible degeneracy")

        # Check for multiple roots at similar locations
        if dimension == 1 and len(boxes['resolved']) > 1:
            centers = [b['center'][0] for b in boxes['resolved']]
            centers.sort()
            for i in range(len(centers) - 1):
                dist = centers[i+1] - centers[i]
                if dist < 1e-6:
                    print(f"⚠️  WARNING: Multiple roots very close together (distance: {dist:.3e})")

    print("=" * 60 + "\n")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Visualize root finding summary from result dump')
    parser.add_argument('dump_file', help='Path to result dump file')
    parser.add_argument('--output', '-o', help='Output image file (default: show plot)')
    parser.add_argument('--expected-roots', nargs='+', type=float,
                       help='Expected root locations (1D: x1 x2 ..., 2D: x1 y1 x2 y2 ...)')

    args = parser.parse_args()

    # Parse dump file
    dimension, boxes, degeneracy_detected = parse_result_dump(args.dump_file)

    # Print summary
    print_summary(dimension, boxes, degeneracy_detected)

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

