#!/usr/bin/env python3
"""Visualize circle test results in 2D with boxes overlaid."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
import os

def parse_boxes_from_result(result_file):
    """Parse boxes from result file (unresolved boxes only)."""
    boxes = []

    with open(result_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            # Format: x_min x_max y_min y_max
            coords = [float(x) for x in line.split()]
            if len(coords) == 4:
                boxes.append(coords)

    return boxes

def visualize_circle_boxes_2d(result_file, output_file='circle_boxes_2d.png'):
    """Create 2D visualization of circle with boxes."""
    boxes = parse_boxes_from_result(result_file)

    print(f"Parsed {len(boxes)} unresolved boxes from {result_file}")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # Draw unit circle
    theta = np.linspace(0, 2*np.pi, 1000)
    x_circle = np.cos(theta)
    y_circle = np.sin(theta)
    ax.plot(x_circle, y_circle, 'b-', linewidth=3, label='Unit circle: x² + y² = 1', zorder=10)

    # Draw unresolved boxes (all kept boxes)
    for box in boxes:
        x_min, x_max, y_min, y_max = box
        width = x_max - x_min
        height = y_max - y_min
        rect = Rectangle((x_min, y_min), width, height,
                        linewidth=1.5, edgecolor='red', facecolor='pink',
                        alpha=0.4, zorder=2)
        ax.add_patch(rect)

    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [
        plt.Line2D([0], [0], color='b', linewidth=3, label='Unit circle'),
        Patch(facecolor='pink', edgecolor='red', alpha=0.4, label=f'Unresolved boxes ({len(boxes)})')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12)
    
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('x', fontsize=14)
    ax.set_ylabel('y', fontsize=14)
    ax.set_title('Circle Test: x² + y² - 1 = 0\nBox Coverage Analysis', fontsize=16)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n✅ Saved: {output_file}")
    plt.close()

if __name__ == '__main__':
    result_file = 'dumps/circle_result_boxes.txt'
    output_file = 'dumps/circle_boxes_2d.png'

    if len(sys.argv) > 1:
        result_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]

    visualize_circle_boxes_2d(result_file, output_file)

