#!/usr/bin/env python3
"""Visualize the graph hull computation step by step."""

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def visualize_1d():
    """Visualize 1D case: p(t) = t - 0.5"""
    df = pd.read_csv('build/graph_hull_1d_debug.csv')
    
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    # Step 0: Control points
    ax = axes[0]
    control = df[df['type'] == 'control_point']
    ax.plot(control['coord0'], control['coord1'], 'ro', markersize=10, label='Control points')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3, label='y=0 (root line)')
    ax.set_xlabel('t')
    ax.set_ylabel('p(t)')
    ax.set_title('Step 0: Control Points')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Step 1: Convex hull
    ax = axes[1]
    hull = df[df['type'] == 'hull_vertex']
    ax.plot(control['coord0'], control['coord1'], 'ro', markersize=10, label='Control points')
    if len(hull) > 0:
        hull_sorted = hull.sort_values('coord0')
        ax.plot(hull_sorted['coord0'], hull_sorted['coord1'], 'b-', linewidth=2, label='Convex hull')
        ax.plot(hull['coord0'], hull['coord1'], 'bs', markersize=8)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3, label='y=0')
    ax.set_xlabel('t')
    ax.set_ylabel('p(t)')
    ax.set_title('Step 1: Convex Hull')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Step 2: Hyperplane intersection
    ax = axes[2]
    ax.plot(control['coord0'], control['coord1'], 'ro', markersize=10, alpha=0.3, label='Control points')
    if len(hull) > 0:
        hull_sorted = hull.sort_values('coord0')
        ax.plot(hull_sorted['coord0'], hull_sorted['coord1'], 'b-', linewidth=2, alpha=0.3, label='Convex hull')
    hyperplane = df[df['type'] == 'hyperplane_intersection']
    if len(hyperplane) > 0:
        ax.plot(hyperplane['coord0'], hyperplane['coord1'], 'go', markersize=15, label='Intersection with y=0')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('t')
    ax.set_ylabel('p(t)')
    ax.set_title('Step 2: Intersection with y=0')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Step 3: Bounding box
    ax = axes[3]
    bbox_min = df[df['type'] == 'bbox_min']
    bbox_max = df[df['type'] == 'bbox_max']
    if len(bbox_min) > 0 and len(bbox_max) > 0:
        t_min = bbox_min['coord0'].values[0]
        t_max = bbox_max['coord0'].values[0]
        ax.axvspan(t_min, t_max, alpha=0.3, color='green', label=f'Bounding box: [{t_min:.3f}, {t_max:.3f}]')
    if len(hyperplane) > 0:
        ax.plot(hyperplane['coord0'], hyperplane['coord1'], 'go', markersize=15, label='Root')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(x=0.5, color='r', linestyle=':', alpha=0.5, label='True root: t=0.5')
    ax.set_xlabel('t')
    ax.set_ylabel('p(t)')
    ax.set_title('Step 3: Bounding Box')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 0.1)
    
    plt.tight_layout()
    plt.savefig('graph_hull_1d_visualization.png', dpi=150)
    print("1D visualization saved to graph_hull_1d_visualization.png")
    plt.show()

def visualize_2d():
    """Visualize 2D case: f1(x,y) = x - 0.5, f2(x,y) = y - 0.3"""
    df = pd.read_csv('build/graph_hull_2d_debug.csv')
    
    fig = plt.figure(figsize=(18, 12))
    
    # Equation 1: f1(x, y) = x - 0.5
    # Step 0: Control points for equation 1
    ax = fig.add_subplot(3, 3, 1, projection='3d')
    eq1_control = df[(df['equation'] == 1) & (df['type'] == 'control_point')]
    ax.scatter(eq1_control['coord0'], eq1_control['coord1'], eq1_control['coord2'], 
               c='red', s=100, marker='o', label='Control points')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('f1(x,y)')
    ax.set_title('Eq1: Control Points')
    ax.legend()
    
    # Step 1: Hull for equation 1
    ax = fig.add_subplot(3, 3, 2, projection='3d')
    eq1_hull = df[(df['equation'] == 1) & (df['type'] == 'hull_vertex')]
    ax.scatter(eq1_control['coord0'], eq1_control['coord1'], eq1_control['coord2'], 
               c='red', s=100, marker='o', alpha=0.3, label='Control points')
    if len(eq1_hull) > 0:
        ax.plot(eq1_hull['coord0'], eq1_hull['coord1'], eq1_hull['coord2'], 
                'b-', linewidth=2, label='Convex hull')
        ax.scatter(eq1_hull['coord0'], eq1_hull['coord1'], eq1_hull['coord2'], 
                   c='blue', s=80, marker='s')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('f1(x,y)')
    ax.set_title('Eq1: Convex Hull')
    ax.legend()
    
    # Equation 2: f2(x, y) = y - 0.3
    # Step 0: Control points for equation 2
    ax = fig.add_subplot(3, 3, 4, projection='3d')
    eq2_control = df[(df['equation'] == 2) & (df['type'] == 'control_point')]
    ax.scatter(eq2_control['coord0'], eq2_control['coord1'], eq2_control['coord2'], 
               c='red', s=100, marker='o', label='Control points')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('f2(x,y)')
    ax.set_title('Eq2: Control Points')
    ax.legend()
    
    # Step 1: Hull for equation 2
    ax = fig.add_subplot(3, 3, 5, projection='3d')
    eq2_hull = df[(df['equation'] == 2) & (df['type'] == 'hull_vertex')]
    ax.scatter(eq2_control['coord0'], eq2_control['coord1'], eq2_control['coord2'], 
               c='red', s=100, marker='o', alpha=0.3, label='Control points')
    if len(eq2_hull) > 0:
        ax.plot(eq2_hull['coord0'], eq2_hull['coord1'], eq2_hull['coord2'], 
                'b-', linewidth=2, label='Convex hull')
        ax.scatter(eq2_hull['coord0'], eq2_hull['coord1'], eq2_hull['coord2'], 
                   c='blue', s=80, marker='s')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('f2(x,y)')
    ax.set_title('Eq2: Convex Hull')
    ax.legend()

    # Step 2: Intersection of hulls
    ax = fig.add_subplot(3, 3, 7, projection='3d')
    if len(eq1_hull) > 0:
        ax.plot(eq1_hull['coord0'], eq1_hull['coord1'], eq1_hull['coord2'],
                'b-', linewidth=1, alpha=0.3, label='Hull 1')
    if len(eq2_hull) > 0:
        ax.plot(eq2_hull['coord0'], eq2_hull['coord1'], eq2_hull['coord2'],
                'g-', linewidth=1, alpha=0.3, label='Hull 2')
    hull_inter = df[df['type'] == 'hull_intersection']
    if len(hull_inter) > 0:
        ax.scatter(hull_inter['coord0'], hull_inter['coord1'], hull_inter['coord2'],
                   c='purple', s=100, marker='D', label='Hull intersection')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Step 2: Hull Intersection')
    ax.legend()

    # Step 3: Hyperplane intersection
    ax = fig.add_subplot(3, 3, 8, projection='3d')
    if len(hull_inter) > 0:
        ax.scatter(hull_inter['coord0'], hull_inter['coord1'], hull_inter['coord2'],
                   c='purple', s=100, marker='D', alpha=0.3, label='Hull intersection')
    hyperplane = df[df['type'] == 'hyperplane_intersection']
    if len(hyperplane) > 0:
        ax.scatter(hyperplane['coord0'], hyperplane['coord1'], hyperplane['coord2'],
                   c='green', s=150, marker='*', label='Intersection with z=0')
    # Draw z=0 plane
    xx, yy = np.meshgrid(np.linspace(-0.1, 1.1, 10), np.linspace(-0.1, 1.1, 10))
    zz = np.zeros_like(xx)
    ax.plot_surface(xx, yy, zz, alpha=0.2, color='gray')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Step 3: Intersection with z=0')
    ax.legend()

    # Step 4: Bounding box (2D projection)
    ax = fig.add_subplot(3, 3, 9)
    bbox_min = df[df['type'] == 'bbox_min']
    bbox_max = df[df['type'] == 'bbox_max']
    if len(bbox_min) > 0 and len(bbox_max) > 0:
        x_min = bbox_min['coord0'].values[0]
        y_min = bbox_min['coord1'].values[0]
        x_max = bbox_max['coord0'].values[0]
        y_max = bbox_max['coord1'].values[0]
        rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                              fill=True, alpha=0.3, color='green',
                              label=f'Bbox: [{x_min:.3f},{x_max:.3f}]×[{y_min:.3f},{y_max:.3f}]')
        ax.add_patch(rect)
    if len(hyperplane) > 0:
        ax.plot(hyperplane['coord0'], hyperplane['coord1'], 'g*', markersize=20, label='Root')
    ax.plot(0.5, 0.3, 'r+', markersize=20, markeredgewidth=3, label='True root: (0.5, 0.3)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Step 4: Bounding Box (2D projection)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig('graph_hull_2d_visualization.png', dpi=150)
    print("2D visualization saved to graph_hull_2d_visualization.png")
    plt.show()

if __name__ == '__main__':
    print("Visualizing 1D case...")
    visualize_1d()
    print("\nVisualizing 2D case...")
    visualize_2d()

