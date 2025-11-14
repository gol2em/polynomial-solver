#!/usr/bin/env python3
"""Visualize polynomial graph test data produced by the C++ tests.

This script reads CSV files written by tests/test_polynomial_graph.cpp into the
build/ directory and produces PNGs showing

- 1D graph and control points, and
- 2D surface and 3D control net.

Run the C++ tests first (from the build directory):

    cmake .. && make -j4 && ctest --output-on-failure

Then, from the repository root, run:

    python python/visualize_graph_from_test.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parent.parent
BUILD_DIR = BASE_DIR / "build"
OUT_DIR = Path(__file__).resolve().parent


def _load_csv(path: Path) -> np.ndarray:
    """Load a simple comma-separated file with a single header line."""
    return np.loadtxt(path, delimiter=",", skiprows=1)


def plot_1d() -> None:
    """Plot the 1D graph and its control points from CSV files."""

    curve_path = (BUILD_DIR / "tests") / "graph_1d_curve.csv"
    ctrl_path = (BUILD_DIR / "tests") / "graph_1d_control_points.csv"

    if not curve_path.exists() or not ctrl_path.exists():
        print("[1D] CSV files not found. Make sure PolynomialGraphTest has run.")
        print(f"Expected: {curve_path} and {ctrl_path}")
        return

    curve = _load_csv(curve_path)
    ctrl = _load_csv(ctrl_path)

    t, f = curve[:, 0], curve[:, 1]
    cx, cy = ctrl[:, 0], ctrl[:, 1]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, f, label="graph f(t)")
    ax.plot(cx, cy, "o--", label="control points")

    ax.set_xlabel("t")
    ax.set_ylabel("f(t)")
    ax.set_title("1D polynomial graph and control points (from C++ test)")
    ax.legend()

    fig.tight_layout()
    out_path = OUT_DIR / "graph_1d_control_points.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    print(f"[1D] Wrote {out_path}")


def plot_2d() -> None:
    """Plot the 2D surface and its control net from CSV files."""

    surface_path = (BUILD_DIR / "tests") / "graph_2d_surface.csv"
    ctrl_path = (BUILD_DIR / "tests") / "graph_2d_control_points.csv"

    if not surface_path.exists() or not ctrl_path.exists():
        print("[2D] CSV files not found. Make sure PolynomialGraphTest has run.")
        print(f"Expected: {surface_path} and {ctrl_path}")
        return

    surface = _load_csv(surface_path)
    ctrl = _load_csv(ctrl_path)

    # Surface samples: x, y, z on a regular grid.
    x, y, z = surface[:, 0], surface[:, 1], surface[:, 2]
    xs = np.unique(x)
    ys = np.unique(y)
    nx, ny = len(xs), len(ys)

    X = xs.reshape(nx, 1).repeat(ny, axis=1)
    Y = ys.reshape(1, ny).repeat(nx, axis=0)
    Z = z.reshape(nx, ny)

    # Plot the surface.
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, cmap="viridis", alpha=0.7, linewidth=0)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x, y)")
    ax.set_title("2D polynomial surface (from C++ test)")

    fig.tight_layout()
    out_surface = OUT_DIR / "graph_2d_surface.png"
    fig.savefig(out_surface, dpi=150)
    plt.close(fig)

    print(f"[2D] Wrote {out_surface}")

    # Control net: x, y, z for each Bernstein coefficient.
    cx, cy, cz = ctrl[:, 0], ctrl[:, 1], ctrl[:, 2]
    ux = np.unique(cx)
    uy = np.unique(cy)
    nx_ctrl, ny_ctrl = len(ux), len(uy)

    CX = cx.reshape(nx_ctrl, ny_ctrl)
    CY = cy.reshape(nx_ctrl, ny_ctrl)
    CZ = cz.reshape(nx_ctrl, ny_ctrl)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")

    # Draw the grid of control points.
    for i in range(nx_ctrl):
        ax.plot(CX[i, :], CY[i, :], CZ[i, :], "k-", alpha=0.6)
    for j in range(ny_ctrl):
        ax.plot(CX[:, j], CY[:, j], CZ[:, j], "k-", alpha=0.6)

    ax.scatter(CX, CY, CZ, c="red", s=40)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x, y)")
    ax.set_title("2D control net (from C++ test)")

    fig.tight_layout()
    out_control = OUT_DIR / "graph_2d_control_net.png"
    fig.savefig(out_control, dpi=150)
    plt.close(fig)

    print(f"[2D] Wrote {out_control}")


def main() -> None:
    plot_1d()
    plot_2d()


if __name__ == "__main__":
    main()

