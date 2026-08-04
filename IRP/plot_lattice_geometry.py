#!/usr/bin/env python3
# ============================================================================
# plot_lattice_geometry.py
#
# Visualises the output of generate_lattice_geometry.cpp. Reads
# lattice_surfaces.csv (airfoil outlines) and lattice_centers.csv (pivot +
# angle per unit) from the current directory and plots the lattice, with
# the derived domain box and margins drawn around it.
#
# Usage:
#   python3 plot_lattice_geometry.py
#   python3 plot_lattice_geometry.py --outdir /path/to/csvs
#
# Requires: matplotlib, pandas (pip3 install matplotlib pandas)
# ============================================================================

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".",
                     help="folder containing lattice_centers.csv / lattice_surfaces.csv")
    ap.add_argument("--save", default="lattice_geometry.png",
                     help="output image filename (set to '' to skip saving)")
    ap.add_argument("--N_chord", type=float, default=None,
                     help="chord length in lattice units, for chord-unit axis labels "
                          "(optional, defaults to reading from lattice_centers.csv spacing if omitted)")
    args = ap.parse_args()

    centers_path = os.path.join(args.outdir, "lattice_centers.csv")
    surfaces_path = os.path.join(args.outdir, "lattice_surfaces.csv")

    centers = pd.read_csv(centers_path)
    surfaces = pd.read_csv(surfaces_path)

    fig, ax = plt.subplots(figsize=(12, 6))

    # --- draw every airfoil outline (upper surface out, lower surface back) ---
    for (i, j), grp in surfaces.groupby(["i", "j"]):
        grp = grp.reset_index(drop=True)
        x = list(grp["x_upper"]) + list(grp["x_lower"])[::-1]
        y = list(grp["y_upper"]) + list(grp["y_lower"])[::-1]
        ax.fill(x, y, color="teal", alpha=0.6, linewidth=0.6, edgecolor="black")

    # --- mark pivots and label angle for a spot check ---
    ax.scatter(centers["cx"], centers["cy"], s=6, color="black", zorder=5)
    for _, row in centers.iterrows():
        ax.annotate(f"{row.theta_deg:.1f}°", (row.cx, row.cy),
                     textcoords="offset points", xytext=(0, 8),
                     fontsize=7, ha="center", color="dimgray")

    # --- domain box from the data extent (matches what buildLattice derived) ---
    xmin, xmax = 0, surfaces[["x_upper", "x_lower"]].values.max() * 1.0
    ymin, ymax = 0, surfaces[["y_upper", "y_lower"]].values.max() * 1.0
    # Better: infer full domain from console output if you want exact Lx/Ly;
    # here we just pad past the outermost geometry so margins are visible.
    pad_x = 0.15 * (xmax - xmin)
    pad_y = 0.4 * (ymax - ymin)

    ax.set_xlim(xmin - pad_x * 0.3, xmax + pad_x)
    ax.set_ylim(ymin - pad_y * 0.3, ymax + pad_y)
    ax.set_aspect("equal")
    ax.set_title("Airfoil lattice geometry (angle in degrees, black dot = pivot)")
    ax.set_xlabel("x [lattice units]")
    ax.set_ylabel("y [lattice units]")

    plt.tight_layout()
    if args.save:
        plt.savefig(os.path.join(args.outdir, args.save), dpi=150)
        print(f"Saved plot to {os.path.join(args.outdir, args.save)}")
    plt.show()


if __name__ == "__main__":
    main()