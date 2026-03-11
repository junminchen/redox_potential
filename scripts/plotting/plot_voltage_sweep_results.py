#!/usr/bin/env python
"""
Plot charge-voltage trends from a voltage sweep result CSV.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Plot voltage sweep results")
    parser.add_argument("csv", help="Path to voltage_sweep_results.csv")
    parser.add_argument(
        "--output",
        help="Output image path. Defaults to <csv_dir>/voltage_sweep_summary.png",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv).resolve()
    out_path = (
        Path(args.output).resolve()
        if args.output
        else csv_path.parent / "voltage_sweep_summary.png"
    )

    df = pd.read_csv(csv_path)
    dqdv = df["q_cathode_mean"].diff() / df["voltage_v"].diff()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].plot(df["voltage_v"], df["q_cathode_mean"], "o-", label="Cathode")
    axes[0].plot(df["voltage_v"], -df["q_anode_mean"], "s-", label="-Anode")
    axes[0].set_xlabel("Total Voltage (V vs vacuum)")
    axes[0].set_ylabel("Electrode Charge (e)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(df["voltage_v"].iloc[1:], dqdv.iloc[1:], "o-", color="tab:red")
    axes[1].set_xlabel("Total Voltage (V vs vacuum)")
    axes[1].set_ylabel("dQ/dV (e/V)")
    axes[1].grid(True, alpha=0.3)

    fig.suptitle("Voltage Sweep Summary")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(out_path)


if __name__ == "__main__":
    main()
