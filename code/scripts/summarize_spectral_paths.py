#!/usr/bin/env python3
"""Create the combined, interruption-safe spectral-path report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=Path("results"))
    args = parser.parse_args()
    endpoints = {
        row["topology"]: row
        for row in read_rows(args.results / "experiment2_spectral_endpoints.csv")
    }
    cases = [
        (
            "cross-channel",
            args.results / "experiment2_cross_channel_spectral_path.csv",
        ),
        (
            "winding-ring",
            args.results / "experiment2_winding_ring_spectral_path.csv",
        ),
    ]
    summaries = []
    for topology, path in cases:
        rows = read_rows(path)
        spectral = min(rows, key=lambda row: float(row["rho_TG"]))
        effective = min(rows, key=lambda row: float(row["rho_eff"]))
        summaries.append(
            (
                topology,
                spectral["m"], float(spectral["rho_TG"]),
                effective["m"], float(effective["rho_eff"]),
                float(endpoints[topology]["rho_TG"]),
                float(endpoints[topology]["rho_eff"]),
                max(float(row["spectral_stage_difference"]) for row in rows),
            )
        )

    output = args.results / "experiment2_spectral_path.txt"
    with output.open("w", encoding="utf-8") as stream:
        stream.write("True two-grid spectral-radius paths along finite PCG\n")
        stream.write("====================================================\n\n")
        stream.write("Version                   : 9.1.0\n")
        stream.write("Problems                  : 128/16, contrast 1e4\n")
        stream.write("Scanned interval          : m=1,...,128\n")
        stream.write("Lanczos iterations        : 200\n")
        stream.write("RHS solve tolerance       : 1e-6\n")
        stream.write("Endpoint column tolerance : 1e-10\n\n")
        stream.write(
            "Topology          Min-rho m   Min rho_TG   Min-eff m  "
            "Min rho_eff  Endpoint rho_TG  Endpoint rho_eff  Max stage diff\n"
        )
        stream.write("-" * 119 + "\n")
        for row in summaries:
            stream.write(
                f"{row[0]:15s} {row[1]:>9s}   {row[2]:.9f}  "
                f"{row[3]:>10s}  {row[4]:.9f}     {row[5]:.9f}       "
                f"{row[6]:.9f}       {row[7]:.3e}\n"
            )
        stream.write(
            "\nNote: rho_TG is a cold-start matrix-free Lanczos/Ritz estimate "
            "in the A-energy inner product. The stage difference compares "
            "the 100- and 200-step Ritz values. rho_eff is the independent "
            "constant-RHS residual-history factor. No energy curve, online "
            "selector, or oracle is used.\n"
        )


if __name__ == "__main__":
    main()
