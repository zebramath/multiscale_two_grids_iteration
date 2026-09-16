#!/usr/bin/env python3
import csv
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
with open(sys.argv[1], newline="", encoding="utf-8") as stream:
    data = list(csv.DictReader(stream))
steps = [int(row["m"]) for row in data]
spectral = [float(row["rho_TG"]) for row in data]
effective = [float(row["rho_eff"]) for row in data]
spectral_best = min(range(len(spectral)), key=spectral.__getitem__)
effective_best = min(range(len(effective)), key=effective.__getitem__)
figure, axis = plt.subplots(figsize=(7.2, 4.5))
axis.plot(
    steps, spectral, color="#174a7e", linewidth=1.7,
    label=r"true $\rho_{\mathrm{TG}}$",
)
axis.plot(
    steps, effective, color="#bc5a2a", linewidth=1.35, linestyle="--",
    label=r"constant-RHS $\rho_{\mathrm{eff}}$",
)
axis.scatter(
    [steps[spectral_best]], [spectral[spectral_best]], marker="*", s=115,
    color="#174a7e", zorder=4, label=r"minimum $\rho_{\mathrm{TG}}$",
)
axis.scatter(
    [steps[effective_best]], [effective[effective_best]], marker="o", s=40,
    facecolors="white", edgecolors="#bc5a2a", linewidths=1.4, zorder=4,
    label=r"minimum $\rho_{\mathrm{eff}}$",
)
axis.set_xlabel("Finite-PCG steps, m")
axis.set_ylabel("Contraction factor")
axis.set_title(sys.argv[3])
axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.8)
axis.spines["top"].set_visible(False)
axis.spines["right"].set_visible(False)
axis.legend(frameon=False, ncol=2)
figure.tight_layout()
figure.savefig(sys.argv[2], dpi=220, bbox_inches="tight")
