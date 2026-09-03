import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


with open(sys.argv[1], newline="", encoding="utf-8") as stream:
    data = list(csv.DictReader(stream))

steps = [int(row["m"]) for row in data]
factors = [float(row["Effective factor"]) for row in data]
best_index = min(range(len(factors)), key=factors.__getitem__)
adaptive_step = 43
adaptive_index = steps.index(adaptive_step)

figure, axis = plt.subplots(figsize=(7.4, 4.6))
axis.plot(
    steps, factors, color="#27647b", marker="o", markersize=2.8,
    linewidth=1.45, label="finite-PCG path")
axis.scatter(
    [steps[best_index]], [factors[best_index]], color="#c96f00",
    marker="*", s=125, zorder=4, label="minimum factor")
axis.scatter(
    [steps[adaptive_index]], [factors[adaptive_index]],
    facecolors="white", edgecolors="#8b3f73", linewidths=1.6,
    s=58, zorder=4, label="adaptive choice")
axis.axvline(
    steps[best_index], color="#c96f00", linestyle="--",
    linewidth=0.9, alpha=0.75)
axis.set_xlabel("Finite-PCG steps, m")
axis.set_ylabel("Effective convergence factor")
axis.set_title("Central 128/16 cross-channel problem")
axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.8)
axis.spines["top"].set_visible(False)
axis.spines["right"].set_visible(False)
axis.legend(frameon=False)
figure.tight_layout()
figure.savefig(sys.argv[2], dpi=220, bbox_inches="tight")
