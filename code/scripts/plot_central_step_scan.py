import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


source = sys.argv[1]
target = sys.argv[2]
with open(source, newline="", encoding="utf-8") as stream:
    data = list(csv.DictReader(stream))

steps = [int(row["m"]) for row in data]
cycles = [int(row["Cycles"]) for row in data]
factors = [float(row["Effective factor"]) for row in data]
converged = [row["Status"] == "converged" for row in data]
best_index = min(
    (index for index, value in enumerate(converged) if value),
    key=lambda index: cycles[index])

figure, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)
axes[0].plot(steps, cycles, color="#1f4e79", linewidth=1.7)
axes[0].scatter(
    [steps[index] for index, value in enumerate(converged) if value],
    [cycles[index] for index, value in enumerate(converged) if value],
    color="#1f4e79", s=18, zorder=3)
axes[0].scatter(
    [steps[index] for index, value in enumerate(converged) if not value],
    [cycles[index] for index, value in enumerate(converged) if not value],
    facecolors="white", edgecolors="#b33a3a", s=28, zorder=4,
    label="cycle cap")
axes[0].scatter(
    [steps[best_index]], [cycles[best_index]], color="#d17a00",
    marker="*", s=110, zorder=5, label="minimum")
axes[0].set_ylabel("Two-grid cycles")
axes[0].legend(frameon=False, ncol=2)

axes[1].plot(
    steps, factors, color="#2a7f62", marker="o", markersize=3,
    linewidth=1.5)
axes[1].scatter(
    [steps[best_index]], [factors[best_index]], color="#d17a00",
    marker="*", s=110, zorder=5)
axes[1].set_xlabel("Finite-PCG steps, m")
axes[1].set_ylabel("Effective convergence factor")

for axis in axes:
    axis.axvline(
        steps[best_index], color="#d17a00", linestyle="--",
        linewidth=1.0, alpha=0.8)
    axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

figure.suptitle("Central 128/16 cross-channel problem")
figure.tight_layout()
figure.savefig(target, dpi=220, bbox_inches="tight")
