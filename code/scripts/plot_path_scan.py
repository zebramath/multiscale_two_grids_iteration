import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


with open(sys.argv[1], newline="", encoding="utf-8") as stream:
    data = list(csv.DictReader(stream))

steps = [int(row["m"]) for row in data]
energies = [float(row["Normalized energy excess"]) for row in data]
factors = [float(row["Effective factor"]) for row in data]
best_index = min(range(len(factors)), key=factors.__getitem__)
adaptive_step = int(sys.argv[4])
adaptive_index = steps.index(adaptive_step)

figure, (energy_axis, factor_axis) = plt.subplots(
    2, 1, figsize=(7.5, 7.4), sharex=True)
energy_axis.semilogy(
    steps, energies, color="#27647b", linewidth=1.55)
energy_axis.set_ylabel("Normalized energy excess")
energy_axis.set_title(sys.argv[3])
energy_axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.8)

factor_axis.plot(
    steps, factors, color="#27647b", linewidth=1.45,
    label=r"finite-PCG path, $\rho_{\mathrm{eff}}$")
factor_axis.scatter(
    [steps[best_index]], [factors[best_index]], color="#c96f00",
    marker="*", s=125, zorder=4, label="minimum in scanned interval")
factor_axis.scatter(
    [steps[adaptive_index]], [factors[adaptive_index]],
    facecolors="white", edgecolors="#8b3f73", linewidths=1.6,
    s=58, zorder=4, label="adaptive choice")
factor_axis.axvline(
    steps[best_index], color="#c96f00", linestyle="--",
    linewidth=0.9, alpha=0.75)
factor_axis.set_xlabel("Finite-PCG steps, m")
factor_axis.set_ylabel(r"Effective factor, $\rho_{\mathrm{eff}}$")
factor_axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.8)
factor_axis.legend(frameon=False)
for axis in (energy_axis, factor_axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
figure.tight_layout()
figure.savefig(sys.argv[2], dpi=220, bbox_inches="tight")
