#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la

def symmetric_power(matrix: np.ndarray, exponent: float) -> np.ndarray:
    values, vectors = la.eigh(matrix, check_finite=True)
    scale = max(1.0, float(np.max(np.abs(values))))
    if float(np.min(values)) <= 100.0 * np.finfo(float).eps * scale:
        raise RuntimeError("expected a numerically positive-definite matrix")
    return (vectors * values**exponent) @ vectors.T

def load_export(path: Path) -> tuple[np.ndarray, dict[int, np.ndarray], np.ndarray]:
    records: list[tuple[str, int, int, int, float]] = []
    maximum_a = -1
    maximum_p_row = -1
    maximum_p_col = -1
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            record = (
                row["kind"], int(row["m"]), int(row["row"]),
                int(row["col"]), float(row["value"]),
            )
            records.append(record)
            if record[0] == "A":
                maximum_a = max(maximum_a, record[2], record[3])
            else:
                maximum_p_row = max(maximum_p_row, record[2])
                maximum_p_col = max(maximum_p_col, record[3])
    size = maximum_a + 1
    coarse_size = maximum_p_col + 1
    if maximum_p_row + 1 != size:
        raise RuntimeError("A and P dimensions in export are inconsistent")
    a = np.zeros((size, size))
    finite: dict[int, np.ndarray] = {}
    endpoint = np.zeros((size, coarse_size))
    for kind, step, row, col, value in records:
        if kind == "A":
            a[row, col] = value
        elif kind == "P":
            finite.setdefault(step, np.zeros((size, coarse_size)))[row, col] = value
        elif kind == "P_endpoint":
            endpoint[row, col] = value
        else:
            raise RuntimeError(f"unknown record kind {kind!r}")
    return a, finite, endpoint

def grid_partition(intervals: int, ratio: int) -> tuple[np.ndarray, np.ndarray]:
    side = intervals - 1
    coarse: list[int] = []
    fine: list[int] = []
    for row in range(side):
        for col in range(side):
            node = row * side + col
            if (col + 1) % ratio == 0 and (row + 1) % ratio == 0:
                coarse.append(node)
            else:
                fine.append(node)
    return np.asarray(fine, dtype=int), np.asarray(coarse, dtype=int)

def graph_factor(
    z: np.ndarray, t_f: np.ndarray, t_star: np.ndarray
) -> float:
    left = symmetric_power(np.eye(z.shape[0]) + z @ z.T, -0.5)
    return float(la.svdvals(left @ (t_f - z @ t_star))[0])

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--intervals", type=int, default=16)
    parser.add_argument("--coarse-intervals", type=int, default=4)
    args = parser.parse_args()
    source = args.results / "experiment4_local_diagnostic_matrices.csv"
    a, path, p_endpoint = load_export(source)
    ratio = args.intervals // args.coarse_intervals
    f_index, c_index = grid_partition(args.intervals, ratio)
    b = a[np.ix_(f_index, f_index)]
    c = a[np.ix_(f_index, c_index)]
    w_star = -la.solve(b, c, assume_a="pos")
    endpoint_error = float(la.norm(p_endpoint[f_index, :] - w_star, ord="fro"))

    s = a[np.ix_(c_index, c_index)] - c.T @ la.solve(b, c, assume_a="pos")
    b_half = symmetric_power(b, 0.5)
    b_inverse_half = symmetric_power(b, -0.5)
    s_inverse_half = symmetric_power(s, -0.5)
    a_half = symmetric_power(a, 0.5)
    a_inverse_half = symmetric_power(a, -0.5)

    lower = np.tril(a)
    g = np.eye(a.shape[0]) - la.solve_triangular(
        lower, a, lower=True, check_finite=True
    )
    t = a_half @ g @ a_inverse_half

    q_f = np.zeros((a.shape[0], f_index.size))
    q_f[f_index, :] = b_inverse_half
    q_star = p_endpoint @ s_inverse_half
    q_f_tilde = a_half @ q_f
    q_star_tilde = a_half @ q_star
    orthogonality_defect = float(
        la.norm(
            np.block([q_f_tilde, q_star_tilde]).T
            @ np.block([q_f_tilde, q_star_tilde])
            - np.eye(a.shape[0]),
            ord=2,
        )
    )
    t_f = q_f_tilde.T @ t
    t_star = q_star_tilde.T @ t
    u, singular_values, vh = la.svd(t_f, full_matrices=False)
    sigma_1 = float(singular_values[0])
    sigma_2 = float(singular_values[1])
    gap = sigma_1 - sigma_2
    u_f = u[:, 0]
    v = vh[0, :]
    y = t_star @ v
    tau = float(la.svdvals(t_star)[0])

    steps = sorted(path)
    z_path: dict[int, np.ndarray] = {}
    f_path: dict[int, float] = {}
    rho_direct: dict[int, float] = {}
    angles: dict[int, float] = {}
    z_norms: dict[int, float] = {}
    local_flags: dict[int, bool] = {}
    formula_defects: list[float] = []
    for step in steps:
        w = path[step][f_index, :]
        z = b_half @ (w - w_star) @ s_inverse_half
        z_path[step] = z
        z_norm = float(la.svdvals(z)[0])
        z_norms[step] = z_norm
        angles[step] = float(np.arctan(z_norm))
        epsilon = sigma_1 * (1.0 - 1.0 / np.sqrt(1.0 + z_norm**2))
        epsilon += tau * z_norm / np.sqrt(1.0 + z_norm**2)
        local_flags[step] = epsilon < gap / 2.0
        f_value = graph_factor(z, t_f, t_star)
        f_path[step] = f_value

        graph_basis = np.vstack([z, np.eye(z.shape[1])])
        projector_coordinates = graph_basis @ la.solve(
            graph_basis.T @ graph_basis, graph_basis.T, assume_a="pos"
        )
        q = np.block([q_f_tilde, q_star_tilde])
        projector = q @ projector_coordinates @ q.T
        energy_operator = t.T @ (np.eye(a.shape[0]) - projector) @ t
        direct = float(la.eigvalsh(energy_operator)[-1])
        rho_direct[step] = direct
        formula_defects.append(abs(f_value**2 - direct))

    output_rows: list[dict[str, object]] = []
    sign_matches = 0
    informative = 0
    local_sign_matches = 0
    local_informative = 0
    for left, right in zip(steps[:-1], steps[1:]):
        if right != left + 1:
            continue
        d = z_path[left] - z_path[right]
        prediction = float(u_f @ d @ y)
        actual = f_path[right] - f_path[left]
        scale = max(abs(actual), abs(prediction), 1.0e-15)
        relative_remainder = abs(actual - prediction) / scale
        is_informative = abs(actual) > 1.0e-12 and abs(prediction) > 1.0e-12
        sign_match = bool(is_informative and np.sign(actual) == np.sign(prediction))
        if is_informative:
            informative += 1
            sign_matches += int(sign_match)
        is_local = local_flags[left] and local_flags[right]
        if is_local and is_informative:
            local_informative += 1
            local_sign_matches += int(sign_match)
        output_rows.append(
            {
                "m": left,
                "rho_TG": rho_direct[left],
                "f": f_path[left],
                "theta_max_rad": angles[left],
                "Z_spectral_norm": z_norms[left],
                "local_gap_condition": "yes" if local_flags[left] else "no",
                "actual_delta_f": actual,
                "first_order_prediction": prediction,
                "sign_match": "yes" if sign_match else "no",
                "relative_remainder": relative_remainder,
                "D_frobenius_norm": float(la.norm(d, ord="fro")),
            }
        )

    csv_path = args.results / "experiment4_local_gap_direction.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    local_steps = [step for step in steps if local_flags[step]]
    summary_path = args.results / "experiment4_local_gap_direction.txt"
    with summary_path.open("w", encoding="utf-8") as stream:
        stream.write("Small-scale singular-gap and first-order diagnostic\n")
        stream.write("===================================================\n\n")
        stream.write(f"Grid                         : {args.intervals}/{args.coarse_intervals}\n")
        stream.write("Topology                     : cross-channel\n")
        stream.write("Contrast                     : 1e4\n")
        stream.write(f"sigma_1(T_F)                 : {sigma_1:.12e}\n")
        stream.write(f"sigma_2(T_F)                 : {sigma_2:.12e}\n")
        stream.write(f"simple singular-value gap    : {gap:.12e}\n")
        stream.write(f"relative gap                 : {gap / sigma_1:.12e}\n")
        stream.write(f"||T_*||_2                    : {tau:.12e}\n")
        stream.write(f"endpoint W consistency       : {endpoint_error:.12e}\n")
        stream.write(f"basis orthogonality defect   : {orthogonality_defect:.12e}\n")
        stream.write(f"max formula/direct rho defect: {max(formula_defects):.12e}\n")
        stream.write(
            "local-condition steps         : "
            + (f"{min(local_steps)}..{max(local_steps)}"
               if local_steps else "none")
            + "\n"
        )
        stream.write(
            f"all informative sign matches  : {sign_matches}/{informative}\n"
        )
        stream.write(
            "local informative sign matches: "
            f"{local_sign_matches}/{local_informative}\n\n"
        )
        stream.write(
            "The diagnostic tests the theorem's actual objects.  The gap is "
            "that of T_F at the energy endpoint; theta_max is the largest "
            "A-principal angle via tan(theta_max)=||Z||_2; and the one-step "
            "prediction is u_F^T (Z_m-Z_{m+1}) y for "
            "f=sqrt(rho_TG). The quantities diagnose the local expansion.\n"
        )

    figure, axes = plt.subplots(2, 1, figsize=(7.0, 6.8), sharex=True)
    plot_steps = [int(row["m"]) for row in output_rows]
    axes[0].plot(plot_steps, [row["rho_TG"] for row in output_rows], lw=1.8)
    axes[0].set_ylabel(r"$\rho_{\mathrm{TG}}$")
    axes[0].grid(alpha=0.25)
    axes[1].plot(
        plot_steps, [row["actual_delta_f"] for row in output_rows],
        label="actual", lw=1.6,
    )
    axes[1].plot(
        plot_steps, [row["first_order_prediction"] for row in output_rows],
        label="first order", lw=1.3, ls="--",
    )
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set_xlabel("PCG step m")
    axes[1].set_ylabel(r"$f_{m+1}-f_m$")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(
        args.results / "experiment4_local_gap_direction.png", dpi=220
    )
    plt.close(figure)

if __name__ == "__main__":
    main()
