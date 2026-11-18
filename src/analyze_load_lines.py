"""Load-line analysis of a carbon-fibre tibial prosthesis blade.

For each compression test (load applied at a different distance from the
blade edge, i.e. a different blade angle), the script computes:
  * the peak force and displacement,
  * the energy stored while loading and the energy returned while unloading,
  * the hysteresis loss and the energy-return ratio,
  * an average stiffness,
  * a quadratic least-squares fit of the last points of the raw unloading
    curve, used to correct the zero offset of the test bench.
"""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

TESTS = [
    # folder, angle (deg), load line distance from blade edge (mm), colour
    ("test_1_20deg_10mm", 20, 10, "deepskyblue"),
    ("test_2_17deg_8mm", 17, 8, "darkviolet"),
    ("test_3_14deg_6mm", 14, 6, "darkorange"),
]
N_FIT = 20  # number of last raw points used for the quadratic fit


def read_columns(path, ncols, skip=2):
    """Read whitespace-separated numeric columns, skipping the header lines."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    rows = [l.split() for l in lines[skip:] if l.strip()]
    return np.array([[float(x.replace(",", ".")) for x in r[:ncols]] for r in rows])


def energies(disp_mm, force_n):
    """Energy stored (loading) and returned (unloading), in joules."""
    i_max = int(np.argmax(force_n))
    trapz = getattr(np, "trapezoid", None) or np.trapz
    stored = trapz(force_n[: i_max + 1], disp_mm[: i_max + 1]) * 1e-3
    returned = -trapz(force_n[i_max:], disp_mm[i_max:]) * 1e-3
    return stored, returned


def quadratic_fit(raw_disp, raw_force):
    """Least squares fit F = a x^2 + b x + c on the last N_FIT raw points."""
    x, y = raw_disp[-N_FIT:], raw_force[-N_FIT:]
    design = np.vstack([x**2, x, np.ones_like(x)]).T
    model = LinearRegression(fit_intercept=False).fit(design, y)
    return model.coef_, model.score(design, y)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    ap.add_argument("--show", action="store_true", help="display the figures")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = []
    fig_all, ax_all = plt.subplots()

    for folder, angle, mm, colour in TESTS:
        d = Path(args.data) / folder
        raw = read_columns(d / "raw.txt", 4)
        cor = read_columns(d / "corrected.txt", 2)
        raw_disp, raw_force = -raw[:, 2], -raw[:, 3]
        disp, force = -cor[:, 0], -cor[:, 1]

        stored, returned = energies(disp, force)
        coeffs, r2 = quadratic_fit(raw_disp, raw_force)
        summary.append({
            "angle_deg": angle,
            "load_line_mm": mm,
            "peak_force_N": round(float(force.max()), 1),
            "max_displacement_mm": round(float(disp.max()), 2),
            "stiffness_N_per_m": round(float(force.max() / disp.max() * 1e3)),
            "energy_stored_J": round(float(stored), 2),
            "energy_returned_J": round(float(returned), 2),
            "energy_lost_J": round(float(stored - returned), 2),
            "energy_return_pct": round(float(100 * returned / stored), 1),
            "fit_a": round(float(coeffs[0]), 4),
            "fit_b": round(float(coeffs[1]), 4),
            "fit_c": round(float(coeffs[2]), 4),
            "fit_r2": round(float(r2), 5),
        })

        fig, ax = plt.subplots()
        ax.plot(disp, force, color=colour, label=f"Load line {mm} mm ({angle}°)")
        ax.set(title="Compression test", xlabel="Displacement (mm)", ylabel="Force (N)")
        ax.grid(True); ax.legend(loc="best")
        fig.savefig(out / f"load_line_{mm}mm.png", dpi=150, bbox_inches="tight")
        ax_all.plot(disp, force, color=colour, label=f"{mm} mm ({angle}°)")
        if not args.show:
            plt.close(fig)

    ax_all.set(title="Load lines compared", xlabel="Displacement (mm)", ylabel="Force (N)")
    ax_all.grid(True); ax_all.legend(loc="best")
    fig_all.savefig(out / "load_lines_comparison.png", dpi=150, bbox_inches="tight")

    with open(out / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]))
        w.writeheader(); w.writerows(summary)
    for row in summary:
        print(row)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
