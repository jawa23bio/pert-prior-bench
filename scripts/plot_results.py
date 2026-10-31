"""One figure per experiment: per-perturbation accuracy for each model.

    python scripts/plot_results.py --results results

Each dot is one held-out perturbation; the dark bar is the model's mean.
Single series, so one neutral hue; models are named on the axis, not by color.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

INK, MUTED, DOT = "#1f2328", "#6b7280", "#5b7fa6"


def strip(df: pd.DataFrame, metric: str, title: str, path: Path):
    models = list(dict.fromkeys(df["model"]))
    fig, ax = plt.subplots(figsize=(7, 0.75 * len(models) + 1.4), dpi=160)
    rng = np.random.default_rng(0)
    for i, m in enumerate(models):
        v = df.loc[df["model"] == m, metric].dropna().values
        y = i + rng.uniform(-0.18, 0.18, len(v))
        ax.scatter(v, y, s=12, color=DOT, alpha=0.45, linewidths=0)
        ax.plot([v.mean()] * 2, [i - 0.3, i + 0.3], color=INK, lw=2, solid_capstyle="round")
        ax.text(v.mean(), i - 0.34, f"mean {v.mean():.2f}", ha="center", va="bottom",
                fontsize=8, color=INK, bbox=dict(facecolor="white", edgecolor="none", pad=1))
    ax.set_yticks(range(len(models)), models, fontsize=9, color=INK)
    ax.invert_yaxis()
    ax.set_xlabel(f"{metric}  (higher is better)", fontsize=9, color=MUTED)
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(axis="x", colors=MUTED, labelsize=8)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#e5e7eb", lw=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--metric", default="pearson_top20")
    a = ap.parse_args()
    r = Path(a.results)
    for fname, title in [("singles_per_perturbation.csv", "Unseen single perturbations"),
                         ("doubles_per_perturbation.csv", "Unseen double perturbations")]:
        if (r / fname).exists():
            out = r / fname.replace("_per_perturbation.csv", f"_{a.metric}.png")
            strip(pd.read_csv(r / fname), a.metric, title, out)
            print("wrote", out)


if __name__ == "__main__":
    main()
