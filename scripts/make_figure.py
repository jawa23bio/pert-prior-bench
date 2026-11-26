"""Overview figure for the README: does the network help, and could we have seen it?

    python scripts/make_figure.py --out docs/figures/overview.png

Reads only committed result files:
  results_v2/singles_mean.csv          per-perturbation scores (for the intervals)
  results_v2/null_draws.csv            score of each permuted-network draw
  results_v2_direct/null_draws.csv     the 99 draws for direct[prior]
  results_v2/diagnostics.csv           which TFs have a measured target
  results_noise_ceiling/noise_ceiling_per_perturbation.csv
  results_spike_in/spike_in.csv        detection limit
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pertbench.protocol import boot_ci  # noqa: E402

# colors: categorical slots 1-2 (validated light-mode pair), neutrals for baselines and nulls
BLUE, ORANGE = "#257de0", "#f47230"
NULL, BASE = "#b3b3b3", "#8c8c8c"
TEXT, TEXT2, GRID, SURFACE = "#191919", "#666666", "#e6e6e6", "#ffffff"
METRIC = "pearson_ctr_top20"


def load(root: Path):
    sm = pd.read_csv(root / "results_v2/singles_mean.csv")
    cov = set(pd.read_csv(root / "results_v2/diagnostics.csv")["perturbed_tf"])
    sm = sm[sm["perturbation"].isin(cov)]

    def ci(model, network):
        x = sm[(sm["model"] == model) & (sm["network"] == network)][METRIC].values
        return boot_ci(x)

    draws = pd.read_csv(root / "results_v2/null_draws.csv")
    direct = pd.read_csv(root / "results_v2_direct/null_draws.csv")
    nc = pd.read_csv(root / "results_noise_ceiling/noise_ceiling_per_perturbation.csv")
    ceiling = nc[nc["perturbation"].isin(cov)][METRIC].mean()
    spike = pd.read_csv(root / "results_spike_in/spike_in.csv")
    tests = pd.read_csv(root / "results_v2_direct/null_tests.csv")
    p_direct = tests[(tests["null"] == "regulator-permuted") & (tests["metric"] == METRIC)]["p"].iloc[0]
    return dict(n=len(cov), expr=ci("linear[expression]", "none"), real=ci("linear[both]", "real"),
                perm=draws[(draws["model"] == "linear[both]") & (draws["network"] == "regulator-permuted")][METRIC].values,
                direct=ci("direct[prior]", "real")[0],
                direct_null=direct[(direct["model"] == "direct[prior]")
                                   & (direct["network"] == "regulator-permuted")][METRIC].values,
                p_direct=p_direct, ceiling=ceiling, spike=spike)


def style(ax):
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=TEXT2, labelsize=8, length=0)
    ax.set_facecolor(SURFACE)


def panel_models(ax, d):
    rows = ["Mean of training\nperturbations", "No network", "With CollecTRI\nnetwork",
            "With permuted\nnetworks (20)", "Network edges\nread directly"]
    y = np.arange(len(rows))[::-1]
    rng = np.random.default_rng(1)
    ax.set_xlim(-0.3, 1.05)
    for x in [0, 0.25, 0.5, 0.75, 1.0]:
        ax.axvline(x, color="#cfcfcf" if x == 0 else GRID, lw=0.8, zorder=0)
    # replicate ceiling
    ax.axvline(d["ceiling"], color=TEXT2, lw=1.2, zorder=1)
    ax.text(d["ceiling"] - 0.02, y[0] + 0.55, f"two halves of the\nsame cells agree: {d['ceiling']:.2f}",
            ha="right", va="top", fontsize=7.5, color=TEXT2)

    def point(yy, m, lo, hi, color, label):
        if lo is not None:
            ax.plot([lo, hi], [yy, yy], color=color, lw=2, solid_capstyle="round", zorder=2)
        ax.scatter([m], [yy], s=46, color=color, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.text((hi if hi is not None else m) + 0.03, yy, label, va="center", fontsize=8, color=TEXT)

    point(y[0], 0.0, None, None, BASE, "0 (by definition)")
    m, lo, hi = d["expr"]
    point(y[1], m, lo, hi, BLUE, f"{m:.2f}")
    m, lo, hi = d["real"]
    point(y[2], m, lo, hi, ORANGE, f"{m:.2f}")
    ax.scatter(d["perm"], y[3] + rng.uniform(-0.13, 0.13, len(d["perm"])), s=10, color=NULL, zorder=2,
               linewidth=0)
    ax.text(d["perm"].max() + 0.03, y[3], f"mean {d['perm'].mean():.2f}", va="center", fontsize=8, color=TEXT2)
    ax.scatter(d["direct_null"], y[4] + rng.uniform(-0.13, 0.13, len(d["direct_null"])), s=6, color=NULL,
               zorder=2, linewidth=0)
    ax.scatter([d["direct"]], [y[4]], s=46, color=ORANGE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.text(d["direct_null"].max() + 0.03, y[4], f"\u2248 0  (p = {d['p_direct']:.2f} vs 99 permuted)",
            va="center", fontsize=8, color=TEXT)
    ax.set_yticks(y, rows, fontsize=8, color=TEXT)
    ax.set_ylim(y[-1] - 0.6, y[0] + 0.7)
    ax.set_xlabel(f"Correlation with each perturbation's specific effect ({d['n']} TFs, 95% CI)",
                  fontsize=8, color=TEXT2)
    ax.set_title("Adding the network does not help", loc="left", fontsize=10, color=TEXT,
                 fontweight="bold", pad=8)
    style(ax)


def panel_detection(ax, d):
    s = d["spike"]
    b = s["beta"].values
    real, null = s["v2 direct[prior] real"].values, s["v2 direct[prior] null"].values
    hit = s["v2 direct[prior] detects"].astype(bool).values
    ax.set_xlim(-0.05, 1.08)
    ax.set_ylim(-0.14, 1.0)
    for v in [0, 0.25, 0.5, 0.75]:
        ax.axhline(v, color="#cfcfcf" if v == 0 else GRID, lw=0.8, zorder=0)
    ax.plot(b, null, color=NULL, lw=2, zorder=1)
    ax.scatter(b, null, s=30, color=NULL, edgecolor=SURFACE, linewidth=1.5, zorder=2)
    ax.plot(b, real, color=ORANGE, lw=2, zorder=1)
    ax.scatter(b[hit], real[hit], s=46, color=ORANGE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.scatter(b[~hit], real[~hit], s=46, facecolor=SURFACE, edgecolor=ORANGE, linewidth=1.8, zorder=3)
    ax.text(b[-1] + 0.02, real[-1], "real\nnetwork", va="center", fontsize=8, color=TEXT)
    ax.text(b[-1] + 0.02, null[-1], "permuted\n(19)", va="center", fontsize=8, color=TEXT2)
    first = b[hit & (b > 0)].min()
    i = int(np.where(b == first)[0][0])
    ax.annotate(f"detected at {first:g}", xy=(first, real[i]), xytext=(first + 0.13, real[i] - 0.02),
                fontsize=8, color=TEXT, va="center",
                arrowprops=dict(arrowstyle="-", color=TEXT2, lw=0.8, shrinkA=2, shrinkB=5))
    ax.text(0, -0.085, "the actual data", ha="left", fontsize=7.5, color=TEXT2)
    v1 = s["beta"][s["v1 detects"].astype(bool) & (s["beta"] > 0)]
    if len(v1):
        ax.text(0.02, 0.97, f"filled = detected (p \u2264 0.05)\nthe original analysis first\ndetected it at {v1.min():g}",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.5, color=TEXT2, linespacing=1.4)
    ax.set_xticks(b, [f"{x:g}" for x in b], fontsize=8)
    ax.set_xlabel("Planted network effect, relative to each TF's own effect", fontsize=8, color=TEXT2)
    ax.set_ylabel("Direct network read-out, same score", fontsize=8, color=TEXT2)
    ax.set_title("An effect a tenth that size would have been caught", loc="left", fontsize=10, color=TEXT,
                 fontweight="bold", pad=8)
    style(ax)
    ax.spines["left"].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--out", default=str(ROOT / "docs/figures/overview.png"))
    a = ap.parse_args()
    d = load(Path(a.root))
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 3.9), gridspec_kw=dict(width_ratios=[1.45, 1]),
                                      facecolor=SURFACE)
    panel_models(left, d)
    panel_detection(right, d)
    fig.text(0.01, 0.985, "Does a curated TF network (CollecTRI) help predict unseen perturbations?  "
             "Norman et al. 2019 CRISPRa Perturb-seq, K562", fontsize=8, color=TEXT2, va="top")
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=3)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=200, facecolor=SURFACE)
    fig.savefig(str(Path(a.out).with_suffix(".svg")), facecolor=SURFACE)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
