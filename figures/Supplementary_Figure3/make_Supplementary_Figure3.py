import json
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import plot_style as style
from plot_style import C, panel, ptitle

DATA = os.path.join(HERE, "data")
PROTECTION_DATA = os.path.join(os.path.dirname(HERE), "Figure5", "data")
style.set_style()
WS = [pe.withStroke(linewidth=2.0, foreground="white")]
INC_C = {"HIC": "#264F73", "UMC": "#5E9BC8", "LMC": "#E39A3B", "LIC": "#A63B2A"}
FW = 7.2


def ax_at(fig, x, y, w, h, HF):
    return fig.add_axes([x / FW, y / HF, w / FW, h / HF])


S = json.load(open(os.path.join(DATA, "protection_sensitivity.json")))
rps = ["RP50", "RP100", "RP200", "RP500"]
d_cty = pd.read_csv(os.path.join(PROTECTION_DATA, "protection_by_country.csv")).head(15)


def render_a(ax):
    x = np.arange(4)
    ax.bar(x - 0.2, [S[r]["nominal"] / 1000 for r in rps], 0.38, color=C["lgray"], label="nominally exposed")
    ax.bar(
        x + 0.2, [S[r]["effective"] / 1000 for r in rps], 0.38, color=C["risk"], label="effectively exposed"
    )
    for xi, r in zip(x, rps):
        ax.text(
            xi - 0.2,
            S[r]["nominal"] / 1000 + 4,
            f"{S[r]['nominal'] / 1000:.0f}k",
            ha="center",
            fontsize=6.1,
            color=C["spine"],
        )
        ax.text(
            xi + 0.2,
            S[r]["effective"] / 1000 + 4,
            f"{S[r]['effective'] / 1000:.0f}k",
            ha="center",
            fontsize=6.1,
            color=C["risk_dark"],
            fontweight="bold",
        )
        ax.text(
            xi,
            190,
            f"{S[r]['pct_protected']:.1f}%",
            ha="center",
            fontsize=6.6,
            color=C["ink"],
            fontweight="bold",
        )
    ax.text(
        1.5,
        208,
        "share of nominal exposure behind adequate protection",
        ha="center",
        fontsize=6.1,
        color=C["ink"],
    )
    ax.set_xticks(x)
    ax.set_xticklabels(["50-yr", "100-yr (main)", "200-yr", "500-yr"], fontsize=6.8)
    ax.set_ylim(0, 222)
    ax.set_ylabel("Chargers (thousands)", fontsize=6.8)
    ax.legend(frameon=False, fontsize=6.1, loc="upper left", bbox_to_anchor=(0.0, 0.88), labelspacing=0.3)
    ptitle(ax, "Protection vs event magnitude", pad=4)


def render_b(ax):
    d = d_cty
    ax.plot([500, 150000.0], [500, 150000.0], color=C["lgray"], lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.scatter(
        d.nominal,
        d.effective.clip(lower=0.3),
        s=26,
        c=[INC_C.get(g, "#999999") for g in d.inc_grp],
        edgecolor="white",
        lw=0.5,
        zorder=3,
    )
    sh = {"United States of America": "United States"}
    off = {
        "China": (-6, 5),
        "Netherlands": (6, 7),
        "Japan": (7, -2),
        "United States of America": (-50, -13),
        "India": (-22, -3),
        "Germany": (5, 4),
        "Thailand": (4, -3),
        "France": (4, -9),
    }
    for _, r in d.iterrows():
        if r.nominal >= 3000 or r.effective == 0:
            ax.annotate(
                sh.get(r.country, r.country),
                (r.nominal, max(r.effective, 0.3)),
                textcoords="offset points",
                xytext=off.get(r.country, (4, 3)),
                fontsize=5.9,
                color=C["ink"],
                path_effects=WS,
                zorder=5,
                arrowprops=(
                    dict(arrowstyle="-", lw=0.4, color=C["gray"], shrinkA=0, shrinkB=2)
                    if r.country.startswith("United States")
                    else None
                ),
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(700, 160000.0)
    ax.set_ylim(0.2, 160000.0)
    ax.set_yticks([0.3, 1, 10, 100, 1000.0, 10000.0, 100000.0])
    ax.set_yticklabels(["0", "1", "10", "100", "1k", "10k", "100k"], fontsize=6.6)
    ax.minorticks_off()
    ax.set_xticks([1000.0, 10000.0, 100000.0])
    ax.set_xticklabels(["1k", "10k", "100k"], fontsize=6.6)
    ax.set_xlabel("Nominally exposed chargers (log)", fontsize=6.8)
    ax.set_ylabel("Effectively exposed (log)", fontsize=6.8)
    hs = [
        plt.Line2D([], [], marker="o", ls="", mfc=INC_C[k], mec="white", ms=4.6, label=l)
        for k, l in [("HIC", "High income"), ("UMC", "Upper-middle"), ("LMC", "Lower-middle")]
    ]
    ax.legend(
        handles=hs,
        frameon=False,
        fontsize=6.0,
        loc="upper left",
        handletextpad=0.2,
        borderaxespad=0.05,
        labelspacing=0.3,
    )
    ptitle(ax, "Country-level adjustment (top 15)", pad=4)


def composite():
    HF = 2.5
    fig = plt.figure(figsize=(FW, HF))
    a = ax_at(fig, 0.55, 0.44, 2.6, 1.55, HF)
    render_a(a)
    panel(a, "a", dx=-0.17, dy=1.15)
    b = ax_at(fig, 4.1, 0.44, 2.5, 1.55, HF)
    render_b(b)
    panel(b, "b", dx=-0.26, dy=1.15)
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Supplementary_Figure3"))


if __name__ == "__main__":
    main()
