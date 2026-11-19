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
style.set_style()
WS = [pe.withStroke(linewidth=2.0, foreground="white")]
FW = 7.2


def ax_at(fig, x, y, w, h, HF):
    return fig.add_axes([x / FW, y / HF, w / FW, h / HF])


d_iea = pd.read_csv(os.path.join(DATA, "crosscheck_stations_vs_iea.csv")).sort_values("points_per_station")
d_cov = pd.read_csv(os.path.join(DATA, "coverage_table.csv"))
_c1 = pd.read_csv(os.path.join(DATA, "cdf_charger_to_power_line.csv"))
_c2 = pd.read_csv(os.path.join(DATA, "cdf_charger_to_feeder_site.csv"))
_cs = json.load(open(os.path.join(DATA, "coupling_summary.json")))
att, cdf_a = (_c1.distance_km.values, _c1.cumulative_share.values)
km, cdf = (_c2.distance_km.values, _c2.cumulative_share.values)
med_a = _cs["median_attach_km"]
med = _cs["median_network_km_within_10km"]
within10 = _cs["within_10km_pct"]


def render_a(ax):
    d = d_iea
    sh = {"United States of America": "United States"}
    ax.barh(np.arange(len(d)), d.points_per_station, 0.62, color=C["ev"])
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels([sh.get(c, c) for c in d.country], fontsize=6.8)
    for i, v in enumerate(d.points_per_station):
        ax.text(v + 0.35, i, f"{v:.1f}", va="center", fontsize=6.6, color=C["ink"])
    ax.set_xlim(0, 19.5)
    ax.set_xlabel("IEA charging points per\nmapped station (2025)", fontsize=6.8)
    ptitle(ax, "IEA cross-audit", pad=4)


def render_b(ax):
    d = d_cov
    ok = d[d.qa_pass]
    ko = d[~d.qa_pass]
    ax.scatter(
        ok.max_year_share * 100,
        ok.n,
        s=16,
        facecolor=C["ev"],
        edgecolor="white",
        lw=0.4,
        zorder=4,
        label=f"pass ({len(ok)})",
    )
    ax.scatter(
        ko.max_year_share * 100,
        ko.n,
        s=14,
        facecolor=C["gray"],
        edgecolor="white",
        lw=0.4,
        zorder=3,
        label="fail",
    )
    ax.axvline(55, color=C["risk"], lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.axhline(1200, color=C["risk"], lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.text(58, 155000.0, "single-year ≤55%", fontsize=5.9, color=C["risk_dark"])
    ax.text(2.5, 1450, "n≥1,200", fontsize=5.9, color=C["risk_dark"], ha="left", va="bottom")
    for nm, lab, xy in [
        ("Netherlands", "Netherlands", (-2, -11)),
        ("India", "India", (3, 4)),
        ("United States of America", "United States", (-46, -12)),
        ("United Kingdom", "UK", (5, 6)),
    ]:
        r = d[d.country == nm]
        if len(r):
            r = r.iloc[0]
            ax.annotate(
                lab,
                (r.max_year_share * 100, r.n),
                textcoords="offset points",
                xytext=xy,
                fontsize=6.0,
                color=C["ink"],
                path_effects=WS,
                zorder=6,
            )
    ax.set_yscale("log")
    ax.set_xlim(0, 104)
    ax.set_ylim(80, 360000.0)
    ax.set_xlabel("Largest single registration year (%)", fontsize=6.8)
    ax.set_ylabel("Dated stations (log)", fontsize=6.8)
    ax.legend(
        frameon=False, fontsize=6.0, loc="upper left", handletextpad=0.2, borderaxespad=0.02, labelspacing=0.3
    )
    ptitle(ax, "OCM coverage gates", pad=4)


def render_c(ax):
    ax.plot(att, cdf_a * 100, color=C["ev"], lw=1.3)
    ax.plot(km, cdf * 100, color=C["ev_dark"], lw=1.5)
    ax.set_xscale("log")
    ax.set_xlim(0.02, 60)
    ax.set_ylim(0, 103)
    ax.axvline(10, color=C["risk"], lw=0.8, ls=(0, (3, 2)))
    ax.text(
        10 * 1.2,
        34,
        f"10 km cap\n{within10:.1f}% network-\nassigned",
        fontsize=6.0,
        color=C["risk_dark"],
        path_effects=WS,
    )
    ax.text(0.028, 84, "to power line", fontsize=6.2, color=C["ev"], fontweight="bold", path_effects=WS)
    ax.text(0.028, 73, f"median {med_a * 1000:.0f} m", fontsize=5.8, color=C["ev"], path_effects=WS)
    ax.text(
        2.6,
        16.5,
        "to feeder site\n(network)",
        fontsize=6.2,
        color=C["ev_dark"],
        fontweight="bold",
        path_effects=WS,
        linespacing=1.15,
    )
    ax.text(2.6, 5, f"median {med:.2f} km", fontsize=5.8, color=C["ev_dark"], path_effects=WS)
    ax.set_xlabel("Distance (km, log)", fontsize=6.8)
    ax.set_ylabel("Cumulative share of chargers (%)", fontsize=6.8)
    ptitle(ax, "Charger–grid coupling", pad=4)


def composite():
    HF = 2.45
    fig = plt.figure(figsize=(FW, HF))
    a = ax_at(fig, 0.72, 0.48, 1.44, 1.55, HF)
    render_a(a)
    panel(a, "a", dx=-0.42, dy=1.16)
    b = ax_at(fig, 2.92, 0.48, 1.52, 1.55, HF)
    render_b(b)
    panel(b, "b", dx=-0.34, dy=1.16)
    c = ax_at(fig, 5.42, 0.48, 1.5, 1.55, HF)
    render_c(c)
    panel(c, "c", dx=-0.36, dy=1.16)
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Supplementary_Figure1"))


if __name__ == "__main__":
    main()
