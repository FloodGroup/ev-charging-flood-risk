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
CASCADE_DATA = os.path.join(os.path.dirname(HERE), "Figure4", "data")
style.set_style()
FW = 7.2
GRID_C = "#6A51A3"
ISG = json.load(open(f"{DATA}/islanding_global.json"))
ISC = pd.read_csv(f"{DATA}/islanding_country.csv", index_col=0)
rp = pd.read_csv(f"{CASCADE_DATA}/cascade_rp_curve.csv")
mtn = pd.read_csv(f"{CASCADE_DATA}/cascade_metro_named.csv")
r100 = rp.loc[rp.rp == 100].iloc[0]


def ax_at(fig, x, y, w, h, HF):
    return fig.add_axes([x / FW, y / HF, w / FW, h / HF])


def render_a(ax):
    rows = [
        ("Deterministic thresholds", float(r100.A_threshold), None, "#9A9A9A"),
        (
            "Fragility curves",
            float(r100.A_fragility),
            (
                float(min(r100.A_hard, r100.A_beta04, r100.A_beta08)),
                float(max(r100.A_soft, r100.A_beta04, r100.A_beta08)),
            ),
            GRID_C,
        ),
        ("+ backup feeder", float(r100.A_redundant), None, "#0072B2"),
    ]
    y = np.arange(len(rows))[::-1]
    for yi, (lab, v, band, col) in zip(y, rows):
        ax.barh(yi, v - 1.0, left=1.0, height=0.58, color=col, alpha=0.92)
        if band:
            ax.plot(band, [yi, yi], color="#1A1A1A", lw=1.0, zorder=5)
            for b in band:
                ax.plot([b, b], [yi - 0.16, yi + 0.16], color="#1A1A1A", lw=1.0, zorder=5)
        ax.text(
            max(v, band[1] if band else v) + 0.006,
            yi,
            f"×{v:.3f}",
            fontsize=6.6,
            va="center",
            color="#1A1A1A",
        )
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.0)
    ax.set_xlim(1.0, 1.3)
    ax.set_xlabel("Grid amplification ×A (RP100)", fontsize=6.8)
    ax.tick_params(axis="y", length=0)
    ptitle(ax, "Failure-model comparison", pad=4)


def render_b(ax):
    m = ISC[(ISC.n >= 1000) & (ISC.direct >= 100)].copy()
    m = m.nlargest(12, "n").sort_values("A_service")
    short = {"United States of America": "USA", "United Kingdom": "UK"}
    names = [short.get(i, i) for i in m.index]
    y = np.arange(len(m))
    for yi, s_, c_ in zip(y, m.A_service, m.A_conn):
        ax.plot([s_, c_], [yi, yi], color=C["lgray"], lw=1.2, zorder=2)
    ax.scatter(
        m.A_service,
        y,
        s=17,
        facecolor="white",
        edgecolor=C["ev_dark"],
        lw=0.9,
        zorder=3,
        label="feeder cascade",
    )
    ax.scatter(
        m.A_conn,
        y,
        s=19,
        facecolor=C["risk"],
        edgecolor="white",
        lw=0.5,
        zorder=4,
        label="+ source connectivity",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.2)
    xmax = float(max(1.65, m.A_conn.max() + 0.08))
    ax.set_xlim(0.98, xmax)
    ax.set_xlabel("×A at P_fail ≥ 0.5 (RP100)", fontsize=6.8)
    ax.legend(
        frameon=False,
        fontsize=5.8,
        loc="lower right",
        borderaxespad=0.05,
        handletextpad=0.2,
        labelspacing=0.3,
    )
    dA = float(ISG["delta_A"])
    dA_txt = f"{dA:.4f}" if dA < 0.01 else f"{dA:.3f}"
    ax.text(
        0.03,
        0.985,
        f"global: +{int(ISG['islanded_add']):,}\nΔA = {dA_txt}",
        transform=ax.transAxes,
        fontsize=5.9,
        va="top",
        color="#3C3C3C",
    )
    ptitle(ax, "Islanding, global model", pad=4)


def render_c(ax):
    d_met = mtn[(mtn.direct >= 50) & mtn.A_total.notna()]
    med2, p95 = (d_met.A_total.median(), d_met.A_total.quantile(0.95))
    ax.hist(
        d_met.A_total.clip(upper=2.6),
        bins=np.arange(1.0, 2.65, 0.05),
        color=C["ev"],
        edgecolor="white",
        lw=0.35,
    )
    yl = ax.get_ylim()[1]
    ax.axvline(med2, color=C["ink"], lw=0.8, ls=(0, (4, 2)))
    ax.axvline(p95, color=C["risk_dark"], lw=0.8, ls=(0, (4, 2)))
    ax.text(med2 + 0.03, yl * 0.66, f"median ×{med2:.2f}", fontsize=6.3, color=C["ink"])
    ax.text(p95 + 0.03, yl * 0.46, f"P95 ×{p95:.2f}", fontsize=6.3, color=C["risk_dark"])
    ax.set_xlim(1.0, 2.62)
    ax.set_xlabel("Amplification ×A (RP100)", fontsize=6.8)
    ax.set_ylabel("Urban centres", fontsize=6.8)
    ax.text(
        0.97,
        0.97,
        f"{len(d_met)} centres\n≥50 flooded",
        transform=ax.transAxes,
        fontsize=5.9,
        va="top",
        ha="right",
        color="#3C3C3C",
    )
    ptitle(ax, "Urban-centre spread", pad=4)
    return (med2, p95, len(d_met))


def composite():
    HF = 2.55
    fig = plt.figure(figsize=(FW, HF))
    a = ax_at(fig, 1.1, 0.46, 1.55, 1.66, HF)
    render_a(a)
    panel(a, "a", dx=-0.66, dy=1.13)
    b = ax_at(fig, 3.5, 0.46, 1.55, 1.66, HF)
    render_b(b)
    panel(b, "b", dx=-0.44, dy=1.13)
    c = ax_at(fig, 5.62, 0.46, 1.43, 1.66, HF)
    stats_c = render_c(c)
    panel(c, "c", dx=-0.26, dy=1.13)
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Supplementary_Figure2"))


if __name__ == "__main__":
    main()
