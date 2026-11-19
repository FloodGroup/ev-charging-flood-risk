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
from plot_style import C, draw_world, panel, ptitle, set_style, world

DATA = os.path.join(HERE, "data")
style.set_style()
plt.rcParams.update(
    {
        "font.size": 8,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.0,
        "axes.titlesize": 8.4,
        "axes.labelsize": 8,
    }
)
GAS = "#B3542F"
INC_C = {"HIC": "#264F73", "UMC": "#5E9BC8", "LMC": "#E39A3B", "LIC": "#A63B2A"}
INC_LAB = {"HIC": "High income", "UMC": "Upper-middle", "LMC": "Lower-middle", "LIC": "Low income"}
RAMP = ["#F7ECE1", "#F0C9A0", "#E08A45", "#BC4B12", "#7E2708"]
PROT_C = {
    "<10": "#7E2708",
    "10-50": "#D55E00",
    "50-100": "#F2B387",
    "100-1000": "#9FBFDD",
    ">=1000": "#3D74A6",
}
WS = [pe.withStroke(linewidth=2.2, foreground="white")]
XMAX = 17005000.0
XL, XR = (-XMAX * (170 / 180), XMAX)
inc = pd.read_csv(os.path.join(DATA, "income_group_exposure.csv")).set_index("inc_grp")
cc = pd.read_csv(os.path.join(DATA, "concentration_curve.csv"))
ci = json.load(open(os.path.join(DATA, "concentration_indices.json")))
pbc = pd.read_csv(os.path.join(DATA, "protection_by_country.csv"))
pbn = pd.read_csv(os.path.join(DATA, "protection_bins_income.csv")).set_index("inc_grp")
qw = pd.read_csv(os.path.join(DATA, "quintile_within_country.csv"))
qp = pd.read_csv(os.path.join(DATA, "quintile_pooled.csv"))
inh = pd.read_csv(os.path.join(DATA, "inheritance_country.csv"))
S = json.load(open(os.path.join(DATA, "fig05_stats.json")))


def render_a(ax):
    w = world().copy()
    w["v"] = w["NAME"].map(inh.set_index("country")["gas_eff_pct"])
    draw_world(ax, graticule=True)
    eff = inh.set_index("country")["gas_eff_pct"]
    pc = pbc.set_index("country")
    bins = [0, 2, 5, 10, 20, 1000]
    for i in range(5):
        sel = w[(w["v"] >= bins[i]) & (w["v"] < bins[i + 1])]
        if len(sel):
            sel.plot(ax=ax, facecolor=RAMP[i], edgecolor="white", linewidth=0.25, zorder=3)
    pts = cc[cc["eff_ev"] > 100].copy()
    geo = w.set_index("NAME").geometry
    pts["g"] = pts["country"].map(geo)
    pts = pts.dropna(subset=["g"])
    pts["x"] = pts["g"].apply(lambda p: p.representative_point().x)
    pts["y"] = pts["g"].apply(lambda p: p.representative_point().y)
    smax = float(pts["eff_ev"].max())
    ax.scatter(
        pts["x"],
        pts["y"],
        s=4 + 260 * pts["eff_ev"] / smax,
        facecolor=C["ev"],
        edgecolor="white",
        linewidth=0.5,
        alpha=0.88,
        zorder=6,
    )
    ax.set_xlim(XL, XR)
    ax.set_ylim(-6350000.0, 7850000.0)

    def note(txt, xy, xt, fs=6.8):
        ax.annotate(
            txt,
            xy=xy,
            xytext=xt,
            fontsize=fs,
            color=C["ink"],
            path_effects=WS,
            arrowprops=dict(arrowstyle="-", lw=0.55, color=C["gray"], shrinkA=0, shrinkB=0),
            zorder=8,
            ha="left",
            va="center",
        )

    gx = lambda n: (geo[n].representative_point().x, geo[n].representative_point().y)
    note(f"China {int(pc.loc['China', 'effective']):,}", gx("China"), (13600000.0, 6800000.0))
    note(
        f"Netherlands {int(pc.loc['Netherlands', 'nominal']):,} → {int(pc.loc['Netherlands', 'effective']):,}",
        gx("Netherlands"),
        (-8900000.0, 7000000.0),
    )
    note(f"Vietnam {eff['Vietnam']:.0f}%", gx("Vietnam"), (11900000.0, -1700000.0))
    note(f"Bangladesh {eff['Bangladesh']:.0f}%", gx("Bangladesh"), (9600000.0, -3300000.0))
    note(f"Pakistan {eff['Pakistan']:.0f}%", gx("Pakistan"), (3900000.0, -2500000.0))
    note(f"Egypt {eff['Egypt']:.0f}%", gx("Egypt"), (-3300000.0, 1100000.0))
    x0, yl = (XL + 900000.0, -5050000.0)
    ax.text(
        x0,
        yl + 950000.0,
        "Fuel network eff. exposed (%)",
        fontsize=6.1,
        color=C["ink"],
        zorder=7,
        path_effects=WS,
    )
    for i, lb in enumerate(["0–2", "2–5", "5–10", "10–20", ">20", "n/a"]):
        fc = RAMP[i] if i < 5 else "#E7E7E2"
        ax.add_patch(
            plt.Rectangle((x0 + i * 1840000.0, yl), 1760000.0, 520000.0, fc=fc, ec="white", lw=0.3, zorder=7)
        )
        ax.text(
            x0 + i * 1840000.0 + 880000.0,
            yl - 300000.0,
            lb,
            ha="center",
            va="top",
            fontsize=5.8,
            color=C["ink"],
            zorder=7,
        )
    bx = x0 + 13200000.0
    ax.text(
        bx - 500000.0,
        yl + 1150000.0,
        "Chargers eff. exposed",
        fontsize=6.1,
        color=C["ink"],
        zorder=7,
        path_effects=WS,
    )
    for v, lab, dx in [
        (1000, "1k", 0),
        (10000, "10k", 1500000.0),
        (int(smax), f"{smax / 1000:.0f}k", 3400000.0),
    ]:
        ax.scatter(
            [bx + dx],
            [yl + 250000.0],
            s=4 + 260 * v / smax,
            facecolor=C["ev"],
            edgecolor="white",
            lw=0.5,
            zorder=7,
        )
        ax.text(bx + dx, yl - 950000.0, lab, ha="center", va="top", fontsize=5.8, color=C["ink"], zorder=7)


def render_b(ax):
    g = inc.loc[["HIC", "UMC", "LMC", "LIC"]]
    x = np.arange(4)
    ax.bar(x - 0.19, g["ev_pct"], 0.36, color=C["ev"], label="Charging")
    ax.bar(x + 0.19, g["gas_pct"], 0.36, color=GAS, label="Fuel")
    ax.errorbar(
        x - 0.19,
        g["ev_pct"],
        yerr=[g["ev_pct"] - g["ev_pct_lo"], g["ev_pct_hi"] - g["ev_pct"]],
        fmt="none",
        ecolor=C["ink"],
        elinewidth=0.6,
        capsize=1.4,
    )
    ax.errorbar(
        x + 0.19,
        g["gas_pct"],
        yerr=[g["gas_pct"] - g["gas_pct_lo"], g["gas_pct_hi"] - g["gas_pct"]],
        fmt="none",
        ecolor=C["ink"],
        elinewidth=0.6,
        capsize=1.4,
    )
    for i, r in enumerate(g["ratio"]):
        yv = max(g["ev_pct_hi"].iloc[i], g["gas_pct_hi"].iloc[i]) + 0.8
        ax.text(
            i,
            yv,
            f"×{r:.2f}",
            ha="center",
            fontsize=7.2,
            fontweight="bold",
            color=C["risk"] if r > 1 else C["ink"],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(["High\nincome", "Upper-\nmiddle", "Lower-\nmiddle", "Low\nincome"], fontsize=7.0)
    ax.set_ylabel("Share in 100-yr floodplain (%)")
    ax.set_ylim(0, 28.5)
    ax.legend(frameon=False, loc="upper right", borderaxespad=0.1, handlelength=1.1)
    ptitle(ax, "Exposure by income group", pad=4)


def render_c(ax):
    x = np.r_[0, cc["cum_pop_nat"].values]
    ax.plot([0, 1], [0, 1], color=C["lgray"], lw=0.8, ls=(0, (3, 2)))
    for col, cl, ls, lw, lab in [
        ("cum_n_ev", C["ev"], (0, (1.5, 1.5)), 1.25, f"charging stock  C={ci['ev_stock']:.2f}"),
        ("cum_exp_ev", C["ev"], "-", 1.8, f"charging exposed  C={ci['ev_exposed']:.2f}"),
        ("cum_eff_ev", C["risk"], "-", 1.8, f"charging effective  C={ci['ev_effective']:.2f}"),
        ("cum_n_gas", GAS, (0, (1.5, 1.5)), 1.25, f"fuel stock  C={ci['gas_stock']:.2f}"),
        ("cum_exp_gas", GAS, (0, (4, 2)), 1.4, f"fuel exposed  C={ci['gas_exposed']:.2f}"),
    ]:
        ax.plot(x, np.r_[0, cc[col].values], color=cl, ls=ls, lw=lw, label=lab)
    ax.set_xlabel("Cumulative population share (poorest " + "$\\rightarrow$" + " richest)")
    ax.set_ylabel("Cumulative share of burden")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.legend(
        frameon=False, loc="upper left", fontsize=6.0, handlelength=1.7, borderaxespad=0.02, labelspacing=0.32
    )
    ptitle(ax, "Concentration of burden", pad=4)


def render_d(ax):
    d = pbc.head(10).iloc[::-1]
    y = np.arange(len(d))
    xeff = d["effective"].clip(lower=0.9)
    for yi, (n0, e0) in enumerate(zip(d["nominal"], xeff)):
        ax.plot([e0, n0], [yi, yi], color=C["lgray"], lw=1.3, zorder=2)
    ax.scatter(d["nominal"], y, s=26, facecolor="white", edgecolor=C["ink"], lw=0.9, zorder=3)
    ax.scatter(xeff, y, s=30, facecolor=[INC_C[i] for i in d["inc_grp"]], edgecolor="white", lw=0.6, zorder=4)
    ax.set_xscale("log")
    ax.set_xlim(0.55, 240000.0)
    ax.set_yticks(y)
    ax.set_yticklabels(d["country"].replace({"United States of America": "United States"}), fontsize=6.8)
    ax.set_xticks([1, 10, 100, 1000.0, 10000.0, 100000.0])
    ax.set_xticklabels(["1", "10", "100", "1k", "10k", "100k"])
    ax.set_xlabel("Exposed chargers,  nominal " + "$\\circ\\;\\rightarrow$" + " effective " + "$\\bullet$")
    for yi, e0, xv in zip(y, d["effective"].to_numpy(), xeff.to_numpy()):
        if e0 == 0:
            ax.text(
                float(xv),
                yi - 0.25,
                "0",
                fontsize=6.2,
                ha="center",
                va="top",
                color=C["ink"],
                zorder=10,
                clip_on=False,
            )
    hs = [
        plt.Line2D([], [], marker="o", ls="", mfc=INC_C[k], mec="white", ms=5.5, label=INC_LAB[k])
        for k in ["HIC", "UMC", "LMC"]
    ]
    ax.legend(
        handles=hs,
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.27, 0.99),
        ncol=3,
        fontsize=6.0,
        borderaxespad=0.0,
        handletextpad=0.05,
        columnspacing=0.7,
    )
    ptitle(ax, "Protection adjustment: top 10", pad=15)


def render_e(ax):
    rows = ["HIC", "UMC", "LMC"]
    labs = ["<10", "10-50", "50-100", "100-1000", ">=1000"]
    sh = S["protection"]["share_effective_by_income"]
    eff_lab = {"HIC": f"{sh['HIC']:.0f}%", "UMC": f"{sh['UMC']:.0f}%", "LMC": f"{sh['LMC']:.1f}%"}
    y = np.arange(len(rows))[::-1]
    left = np.zeros(len(rows))
    short = {"<10": "<10", "10-50": "10–50", "50-100": "50–100", "100-1000": "100–1k", ">=1000": "≥1k"}
    for lb in labs:
        v = pbn.loc[rows, f"p_{lb}"].values
        ax.barh(y, v, left=left, height=0.58, color=PROT_C[lb], label=short[lb], ec="white", lw=0.4)
        left += v
    for yi, r in zip(y, rows):
        ax.text(-2.5, yi, INC_LAB[r], fontsize=7.0, va="center", ha="right")
        ax.text(
            98.5,
            yi + 0.4,
            eff_lab[r],
            fontsize=7.0,
            va="center",
            ha="right",
            color=C["risk_dark"],
            fontweight="bold",
            path_effects=WS,
        )
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.55, 2.85)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Share of exposed chargers (%),\nby protection standard (RP, yr)", fontsize=7.4)
    ax.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.42, 0.92),
        ncol=3,
        fontsize=7.2,
        handlelength=1.5,
        handleheight=1.1,
        handletextpad=0.4,
        columnspacing=1.0,
        labelspacing=0.42,
    )
    ptitle(ax, "Protection standards", pad=18)


def render_f(ax):
    for iso, g in qw.groupby("iso3"):
        g = g.sort_values("q")
        ax.plot(g["q"], g["ev_pct"], color="#B9C6CE", lw=0.55, alpha=0.6, zorder=2)
    for kind, cl, lab, mk in [("ev", C["ev"], "charging (pooled)", "o"), ("gas", GAS, "fuel (pooled)", "s")]:
        g = qp[qp["kind"] == kind].sort_values("q")
        ax.plot(
            g["q"], g["pct"], color=cl, lw=2.0, marker=mk, ms=3.6, mec="white", mew=0.5, label=lab, zorder=4
        )
        ax.text(
            5.14,
            g["pct"].iloc[-1],
            f"{g['pct'].iloc[-1]:.0f}%",
            color=cl,
            fontsize=7.0,
            va="center",
            fontweight="bold",
        )
    chn = qw[qw["iso3"] == "CHN"].sort_values("q")
    ax.plot(chn["q"], chn["ev_pct"], color=C["ev_dark"], lw=1.1, zorder=3)
    ax.text(
        4.62, 45.3, f"China {chn['ev_pct'].iloc[-1]:.0f}%", fontsize=6.6, color=C["ev_dark"], path_effects=WS
    )
    ax.text(4.05, 2.0, "India", fontsize=6.4, color=C["ink"], path_effects=WS)
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xlabel("Subnational GDP p.c. quintile within country\n(1 = poorest fifth, 5 = richest fifth)")
    ax.set_ylabel("Share in 100-yr floodplain (%)")
    ax.set_xlim(0.8, 5.55)
    ax.set_ylim(0, 48)
    ax.legend(frameon=False, loc="upper left", fontsize=6.4, borderaxespad=0.05)
    ptitle(ax, "Within-country income gradient", pad=4)


def render_g(ax):
    d = inh.dropna(subset=["ev_per_gas", "gas_pct"]).copy()
    d["x"] = d["ev_per_gas"].clip(lower=0.0016)
    pbins = [0, 10, 50, 100, 1000000000.0]
    pcols = ["#7E2708", "#D55E00", "#F2B387", "#3D74A6"]
    col = [
        "#B9B9B9" if np.isnan(v) else pcols[int(np.digitize(v, pbins)) - 1]
        for v in d["prot_median_gas_exposed"]
    ]
    smax = d["n_gas"].max()
    ax.axvspan(0.0011, 0.1, color=C["risk"], alpha=0.055, zorder=1)
    ax.text(
        0.00135,
        55.3,
        "pre-transition markets",
        fontsize=6.9,
        color=C["risk_dark"],
        va="top",
        zorder=6,
        fontweight="bold",
    )
    ax.scatter(
        d["x"],
        d["gas_pct"],
        s=8 + 170 * np.sqrt(d["n_gas"] / smax),
        c=col,
        alpha=0.82,
        edgecolor="white",
        lw=0.45,
        zorder=5,
    )
    names = {
        "Pakistan": 2.6,
        "Egypt": -3.6,
        "Bangladesh": 2.6,
        "Vietnam": 2.8,
        "Thailand": 2.4,
        "Indonesia": -3.4,
        "India": 2.4,
        "Nigeria": 2.4,
        "China": 2.8,
        "Netherlands": 2.8,
        "Japan": 2.4,
        "United States of America": 2.4,
        "Norway": 2.4,
    }
    for nm, dy in names.items():
        r = d[d["country"] == nm]
        if not len(r):
            continue
        r = r.iloc[0]
        lab = "United States" if nm.startswith("United States") else nm
        ax.text(
            r["x"],
            r["gas_pct"] + dy,
            lab,
            fontsize=6.3,
            ha="center",
            color=C["ink"],
            path_effects=WS,
            zorder=7,
        )
    mya = d[d["country"] == "Myanmar"]
    if len(mya):
        ax.annotate(
            "Myanmar",
            (mya["x"].iloc[0], mya["gas_pct"].iloc[0]),
            xytext=(0.0118, 21.0),
            fontsize=6.3,
            ha="center",
            va="center",
            color=C["ink"],
            path_effects=WS,
            zorder=7,
            arrowprops=dict(arrowstyle="-", lw=0.5, color=C["gray"], shrinkA=1.5, shrinkB=2.5),
        )
    ax.set_xscale("log")
    ax.set_xlim(0.0011, 30)
    ax.set_ylim(-1.5, 57)
    ax.set_xticks([0.01, 0.1, 1, 10])
    ax.set_xticklabels(["0.01", "0.1", "1", "10"])
    ax.set_xlabel("Charging maturity: chargers per fuel station (log)")
    ax.set_ylabel("Fuel network in 100-yr floodplain (%)")
    hs = [
        plt.Line2D([], [], marker="o", ls="", mfc=c_, mec="white", ms=5, label=l_)
        for c_, l_ in zip(pcols, ["<10", "10–50", "50–100", "≥100"])
    ]
    lg = ax.legend(
        handles=hs,
        frameon=False,
        loc="upper right",
        fontsize=6.1,
        title="median protection (RP, yr)",
        title_fontsize=6.3,
        handletextpad=0.25,
        labelspacing=0.3,
        borderaxespad=0.02,
    )
    lg._legend_box.align = "left"
    ax.text(0.44, 52.4, "fuel stations", fontsize=5.9, color=C["ink"], ha="center")
    for v, lab, xv in [(10000, "10k", 0.3), (130000, "130k", 0.66)]:
        ax.scatter([xv], [47.2], s=8 + 170 * np.sqrt(v / smax), facecolor="none", edgecolor=C["gray"], lw=0.6)
        ax.text(xv, 41.4, lab, fontsize=5.6, color=C["ink"], ha="center")
    ptitle(ax, "Fuel exposure vs charging maturity", pad=4)


H1, H2, H3 = (1.874, 1.66, 1.66)
GAP1, GAP2 = (0.85, 0.8)
MT, MB = (0.12, 0.52)
HF = MT + H1 + GAP1 + H2 + GAP2 + H3 + MB
FW = 7.2


def ax_at(fig, x, y, w, h):
    return fig.add_axes([x / FW, y / HF, w / FW, h / HF])


def composite():
    fig = plt.figure(figsize=(FW, HF))
    y1 = HF - MT - H1
    a = ax_at(fig, 0.02, y1, H1 * 2.329, H1)
    render_a(a)
    a.text(0.012, 0.97, "a", transform=a.transAxes, fontsize=13, fontweight="bold", va="top")
    b = ax_at(fig, 5.06, y1, 2.02, H1 - 0.2)
    render_b(b)
    panel(b, "b", dx=-0.23, dy=1.115, fs=13)
    y2 = y1 - GAP1 - H2
    c = ax_at(fig, 0.5, y2, 1.72, H2)
    render_c(c)
    panel(c, "c", dx=-0.27, dy=1.22, fs=13)
    d = ax_at(fig, 3.32, y2, 1.5, H2)
    render_d(d)
    panel(d, "d", dx=-0.6, dy=1.22, fs=13)
    e = ax_at(fig, 5.68, y2, 1.4, H2)
    render_e(e)
    panel(e, "e", dx=-0.51, dy=1.22, fs=13)
    y3 = y2 - GAP2 - H3
    f = ax_at(fig, 0.5, y3, 2.1, H3)
    render_f(f)
    panel(f, "f", dx=-0.22, dy=1.13, fs=13)
    g = ax_at(fig, 3.45, y3, 3.63, H3)
    render_g(g)
    panel(g, "g", dx=-0.15, dy=1.13, fs=13)
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Figure5"))


if __name__ == "__main__":
    main()
