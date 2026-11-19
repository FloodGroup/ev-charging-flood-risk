import json
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import plot_style as style
from plot_style import C

DATA = os.path.join(HERE, "data")
style.set_style()
plt.rcParams.update(
    {
        "font.size": 8,
        "axes.labelsize": 8.2,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.2,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
    }
)
S = json.load(open(f"{DATA}/fig02_stats.json"))
S2 = json.load(open(f"{DATA}/exposure_decomposition.json"))
rp = pd.read_csv(f"{DATA}/rp_curve.csv")
dc = pd.read_csv(f"{DATA}/depth_curve.csv")
mr = pd.read_csv(f"{DATA}/metro_ratios.csv")
mrv = mr.dropna(subset=["ratio"])
mrv = mrv[(mrv.exp_ev > 0) & (mrv.exp_gas > 0)]
cm_ev = pd.read_csv(f"{DATA}/city_map_exposure.csv")
cm_gas = pd.read_csv(f"{DATA}/city_map_exposure_gas.csv")
cs = pd.read_csv(f"{DATA}/country_scatter.csv")
stock = pd.read_csv(f"{DATA}/exposed_stock.csv")
XMAX = 17005000.0
XL, XR = (-XMAX * (170 / 180), XMAX)
BINS = [0, 1e-09, 5, 15, 30, 100.001]
BCOL = ["#9EC4BC", "#F2D488", "#E89143", "#C4441C", "#7C120A"]
BLAB = ["none", "0–5", "5–15", "15–30", ">30"]


def map_legend(ax, word):
    lx = ax.inset_axes([0.335, 0.004, 0.59, 0.175])
    lx.set_xlim(0, 1)
    lx.set_ylim(0, 1)
    lx.axis("off")
    lx.add_patch(Rectangle((0, 0), 1, 1, fc="none", ec="#D8D8D4", lw=0.5, transform=lx.transAxes))
    lx.text(
        0.025,
        0.85,
        f"dot colour — % of city's {word} in 100-yr floodplain",
        fontsize=6.1,
        va="center",
        color="#2E2E2E",
    )
    x0, bw = (0.025, 0.135)
    for i, (cc, lb) in enumerate(zip(BCOL, BLAB)):
        lx.add_patch(Rectangle((x0 + i * bw, 0.36), bw - 0.006, 0.28, fc=cc, ec="none"))
        lx.text(
            x0 + i * bw + (bw - 0.006) / 2,
            0.5,
            lb,
            fontsize=5.8,
            ha="center",
            va="center",
            color="white" if i >= 3 else "#2E2E2E",
        )
    lx.text(0.025, 0.13, "dot area — stations per city", fontsize=6.1, va="center", color="#2E2E2E")
    for nn, xx in [(50, 0.58), (500, 0.7), (5000, 0.84)]:
        lx.scatter(
            [xx],
            [0.15],
            s=np.clip(0.42 * np.sqrt(nn), 2, 55),
            c="#9A9A9A",
            lw=0.3,
            edgecolor="white",
            transform=lx.transAxes,
            clip_on=False,
        )
        lx.text(xx + 0.022, 0.15, f"{nn:,}", fontsize=5.8, va="center", color="#2E2E2E")


def city_map(ax, df, pct_col, n_col, headline, subline, word):
    style.draw_world(ax)
    d = df.dropna(subset=["GCPNT_LAT"]).copy()
    x, y = style.to_rob(d.GCPNT_LON.values, d.GCPNT_LAT.values)
    v = d[pct_col].values
    n = d[n_col].values
    s = np.clip(0.42 * np.sqrt(n), 2.0, 55.0)
    ci = np.clip(np.digitize(v, BINS) - 1, 0, 4)
    order = np.argsort(-s)
    ax.scatter(
        x[order],
        y[order],
        s=s[order],
        c=np.array(BCOL)[ci[order]],
        alpha=0.88,
        lw=0.18,
        edgecolor="white",
        zorder=4,
        rasterized=True,
    )
    ax.set_xlim(XL, XR)
    ax.set_ylim(-6350000.0, 7850000.0)
    ax.text(
        0.992,
        0.965,
        headline,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
        color="#2E2E2E",
        zorder=30,
        bbox=dict(fc="white", ec="none", alpha=0.92, pad=1.8),
    )
    ax.text(
        0.992,
        0.865,
        subline,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.6,
        color="#6E6E6E",
        zorder=30,
        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2),
    )
    map_legend(ax, word)


def note_city(ax, lon, lat, txt, dxy, ha="left"):
    x, y = style.to_rob(np.array([lon]), np.array([lat]))
    ax.annotate(
        txt,
        (x[0], y[0]),
        xytext=(x[0] + dxy[0], y[0] + dxy[1]),
        fontsize=6.0,
        color="#1A1A1A",
        ha=ha,
        va="center",
        zorder=31,
        arrowprops=dict(arrowstyle="-", lw=0.55, color="#6E6E6E", shrinkA=0, shrinkB=1.5),
        bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.9),
    )


def _city_pct(name, country):
    r = cm_ev[(cm_ev.UC_NM_MN == name) & (cm_ev.CTR_MN_NM == country)]
    return float(r.ev_pct.iloc[0])


def panel_a(ax, letter=None):
    city_map(
        ax,
        cm_ev,
        "ev_pct",
        "n_ev",
        f"charging: {S['rp100_ev']:.1f}% exposed (n = {S['n_ev']:,})",
        f"{S2['n_cities_a']:,} urban centres shown",
        "chargers",
    )
    note_city(ax, 121.4, 31.2, f"Shanghai {_city_pct('Shanghai', 'China'):.0f}%", (2200000.0, 1100000.0))
    note_city(ax, 100.5, 13.8, f"Bangkok {_city_pct('Bangkok', 'Thailand'):.0f}%", (600000.0, -2300000.0))
    note_city(
        ax,
        4.4,
        52.0,
        f"Rotterdam–The Hague {_city_pct('Rotterdam [The Hague]', 'Netherlands'):.0f}%",
        (-3100000.0, 1900000.0),
        ha="center",
    )
    if letter:
        ax.text(
            0.011, 0.972, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", zorder=10
        )


def panel_c(ax, letter=None):
    city_map(
        ax,
        cm_gas,
        "gas_pct",
        "n_gas",
        f"fuel: {S['rp100_gas']:.1f}% exposed (n = {S['n_gas']:,})",
        f"{len(cm_gas):,} urban centres shown",
        "fuel stations",
    )
    if letter:
        ax.text(
            0.011, 0.972, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", zorder=10
        )


def panel_b(ax):
    x = rp.rp.values
    ax.fill_between(x, rp.gas_pct, rp.ev_pct, color=C["risk"], alpha=0.12, lw=0)
    ax.plot(x, rp.ev_pct, "-o", color=C["ev_dark"], lw=1.6, ms=3.6, label="Charging")
    ax.plot(x, rp.gas_pct, "-s", color=C["gas"], lw=1.4, ms=3.2, label="Fuel")
    for i, ha, xoff in [(0, "left", 1.07), (4, "center", 1.0), (6, "right", 0.93)]:
        ax.text(
            x[i] * xoff,
            (rp.ev_pct.iloc[i] + rp.gas_pct.iloc[i]) / 2,
            f"×{rp.ratio.iloc[i]:.2f}",
            fontsize=6.6,
            ha=ha,
            va="center",
            color=C["risk_dark"],
            fontweight="bold",
        )
    ax.set_xscale("log")
    ax.set_xticks([10, 20, 50, 100, 200, 500])
    ax.set_xticklabels(["10", "20", "50", "100", "200", "500"])
    ax.minorticks_off()
    ax.set_xlabel("Flood return period (years)")
    ax.set_ylabel("Stations in floodplain (%)")
    ax.set_ylim(7.5, 21.0)
    ax.set_yticks([8, 10, 12, 14, 16, 18, 20])
    ax.legend(loc="lower right", frameon=False, fontsize=6.8)
    ins = ax.inset_axes([0.1, 0.655, 0.42, 0.275], facecolor="none")
    ins.patch.set_facecolor((1, 1, 1, 0))
    ins.patch.set_alpha(0)
    ins.plot(dc.thresh, dc.ratio, "-o", color=C["risk_dark"], lw=1.2, ms=2.6)
    ins.set_xticks([0, 0.5, 1, 1.5])
    ins.set_xticklabels(["0", "0.5", "1", "1.5"], fontsize=5.9)
    ins.set_yticks([1.2, 1.3])
    ins.set_yticklabels(["1.2", "1.3"], fontsize=5.9)
    ins.set_ylim(1.14, 1.33)
    ttl = ins.set_title("Charging:fuel ratio vs depth ≥ x m", fontsize=5.8, pad=2)
    ttl.set_x(0.58)
    for sp in ins.spines.values():
        sp.set_linewidth(0.7)
    ins.tick_params(width=0.7, length=2)
    ins.tick_params(axis="y", pad=1.3)


def panel_d(ax):
    r = mrv.ratio.values
    bins = np.logspace(np.log10(0.15), np.log10(8), 25)
    ax.hist(
        np.clip(r, bins[0] * 1.0001, bins[-1] * 0.9999),
        bins=bins,
        color=C["ev_light"],
        edgecolor="white",
        lw=0.5,
    )
    ax.set_xscale("log")
    med = float(np.median(r))
    ax.axvline(1, color=C["gas"], lw=1.1, ls="--")
    ax.axvline(med, color=C["ink"], lw=1.3)
    ax.text(
        med * 0.92,
        0.985,
        f"median ×{med:.2f}",
        fontsize=6.8,
        ha="right",
        va="top",
        transform=ax.get_xaxis_transform(),
        color=C["ink"],
    )
    ax.text(1.06, 0.985, "parity", fontsize=6.6, va="top", transform=ax.get_xaxis_transform(), color=C["gas"])
    ax.text(
        0.97,
        0.82,
        f"{S['n_metros']} urban centres\n{S['metro_share_gt1']:.0f}% above parity",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
    )
    ax.set_xticks([0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.minorticks_off()
    ax.set_xlabel("Within-city exposure ratio\n(charging % ÷ fuel %, RP100)", labelpad=1)
    ax.xaxis.set_label_coords(0.5, -0.11)
    ax.set_ylabel("Urban centres")


CONT_C = style.CONT_C


def panel_e(ax):
    d = cs.dropna(subset=["continent"]).copy()
    lo, hi = (0.07, 60)
    ax.plot([lo, hi], [lo, hi], color="#B0B0B0", lw=0.9, ls="--", zorder=1)
    ax.text(
        0.04,
        0.9,
        "chargers more exposed",
        transform=ax.transAxes,
        fontsize=6.2,
        color="#8A8A8A",
        style="italic",
    )
    d["s"] = np.clip(0.5 * np.sqrt(d.n_ev), 4, 300)
    d = d.sort_values("s", ascending=False)
    ax.scatter(
        d.gas_pct, d.ev_pct, s=d.s, c=d.continent.map(CONT_C), alpha=0.78, lw=0.4, edgecolor="white", zorder=3
    )
    ax.scatter([S["rp100_gas"]], [S["rp100_ev"]], marker="*", s=95, c="#1A1A1A", zorder=6)
    ax.text(
        S["rp100_gas"] * 0.8,
        S["rp100_ev"],
        f"World ×{S['rp100_ratio']:.2f}",
        ha="right",
        fontsize=6.6,
        fontweight="bold",
        va="center",
        zorder=6,
    )
    short = {"United States of America": "United States"}
    near = {
        "China": (1.5, 0.8, "left"),
        "Netherlands": (0.7, 1.05, "right"),
        "United Kingdom": (0.7, 1.05, "right"),
        "Denmark": (1.55, 1.15, "left"),
        "Brazil": (1.15, 0.52, "left"),
        "Germany": (0.55, 0.88, "right"),
    }
    far = {
        "Japan": (30, 12.5, "left"),
        "India": (18.5, 7.4, "left"),
        "United States of America": (11.5, 5.2, "left"),
        "Vietnam": (28, 54, "center"),
        "Thailand": (11.5, 43, "right"),
    }
    arr = dict(arrowstyle="-", lw=0.5, color="#8A8A8A", shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0")
    for name, (fx, fy, ha) in near.items():
        r = d[d.country == name]
        if len(r):
            gx, ey = (float(r.gas_pct.iloc[0]), float(r.ev_pct.iloc[0]))
            tx, ty = (gx * fx, ey * fy)
            ax.plot([gx, tx], [ey, ty], color="#8A8A8A", lw=0.5, zorder=6, clip_on=False)
            ax.text(tx, ty, short.get(name, name), fontsize=6.2, va="center", ha=ha, zorder=7)
    for name, (tx, ty, ha) in far.items():
        r = d[d.country == name]
        if len(r):
            gx, ey = (float(r.gas_pct.iloc[0]), float(r.ev_pct.iloc[0]))
            ax.plot([gx, tx], [ey, ty], color="#8A8A8A", lw=0.5, zorder=6, clip_on=False)
            ax.text(tx, ty, short.get(name, name), fontsize=6.2, va="center", ha=ha, zorder=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    tks = [0.1, 0.3, 1, 3, 10, 30]
    ax.set_xticks(tks)
    ax.set_xticklabels([str(t) for t in tks])
    ax.set_yticks(tks)
    ax.set_yticklabels([str(t) for t in tks])
    ax.minorticks_off()
    ax.set_xlabel("Fuel stations in floodplain (%)")
    ax.set_ylabel("Chargers in floodplain (%)")
    present = set(d.continent)
    hnd = [
        Line2D([], [], marker="o", ls="none", ms=4.2, mfc=cc, mec="none", label=k)
        for k, cc in CONT_C.items()
        if k in present
    ]
    ax.legend(
        handles=hnd,
        loc="lower right",
        frameon=False,
        fontsize=6.0,
        borderaxespad=0.2,
        handletextpad=0.15,
        labelspacing=0.32,
        bbox_to_anchor=(1.02, 0.0),
    )


def panel_f(ax):
    g, comp, sit, e = (S2["kita_gas"], S2["kita_comp"], S2["kita_sit"], S2["kita_ev"])
    peak = g + comp
    y0 = 11.55
    bw = 0.62
    ax.bar(0, g - y0, bw, bottom=y0, color=C["gas"], alpha=0.92)
    ax.bar(1, comp, bw, bottom=g, color=C["risk"], alpha=0.92)
    ax.bar(2, sit, bw, bottom=peak, color=C["ev"], alpha=0.92)
    ax.bar(3, e - y0, bw, bottom=y0, color=C["ev_dark"], alpha=0.92)
    for xa, xb, lv in [(0, 1, g), (1, 2, peak), (2, 3, e)]:
        ax.plot([xa + bw / 2, xb - bw / 2], [lv, lv], color="#6E6E6E", lw=0.7, ls=":")
    ax.text(0, g + 0.12, f"{g:.1f}%", ha="center", fontsize=7.0, fontweight="bold", color=C["gas"])
    ax.text(
        1, peak + 0.12, f"+{comp:.1f} pp", ha="center", fontsize=7.0, fontweight="bold", color=C["risk_dark"]
    )
    ax.text(
        2,
        peak + 0.12,
        f"−{abs(sit):.1f} pp",
        ha="center",
        fontsize=7.0,
        fontweight="bold",
        color=C["ev_dark"],
    )
    ax.text(3, e + 0.12, f"{e:.1f}%", ha="center", fontsize=7.0, fontweight="bold", color=C["ev_dark"])
    ax.set_xticks(range(4))
    ax.set_xticklabels(["Fuel", "Country\nmix", "Within-\ncountry\nsiting", "Chargers"], fontsize=5.7)
    ax.plot(
        [-0.022, 0.022],
        [0.028, 0.048],
        transform=ax.transAxes,
        color=C["spine"],
        lw=0.9,
        clip_on=False,
        zorder=20,
    )
    ax.plot(
        [-0.022, 0.022],
        [0.05, 0.07],
        transform=ax.transAxes,
        color=C["spine"],
        lw=0.9,
        clip_on=False,
        zorder=20,
    )
    ax.set_ylim(y0, 16.8)
    ax.set_yticks([12, 13, 14, 15, 16])
    ax.set_ylabel("Stations in floodplain (%)")
    ax.text(
        0.02,
        0.99,
        "Kitagawa decomposition, country strata:\ncountry mix explains >100% of the gap;\nwithin-country siting is slightly protective",
        transform=ax.transAxes,
        fontsize=6.0,
        va="top",
        color="#3C3C3C",
    )


def panel_g(ax):
    st = stock[stock.country.notna() & (stock.country != "")]
    top = st.head(8).copy()
    total = int(S["ev_exp100"])
    rest = total - int(top.exp_ev.sum())
    short = {"United States of America": "United States"}
    names = [short.get(n, n) for n in top.country] + ["Rest of world"]
    vals = list(top.exp_ev / 1000) + [rest / 1000]
    shares = list(top.share) + [100 * rest / total]
    cmap = dict(zip(cs.country, cs.continent))
    cols = [CONT_C.get(cmap.get(c), "#B0B0B0") for c in top.country] + ["#B0B0B0"]
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, height=0.68, color=cols, alpha=0.92)
    for i, (yi, v, sh) in enumerate(zip(y, vals, shares)):
        if v > 50:
            ax.text(
                v - 2,
                yi,
                f"{v:.1f}k · {sh:.0f}%",
                ha="right",
                va="center",
                fontsize=6.0,
                color="white",
                fontweight="bold",
            )
        elif i == 1:
            ax.text(v + 2, yi, f"{v:.1f}k · {sh:.0f}%", ha="left", va="center", fontsize=6.0, color="#3C3C3C")
        else:
            ax.text(v + 2, yi, f"{v:.1f}k", ha="left", va="center", fontsize=6.0, color="#3C3C3C")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.3)
    ax.set_xlim(0, 80)
    ax.set_xticks([0, 25, 50, 75])
    ax.set_xlabel("Chargers in 100-yr floodplain\n(thousands, RP100)")
    top2_share = 100.0 * float(st.head(2).exp_ev.sum()) / total
    ax.text(
        0.97,
        0.4,
        f"China + the Netherlands hold\n{top2_share:.0f}% of the world's\nflood-exposed chargers",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=6.2,
        color="#3C3C3C",
    )
    ax.tick_params(axis="y", length=0)


def composite():
    ASP = (XR - XL) / (7850000.0 + 6350000.0)
    H = 1.874
    WM = H * ASP
    L, GUT, R = (0.09, 0.8, 0.07)
    WC = 7.2 - L - WM - GUT - R
    H3 = 2.02
    BOT, GAPR, GAP3, TOP = (0.5, 0.36, 0.44, 0.06)
    HF = BOT + H3 + GAP3 + 2 * H + GAPR + TOP
    fig = plt.figure(figsize=(7.2, HF))

    def ax_at(x, y, w, h):
        return fig.add_axes([x / 7.2, y / HF, w / 7.2, h / HF])

    y3 = BOT
    y2 = BOT + H3 + GAP3
    y1 = y2 + H + GAPR
    ax_a = ax_at(L, y1, WM, H)
    ax_b = ax_at(L + WM + GUT, y1, WC, H)
    ax_c = ax_at(L, y2, WM, H)
    ax_d = ax_at(L + WM + GUT, y2, WC, H)
    ax_e = ax_at(0.42, y3, 2.02, H3)
    ax_f = ax_at(3.0, y3, 1.72, H3)
    ax_g = ax_at(5.4, y3, 1.73, H3)
    panel_a(ax_a, "a")
    panel_b(ax_b)
    panel_c(ax_c, "c")
    panel_d(ax_d)
    panel_e(ax_e)
    panel_f(ax_f)
    panel_g(ax_g)
    ax_d.set_zorder(10)
    ax_g.set_zorder(1)
    ax_e.set_zorder(2)
    ax_f.set_zorder(2)
    ax_d.patch.set_facecolor((1, 1, 1, 0))
    ax_d.set_clip_on(False)
    ax_d.xaxis.get_label().set_zorder(20)
    ax_d.xaxis.get_label().set_clip_on(False)
    xg = (L + WM + 0.055) / 7.2
    for lt, ax_ in [("b", ax_a), ("d", ax_c)]:
        p_ = ax_.get_position()
        fig.text(xg, p_.y1 - 0.028 * p_.height, lt, fontsize=13, fontweight="bold", va="top", ha="left")
    for lt, ax_ in [("e", ax_e), ("f", ax_f), ("g", ax_g)]:
        p_ = ax_.get_position()
        fig.text(p_.x0 - 0.05, p_.y1 + 0.012, lt, fontsize=13, fontweight="bold", va="top", ha="left")
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Figure2"))


if __name__ == "__main__":
    main()
