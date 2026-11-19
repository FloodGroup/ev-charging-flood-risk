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
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
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
        "font.size": 8.6,
        "axes.labelsize": 8.8,
        "xtick.labelsize": 7.9,
        "ytick.labelsize": 7.9,
        "legend.fontsize": 7.8,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
    }
)
S4 = json.load(open(f"{DATA}/cohort_stats.json"))
cc = pd.read_csv(f"{DATA}/cohort_curves.csv")
ss = pd.read_csv(f"{DATA}/country_slopes.csv").sort_values("slope_pp_yr", ascending=False)
mm = pd.read_csv(f"{DATA}/metro_slopes.csv")
SHORT = {"United States of America": "United States", "United Kingdom": "UK"}
CONT = {
    "United States of America": "North America",
    "Canada": "North America",
    "United Kingdom": "Europe",
    "Germany": "Europe",
    "Spain": "Europe",
    "France": "Europe",
    "Italy": "Europe",
    "Sweden": "Europe",
    "Norway": "Europe",
    "Portugal": "Europe",
    "Ireland": "Europe",
    "Austria": "Europe",
    "Belgium": "Europe",
    "Turkey": "Asia",
    "Australia": "Oceania",
}
CONT_C = style.CONT_C
EV = C["ev"]
EVD = C["ev_dark"]
GAS = C["gas"]
GRID_C = "#6A51A3"
halo = [pe.withStroke(linewidth=1.8, foreground="white", alpha=0.95)]


def _metro_slope(metro_id):
    v = float(mm.loc[mm.metro == metro_id, "slope_pp_yr"].iloc[0])
    return f"{v:+.2f}".replace("-", "−")


def mini(ax, cname, first_col, last_row):
    g = cc[(cc.country == cname) & (cc.n >= 30)].sort_values("mid")
    r = ss[ss.country == cname].iloc[0]
    ax.fill_between(g["mid"], g.lo, g.hi, color=EV, alpha=0.18, lw=0)
    ax.plot(g["mid"], g.exp_pct, "-o", color=EVD, lw=1.3, ms=2.8)
    ax.axhline(r.stock_gas_pct, color=GAS, lw=1.0, ls=(0, (3, 2)))
    scol = C["risk_dark"] if r.lo > 0 else EVD if r.hi < 0 else "#3C3C3C"
    ax.text(
        0.05,
        0.985,
        SHORT.get(cname, cname),
        transform=ax.transAxes,
        fontsize=8.0,
        fontweight="bold",
        va="top",
    )
    ax.text(
        0.97,
        0.985,
        f"{r.slope_pp_yr:+.2f}".replace("-", "−"),
        transform=ax.transAxes,
        fontsize=7.4,
        fontweight="bold",
        va="top",
        ha="right",
        color=scol,
    )
    ax.set_xlim(2011.5, 2026)
    ax.set_ylim(0, 24.5)
    ax.set_xticks([2015, 2020, 2025])
    if last_row:
        ax.set_xticklabels(["2015", "2020", "2025"], fontsize=7.0)
    else:
        ax.set_xticklabels([])
    ax.set_yticks([0, 10, 20])
    if first_col:
        ax.set_yticklabels(["0", "10", "20"], fontsize=7.0)
    else:
        ax.set_yticklabels([])
    ax.tick_params(length=2.2, width=0.8)
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)


def panel_b(ax):
    d = ss.iloc[::-1]
    y = np.arange(len(d))
    ax.axvline(0, color="#B0B0B0", lw=0.9)
    ax.hlines(y, d.lo, d.hi, color="#9FB8B3", lw=1.4, zorder=2)
    ax.scatter(d.slope_pp_yr, y, s=26, c=EVD, zorder=4, label="national")
    ok = d.fe_slope_pp_yr.notna()
    ax.scatter(
        d.fe_slope_pp_yr[ok],
        y[ok],
        s=30,
        facecolors="none",
        edgecolors=GRID_C,
        lw=1.3,
        zorder=5,
        label="within-city",
    )
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT.get(c, c) for c in d.country], fontsize=7.0)
    ax.set_xlabel("Cohort trend in floodplain share\n(pp per year, 2010–2025)")
    ax.set_xlim(-1.15, 1.25)
    ax.tick_params(axis="y", length=0)
    ax.legend(
        loc="lower right",
        frameon=False,
        fontsize=7.0,
        handletextpad=0.25,
        borderaxespad=0.02,
        labelspacing=0.3,
        bbox_to_anchor=(1.0, 0.0),
    )
    hdr = f"median {S4['slope_median']:+.2f} pp/yr — {S4['n_countries'] - S4['pos_sig'] - S4['neg_sig']} of {S4['n_countries']} flat".replace(
        "-", "−"
    )
    ax.text(0.02, 1.015, hdr, transform=ax.transAxes, fontsize=7.4, color="#1A1A1A", va="bottom")


def panel_c(ax):
    v = mm.slope_pp_yr.values
    bins = np.arange(-2.4, 1.5, 0.2)
    ax.hist(
        np.clip(v, bins[0] + 0.01, bins[-1] - 0.01), bins=bins, color=C["ev_light"], edgecolor="white", lw=0.6
    )
    ax.axvline(0, color="#B0B0B0", lw=0.9)
    med = float(np.median(v))
    ax.axvline(med, color=C["ink"], lw=1.3)
    ax.text(
        med + 0.1,
        0.9,
        f"median\n{med:+.2f}".replace("-", "−"),
        fontsize=7.4,
        va="top",
        ha="left",
        transform=ax.get_xaxis_transform(),
        color=C["ink"],
        linespacing=1.05,
    )
    for nm, mid, ha in [("Paris", 2046, "center"), ("Málaga", 1677, "left")]:
        tv = float(mm.loc[mm.metro == mid, "slope_pp_yr"].iloc[0])
        xv = float(np.clip(tv, -2.3, 1.3))
        if xv != tv:
            nm = f"{nm} ({tv:.2f})".replace("-", "−")
        ax.annotate(
            nm,
            (xv, 1.5),
            xytext=(xv - (0.12 if ha == "left" else 0), 6.2),
            fontsize=6.6,
            ha=ha,
            color="#1A1A1A",
            arrowprops=dict(
                arrowstyle="-", lw=0.6, color="#3C3C3C", relpos=(0.08 if ha == "left" else 0.5, 0)
            ),
        )
    ax.text(
        0.03,
        0.97,
        f"{S4['n_metros']} centres\n{S4['metro_pos_sig']}↑ · {S4['metro_neg_sig']}↓ sig.",
        transform=ax.transAxes,
        fontsize=6.8,
        va="top",
        color="#1A1A1A",
    )
    ax.set_xlabel("City cohort trend (pp per year)")
    ax.set_ylabel("Urban centres")
    ax.set_xlim(-2.45, 1.4)


def panel_d(ax):
    d = ss.dropna(subset=["newest_pct", "stock_gas_pct"])
    lo, hi = (0, 17)
    ax.plot([lo, hi], [lo, hi], color="#B0B0B0", lw=0.9, ls="--", zorder=1)
    for _, r in d.iterrows():
        cc_ = CONT_C.get(CONT.get(r.country, ""), "#8A8A8A")
        ax.scatter(
            r.stock_gas_pct,
            r.newest_pct,
            s=np.clip(0.5 * np.sqrt(r.n), 8, 120),
            c=cc_,
            alpha=0.85,
            lw=0.4,
            edgecolor="white",
            zorder=3,
        )
    lab = {
        "United Kingdom": (1.6, 7.6, "center"),
        "United States of America": (9.0, 2.9, "left"),
        "Canada": (4.6, 9.3, "center"),
        "Portugal": (0.5, 2.3, "left"),
        "Norway": (11.3, 13.4, "left"),
        "Turkey": (5.2, 1.7, "left"),
        "Australia": (13.4, 4.2, "center"),
        "Austria": (13.4, 8.6, "center"),
    }
    for nm, (tx, ty, ha) in lab.items():
        r = d[d.country == nm]
        if len(r):
            ax.annotate(
                SHORT.get(nm, nm),
                (r.stock_gas_pct.iloc[0], r.newest_pct.iloc[0]),
                xytext=(tx, ty),
                fontsize=7.0,
                ha=ha,
                va="center",
                zorder=6,
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#8A8A8A", shrinkA=1.5, shrinkB=2.5),
            )
    ax.text(
        0.04,
        0.93,
        "newest chargers\nmore exposed",
        transform=ax.transAxes,
        fontsize=6.8,
        color="#555555",
        style="italic",
        va="top",
    )
    ax.text(
        0.96,
        0.05,
        "less exposed",
        transform=ax.transAxes,
        fontsize=6.2,
        color="#555555",
        style="italic",
        ha="right",
    )
    ax.set_xlabel("Fuel stock in floodplain (%)")
    ax.set_ylabel("2024–25 charger cohort (%)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.text(
        0.04,
        0.74,
        f"{S4['newest_below_gas']}/{S4['n_countries']} below fuel\nbenchmark",
        transform=ax.transAxes,
        fontsize=7.2,
        color="#1A1A1A",
        va="top",
    )


GRP_C = ["#4FA79B", "#BCC8C5", "#D55E00"]
RELIEF_CM = LinearSegmentedColormap.from_list(
    "relief", [(0.0, "#9C988E"), (0.62, "#F4F2EC"), (1.0, "#FFFFFF")]
)
ROAD_LW = {0: 1.0, 1: 0.85, 2: 0.7, 3: 0.5, 4: 0.4}


def draw_basemap(ax, mid, scale=1.0):
    rel = np.load(f"{DATA}/private/relief_{mid}.npz")
    ax.imshow(
        rel["shade"],
        extent=tuple(rel["extent"]),
        origin="upper",
        cmap=RELIEF_CM,
        vmin=107,
        vmax=255,
        interpolation="bilinear",
        zorder=0.5,
        aspect="auto",
    )
    bm = np.load(f"{DATA}/private/basemap_{mid}.npz")
    sea = bm["sea"].astype(float)
    sea[sea == 0] = np.nan
    ax.imshow(
        sea,
        extent=tuple(bm["sea_extent"]),
        origin="upper",
        cmap=ListedColormap(["#DFE9F0"]),
        interpolation="nearest",
        zorder=0.6,
        aspect="auto",
    )
    wo, wxy = (bm["water_off"], bm["water_xy"])
    if len(wo) > 1:
        codes = np.full(len(wxy), MplPath.LINETO, dtype=np.uint8)
        codes[wo[:-1]] = MplPath.MOVETO
        codes[wo[1:] - 1] = MplPath.CLOSEPOLY
        ax.add_patch(
            PathPatch(MplPath(wxy, codes), facecolor="#A9CBE6", edgecolor="none", zorder=0.7, rasterized=True)
        )
    ro, rxy, rc = (bm["road_off"], bm["road_xy"], bm["road_class"])
    segs = [rxy[ro[i] : ro[i + 1]] for i in range(len(ro) - 1)]
    if segs:
        lw = np.array([ROAD_LW[int(c)] for c in rc]) * scale
        ax.add_collection(
            LineCollection(
                segs, colors="#CFCBC0", linewidths=lw + 0.45, zorder=0.8, capstyle="round", rasterized=True
            )
        )
        ax.add_collection(
            LineCollection(
                segs, colors="#FFFFFF", linewidths=lw, zorder=0.85, capstyle="round", rasterized=True
            )
        )


def city_map(ax, mid, name, note, sbkm, letter):
    z = np.load(f"{DATA}/private/citymap_{mid}.npz")
    w, s, e, n = z["bbox"]
    ax.set_facecolor("#F4F2EC")
    draw_basemap(ax, mid, 0.8)
    fl = z["flood"].astype(float)
    fl[fl == 0] = np.nan
    ax.imshow(
        fl,
        extent=(w, e, s, n),
        origin="upper",
        cmap=ListedColormap(["#7FB2DD"]),
        alpha=0.55,
        interpolation="nearest",
        zorder=1,
        aspect="auto",
    )
    m = (z["lon"] >= w) & (z["lon"] <= e) & (z["lat"] >= s) & (z["lat"] <= n)
    exp = z["exp"].astype(bool)
    for grp, col, zo in [(2, GRP_C[2], 5), (0, GRP_C[0], 6)]:
        k = m & (z["grp"] == grp)
        ax.scatter(
            z["lon"][k & ~exp],
            z["lat"][k & ~exp],
            s=3.4,
            c=col,
            lw=0,
            alpha=0.42,
            zorder=zo - 3,
            rasterized=True,
        )
        ax.scatter(
            z["lon"][k & exp],
            z["lat"][k & exp],
            s=11,
            c=col,
            lw=0.35,
            edgecolor="white",
            zorder=zo,
            rasterized=True,
        )
    ax.set_xlim(w, e)
    ax.set_ylim(s, n)
    ax.set_aspect(1 / np.cos(np.deg2rad((s + n) / 2)))
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if letter:
        ax.text(
            0.035,
            0.975,
            letter,
            transform=ax.transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
            zorder=30,
            path_effects=halo,
        )
    ax.text(
        0.97,
        0.975,
        name,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.6,
        fontweight="bold",
        color="#1A1A1A",
        zorder=30,
        path_effects=halo,
    )
    ax.text(
        0.97,
        0.895,
        note,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.0,
        color="#1A1A1A",
        zorder=30,
        path_effects=halo,
    )
    lat0 = (s + n) / 2
    dlon = sbkm / (111.32 * np.cos(np.deg2rad(lat0)))
    x0 = e - dlon - (e - w) * 0.05
    y0 = s + (n - s) * 0.045
    ax.plot(
        [x0, x0 + dlon],
        [y0, y0],
        color="#3C3C3C",
        lw=1.4,
        zorder=30,
        path_effects=[pe.withStroke(linewidth=3.0, foreground="white", alpha=0.9)],
    )
    ax.text(
        x0 + dlon / 2,
        y0 + (n - s) * 0.018,
        f"{sbkm} km",
        ha="center",
        fontsize=6.6,
        color="#1A1A1A",
        zorder=30,
        path_effects=[pe.withStroke(linewidth=1.6, foreground="white", alpha=0.95)],
    )


def panel_leg(ax):
    ax.axis("off")
    items = [
        ("dot", GRP_C[0], 26, 1.0, "in floodplain, opened ≤2018"),
        ("dot", GRP_C[2], 26, 1.0, "in floodplain, opened 2022–25"),
        ("pair", None, 9, 0.55, "outside floodplain (same colours)"),
        ("patch", "#B4CFE4", 0, 1.0, "100-yr floodplain"),
    ]
    for i, (kind, col, ssz, al, lab) in enumerate(items):
        yy = 0.95 - i * 0.105
        if kind == "patch":
            ax.add_patch(Rectangle((0.02, yy - 0.035), 0.1, 0.07, fc=col, ec="none", transform=ax.transAxes))
        elif kind == "pair":
            ax.scatter(
                [0.045, 0.095],
                [yy, yy],
                s=ssz,
                c=[GRP_C[0], GRP_C[2]],
                lw=0,
                alpha=al,
                transform=ax.transAxes,
                clip_on=False,
            )
        else:
            ax.scatter(
                [0.07], [yy], s=ssz, c=col, lw=0.3, edgecolor="white", transform=ax.transAxes, clip_on=False
            )
        ax.text(0.16, yy, lab, fontsize=7.6, va="center", transform=ax.transAxes)
    zp = np.load(f"{DATA}/private/citymap_2046.npz")["stats"]
    zv = np.load(f"{DATA}/private/citymap_1677.npz")["stats"]
    ax.text(
        0.02,
        0.535,
        "all dated stations of the\nearliest and latest cohorts",
        fontsize=6.6,
        color="#1A1A1A",
        transform=ax.transAxes,
        va="top",
    )
    ax.text(0.02, 0.395, "cohort floodplain share", fontsize=7.6, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.02,
        0.372,
        "≤2018 → 2022–25 (n, stations)",
        fontsize=6.6,
        color="#1A1A1A",
        transform=ax.transAxes,
        va="top",
    )
    for nm, zz, col, y0 in [("Paris", zp, C["risk_dark"], 0.285), ("Málaga", zv, EVD, 0.135)]:
        ax.text(
            0.02,
            y0,
            f"{nm}  {zz[1]:.1f}% → {zz[3]:.1f}%",
            fontsize=7.5,
            color=col,
            transform=ax.transAxes,
            va="top",
        )
        ax.text(
            0.02,
            y0 - 0.068,
            f"n = {int(zz[0]):,} → {int(zz[2]):,}",
            fontsize=6.6,
            color=col,
            transform=ax.transAxes,
            va="top",
        )


def composite():
    NC, NR = (5, 3)
    MW, MH, MGX, MGY = (1.24, 0.86, 0.1, 0.24)
    BW = NC * MW + (NC - 1) * MGX
    H2, H3 = (1.95, 2.1)
    L0 = 0.52
    TOP, G1, G2, BOT = (0.32, 0.62, 0.62, 0.12)
    HF = BOT + H3 + G2 + H2 + G1 + (NR * MH + (NR - 1) * MGY) + TOP
    fig = plt.figure(figsize=(7.2, HF))

    def ax_at(x, y, w, h):
        return fig.add_axes([x / 7.2, y / HF, w / 7.2, h / HF])

    order = ss.country.tolist()
    y_top = HF - TOP
    for i, cname in enumerate(order):
        r, cix = divmod(i, NC)
        ax = ax_at(L0 + cix * (MW + MGX), y_top - (r + 1) * MH - r * MGY, MW, MH)
        mini(ax, cname, first_col=cix == 0, last_row=r == NR - 1)
    fig.text((L0 - 0.4) / 7.2, (y_top - 0.02) / HF, "a", fontsize=13, fontweight="bold", va="top")
    fig.text(
        (L0 - 0.06) / 7.2,
        (y_top + 0.17) / HF,
        "Floodplain share (%) of charging cohorts by station opening year — dashed line: national fuel benchmark",
        fontsize=8.0,
        color="#1A1A1A",
        va="top",
    )
    yx = y_top - NR * MH - (NR - 1) * MGY - 0.3
    fig.text(
        (L0 + BW / 2) / 7.2,
        yx / HF,
        "Station opening year",
        fontsize=8.2,
        color="#1A1A1A",
        ha="center",
        va="top",
    )
    y2 = BOT + H3 + G2
    ax_b = ax_at(0.95, y2, 1.75, H2)
    ax_c = ax_at(3.42, y2, 1.55, H2)
    ax_d = ax_at(5.54, y2, 1.55, H2)
    panel_b(ax_b)
    panel_c(ax_c)
    panel_d(ax_d)
    for lt, ax_ in [("b", ax_b), ("c", ax_c), ("d", ax_d)]:
        p_ = ax_.get_position()
        fig.text(p_.x0 - 0.075, p_.y1 + 0.012, lt, fontsize=13, fontweight="bold", va="top")
    y3 = BOT
    ax_e = ax_at(0.1, y3, 2.3, H3)
    ax_f = ax_at(2.62, y3, 2.3, H3)
    ax_g = ax_at(5.12, y3, 2.0, H3)
    city_map(ax_e, 2046, "Paris", f"drifting into the floodplain ({_metro_slope(2046)} pp/yr)", 5, "e")
    city_map(ax_f, 1677, "Málaga", f"retreating from the floodplain ({_metro_slope(1677)} pp/yr)", 5, "f")
    panel_leg(ax_g)
    ax_e.set_zorder(0)
    for _ax in fig.axes:
        if _ax is not ax_e:
            _ax.set_zorder(max(_ax.get_zorder(), 1))
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Figure3"))


if __name__ == "__main__":
    main()
