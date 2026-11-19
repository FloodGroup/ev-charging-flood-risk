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
from matplotlib.patches import FancyArrowPatch

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
S3 = json.load(open(f"{DATA}/cascade_stats.json"))
rp = pd.read_csv(f"{DATA}/cascade_rp_curve.csv")
cy = pd.read_csv(f"{DATA}/cascade_by_country.csv", index_col=0)
mt = pd.read_csv(f"{DATA}/cascade_metro_named.csv", index_col=0)
cm2 = pd.read_csv(f"{DATA}/urban_centre_coordinates.csv")
mt = mt.join(cm2.set_index("metro"), how="left")
XMAX = 17005000.0
XL, XR = (-XMAX * (170 / 180), XMAX)
GRID_C = "#6A51A3"
RELIEF_CM = LinearSegmentedColormap.from_list(
    "relief", [(0.0, "#9C988E"), (0.62, "#F4F2EC"), (1.0, "#FFFFFF")]
)
RESC_C = "#0072B2"
FAIL_C = "#8C2D12"
ABINS = [1.0, 1.05, 1.15, 1.3, 1.6, 99]
ACOL = ["#C9CDE8", "#A79ED2", "#8470BA", "#59389B", "#2B115C"]
ALAB = ["<1.05", "1.05–1.15", "1.15–1.3", "1.3–1.6", ">1.6"]
CONT_C = {
    "Asia": "#D55E00",
    "Europe": "#0E8074",
    "North America": "#2E6FB7",
    "South America": "#E7B33C",
    "Africa": "#6A51A3",
    "Oceania": "#8A8A8A",
}
CONT_ORDER = ["Asia", "Europe", "North America", "South America", "Oceania", "Africa"]
SHORT = {"United States of America": "United States", "United Kingdom": "UK"}
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
        cmap=ListedColormap(["#C9DCEB"]),
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


def city_A(mid):
    r = mt.loc[mt.index == mid]
    return float(r.A_total.iloc[0]) if len(r) else np.nan


def city_Ared(mid):
    r = mt.loc[mt.index == mid]
    return float(r.A_redundant.iloc[0]) if len(r) else np.nan


CITY_INFO = [(2288, "Düsseldorf"), (12703, "Tainan"), (10912, "Handan"), (1048, "Québec City")]
CITIES = sorted(CITY_INFO, key=lambda t: -city_A(t[0]))


def site_order(z):
    nd = z["site_n_disp"]
    o = np.argsort(-nd)
    return o[nd[o] > 0]


def story_sites(z, cover=0.7, max_sites=3, min_dep=10):
    oo = site_order(z)
    nd = z["site_n_disp"]
    tot = max(int(nd[oo].sum()), 1)
    sel, acc = ([], 0)
    for i in oo:
        if len(sel) >= max_sites or int(nd[i]) < min_dep:
            break
        sel.append(i)
        acc += int(nd[i])
        if acc / tot >= cover:
            break
    return sel


def draw_feeder_view(ax, mid, name, letter):
    z = np.load(f"{DATA}/private/feeder_map_{mid}.npz")
    W0, S0, E0, N0 = z["bbox"]
    sid = z["site_id"]
    sl, sb = (z["site_lon"], z["site_lat"])
    sp = z["site_p"]
    sfail = z["site_fail"].astype(bool)
    pos = {int(s): (float(sl[k]), float(sb[k])) for k, s in enumerate(sid)}
    sel = story_sites(z)
    sel_ids = [int(sid[i]) for i in sel]
    lat0 = float(np.mean(z["ch_lat"]))
    kx = 111.32 * np.cos(np.deg2rad(lat0))
    ky = 110.54
    gonly = z["gridonly"].astype(bool)
    resc = z["rescued"].astype(bool)
    dep_any = np.isin(z["ch_s1"], sel_ids) & gonly
    bk_id = z["ch_bk_draw"] if "ch_bk_draw" in z.files else z["ch_s2"]
    bk_pairs = []
    for s in sel_ids:
        dep = (z["ch_s1"] == s) & gonly & resc
        ub, cn = np.unique(bk_id[dep & (bk_id >= 0)], return_counts=True)
        for b, c in zip(ub, cn):
            if c >= 1 and int(b) in pos:
                bk_pairs.append((s, int(b), int(c)))
    xs = [pos[s][0] for s in sel_ids] + list(z["ch_lon"][dep_any]) + [pos[b][0] for _, b, _ in bk_pairs]
    ys = [pos[s][1] for s in sel_ids] + list(z["ch_lat"][dep_any]) + [pos[b][1] for _, b, _ in bk_pairs]
    cx, cyy = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
    side = max((max(xs) - min(xs)) * kx, (max(ys) - min(ys)) * ky) * 1.24
    side = max(side, 3.0)
    for s in sel_ids:
        px, py = pos[s]
        side = max(side, abs(px - cx) * kx / 0.4, abs(py - cyy) * ky / 0.4)
    w, e = (cx - side / kx / 2, cx + side / kx / 2)
    s_, n = (cyy - side / ky / 2, cyy + side / ky / 2)
    ax.set_facecolor("#F4F2EC")
    draw_basemap(ax, mid, 1.0)
    fl = z["flood"].astype(float)
    fl[fl == 0] = np.nan
    ax.imshow(
        fl,
        extent=(W0, E0, S0, N0),
        origin="upper",
        cmap=ListedColormap(["#7FB2DD"]),
        alpha=0.55,
        interpolation="nearest",
        zorder=1,
        aspect="auto",
    )
    LINE_C = "#57544E"
    for key, lsty, lwd in [("osm_seg", "solid", 1.35), ("gf_seg", (0, (3.0, 2.0)), 1.15)]:
        sg = z[key]
        if len(sg):
            ax.add_collection(
                LineCollection(
                    sg.reshape(-1, 2, 2),
                    colors=LINE_C,
                    linewidths=lwd,
                    linestyles=lsty,
                    zorder=2,
                    rasterized=True,
                )
            )
    if "plant_lon" in z.files and len(z["plant_lon"]):
        ax.scatter(
            z["plant_lon"], z["plant_lat"], s=42, c="#3C3C3C", marker="^", lw=0.7, edgecolor="white", zorder=6
        )
    ok = ~(z["direct"] | gonly)
    ax.scatter(
        z["ch_lon"][ok], z["ch_lat"][ok], s=3.0, c="#1A1A1A", lw=0, zorder=3, alpha=0.8, rasterized=True
    )
    ax.scatter(
        z["ch_lon"][z["direct"]],
        z["ch_lat"][z["direct"]],
        s=7,
        c=C["risk"],
        lw=0.25,
        edgecolor="white",
        zorder=6,
    )
    ctx = sfail & ~np.isin(sid, sel_ids) & (sl >= w) & (sl <= e) & (sb >= s_) & (sb <= n)
    ax.scatter(
        sl[ctx], sb[ctx], s=10, c="#C4907E", marker="s", lw=0.3, edgecolor="white", zorder=5, rasterized=True
    )
    good = ~sfail & (sl >= w) & (sl <= e) & (sb >= s_) & (sb <= n)
    ax.scatter(sl[good], sb[good], s=5, c="#B9B4A8", marker="s", lw=0, zorder=4, rasterized=True)
    routed = z["routed"].astype(bool) if "routed" in z.files else np.zeros(len(z["ch_lon"]), bool)
    rescue_routed = (
        z["rescue_routed"].astype(bool) if "rescue_routed" in z.files else np.zeros(len(z["ch_lon"]), bool)
    )
    all_r = []
    for s in sel_ids:
        rk, nk = (f"rflow_{s}_seg", f"rflow_{s}_n")
        if rk in z.files and len(z[rk]):
            rs_ = np.asarray(z[rk], float).reshape(-1, 4)
            rn_ = z[nk].astype(float)
            all_r.append((rs_, rn_))
    if all_r:
        segs = np.vstack([a[0] for a in all_r])
        wts = np.concatenate([a[1] for a in all_r])
        a = np.minimum(segs[:, 0:2], segs[:, 2:4])
        b = np.maximum(segs[:, 0:2], segs[:, 2:4])
        keys = np.round(np.c_[a, b], 6)
        order = np.lexsort((-wts, keys[:, 3], keys[:, 2], keys[:, 1], keys[:, 0]))
        keys_s, segs_s, wts_s = (keys[order], segs[order], wts[order])
        keep = np.r_[True, np.any(np.diff(keys_s, axis=0) != 0, axis=1)]
        segs_u, wts_u = (segs_s[keep], wts_s[keep])
        wds = np.clip(0.9 + 0.45 * np.sqrt(wts_u), 0.9, 2.2)
        ax.add_collection(
            LineCollection(
                segs_u.reshape(-1, 2, 2), colors="white", linewidths=wds + 1.4, zorder=6.2, capstyle="round"
            )
        )
        ax.add_collection(
            LineCollection(
                segs_u.reshape(-1, 2, 2),
                colors=RESC_C,
                alpha=0.95,
                linewidths=wds,
                zorder=6.25,
                capstyle="round",
            )
        )
    for s in sel_ids:
        dep = (z["ch_s1"] == s) & gonly & resc
        targets = (
            bk_id[dep & (bk_id >= 0)] if "ch_bk_draw" in z.files else z["ch_s2"][dep & (z["ch_s2"] >= 0)]
        )
        for b in np.unique(targets):
            sk = f"bsitecon_{int(b)}"
            if sk in z.files and len(z[sk]):
                ax.add_collection(
                    LineCollection(
                        z[sk].reshape(-1, 2, 2),
                        colors=RESC_C,
                        alpha=0.9,
                        linewidths=1.15,
                        zorder=6.28,
                        capstyle="round",
                    )
                )
    for s in sel_ids:
        fk, nk = (f"flow_{s}_seg", f"flow_{s}_n")
        if fk in z.files and len(z[fk]):
            fs_, fn_ = (z[fk], z[nk].astype(float))
            wds = 0.9 + 0.55 * np.sqrt(fn_)
            ax.add_collection(
                LineCollection(
                    fs_.reshape(-1, 2, 2), colors="white", linewidths=wds + 1.5, zorder=6.4, capstyle="round"
                )
            )
            ax.add_collection(
                LineCollection(
                    fs_.reshape(-1, 2, 2),
                    colors=GRID_C,
                    alpha=0.95,
                    linewidths=wds,
                    zorder=6.5,
                    capstyle="round",
                )
            )
        sk = f"sitecon_{s}"
        if sk in z.files and len(z[sk]):
            ax.add_collection(
                LineCollection(z[sk].reshape(-1, 2, 2), colors=GRID_C, alpha=0.9, linewidths=1.1, zorder=6.5)
            )
    if "con_seg" in z.files and len(z["con_seg"]):
        cs_ = z["con_seg"][np.isin(z["con_site"], sel_ids)]
        if len(cs_):
            ax.add_collection(
                LineCollection(
                    cs_.reshape(-1, 2, 2), colors="#A79ED2", linewidths=0.55, alpha=0.9, zorder=6.55
                )
            )
    for s in sel_ids:
        dep = (z["ch_s1"] == s) & gonly & ~routed
        if dep.any():
            x0, y0 = pos[s]
            seg = np.stack(
                [
                    np.column_stack([np.full(dep.sum(), x0), np.full(dep.sum(), y0)]),
                    np.column_stack([z["ch_lon"][dep], z["ch_lat"][dep]]),
                ],
                axis=1,
            )
            ax.add_collection(
                LineCollection(
                    seg, colors=GRID_C, linewidths=0.5, alpha=0.5, linestyles=(0, (2, 2)), zorder=6.6
                )
            )
        miss = (z["ch_s1"] == s) & gonly & resc & ~rescue_routed & (z["ch_s2"] >= 0)
        if miss.any():
            for ci in np.where(miss)[0]:
                b = int(z["ch_s2"][ci])
                if b not in pos:
                    continue
                x1, y1 = pos[b]
                seg = np.array([[[z["ch_lon"][ci], z["ch_lat"][ci]], [x1, y1]]], dtype=float)
                ax.add_collection(
                    LineCollection(
                        seg, colors=RESC_C, linewidths=0.55, alpha=0.55, linestyles=(0, (2, 2)), zorder=6.7
                    )
                )
    rest = gonly & ~dep_any
    ax.scatter(
        z["ch_lon"][rest & ~resc],
        z["ch_lat"][rest & ~resc],
        s=6,
        c=GRID_C,
        lw=0.25,
        edgecolor="white",
        alpha=0.85,
        zorder=8,
    )
    ax.scatter(
        z["ch_lon"][rest & resc],
        z["ch_lat"][rest & resc],
        s=6,
        facecolor="white",
        edgecolor=RESC_C,
        lw=0.55,
        alpha=0.85,
        zorder=8,
    )
    off = dep_any & ~resc
    rec = dep_any & resc
    ax.scatter(z["ch_lon"][off], z["ch_lat"][off], s=11, c=GRID_C, lw=0.35, edgecolor="white", zorder=9)
    ax.scatter(
        z["ch_lon"][rec], z["ch_lat"][rec], s=11, facecolor="white", edgecolor=RESC_C, lw=0.75, zorder=9
    )
    drawn_bk = set()
    for s, b, cnum in bk_pairs:
        if b not in drawn_bk:
            drawn_bk.add(b)
            j = np.where(sid == b)[0][0]
            x1, y1 = pos[b]
            if sp[j] >= 0.5:
                ax.scatter([x1], [y1], s=52, c=FAIL_C, marker="s", lw=0.7, edgecolor="white", zorder=11)
            else:
                ax.scatter(
                    [x1], [y1], s=52, facecolor="white", marker="s", edgecolor=RESC_C, lw=1.3, zorder=11
                )
    for s in sel_ids:
        x0, y0 = pos[s]
        ax.scatter([x0], [y0], s=64, c=FAIL_C, marker="s", lw=0.8, edgecolor="white", zorder=12)
    ax.set_xlim(w, e)
    ax.set_ylim(s_, n)
    ax.set_aspect(1 / np.cos(np.deg2rad((s_ + n) / 2)))
    ax.set_xticks([])
    ax.set_yticks([])
    for spn in ax.spines.values():
        spn.set_visible(True)
        spn.set_edgecolor("#D8D5CC")
        spn.set_linewidth(0.6)
    ax.text(
        0.0,
        1.028,
        letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=C["ink"],
    )
    ax.text(
        0.115,
        1.04,
        name,
        transform=ax.transAxes,
        fontsize=8.4,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=C["ink"],
    )
    ar_ = city_Ared(mid)
    t_teal = f"→ ×{ar_:.2f}"
    ax.text(
        1.0,
        1.04,
        t_teal,
        transform=ax.transAxes,
        fontsize=8.0,
        va="bottom",
        ha="right",
        color=RESC_C,
        fontweight="bold",
    )
    ax.text(
        1.0 - 0.036 * len(t_teal) - 0.014,
        1.04,
        f"×{city_A(mid):.2f}",
        transform=ax.transAxes,
        fontsize=8.0,
        va="bottom",
        ha="right",
        color=GRID_C,
        fontweight="bold",
    )
    ax.annotate(
        "N",
        xy=(0.974, 0.985),
        xytext=(0.974, 0.892),
        xycoords="axes fraction",
        textcoords="axes fraction",
        ha="center",
        va="top",
        fontsize=7.4,
        fontweight="bold",
        color="#1A1A1A",
        zorder=31,
        path_effects=[pe.withStroke(linewidth=1.8, foreground="white", alpha=0.95)],
        arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#1A1A1A", shrinkA=1, shrinkB=0),
    )
    skm = max(1, int(round(side / 5)))
    dlon = skm / kx
    x0 = e - dlon - (e - w) * 0.05
    y0 = s_ + (n - s_) * 0.05
    ax.plot(
        [x0, x0 + dlon],
        [y0, y0],
        color="#3C3C3C",
        lw=1.3,
        zorder=30,
        path_effects=[pe.withStroke(linewidth=2.8, foreground="white", alpha=0.9)],
    )
    ax.add_patch(
        Rectangle(
            (x0 - (e - w) * 0.025, s_ + (n - s_) * 0.018),
            dlon + (e - w) * 0.05,
            (n - s_) * 0.105,
            fc="white",
            ec="none",
            alpha=0.88,
            zorder=29,
        )
    )
    ax.text(
        x0 + dlon / 2,
        y0 + (n - s_) * 0.022,
        f"{skm} km",
        ha="center",
        fontsize=6.6,
        color="#1A1A1A",
        zorder=30,
    )


def legend_strip(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    GL = 0.045
    ROWS = [0.8, 0.5, 0.2]
    COLS = {
        0.015: [
            ("sq", FAIL_C, 26, None, "failed feeder site (P≥0.5)"),
            ("sqo", RESC_C, 26, None, "backup feeder site, dry"),
            ("tri", "#3C3C3C", 24, None, "power plant"),
        ],
        0.265: [
            ("line", "#57544E", "solid", 1.6, "power line (OSM)"),
            ("dline2", "#57544E", (0, (2.4, 1.7)), 1.6, "power line (gridfinder)"),
            ("patch", "#B4CFE4", None, 0, "100-yr floodplain"),
        ],
        0.51: [
            ("line", GRID_C, "solid", 2.6, "supply path along the lines"),
            ("line", RESC_C, "solid", 2.0, "rescue path along the lines"),
            ("dot", "#1A1A1A", 6, None, "unaffected charger"),
        ],
        0.77: [
            ("dot", GRID_C, 12, None, "charger offline via grid"),
            ("doto", RESC_C, 12, None, "recoverable after rerouting"),
            ("dot", C["risk"], 8, None, "flooded charger"),
        ],
    }
    for x0, items in COLS.items():
        for (kind, col, aa, bb, lab), ys in zip(items, ROWS):
            if kind in ("line", "dline", "dline2"):
                ax.plot(
                    [x0, x0 + GL],
                    [ys, ys],
                    color=col,
                    lw=bb or 1.2,
                    ls="solid" if kind == "line" else aa,
                    solid_capstyle="round",
                )
                if kind == "dline":
                    ar = FancyArrowPatch(
                        (x0 + GL * 0.78, ys),
                        (x0 + GL, ys),
                        arrowstyle="-|>",
                        mutation_scale=7,
                        lw=0.0,
                        color=col,
                    )
                    ax.add_patch(ar)
            elif kind == "patch":
                ax.add_patch(Rectangle((x0, ys - 0.12), GL, 0.24, fc=col, ec="none"))
            elif kind == "sqo":
                ax.scatter(
                    [x0 + GL / 2],
                    [ys],
                    marker="s",
                    s=aa,
                    facecolor="white",
                    edgecolor=col,
                    lw=1.2,
                    clip_on=False,
                )
            elif kind == "doto":
                ax.scatter(
                    [x0 + GL / 2],
                    [ys],
                    marker="o",
                    s=aa,
                    facecolor="white",
                    edgecolor=col,
                    lw=0.9,
                    clip_on=False,
                )
            elif kind == "tri":
                ax.scatter(
                    [x0 + GL / 2], [ys], marker="^", s=aa, c=col, lw=0.6, edgecolor="white", clip_on=False
                )
            else:
                mk = "s" if kind == "sq" else "o"
                ax.scatter(
                    [x0 + GL / 2], [ys], marker=mk, s=aa, c=col, lw=0.3, edgecolor="white", clip_on=False
                )
            ax.text(x0 + GL + 0.009, ys, lab, fontsize=6.6, va="center", color="#000000")


LABELS_E = {
    2288: ("Düsseldorf", -2750000.0, -1950000.0, "right", 0.24),
    1048: ("Québec City", 950000.0, 1200000.0, "right", 0.15),
    12703: ("Tainan", 1900000.0, -750000.0, "left", 0.15),
    10912: ("Handan", 2700000.0, 950000.0, "left", -0.15),
}


def panel_e(ax, letter=None):
    style.draw_world(ax)
    d = mt[(mt.direct >= 20) & (mt.n >= 100)].dropna(subset=["GCPNT_LAT", "A_total"])
    x, y = style.to_rob(d.GCPNT_LON.values, d.GCPNT_LAT.values)
    ci = np.clip(np.digitize(d.A_total.values, ABINS) - 1, 0, 4)
    s = np.clip(1.1 * np.sqrt(d.total.values), 3.5, 121)
    o = np.argsort(-s)
    ax.scatter(
        x[o],
        y[o],
        s=s[o],
        c=np.array(ACOL)[ci[o]],
        alpha=0.92,
        lw=0.3,
        edgecolor="white",
        zorder=4,
        rasterized=True,
    )
    ax.set_xlim(XL, XR)
    ax.set_ylim(-6350000.0, 7850000.0)
    halo_w = [pe.withStroke(linewidth=2.2, foreground="white", alpha=0.9)]
    ax.text(
        0.992,
        0.975,
        f"RP100: {S3['direct_rp100']:,} flooded → {S3['total_rp100']:,.0f} expected offline (×{S3['A_rp100']:.2f})",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
        color="#2E2E2E",
        zorder=30,
        path_effects=halo_w,
    )
    pos = {}
    for mid in LABELS_E:
        r = d[d.index == mid]
        if len(r):
            xx, yy = style.to_rob(np.array([r.GCPNT_LON.iloc[0]]), np.array([r.GCPNT_LAT.iloc[0]]))
            pos[mid] = (xx[0], yy[0], float(r.A_total.iloc[0]))
    for mid, (nm, dx, dy, ha, rad) in LABELS_E.items():
        if mid not in pos:
            continue
        xx, yy, aa = pos[mid]
        ax.annotate(
            f"{nm} ×{aa:.2f}",
            (xx, yy),
            xytext=(xx + dx, yy + dy),
            fontsize=6.6,
            ha=ha,
            va="center",
            zorder=31,
            color="#1A1A1A",
            path_effects=[pe.withStroke(linewidth=1.6, foreground="white", alpha=0.9)],
            arrowprops=dict(
                arrowstyle="-",
                lw=0.55,
                color="#6E6E6E",
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=0,
                shrinkB=1.5,
            ),
        )
    lx = ax.inset_axes([0.405, 0.004, 0.575, 0.19])
    lx.set_xlim(0, 1)
    lx.set_ylim(0, 1)
    lx.axis("off")
    lx.text(0.01, 0.88, "dot colour — amplification ×A", fontsize=6.3, va="center", color="#000000")
    x0, bw, y0, y1 = (0.01, 0.196, 0.42, 0.72)
    for i, (cc, lb) in enumerate(zip(ACOL, ALAB)):
        lx.add_patch(Rectangle((x0 + i * bw, y0), bw, y1 - y0, fc=cc, ec="white", lw=0.5))
        lx.text(
            x0 + (i + 0.5) * bw,
            (y0 + y1) / 2,
            lb,
            fontsize=6.0,
            ha="center",
            va="center",
            color="white" if i >= 2 else "#1A1A1A",
        )
    lx.text(0.01, 0.1, "dot area — expected offline", fontsize=6.3, va="center", color="#000000")
    inset_w_pt = 0.575 * WM * 72.0
    xc = 0.525
    for nn in (100, 1000, 10000):
        s_pt2 = float(np.clip(1.1 * np.sqrt(nn), 3.5, 121))
        r_frac = np.sqrt(s_pt2 / np.pi) / inset_w_pt
        lx.scatter([xc], [0.1], s=s_pt2, c="#9A9A9A", lw=0.3, edgecolor="white", clip_on=False)
        lx.text(xc + r_frac + 0.012, 0.1, f"{nn:,}", fontsize=6.2, va="center", ha="left", color="#1A1A1A")
        xc += r_frac + 0.012 + 0.026 * len(f"{nn:,}") + 0.03
    if letter:
        ax.text(
            0.011, 0.972, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", zorder=10
        )


def panel_f(ax):
    x = rp.rp.values
    lo = np.minimum.reduce([rp.A_hard.values, rp.A_beta04.values, rp.A_beta08.values])
    hi = np.maximum.reduce([rp.A_soft.values, rp.A_beta04.values, rp.A_beta08.values])
    ax.fill_between(x, lo, hi, color=GRID_C, alpha=0.13, lw=0, label="fragility parameter band")
    ax.plot(x, rp.A_fragility, "-o", color=GRID_C, lw=1.6, ms=3.6, label="fragility (θ 0.5/2.0 m, β 0.6)")
    ax.plot(x, rp.A_threshold, ":", color="#8A8A8A", lw=1.2, label="deterministic thresholds")
    ax.plot(x, rp.A_redundant, "--s", color=RESC_C, lw=1.3, ms=3.0, label="with backup feeder")
    y_pri = float(rp.loc[rp.rp == 100, "A_fragility"].iloc[0])
    y_bak = float(rp.loc[rp.rp == 100, "A_redundant"].iloc[0])
    ax.annotate(
        f"×{y_pri:.3f}",
        (100, y_pri),
        xytext=(100, y_pri + 0.055),
        ha="center",
        va="bottom",
        fontsize=7.6,
        fontweight="bold",
        color=GRID_C,
        arrowprops=dict(arrowstyle="-", lw=0.6, color=GRID_C),
    )
    ax.annotate(
        f"×{y_bak:.3f}",
        (100, y_bak),
        xytext=(100, y_bak - 0.055),
        ha="center",
        va="top",
        fontsize=7.4,
        fontweight="bold",
        color=RESC_C,
        arrowprops=dict(arrowstyle="-", lw=0.6, color=RESC_C),
    )
    ax.set_xscale("log")
    ax.set_xticks([10, 20, 50, 100, 200, 500])
    ax.set_xticklabels(["10", "20", "50", "100", "200", "500"])
    ax.minorticks_off()
    ax.set_xlabel("Flood return period (years)")
    ax.set_ylabel("Grid amplification (×A)")
    ax.set_ylim(1.0, 1.42)
    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=6.8,
        handlelength=1.6,
        labelspacing=0.35,
        borderaxespad=0.25,
    )


def panel_g(ax):
    m = cy[(cy.n >= 1000) & (cy.direct >= 100)].dropna(subset=["A_total"]).copy()
    ax.axhline(1, color="#B0B0B0", lw=0.8, ls="--")
    for cont in CONT_ORDER:
        if cont not in CONT_C:
            continue
        sub = m[m.continent == cont]
        if not len(sub):
            continue
        ax.scatter(
            sub.direct_pct,
            sub.A_total,
            s=np.clip(0.55 * np.sqrt(sub.n), 5, 200),
            c=CONT_C[cont],
            alpha=0.8,
            lw=0.4,
            edgecolor="white",
            zorder=3,
        )
    b = np.polyfit(np.log(m.direct_pct), m.A_total, 1)
    xs = np.geomspace(m.direct_pct.min() * 0.9, m.direct_pct.max() * 1.1, 50)
    ax.plot(xs, b[0] * np.log(xs) + b[1], color="#6E6E6E", lw=1.0, ls=":", zorder=2)
    ax.text(
        0.97,
        0.665,
        f"r = {S3['mech_corr_logdirect_A']:.2f} ({S3['n_countries_mech']} countries)".replace("-", "−"),
        transform=ax.transAxes,
        fontsize=6.8,
        color="#3C3C3C",
        ha="right",
    )
    lab = {
        "Switzerland": (9.2, 1.8, "left"),
        "Brazil": (3.5, 1.7, "center"),
        "United States of America": (14.0, 1.56, "left"),
        "Spain": (2.95, 1.61, "left"),
        "Germany": (2.95, 1.235, "left"),
        "United Kingdom": (4.3, 1.135, "center"),
        "Japan": (30.0, 1.33, "center"),
        "Netherlands": (9.6, 1.055, "center"),
        "China": (18.0, 1.03, "center"),
        "Thailand": (33.0, 1.02, "center"),
    }
    for nm, (tx, ty, ha) in lab.items():
        r = m[m.index == nm]
        if len(r):
            ax.annotate(
                SHORT.get(nm, nm),
                (r.direct_pct.iloc[0], r.A_total.iloc[0]),
                xytext=(tx, ty),
                fontsize=6.6,
                va="center",
                ha=ha,
                zorder=7,
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#8A8A8A", shrinkA=1.0, shrinkB=2.0),
            )
    ax.set_xscale("log")
    ax.set_xticks([3, 5, 10, 20, 50])
    ax.set_xticklabels(["3", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlim(2.8, 62)
    ax.set_ylim(1.0, 2.06)
    ax.set_yticks([1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    ax.set_xlabel("Chargers directly flooded (%)")
    ax.set_ylabel("Grid amplification (×A)")
    CSHORT = {"North America": "N. America", "South America": "S. America"}
    present = [k for k in CONT_ORDER if k in set(m.continent)]
    hd = {
        k: Line2D([], [], marker="o", ls="none", ms=4.6, mfc=CONT_C[k], mec="none", label=CSHORT.get(k, k))
        for k in present
    }
    blank = Line2D([], [], ls="none", label=" ")
    left = [hd[k] for k in ("Asia", "Europe") if k in hd] + [blank]
    right = [hd[k] for k in ("North America", "South America", "Oceania") if k in hd]
    hnd = left + right
    ax.legend(
        handles=hnd,
        loc="upper right",
        frameon=False,
        fontsize=7.2,
        ncol=2,
        columnspacing=0.9,
        handletextpad=0.25,
        labelspacing=0.55,
        borderaxespad=0.3,
    )


def panel_h(ax):
    m = cy[(cy.n >= 1000) & (cy.direct >= 100)].nlargest(10, "grid_add")
    m = m.iloc[::-1]
    y = np.arange(len(m))
    ax.hlines(y, m.A_redundant, m.A_total, color="#C9C9C9", lw=1.6, zorder=2)
    ax.scatter(m.A_total, y, s=26, c=GRID_C, zorder=3)
    ax.scatter(m.A_redundant, y, s=26, facecolors="white", edgecolors=RESC_C, lw=1.3, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT.get(i, i) for i in m.index], fontsize=6.9)
    ax.set_xlabel("Grid amplification (×A, RP100)")
    ax.set_xlim(1.0, max(1.55, float(m.A_total.max()) + 0.05))
    ax.axvline(1, color="#B0B0B0", lw=0.8)
    ax.tick_params(axis="y", length=0)
    ax.text(
        0.02,
        1.02,
        f"global: ×{S3['A_rp100']:.2f} → ×{S3['A_redundant_rp100']:.2f} with backup",
        transform=ax.transAxes,
        fontsize=7.0,
        color="#3C3C3C",
        va="bottom",
    )


def panel_i(ax):
    m = cy.dropna(subset=["A_total"]).nlargest(8, "grid_add")
    rest = float(S3["grid_add_rp100"] - m.grid_add.sum())
    names = [SHORT.get(i, i) for i in m.index] + ["Rest of world"]
    vals = list(m.grid_add / 1000) + [rest / 1000]
    cols = [CONT_C.get(c, "#B0B0B0") for c in m.continent] + ["#B0B0B0"]
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, height=0.68, color=cols, alpha=0.92)
    shares = list(100 * m.grid_add / S3["grid_add_rp100"]) + [100 * rest / S3["grid_add_rp100"]]
    for i, (yi, v, sh) in enumerate(zip(y, vals, shares)):
        if i == 0:
            ax.text(
                v - 0.12,
                yi,
                f"{v:.1f}k · {sh:.0f}%",
                ha="right",
                va="center",
                fontsize=7.4,
                color="white",
                fontweight="bold",
            )
        else:
            ax.text(v + 0.12, yi, f"{v:.1f}k", ha="left", va="center", fontsize=7.4, color="#333333")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.9)
    ax.set_xlim(0, max(11.4, max(vals) * 1.25))
    ax.set_xlabel("Chargers offline via grid only\n(expected, thousands, RP100)", fontsize=8.0)
    ax.tick_params(axis="y", length=0)


ASP = (XR - XL) / (7850000.0 + 6350000.0)
H = 1.874
WM = H * ASP


def composite():
    ASP = (XR - XL) / (7850000.0 + 6350000.0)
    H = 1.874
    WM = H * ASP
    L, GUT, R = (0.09, 0.8, 0.07)
    WC = 7.2 - L - WM - GUT - R
    VW = 1.72
    VMH = 1.66
    VT = 0.17
    VGAP = (7.2 - 0.1 - 4 * VW) / 3
    H3 = 2.02
    BOT, GAPR, TOPV, LEG = (0.5, 0.62, 0.04, 0.42)
    AIR1 = 0.1
    AIR2 = 0.16
    HF = BOT + H3 + GAPR + H + AIR2 + LEG + AIR1 + VMH + VT + TOPV
    fig = plt.figure(figsize=(7.2, HF))

    def ax_at(x, y, w, h):
        return fig.add_axes([x / 7.2, y / HF, w / 7.2, h / HF])

    yv = HF - TOPV - VT - VMH
    yleg = yv - AIR1 - LEG
    y2 = yleg - AIR2 - H
    y3 = BOT
    vx = [0.05 + i * (VW + VGAP) for i in range(4)]
    ax_v = [ax_at(x, yv, VW, VMH) for x in vx]
    ax_leg = ax_at(0.05, yleg + 0.01, 7.1, LEG - 0.02)
    ax_e = ax_at(L, y2, WM, H)
    ax_f = ax_at(L + WM + GUT, y2, WC, H)
    ax_g = ax_at(0.62, y3, 1.8, H3)
    ax_h = ax_at(3.1, y3, 1.6, H3)
    ax_i = ax_at(5.55, y3, 1.58, H3)
    for axv, (mid, nm), lt in zip(ax_v, CITIES, "abcd"):
        draw_feeder_view(axv, mid, nm, lt)
    legend_strip(ax_leg)
    panel_e(ax_e, "e")
    panel_f(ax_f)
    panel_g(ax_g)
    panel_h(ax_h)
    panel_i(ax_i)
    p_ = ax_f.get_position()
    fig.text(p_.x0 - 0.062, p_.y1 - 0.004, "f", fontsize=13, fontweight="bold", va="top", ha="left")
    for lt, ax_ in [("g", ax_g), ("h", ax_h), ("i", ax_i)]:
        p_ = ax_.get_position()
        fig.text(p_.x0 - 0.062, p_.y1 + 0.01, lt, fontsize=13, fontweight="bold", va="top", ha="left")
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Figure4"))


if __name__ == "__main__":
    main()
