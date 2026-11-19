import json
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm

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
G = json.load(open(os.path.join(DATA, "global_summary.json")))
ct = pd.read_csv(os.path.join(DATA, "country_stats.csv"))
iea_w = pd.read_csv(os.path.join(DATA, "iea_world_by_speed.csv")).set_index("year")
hex_ev = pd.read_csv(os.path.join(DATA, "hexbin_charging_stations.csv"))
hex_gs = pd.read_csv(os.path.join(DATA, "hexbin_fuel_stations.csv"))
poly_ev = np.array(json.load(open(os.path.join(DATA, "hexagon_charging.json"))))
poly_gs = np.array(json.load(open(os.path.join(DATA, "hexagon_fuel.json"))))
cm_ev = LinearSegmentedColormap.from_list("ev", ["#EAF4F2", "#9CCFC8", C["ev"], "#053B35"])
cm_gs = LinearSegmentedColormap.from_list("coral", ["#FDF3EE", "#F6C3A5", "#E4795B", "#7F2A11"])
NAME = {"United Kingdom": "UK", "United States of America": "United States"}
X180 = 17005000.0
XL, XR = (-X180 * (170 / 180), X180)


def hexmap(ax, cells, poly, cmap, n_label, letter=None):
    style.draw_world(ax)
    norm = LogNorm(vmin=1, vmax=1000)
    verts = cells[["x", "y"]].values[:, None, :] + poly[None, :, :]
    pc = PolyCollection(
        verts,
        array=cells["stations"].values,
        cmap=cmap,
        norm=norm,
        edgecolors="face",
        linewidths=0.12,
        zorder=4,
        rasterized=True,
    )
    ax.add_collection(pc)
    ax.set_xlim(XL, XR)
    ax.set_ylim(-6350000.0, 7850000.0)
    cax = ax.inset_axes([0.385, 0.09, 0.4, 0.054])
    cb = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([1, 10, 100, 1000])
    cb.set_ticklabels(["1", "10", "100", "1,000"])
    cb.ax.tick_params(labelsize=6.6, length=2, width=0.7, pad=1.5)
    cb.outline.set_linewidth(0.6)
    cax.set_title("Stations per cell", fontsize=7.2, pad=2)
    ax.text(
        0.992,
        0.962,
        n_label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.4,
        color="#2E2E2E",
        zorder=30,
        bbox=dict(fc="white", ec="none", alpha=0.9, pad=1.8),
    )
    if letter:
        ax.text(
            0.011, 0.972, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", zorder=10
        )


def panel_a(ax, letter=None):
    hexmap(ax, hex_ev, poly_ev, cm_ev, f"n = {G['n_ev']:,} charging stations", letter)


def panel_c(ax, letter=None):
    hexmap(ax, hex_gs, poly_gs, cm_gs, f"n = {G['n_gas']:,} fuel stations", letter)


def panel_b(ax):
    yrs = iea_w.index.values
    m = yrs >= 2015
    slow = iea_w["slow"].values[m] / 1000000.0
    fast = iea_w["fast"].values[m] / 1000000.0
    ultra = iea_w["ultra"].values[m] / 1000000.0
    ax.stackplot(
        yrs[m],
        slow,
        fast,
        ultra,
        colors=["#8CC9C0", "#2E8B7F", "#0B4F47"],
        labels=["Slow (≤22 kW)", "Fast (>22–150 kW)", "Ultra-fast (>150 kW)"],
        lw=0,
    )
    tot = slow + fast + ultra
    ax.plot(yrs[m], tot, color=C["ink"], lw=1.1)
    ax.text(2024.7, tot[-1] + 0.12, f"{tot[-1]:.1f} M", ha="right", fontsize=8.6, fontweight="bold")
    f15 = 100 * (fast[0] + ultra[0]) / tot[0]
    f25 = 100 * (fast[-1] + ultra[-1]) / tot[-1]
    ax.text(
        0.07,
        0.55,
        f"fast + ultra-fast\n{f15:.0f}% → {f25:.0f}%",
        transform=ax.transAxes,
        fontsize=7.4,
        color="#0B4F47",
        va="top",
    )
    ax.set_xlim(2015, 2025.4)
    ax.set_ylim(0, 7.6)
    ax.set_xticks([2015, 2020, 2025])
    ax.set_ylabel("Public charging points (millions)")
    ax.legend(
        loc="upper left",
        frameon=False,
        bbox_to_anchor=(0.0, 1.005),
        handlelength=1.1,
        labelspacing=0.3,
        fontsize=6.8,
    )


def panel_d(ax):
    t = ct[(ct.gas_n >= 3000) & (ct.ev_n >= 1000)].nlargest(12, "ev_per_gas").copy()
    t["short"] = t.country.map(lambda c: NAME.get(c, c))
    yy = np.arange(len(t))[::-1]
    glob = G["n_ev"] / G["n_gas"]
    ax.hlines(yy, 0.02, t.ev_per_gas, color="#F3C2A9", lw=2.2)
    ax.scatter(t.ev_per_gas, yy, s=36, c=C["ev"], zorder=3)
    ax.axvline(1, color=C["gas"], lw=1.0, ls="--")
    ax.axvline(glob, color=C["gray"], lw=0.9, ls=":")
    ax.text(1.07, -1.02, "parity", fontsize=7, color=C["gas"], ha="left", va="center")
    ax.text(glob * 0.93, -1.02, f"{glob:.2f}", fontsize=7, color=C["gray"], ha="right", va="center")
    for yi, v in zip(yy, t.ev_per_gas):
        ax.text(v * 1.2, yi, f"{v:.1f}", fontsize=7.4, va="center", color=C["ev_dark"])
    ax.set_xscale("log")
    ax.set_xlim(0.13, 46)
    ax.set_xticks([0.3, 1, 3, 10, 30])
    ax.set_xticklabels(["0.3", "1", "3", "10", "30"])
    ax.minorticks_off()
    ax.set_yticks(yy)
    ax.set_yticklabels(t.short, fontsize=7.2)
    ax.set_ylim(-1.55, 12.35)
    ax.set_xlabel("Charging stations per fuel station")


def composite():
    ASP = (XR - XL) / (7850000.0 + 6350000.0)
    H = 1.874
    WM = H * ASP
    L, GUT, R = (0.09, 0.8, 0.07)
    WC = 7.2 - L - WM - GUT - R
    BOT, GAPR, TOP = (0.5, 0.28, 0.06)
    HF = BOT + 2 * H + GAPR + TOP
    fig = plt.figure(figsize=(7.2, HF))

    def ax_at(x_in, y_in, w_in, h_in):
        return fig.add_axes([x_in / 7.2, y_in / HF, w_in / 7.2, h_in / HF])

    y2 = BOT
    y1 = BOT + H + GAPR
    ax_a = ax_at(L, y1, WM, H)
    ax_b = ax_at(L + WM + GUT, y1, WC, H)
    ax_c = ax_at(L, y2, WM, H)
    ax_d = ax_at(L + WM + GUT, y2, WC, H)
    panel_a(ax_a, "a")
    panel_b(ax_b)
    panel_c(ax_c, "c")
    panel_d(ax_d)
    xg = (L + WM + 0.03) / 7.2
    for lt, ax_ in [("b", ax_a), ("d", ax_c)]:
        p_ = ax_.get_position()
        fig.text(xg, p_.y1 - 0.028 * p_.height, lt, fontsize=13, fontweight="bold", va="top", ha="left")
    return fig


def main():
    fig = composite()
    style.export(fig, os.path.join(HERE, "Figure1"))


if __name__ == "__main__":
    main()
