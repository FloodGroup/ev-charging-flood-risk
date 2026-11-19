#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_style.py -- shared plotting conventions for every figure in the package.

Layout language (Nature double column): full figure width 7.2 in (183 mm);
world-map row height H = 1.874 in; Robinson viewport aspect ratio ~2.329;
base font 8 pt, tick labels 7.4 pt, panel titles 7.6 pt bold left-aligned;
PNG at 400 dpi plus vector PDF (TrueType fonts embedded, fonttype 42).
Fonts: Arial on Windows, Liberation Sans on Linux -- the two are metrically
compatible, so the same script produces near-identical output on both.


Project: Flood exposure of the global electric vehicle charging network and
         its grid dependence (reproducibility package)
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Locations (derived from this file; see shared/paths.py for the full scheme)
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))          # <package>/shared
PACKAGE_ROOT = os.path.dirname(_HERE)                       # <package>
FIGDIR = os.path.join(PACKAGE_ROOT, "shared", "_draft_figures")

_SOURCE_ROOT = os.environ.get("EVFR_SOURCE_ROOT") or next(
    (p for p in (os.path.join(PACKAGE_ROOT, "source_data"),
                 os.path.join(os.path.dirname(PACKAGE_ROOT), "02_input_data"))
     if os.path.isdir(p)),
    os.path.join(PACKAGE_ROOT, "source_data"))
_NE110 = os.environ.get(
    "EVFR_NE110M",
    os.path.join(_SOURCE_ROOT, "boundaries_natural_earth",
                 "ne_110m_countries_simplified.geojson"))

# ---------------------------------------------------------------------------
# Palette (identical to the published figures)
# ---------------------------------------------------------------------------
C = {
    "ev":        "#0E8074",   # charging network -- teal
    "ev_dark":   "#075E55",   # charging, dark
    "ev_light":  "#7FB3AC",   # charging, light
    "gas":       "#B3542F",   # fuel network -- rust brown
    "gas_dark":  "#8C2D12",   # fuel, dark
    "risk":      "#D55E00",   # flood / risk -- orange
    "risk_dark": "#A34700",   # risk, dark
    "ink":       "#1A1A1A",   # body text / annotations
    "gray":      "#6E6E6E",   # secondary annotations (arrows, size keys)
    "lgray":     "#B0B0B0",   # light reference lines (diagonals, dumbbells)
    "spine":     "#3C3C3C",   # axis lines and ticks
    "ocean":     "#F2F6F9",   # world map ocean fill
    "land":      "#E7E7E2",   # world map land fill (= legend 'n/a' colour)
}

# Continent colours (single source of truth; identical across Figures 2-4).
# These are the published values: South America yellow, Africa purple,
# Oceania neutral grey (grey keeps the CVD-safe warm/cool axis free for the
# larger continents; Africa never enters a gated panel in the current data).
CONT_C = {"Asia": "#D55E00", "Europe": "#0E8074", "North America": "#2E6FB7",
          "South America": "#E7B33C", "Africa": "#6A51A3", "Oceania": "#8A8A8A"}

# ---------------------------------------------------------------------------
# Robinson viewport (locked values -- do not change)
# ---------------------------------------------------------------------------
XMAX = 17.005e6
XL, XR = -XMAX * (170.0 / 180.0), XMAX
YB, YT = -6.35e6, 7.85e6
MAP_ASPECT = (XR - XL) / (YT - YB)        # ~2.329 (map panel width = H * 2.329)

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
def set_style():
    """Project-wide rcParams. Scripts may override font sizes afterwards."""
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica",
                            "DejaVu Sans", "sans-serif"],
        "font.size": 8.0,
        "axes.labelsize": 8.0,
        "axes.titlesize": 8.4,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.0,
        "text.color": C["ink"],
        "axes.labelcolor": C["ink"],
        "xtick.labelcolor": C["spine"],   # tick labels follow the spine colour
        "ytick.labelcolor": C["spine"],
        "axes.edgecolor": C["spine"],
        "xtick.color": C["spine"],
        "ytick.color": C["spine"],
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "xtick.major.size": 2.6, "ytick.major.size": 2.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 400,
        "figure.dpi": 100,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })

# ---------------------------------------------------------------------------
# Robinson projection helpers
# ---------------------------------------------------------------------------
_TRANS = None

def _transformer():
    global _TRANS
    if _TRANS is None:
        from pyproj import Transformer
        _TRANS = Transformer.from_crs("EPSG:4326", "ESRI:54030", always_xy=True)
    return _TRANS

def to_rob(lon, lat):
    """Longitude/latitude (WGS-84) -> Robinson metric coordinates (x, y)."""
    x, y = _transformer().transform(np.asarray(lon, dtype=float),
                                    np.asarray(lat, dtype=float))
    return np.asarray(x), np.asarray(y)

# ---------------------------------------------------------------------------
# World base map
# ---------------------------------------------------------------------------
_WORLD = None

def world():
    """Natural Earth 110 m country outlines (Robinson CRS); module cache."""
    global _WORLD
    if _WORLD is None:
        import geopandas as gpd
        _WORLD = gpd.read_file(_NE110).to_crs("ESRI:54030")
    return _WORLD

def draw_world(ax, graticule=False):
    """Standard world base: ocean fill + land silhouettes (white outline).

    No title; axes off; x/y limits set to the locked viewport; equal aspect.
    Callers overlay thematic layers with zorder >= 3.
    """
    ax.add_patch(plt.Rectangle((XL, YB), XR - XL, YT - YB,
                               facecolor=C["ocean"], edgecolor="none", zorder=0))
    if graticule:
        phi = np.linspace(-90, 90, 181)
        lam = np.linspace(-180, 180, 361)
        for lo in range(-150, 181, 30):
            gx, gy = to_rob(np.full(181, lo, dtype=float), phi)
            ax.plot(gx, gy, color="white", lw=0.4, zorder=0.5)
        for la in range(-60, 91, 30):
            gx, gy = to_rob(lam, np.full(361, la, dtype=float))
            ax.plot(gx, gy, color="white", lw=0.4, zorder=0.5)
    world().plot(ax=ax, facecolor=C["land"], edgecolor="white",
                 linewidth=0.25, zorder=1)
    ax.set_xlim(XL, XR)
    ax.set_ylim(YB, YT)
    ax.set_aspect("equal")
    ax.axis("off")

# ---------------------------------------------------------------------------
# Panel utilities
# ---------------------------------------------------------------------------
def panel(ax, label, dx=0.0, dy=1.0, fs=13):
    """Panel letter (a/b/c...): axes-fraction coordinates, bold, va='top'."""
    ax.text(dx, dy, label, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", va="top", ha="left", color=C["ink"],
            clip_on=False)

def ptitle(ax, title, pad=3, fs=7.6):
    """Panel title: left-aligned bold (Nature style), fixed 7.6 pt."""
    ax.set_title(title, loc="left", pad=pad, fontweight="bold",
                 fontsize=fs, color=C["ink"])

def despine_all(ax):
    """Remove every spine and tick (for maps and schematic panels)."""
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

def save(fig, name, formats=("png", "pdf"), dpi=400):
    """Draft output to shared/_draft_figures/<name>.{png,pdf}.

    Final figures are written by each module's drawing script into the
    module's own folder; this helper only serves ad-hoc drafts.
    """
    os.makedirs(FIGDIR, exist_ok=True)
    paths = []
    for ext in formats:
        p = os.path.join(FIGDIR, f"{name}.{ext}")
        fig.savefig(p, dpi=dpi, facecolor="white", bbox_inches=None)
        paths.append(p)
    return paths
