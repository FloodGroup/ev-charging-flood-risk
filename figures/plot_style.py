import os
import shutil
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_NE110 = os.path.join(_HERE, "world_countries_110m.geojson")
C = {
    "ev": "#0E8074",
    "ev_dark": "#075E55",
    "ev_light": "#7FB3AC",
    "gas": "#B3542F",
    "gas_dark": "#8C2D12",
    "risk": "#D55E00",
    "risk_dark": "#A34700",
    "ink": "#1A1A1A",
    "gray": "#6E6E6E",
    "lgray": "#B0B0B0",
    "spine": "#3C3C3C",
    "ocean": "#F2F6F9",
    "land": "#E7E7E2",
}
CONT_C = {
    "Asia": "#D55E00",
    "Europe": "#0E8074",
    "North America": "#2E6FB7",
    "South America": "#E7B33C",
    "Africa": "#6A51A3",
    "Oceania": "#8A8A8A",
}
XMAX = 17005000.0
XL, XR = (-XMAX * (170.0 / 180.0), XMAX)
YB, YT = (-6350000.0, 7850000.0)
MAP_ASPECT = (XR - XL) / (YT - YB)


def set_style():
    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.4,
            "xtick.labelsize": 7.4,
            "ytick.labelsize": 7.4,
            "legend.fontsize": 7.0,
            "text.color": C["ink"],
            "axes.labelcolor": C["ink"],
            "xtick.labelcolor": C["spine"],
            "ytick.labelcolor": C["spine"],
            "axes.edgecolor": C["spine"],
            "xtick.color": C["spine"],
            "ytick.color": C["spine"],
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.5,
            "ytick.minor.width": 0.5,
            "xtick.major.size": 2.6,
            "ytick.major.size": 2.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 600,
            "figure.dpi": 100,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


_TRANS = None


def to_rob(lon, lat):
    global _TRANS
    if _TRANS is None:
        from pyproj import Transformer

        _TRANS = Transformer.from_crs("EPSG:4326", "ESRI:54030", always_xy=True)
    x, y = _TRANS.transform(np.asarray(lon, dtype=float), np.asarray(lat, dtype=float))
    return (np.asarray(x), np.asarray(y))


_WORLD = None


def world():
    global _WORLD
    if _WORLD is None:
        import geopandas as gpd

        _WORLD = gpd.read_file(_NE110).to_crs("ESRI:54030")
    return _WORLD


def draw_world(ax, graticule=False):
    ax.add_patch(plt.Rectangle((XL, YB), XR - XL, YT - YB, facecolor=C["ocean"], edgecolor="none", zorder=0))
    if graticule:
        phi = np.linspace(-90, 90, 181)
        lam = np.linspace(-180, 180, 361)
        for lo in range(-150, 181, 30):
            gx, gy = to_rob(np.full(181, lo, dtype=float), phi)
            ax.plot(gx, gy, color="white", lw=0.4, zorder=0.5)
        for la in range(-60, 91, 30):
            gx, gy = to_rob(lam, np.full(361, la, dtype=float))
            ax.plot(gx, gy, color="white", lw=0.4, zorder=0.5)
    world().plot(ax=ax, facecolor=C["land"], edgecolor="white", linewidth=0.25, zorder=1)
    ax.set_xlim(XL, XR)
    ax.set_ylim(YB, YT)
    ax.set_aspect("equal")
    ax.axis("off")


def panel(ax, label, dx=0.0, dy=1.0, fs=13):
    ax.text(
        dx,
        dy,
        label,
        transform=ax.transAxes,
        fontsize=fs,
        fontweight="bold",
        va="top",
        ha="left",
        color=C["ink"],
        clip_on=False,
    )


def ptitle(ax, title, pad=3, fs=7.6):
    ax.set_title(title, loc="left", pad=pad, fontweight="bold", fontsize=fs, color=C["ink"])


def despine_all(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def export(fig, base, tight=False, dpi=600):
    kw = dict(facecolor="white")
    if tight:
        kw["bbox_inches"] = "tight"
    fig.savefig(base + ".pdf", **kw)
    fig.savefig(base + ".svg", **kw)
    fig.savefig(base + ".png", dpi=dpi, **kw)
    shutil.copyfile(base + ".pdf", base + ".ai")
    plt.close(fig)
