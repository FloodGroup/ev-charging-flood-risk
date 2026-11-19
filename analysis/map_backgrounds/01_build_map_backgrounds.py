"""
Build the cartographic backgrounds of the city maps (Fig. 3e,f and Fig. 4a-d).

The backgrounds only help the reader to orient; they enter no calculation.
For every window listed in analysis/map_backgrounds/map_windows.csv the script writes

  relief_<id>.npz   multi-directional hillshade of SRTM 1 arc-second elevation
                    (shade: uint8 image; extent: west, east, south, north)
  basemap_<id>.npz  main roads and water bodies from OpenStreetMap, and a sea
                    mask (road_xy / road_off / road_class, water_xy /
                    water_off, sea / sea_extent). The sea mask is either the
                    area at or below sea level in the elevation model
                    (sea = elevation) or the seaward side of the OpenStreetMap
                    coastline (sea = coastline; used where steep shores make
                    the elevation rule inaccurate)

into <figures_directory>/<figure>/data/private/.

Inputs are downloaded on first use and cached under source_data/
map_backgrounds_srtm_osm/ (override with EVFR_MAP_BG_CACHE):

  SRTM elevation tiles (public domain), Terrain Tiles open dataset
      https://registry.opendata.aws/terrain-tiles
  OpenStreetMap roads and water bodies (ODbL 1.0), through an Overpass API
      endpoint. The default endpoint is https://overpass-api.de/api/interpreter;
      set OVERPASS_URL to use another one.

Usage:
  python analysis/map_backgrounds/01_build_map_backgrounds.py <figures_directory> [id ...]
"""
import csv
import gzip
import json
import os
import shutil
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
from rasterio.features import rasterize
from rasterio.transform import Affine
from scipy import ndimage
from shapely.geometry import LineString, Polygon, box
from shapely.geometry.polygon import orient
from shapely.ops import linemerge, polygonize, unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "shared"))
import paths as P  # noqa: E402

CACHE = P.ensure_dir(P.get("MAP_BG_CACHE"))
WINDOWS = os.path.join(HERE, "map_windows.csv")
SRTM_URL = "https://s3.amazonaws.com/elevation-tiles-prod/skadi/{band}/{name}.hgt.gz"
OVERPASS_URL = os.environ.get("OVERPASS_URL", "https://overpass-api.de/api/interpreter")
ROAD_RANK = {"motorway": 0, "trunk": 1, "primary": 2, "secondary": 3, "tertiary": 4}
ROAD_SETS = {"main": ["motorway", "trunk", "primary"],
             "local": ["motorway", "trunk", "primary", "secondary", "tertiary"]}
EDGE = 12  # cells trimmed from each side after smoothing


def download(url, dst, data=None, timeout=180, tries=5):
    """Fetch *url* into *dst* (atomic write, a few retries)."""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "ev-charging-flood-risk/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(dst + ".part", "wb") as f:
                shutil.copyfileobj(r, f)
            os.replace(dst + ".part", dst)
            return dst
        except Exception as ex:  # network errors are retried
            print(f"  retry {k + 1}/{tries}: {str(ex)[:80]}", flush=True)
            time.sleep(6 + 6 * k)
    raise RuntimeError(f"could not download {url}")


# ---------------------------------------------------------------------------
# Elevation and shaded relief
# ---------------------------------------------------------------------------
def srtm_tile(lat, lon):
    """One-degree SRTM tile (3601 x 3601 samples, 1 arc-second) as float32."""
    name = f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}{'E' if lon >= 0 else 'W'}{abs(lon):03d}"
    f = os.path.join(CACHE, f"{name}.hgt.gz")
    if not os.path.exists(f) or os.path.getsize(f) < 1000:
        download(SRTM_URL.format(band=name[:3], name=name), f)
    a = np.frombuffer(gzip.open(f).read(), dtype=">i2").reshape(3601, 3601).astype(np.float32)
    a[a < -1000] = 0  # voids
    return a


def dem_window(w, s, e, n, step):
    """Elevation resampled bilinearly to a regular grid of *step* degrees."""
    lons = np.arange(w, e, step)
    lats = np.arange(n, s, -step)
    out = np.zeros((len(lats), len(lons)), dtype=np.float32)
    for la in range(int(np.floor(s)), int(np.floor(n)) + 1):
        for lo in range(int(np.floor(w)), int(np.floor(e)) + 1):
            t = srtm_tile(la, lo)
            mi = (lats >= la) & (lats <= la + 1)
            mj = (lons >= lo) & (lons <= lo + 1)
            if not mi.any() or not mj.any():
                continue
            rr = (la + 1 - lats[mi]) * 3600.0
            cc = (lons[mj] - lo) * 3600.0
            R, C = np.meshgrid(rr, cc, indexing="ij")
            out[np.ix_(mi, mj)] = ndimage.map_coordinates(t, [R, C], order=1, mode="nearest")
    return out, lons, lats


def hillshade(z, dx, dy, ve, azimuths=(315, 270, 225, 360), weights=(0.4, 0.2, 0.2, 0.2), altitude=45):
    """Multi-directional hillshade (0 to 1) with vertical exaggeration *ve*."""
    gy, gx = np.gradient(z * ve, dy, dx)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(gy, -gx)
    zen = np.deg2rad(90 - altitude)
    hs = 0
    for az, wt in zip(azimuths, weights):
        a = np.deg2rad(360 - az + 90)
        hs = hs + wt * (np.cos(zen) * np.cos(slope) + np.sin(zen) * np.sin(slope) * np.cos(a - aspect))
    return np.clip(hs, 0, 1)


def build_relief(w, s, e, n, step_arcsec, out):
    step = step_arcsec / 3600.0
    pad = EDGE * step
    z, lons, lats = dem_window(w - pad, s - pad, e + pad, n + pad, step)
    lat0 = (s + n) / 2
    dx = step * 111320 * np.cos(np.deg2rad(lat0))
    dy = step * 110540
    zs = ndimage.gaussian_filter(z, max(1.0, 45.0 / dx))  # about 45 m smoothing
    relief_range = np.percentile(zs, 99) - np.percentile(zs, 1)
    ve = float(np.clip(260.0 / max(relief_range, 1), 1.6, 9.0))  # stronger exaggeration on flat terrain
    hs = hillshade(zs, dx, dy, ve)
    v = 0.78 + (hs - np.sin(np.deg2rad(45))) * 1.15  # flat ground maps to a light tone
    v = np.clip(v, 0.42, 1.0)[EDGE:-EDGE, EDGE:-EDGE]
    np.savez_compressed(out, shade=np.round(v * 255).astype(np.uint8),
                        extent=np.array([lons[EDGE], lons[-EDGE - 1], lats[-EDGE - 1], lats[EDGE]]))
    print(f"  {os.path.basename(out)} {v.shape} relief {relief_range:.0f} m, exaggeration {ve:.2f}", flush=True)


# ---------------------------------------------------------------------------
# OpenStreetMap roads and water bodies
# ---------------------------------------------------------------------------
def overpass(query, dst):
    if os.path.exists(dst) and os.path.getsize(dst) > 200:
        try:
            return json.load(open(dst, encoding="utf8"))
        except Exception:
            pass
    download(OVERPASS_URL, dst, data=urllib.parse.urlencode({"data": query}).encode())
    return json.load(open(dst, encoding="utf8"))


def fetch_osm(wid, w, s, e, n, detail, water_relations):
    bb = f"({s:.5f},{w:.5f},{n:.5f},{e:.5f})"
    roads = "".join(f'way["highway"="{k}"]{bb};' for k in ROAD_SETS[detail])
    water = f'way["natural"="water"]{bb};way["waterway"="riverbank"]{bb};'
    if water_relations:
        water += f'relation["natural"="water"]{bb};'
    rd = overpass(f"[out:json][timeout:150];({roads});out geom;", os.path.join(CACHE, f"osm_{wid}_roads.json"))
    wt = overpass(f"[out:json][timeout:150];({water});out geom;", os.path.join(CACHE, f"osm_{wid}_water.json"))
    return rd["elements"], wt["elements"]


def sea_from_elevation(z):
    """Large connected areas at or below sea level."""
    sea = ndimage.binary_opening(z <= 0.0, iterations=3)
    lab, nl = ndimage.label(sea)
    keep = np.zeros_like(sea)
    if nl:
        sizes = ndimage.sum(sea, lab, range(1, nl + 1))
        for i, sz in enumerate(sizes, start=1):
            if sz > 0.004 * sea.size:
                keep |= lab == i
    return keep


def sea_from_coastline(wid, w, s, e, n, z, lons, lats):
    """Seaward side of the OpenStreetMap coastline, on the elevation grid.

    The window is split along the coastline; a piece counts as sea when most
    of its elevation samples are at or below sea level."""
    bb = f"({s:.5f},{w:.5f},{n:.5f},{e:.5f})"
    coast = overpass(f'[out:json][timeout:150];(way["natural"="coastline"]{bb};);out geom;',
                     os.path.join(CACHE, f"osm_{wid}_coast.json"))["elements"]
    lines = [LineString([(p["lon"], p["lat"]) for p in el["geometry"]])
             for el in coast if el.get("geometry") and len(el["geometry"]) >= 2]
    if not lines:
        return sea_from_elevation(z)
    window = box(w, s, e, n)
    step = float(lons[1] - lons[0])
    transform = Affine(step, 0, float(lons[0]) - step / 2, 0, -step, float(lats[0]) + step / 2)
    sea = np.zeros(z.shape, dtype=bool)
    for face in polygonize(unary_union(lines + [window.boundary])):
        if not window.contains(face.representative_point()):
            continue
        m = rasterize([(face, 1)], out_shape=z.shape, transform=transform, fill=0, dtype="uint8").astype(bool)
        if m.any() and (z[m] <= 0.0).mean() > 0.5:
            sea |= m
    return sea


def build_basemap(wid, w, s, e, n, detail, water_relations, sea_rule, out):
    roads, water = fetch_osm(wid, w, s, e, n, detail, water_relations)
    window = box(w, s, e, n)
    tol = (e - w) / 1500.0

    lines = {k: [] for k in range(5)}
    for el in roads:
        g = el.get("geometry")
        k = ROAD_RANK.get(el.get("tags", {}).get("highway"))
        if not g or k is None or len(g) < 2:
            continue
        ln = LineString([(p["lon"], p["lat"]) for p in g]).intersection(window)
        for part in getattr(ln, "geoms", [ln]):
            if part.geom_type == "LineString" and not part.is_empty:
                lines[k].append(part.simplify(tol))
    xy, off, cls = [], [0], []
    for k in range(5):
        if not lines[k]:
            continue
        merged = unary_union(lines[k])
        merged = linemerge(merged) if merged.geom_type != "LineString" else merged
        for part in getattr(merged, "geoms", [merged]):
            c = np.asarray(part.coords)
            if len(c) < 2:
                continue
            xy.append(c)
            off.append(off[-1] + len(c))
            cls.append(k)

    polys = []
    for el in water:
        if el["type"] == "way" and el.get("geometry"):
            c = [(p["lon"], p["lat"]) for p in el["geometry"]]
            if len(c) >= 4 and c[0] == c[-1]:
                polys.append(Polygon(c))
        elif el["type"] == "relation":
            outer, inner = [], []
            for m in el.get("members", []):
                g = m.get("geometry")
                if not g or len(g) < 2:
                    continue
                (inner if m.get("role") == "inner" else outer).append(LineString([(p["lon"], p["lat"]) for p in g]))
            if outer:
                po = unary_union(list(polygonize(unary_union(outer))))
                if inner:
                    po = po.difference(unary_union(list(polygonize(unary_union(inner)))))
                polys.append(po)
    wxy, woff = [], [0]
    if polys:
        u = unary_union([p.buffer(0) for p in polys if p.is_valid or p.buffer(0).is_valid]).intersection(window).simplify(tol)
        min_area = ((e - w) * (n - s)) * 2e-5
        for part in getattr(u, "geoms", [u]):
            if part.geom_type != "Polygon" or part.area < min_area:
                continue
            part = orient(part, 1.0)  # exterior counter-clockwise, holes clockwise
            for ring in [part.exterior] + [r for r in part.interiors if Polygon(r).area >= min_area]:
                c = np.asarray(ring.coords)
                wxy.append(c)
                woff.append(woff[-1] + len(c))

    z, lons, lats = dem_window(w, s, e, n, 1.0 / 3600)
    keep = sea_from_coastline(wid, w, s, e, n, z, lons, lats) if sea_rule == "coastline" else sea_from_elevation(z)

    np.savez_compressed(
        out,
        road_xy=np.round(np.vstack(xy), 6) if xy else np.zeros((0, 2)), road_off=np.array(off), road_class=np.array(cls),
        water_xy=np.round(np.vstack(wxy), 6) if wxy else np.zeros((0, 2)), water_off=np.array(woff),
        sea=keep.astype(np.uint8), sea_extent=np.array([lons[0], lons[-1], lats[-1], lats[0]]))
    print(f"  {os.path.basename(out)} roads {len(cls)}, water rings {len(woff) - 1}, sea share {keep.mean():.3f}", flush=True)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    figures, only = sys.argv[1], set(sys.argv[2:])
    for row in csv.DictReader(open(WINDOWS, encoding="utf8")):
        if only and row["id"] not in only:
            continue
        w, s, e, n = (float(row[k]) for k in ("west", "south", "east", "north"))
        d = P.ensure_dir(os.path.join(figures, row["figure"], "data", "private"))
        print(f"{row['city']} ({row['id']})", flush=True)
        build_relief(w, s, e, n, float(row["dem_step_arcsec"]), os.path.join(d, f"relief_{row['id']}.npz"))
        build_basemap(row["id"], w, s, e, n, row["road_detail"], row["water_relations"] == "1", row["sea"],
                      os.path.join(d, f"basemap_{row['id']}.npz"))


if __name__ == "__main__":
    main()
