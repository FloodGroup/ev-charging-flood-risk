#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3, step 4 -- extract the city cohort-map data (panels e and f).

For one urban centre, gathers all archive stations, their cohort group
(early <= 2018 / middle / late >= 2022) and exposure flags, plus the local
RP100 inundation mask inside a square window (side length in km) centred
on the flooded late-cohort cluster.

Usage:
  python 04_extract_city_maps.py 2046 28     (Paris,  28-km window)
  python 04_extract_city_maps.py 1677 24     (Malaga, 24-km window)

Outputs: data/citymap_<metro>.npz
"""
import os, sys
import numpy as np, pandas as pd, rasterio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG3_DATA")
JRC = P.get("JRC_DIR")
METRO = int(sys.argv[1]); SIDE = float(sys.argv[2])

st = pd.read_parquet(f"{DATA}/archive_depths.parquet")
meta = pd.read_parquet(f"{DATA}/archive_stations.parquet", columns=["year","country"])
if len(st) != len(meta):
    raise SystemExit(f"archive tables out of sync: depths {len(st):,} rows vs "
                     f"stations {len(meta):,} rows -- re-run stages 0-2 together")
st = pd.concat([st.reset_index(drop=True), meta.reset_index(drop=True)], axis=1)
st["metro"] = np.load(f"{DATA}/archive_metro.npy")
g = st[st.metro == METRO].copy()
g["exp"] = g.d100 > 0
late = g.year >= 2022; early = g.year <= 2018
src_pts = g[late & g.exp] if (late & g.exp).sum() >= 10 else (g[g.exp] if g.exp.sum() >= 5 else g)
cx, cy = float(src_pts.lon.median()), float(src_pts.lat.median())
lat0 = np.deg2rad(cy); kx = np.cos(lat0)*111.32; ky = 110.54
w, e = cx - SIDE/kx/2, cx + SIDE/kx/2
s, n = cy - SIDE/ky/2, cy + SIDE/ky/2
with rasterio.open(P.jrc_raster("100")) as ds:
    win = rasterio.windows.from_bounds(w, s, e, n, ds.transform)
    arr = ds.read(1, window=win)
    arr = np.where(np.isfinite(arr) & (arr < 1e30) & (arr > 0), 1, 0).astype(np.uint8)
grp = np.where(g.year <= 2018, 0, np.where(g.year <= 2021, 1, 2))
np.savez_compressed(f"{DATA}/citymap_{METRO}.npz", bbox=[w, s, e, n],
    lon=g.lon.values, lat=g.lat.values, grp=grp, exp=g.exp.values,
    year=g.year.values, flood=arr,
    stats=[int(early.sum()), float(100*g[early].exp.mean()),
           int(late.sum()), float(100*g[late].exp.mean())])
print(f"metro {METRO}: n={len(g)} early {early.sum()} ({100*g[early].exp.mean():.1f}% exp) "
      f"late {late.sum()} ({100*g[late].exp.mean():.1f}% exp) | flood {arr.shape}", flush=True)
