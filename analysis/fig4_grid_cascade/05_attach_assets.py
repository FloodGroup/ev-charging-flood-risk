#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, stage 2 -- attach plants, substations and chargers to the grid.

Every asset is joined to its nearest line-network VERTEX (exact nearest
neighbour on the unit sphere via a 3-D KD-tree). The join produces both the
logical coupling (vertex id) and the physical connector geometry
(asset point -> vertex point), exported later to the GIS model.

Caps: assets farther than ATTACH_MAX from any line vertex stay unattached
(flagged); chargers then fall back to the euclidean nearest substation.

Output: work/graph/attach_subs.npz     vid,int32  dist_m,float32
        work/graph/attach_plants.npz
        work/graph/attach_chargers.npz
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

FIG4_DATA = P.get("FIG4_DATA")
FIG2_DATA = P.get("FIG2_DATA")
GRAPH = P.ensure_dir(os.path.join(FIG4_DATA, "checkpoints", "graph"))
OUT = FIG4_DATA
import time, gc
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

ATTACH_MAX = 10_000.0          # m
R = 6_371_000.0

t0 = time.time()
z = np.load(f"{GRAPH}/vertices.npz")
vlon, vlat = z["lon"], z["lat"]
NV = len(vlon)
print(f"vertices {NV:,}", flush=True)

def to_xyz(lon, lat):
    lo, la = np.deg2rad(lon), np.deg2rad(lat)
    cl = np.cos(la)
    return np.c_[cl * np.cos(lo), cl * np.sin(lo), np.sin(la)]

tree = cKDTree(to_xyz(vlon, vlat))
print(f"kd-tree built {time.time()-t0:.0f}s", flush=True)

def attach(lon, lat, name):
    out_v = np.full(len(lon), -1, np.int32)
    out_d = np.full(len(lon), np.nan, np.float32)
    CH = 200_000
    for i in range(0, len(lon), CH):
        j = min(i + CH, len(lon))
        dd, vv = tree.query(to_xyz(lon[i:j], lat[i:j]), k=1, workers=2)
        arc = 2.0 * R * np.arcsin(np.clip(dd / 2.0, 0, 1))   # chord -> great circle
        out_v[i:j] = vv.astype(np.int32)
        out_d[i:j] = arc.astype(np.float32)
        print(f"  {name} {j:,}/{len(lon):,}  {time.time()-t0:.0f}s", flush=True)
    ok = out_d <= ATTACH_MAX
    out_v[~ok] = -1
    print(f"{name}: attached {ok.mean()*100:.2f}%  median {np.nanmedian(out_d):,.0f} m "
          f"p95 {np.nanpercentile(out_d,95):,.0f} m", flush=True)
    return out_v, out_d

subs = pd.read_parquet(f"{FIG4_DATA}/substations_global.parquet",
                       columns=["sub_id", "lon", "lat"])
v, d = attach(subs.lon.values, subs.lat.values, "substations")
np.savez(f"{GRAPH}/attach_subs.npz", vid=v, dist_m=d)

pl = pd.read_parquet(f"{FIG4_DATA}/plants_global.parquet",
                     columns=["plant_id", "lon", "lat"])
v, d = attach(pl.lon.values, pl.lat.values, "plants")
np.savez(f"{GRAPH}/attach_plants.npz", vid=v, dist_m=d)

sd = pd.read_parquet(f"{FIG2_DATA}/station_depths.parquet",
                     columns=["lon", "lat", "kind"])
ev = sd[sd.kind == "ev"].reset_index(drop=True)
del sd; gc.collect()
v, d = attach(ev.lon.values, ev.lat.values, "chargers")
np.savez(f"{GRAPH}/attach_chargers.npz", vid=v, dist_m=d)
print(f"DONE {time.time()-t0:.0f}s", flush=True)
