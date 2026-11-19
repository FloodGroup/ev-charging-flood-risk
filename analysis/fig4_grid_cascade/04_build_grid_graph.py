#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, stage 1 -- parse the global power-line network into a graph.

Input : gridfinder_mv_lines.gpkg  layer 'grid' (the official gridfinder
        distribution: 3,678,243 polylines, WGS-84, columns fid/geom/source;
        'source' marks each line as gridfinder-predicted or openstreetmap)
Output: work/graph/vertices.npz    lon,lat (float64), component (int32), degree (uint8)
        work/graph/segments.npz    ia,ib (int32), length_m (float32), fid (uint32)
        work/graph/features.parquet fid, source, n_vert, length_m, comp

Vertices are snapped by rounding to 1e-6 deg (~0.11 m) -- gridfinder endpoints
coincide exactly, so this is loss-free while guarding against float noise.
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
import numpy as np, pandas as pd, pyogrio, shapely

GPKG   = P.get("GRIDFINDER")
LAYER  = "grid"
SLICE  = 250_000
R      = 6_371_000.0

t0 = time.time()
info = pyogrio.read_info(GPKG, layer=LAYER)
NF = int(info["features"])
print(f"features: {NF:,}", flush=True)

# ---------- pass 1: coordinates per slice ----------
ck = os.path.join(GRAPH, "slices")
os.makedirs(ck, exist_ok=True)
n_slices = (NF + SLICE - 1) // SLICE
for si in range(n_slices):
    f = f"{ck}/s{si:03d}.npz"
    if os.path.exists(f):
        continue
    g = pyogrio.read_dataframe(GPKG, layer=LAYER, skip_features=si * SLICE,
                               max_features=SLICE, columns=["source"])
    geo = g.geometry.values
    # multilinestrings (if any) -> parts as separate chains
    if (shapely.get_type_id(geo) != 1).any():
        geo = shapely.line_merge(geo)          # best effort; parts handled below
    coords = shapely.get_coordinates(geo).astype(np.float64)
    ncoor = shapely.get_num_coordinates(geo).astype(np.int32)
    src = (g["source"].astype(str) == "gridfinder").values.astype(np.uint8)
    np.savez_compressed(f + ".tmp.npz", coords=coords, ncoor=ncoor, src=src)
    os.replace(f + ".tmp.npz", f)
    print(f"slice {si+1}/{n_slices}  coords={len(coords):,}  {time.time()-t0:.0f}s",
          flush=True)
    del g, geo, coords, ncoor
    gc.collect()

# ---------- pass 2: global vertex table + segments ----------
print("assembling…", flush=True)
all_coords, all_ncoor, all_src = [], [], []
for si in range(n_slices):
    z = np.load(f"{ck}/s{si:03d}.npz")
    all_coords.append(z["coords"]); all_ncoor.append(z["ncoor"]); all_src.append(z["src"])
coords = np.concatenate(all_coords); del all_coords
ncoor = np.concatenate(all_ncoor);   del all_ncoor
src   = np.concatenate(all_src);     del all_src
gc.collect()
NC = len(coords)
print(f"total coords {NC:,}  ({time.time()-t0:.0f}s)", flush=True)

# integer snap keys (1e-6 deg)
ilon = np.round((coords[:, 0] + 180.0) * 1e6).astype(np.int64)   # 0 .. 360e6
ilat = np.round((coords[:, 1] +  90.0) * 1e6).astype(np.int64)   # 0 .. 180e6
key = ilon * 200_000_007 + ilat
del ilon, ilat
gc.collect()
uniq, inv = np.unique(key, return_inverse=True)
inv = inv.astype(np.int32)
NV = len(uniq)
del key, uniq
gc.collect()
print(f"vertices {NV:,}  ({time.time()-t0:.0f}s)", flush=True)

vx = np.zeros(NV, np.float64); vy = np.zeros(NV, np.float64)
vx[inv] = coords[:, 0]; vy[inv] = coords[:, 1]

# segments: consecutive coordinate pairs inside each feature
idx = np.r_[0, np.cumsum(ncoor)]
keep = np.ones(NC, bool)
keep[idx[1:] - 1] = False            # last coord of each feature starts no segment
a = np.where(keep)[0]
ia = inv[a]; ib = inv[a + 1]
fid_of_seg = np.repeat(np.arange(len(ncoor), dtype=np.uint32), ncoor - 1)
# haversine length per segment
lon1, lat1 = np.deg2rad(coords[a, 0]), np.deg2rad(coords[a, 1])
lon2, lat2 = np.deg2rad(coords[a + 1, 0]), np.deg2rad(coords[a + 1, 1])
dlat = lat2 - lat1; dlon = lon2 - lon1
h = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
seg_len = (2 * R * np.arcsin(np.sqrt(h))).astype(np.float32)
del coords, lon1, lat1, lon2, lat2, dlat, dlon, h, keep, a, inv
gc.collect()

# drop zero-length segments (duplicate consecutive vertices)
nz = (ia != ib)
ia, ib, seg_len, fid_of_seg = ia[nz], ib[nz], seg_len[nz], fid_of_seg[nz]
NS = len(ia)
print(f"segments {NS:,}  total length {seg_len.sum()/1e9*1000:,.0f} thousand km  "
      f"({time.time()-t0:.0f}s)", flush=True)

# degree + connected components
deg = np.bincount(ia, minlength=NV) + np.bincount(ib, minlength=NV)
deg = np.clip(deg, 0, 255).astype(np.uint8)
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
g = coo_matrix((np.ones(NS, np.int8), (ia, ib)), shape=(NV, NV))
ncomp, comp = connected_components(g, directed=False)
comp = comp.astype(np.int32)
del g
gc.collect()
sizes = np.bincount(comp)
print(f"components {ncomp:,}  largest {sizes.max():,} vertices  "
      f"({time.time()-t0:.0f}s)", flush=True)

np.savez(os.path.join(GRAPH, "vertices.npz"), lon=vx, lat=vy, comp=comp, deg=deg)
np.savez(os.path.join(GRAPH, "segments.npz"), ia=ia, ib=ib,
         length_m=seg_len, fid=fid_of_seg)

# per-feature table
first_seg = np.r_[0, np.cumsum((ncoor - 1))][:-1]      # index of feature's first segment (pre-drop)
feat_len = np.zeros(len(ncoor), np.float64)
np.add.at(feat_len, fid_of_seg, seg_len.astype(np.float64))
feat_comp = np.zeros(len(ncoor), np.int32)
feat_comp[fid_of_seg] = comp[ia]                        # any vertex of the feature
pd.DataFrame(dict(fid=np.arange(len(ncoor), dtype=np.int64),
                  source=np.where(src == 1, "gridfinder", "osm"),
                  n_vert=ncoor, length_m=feat_len.astype(np.float32),
                  comp=feat_comp)).to_parquet(os.path.join(GRAPH, "features.parquet"))
print(f"DONE {time.time()-t0:.0f}s", flush=True)
