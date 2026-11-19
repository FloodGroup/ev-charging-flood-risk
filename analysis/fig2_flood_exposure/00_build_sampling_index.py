#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 0 -- build and cache the tile-sorted sampling index.

Converts the coordinates of all 2,291,903 stations (947,602 charging +
1,344,301 fuel) to row/column indices of the JRC 90-m grid, groups them by
128 x 128 raster tile and stores the ordering. Sampling the seven return-
period rasters then reads each tile window once instead of seeking per
station, which is roughly sixty times faster. Run once; the index is shared
by every sampling pass.

Output: data/sample_index.npz
"""
import numpy as np, pandas as pd, rasterio
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

RDIR = P.get("JRC_DIR")
DATA = P.ensure_dir(P.get("FIG2_DATA"))
BS = 128
ev = pd.read_parquet(P.get("DB_CHARGING"))[["lon","lat"]]
gs = pd.read_parquet(P.get("DB_FUEL"))[["lon","lat"]]
n_ev = len(ev)
lon = np.r_[ev.lon.values, gs.lon.values]; lat = np.r_[ev.lat.values, gs.lat.values]
with rasterio.open(P.jrc_raster("100")) as r0:
    inv = ~r0.transform; W, H = r0.width, r0.height
cc_, rr_ = inv * (lon, lat)
cols = np.floor(cc_).astype(np.int32); rows = np.floor(rr_).astype(np.int32)
inside = (cols>=0)&(cols<W)&(rows>=0)&(rows<H)
tile = (rows.astype(np.int64)//BS)*10**7 + (cols//BS)
order = np.argsort(np.where(inside, tile, -1)).astype(np.int32)
r_s = rows[order]; c_s = cols[order]
t_s = np.where(inside, tile, -1)[order]
n0 = int((~inside).sum())
b_ = np.flatnonzero(np.diff(t_s[n0:])) + 1
starts = (np.r_[0, b_] + n0).astype(np.int64)
ends = (np.r_[b_, len(t_s)-n0] + n0).astype(np.int64)
np.savez_compressed(f"{DATA}/sample_index.npz",
                    order=order, r_s=r_s, c_s=c_s, starts=starts, ends=ends,
                    n_ev=n_ev, N=len(lon), W=W, H=H, n_outside=n0)
print("index cached | N:", len(lon), "| tiles:", len(starts), "| outside:", n0)
