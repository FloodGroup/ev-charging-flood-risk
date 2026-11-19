#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 3 -- assign every station to a GHS-UCDB urban centre.

Checkpointed: each call processes one quarter of the stations against the
13,135 GHS-UCDB R2019A urban-centre polygons (STRtree point-in-polygon) and
writes metro_p{i}.npy (int32 UCDB identifier; -1 = outside any urban
centre). Re-run until it prints ALL METRO PARTS DONE (at most five calls);
run_all.py repeats these steps.

Outputs: data/metro_p{0..3}.npy
"""
import os, sys, time
import numpy as np, pandas as pd
from shapely.strtree import STRtree
import shapely, pyogrio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG2_DATA")
UC = P.get("UCDB_GPKG")
NP_ = 4
sd = pd.read_parquet(f"{DATA}/station_depths.parquet", columns=["lon","lat"])
N = len(sd); cuts = [0, N//4, N//2, 3*N//4, N]
i = next((j for j in range(NP_) if not os.path.exists(f"{DATA}/metro_p{j}.npy")), None)
if i is None:
    print("ALL METRO PARTS DONE"); sys.exit(0)
t0 = time.time()
uc = pyogrio.read_dataframe(UC, columns=["ID_HDC_G0"])
tree = STRtree(uc.geometry.values)
ids = uc["ID_HDC_G0"].values.astype(np.int32)
print(f"UCDB loaded: {len(uc)} urban centres {time.time()-t0:.0f}s", flush=True)
s, e = cuts[i], cuts[i+1]
pts = shapely.points(sd.lon.values[s:e], sd.lat.values[s:e])
qi, ti = tree.query(pts, predicate="intersects")
out = np.full(e - s, -1, dtype=np.int32)
out[qi] = ids[ti]
tmp = DATA + "/_metro.tmp.npy"; np.save(tmp, out)
import shutil; shutil.move(tmp, f"{DATA}/metro_p{i}.npy")
print(f"part {i}: {e-s} pts, in-metro {(out>=0).sum()} ({time.time()-t0:.0f}s)")
