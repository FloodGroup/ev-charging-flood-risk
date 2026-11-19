#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, step 2 -- assign a country and continent to every substation
(resumable).

Chunked STRtree point-in-polygon queries against Natural Earth 10-m
boundaries with a ~32 s work budget per call; re-run until the call prints
DONE. run_all.py repeats these steps.

Output: data/substation_country.parquet
"""
import os, sys, time
import numpy as np, pandas as pd, pyogrio
from shapely.strtree import STRtree
import shapely
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG4_DATA")
CK = P.ensure_dir(os.path.join(DATA, "checkpoints"))
BUDGET = 32.0; t0 = time.time()

subs = pd.read_parquet(os.path.join(DATA, "substations_global.parquet"), columns=["sub_id","lon","lat"])
ne = pyogrio.read_dataframe(P.get("NE10M"),
                            columns=["NAME","CONTINENT"])
tree = STRtree(ne.geometry.values)

PART = f"{CK}/sub_country_part.npy"; PROG = f"{CK}/sub_country_prog.txt"
codes = np.load(PART) if os.path.exists(PART) else np.full(len(subs), -1, np.int32)
i = int(open(PROG).read()) if os.path.exists(PROG) else 0
CH = 60000
while i < len(subs) and time.time()-t0 < BUDGET:
    j = min(i+CH, len(subs))
    pts = shapely.points(subs.lon.values[i:j], subs.lat.values[i:j])
    q = tree.query(pts, predicate="intersects")     # (2,N) [pt_idx, geom_idx]
    codes[i:j] = -1
    codes[i + q[0]] = q[1]                          # last match wins (rare overlaps)
    i = j
np.save(PART+".tmp.npy", codes); os.replace(PART+".tmp.npy", PART)
open(PROG,"w").write(str(i))
if i >= len(subs):
    out = subs.copy()
    out["country"] = np.where(codes>=0, ne.NAME.values[np.maximum(codes,0)], "")
    out["continent"] = np.where(codes>=0, ne.CONTINENT.values[np.maximum(codes,0)], "")
    out[["sub_id","country","continent"]].to_parquet(os.path.join(DATA, "substation_country.parquet"))
    print("DONE", len(out), "| unmatched:", int((codes<0).sum()),
          "| top:", out.country.value_counts().head(5).to_dict(), flush=True)
else:
    print(f"progress {i}/{len(subs)} | %.0fs"%(time.time()-t0), flush=True)
