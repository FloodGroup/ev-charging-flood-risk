#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate the distances to the two nearest substations using a haversine BallTree query. The nearest substation supplies the 8 km fallback when network assignment fails. Rows follow the charging-station order in the Figure 2 depth table. The runner repeats checkpointed batches until DONE."""
import os, sys, time
import numpy as np, pandas as pd
from sklearn.neighbors import BallTree
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG4_DATA")
CK = P.ensure_dir(os.path.join(DATA, "checkpoints"))
FIG2 = P.get("FIG2_DATA")
R = 6371000.0; BUDGET = 30.0; t0 = time.time()

subs = pd.read_parquet(os.path.join(DATA, "substations_global.parquet"), columns=["sub_id","lon","lat"])
ev = pd.read_parquet(os.path.join(FIG2, "station_depths.parquet"),
                     columns=["lon","lat","country","kind"])
ev = ev[ev.kind=="ev"].reset_index(drop=True)

tree = BallTree(np.deg2rad(np.c_[subs.lat.values, subs.lon.values]), metric="haversine")
P1 = f"{CK}/assign_sub1.npy"; P2 = f"{CK}/assign_sub2.npy"
D1 = f"{CK}/assign_d1.npy"; D2 = f"{CK}/assign_d2.npy"; PROG = f"{CK}/assign_prog.txt"
s1 = np.load(P1) if os.path.exists(P1) else np.full(len(ev), -1, np.int64)
s2 = np.load(P2) if os.path.exists(P2) else np.full(len(ev), -1, np.int64)
d1 = np.load(D1) if os.path.exists(D1) else np.full(len(ev), np.nan, np.float32)
d2 = np.load(D2) if os.path.exists(D2) else np.full(len(ev), np.nan, np.float32)
i = int(open(PROG).read()) if os.path.exists(PROG) else 0
CH = 120000
while i < len(ev) and time.time()-t0 < BUDGET:
    j = min(i+CH, len(ev))
    q = np.deg2rad(np.c_[ev.lat.values[i:j], ev.lon.values[i:j]])
    dist, ind = tree.query(q, k=2)
    d1[i:j] = dist[:,0]*R; d2[i:j] = dist[:,1]*R
    s1[i:j] = ind[:,0]; s2[i:j] = ind[:,1]
    i = j
for p, a in [(P1,s1),(P2,s2),(D1,d1),(D2,d2)]:
    np.save(p+".tmp.npy", a); os.replace(p+".tmp.npy", p)
open(PROG,"w").write(str(i))
if i >= len(ev):
    out = pd.DataFrame(dict(charger_idx=np.arange(len(ev)), country=ev.country,
                            sub1=s1, dist1_m=d1, sub2=s2, dist2_m=d2))
    out.to_parquet(os.path.join(DATA, "charger_assignment.parquet"))
    print("DONE", len(out), "| median dist %.0f m | p95 %.0f m | >8km: %.1f%%" % (
        np.median(d1), np.percentile(d1,95), 100*(d1>8000).mean()), flush=True)
else:
    print(f"progress {i}/{len(ev)} | %.0fs"%(time.time()-t0), flush=True)
