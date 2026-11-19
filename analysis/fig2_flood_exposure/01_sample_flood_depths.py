#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 1 -- sample the JRC water-depth rasters at every station.

Checkpointed heavy step. Each call samples ONE QUARTER of the tiles for the
next incomplete return period (roughly 8-60 s per call, depending on disk),
writes an atomic part file and exits; the four parts are then merged into
depth_RP{rp}.npy with a shape check. Re-run the script until it prints
ALL RPS DONE -- run_all.py repeats these steps. Interrupting it never
corrupts state.

JRC nodata semantics: a value of 3.4e38 means 'not inundated at this return
period' and is converted to depth 0; any finite value is water depth in
metres.

Outputs: data/depth_RP{010,020,050,075,100,200,500}.npy
"""
import os, sys, time, shutil
os.environ["GDAL_CACHEMAX"] = "128"
import numpy as np, rasterio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402
RDIR = P.get("JRC_DIR")
DATA = P.get("FIG2_DATA")
RPS = ["010","020","050","075","100","200","500"]
BS = 128; NPARTS = 4
idx = np.load(f"{DATA}/sample_index.npz")
order, r_s, c_s = idx["order"], idx["r_s"], idx["c_s"]
starts, ends = idx["starts"], idx["ends"]
N, W, H = int(idx["N"]), int(idx["W"]), int(idx["H"])
NT = len(starts)
cuts = [0, NT//4, NT//2, 3*NT//4, NT]

def ok(f):
    if not os.path.exists(f): return False
    try:
        a = np.load(f, mmap_mode="r"); good = a.shape == (N,); del a; return good
    except Exception: return False

def run_range(rp, k0, k1):
    out = np.full(N, np.nan, dtype=np.float32)
    with rasterio.open(f"{RDIR}/RP{rp}.tif") as r:
        for k in range(k0, k1):
            s, e = starts[k], ends[k]
            tr0 = int(r_s[s])//BS*BS; tc0 = int(c_s[s])//BS*BS
            win = rasterio.windows.Window(tc0, tr0, min(BS, W-tc0), min(BS, H-tr0))
            block = r.read(1, window=win)
            out[order[s:e]] = block[r_s[s:e]-tr0, c_s[s:e]-tc0]
    return out

for rp in RPS:
    fin = f"{DATA}/depth_RP{rp}.npy"
    if ok(fin): continue
    parts = [f"{DATA}/depth_RP{rp}_p{i}.npy" for i in range(NPARTS)]
    t0 = time.time()
    for i in range(NPARTS):
        if ok(parts[i]): continue
        if os.path.exists(parts[i]): os.remove(parts[i])
        tmp = os.path.join(DATA, "_part.tmp.npy")
        np.save(tmp, run_range(rp, cuts[i], cuts[i+1]))
        shutil.move(tmp, parts[i])
        print(f"RP{rp} part {i} done {time.time()-t0:.0f}s"); sys.exit(0)
    out = np.load(parts[0])
    for i in range(1, NPARTS):
        b = np.load(parts[i]); out = np.where(np.isnan(out), b, out)
    out[out >= 1e30] = 0.0
    tmp = os.path.join(DATA, "_fin.tmp.npy"); np.save(tmp, out); shutil.move(tmp, fin)
    if not ok(fin): print("MERGE WRITE FAILED", rp); sys.exit(1)
    for p in parts: os.remove(p)
    print(f"RP{rp} MERGED | flooded>0m: {int((out>0).sum())}"); sys.exit(0)
print("ALL RPS DONE")
