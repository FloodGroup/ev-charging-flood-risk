#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3, step 1 -- JRC flood depths at every archive station (resumable).

Tile-grouped window reads with a ~34 s work budget per call; progress is
checkpointed in data/checkpoints so interrupted runs resume exactly where
they stopped. Call repeatedly with a return period (RP010 ... RP500) until
the call prints DONE, then once with 'assemble' to write the joined table.
run_all.py repeats these steps.

  python 01_sample_flood_depths.py RP100     (repeat until DONE)
  python 01_sample_flood_depths.py assemble

Outputs: data/checkpoints/archdepth_RP*.npy, data/archive_depths.parquet
"""
import os, sys, time
import numpy as np, pandas as pd, rasterio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG3_DATA")
CK = P.ensure_dir(os.path.join(DATA, "checkpoints"))
JRC = P.get("JRC_DIR")
RPS = ["RP010","RP020","RP050","RP075","RP100","RP200","RP500"]
BUDGET = 34.0
step = sys.argv[1]; t0 = time.time()
subs = pd.read_parquet(os.path.join(DATA, "archive_stations.parquet"), columns=["lon","lat"])

if step == "assemble":
    df = subs.copy().reset_index(drop=True)
    for rp in RPS:
        df["d"+rp[2:]] = np.load(f"{CK}/archdepth_{rp}.npy")
    df.to_parquet(os.path.join(DATA, "archive_depths.parquet"))
    print("saved", len(df), "| % flooded:",
          {rp: round(float((df["d"+rp[2:]]>0).mean()*100),2) for rp in RPS}, flush=True)
    print("DONE")
    sys.exit()

rp = step
IDX = f"{CK}/archidx.npz"
if not os.path.exists(IDX):
    with rasterio.open(P.jrc_raster("100")) as ds:
        inv = ~ds.transform
        cols, rows = inv * (subs.lon.values, subs.lat.values)
        rows = np.floor(rows).astype(np.int64); cols = np.floor(cols).astype(np.int64)
        ok = (rows>=0)&(rows<ds.height)&(cols>=0)&(cols<ds.width)
    TS = 128
    tid = (rows//TS)*1000000 + (cols//TS)
    idx = np.where(ok)[0][np.argsort(tid[ok], kind="stable")]
    tsort = tid[idx]
    starts = np.r_[0, np.where(np.diff(tsort)!=0)[0]+1]
    tmp = IDX + ".tmp.npz"
    np.savez(tmp, idx=idx, rows=rows, cols=cols, starts=starts,
             tiles=tsort[starts], TS=TS)
    os.replace(tmp, IDX)
    print("index built:", len(starts), "tiles | %.0fs"%(time.time()-t0), flush=True)
z = np.load(IDX)
idx, rows, cols, starts, tiles, TS = z["idx"], z["rows"], z["cols"], z["starts"], z["tiles"], int(z["TS"])
PART = f"{CK}/archdepth_{rp}_part.npy"; PROG = f"{CK}/archdepth_{rp}_prog.txt"
out = np.load(PART) if os.path.exists(PART) else np.full(len(subs), np.nan, np.float32)
k = int(open(PROG).read()) if os.path.exists(PROG) else 0
with rasterio.open(f"{JRC}/{rp}.tif") as ds:
    H, W = ds.height, ds.width
    while k < len(starts) and time.time()-t0 < BUDGET:
        i0 = starts[k]; i1 = starts[k+1] if k+1 < len(starts) else len(idx)
        sel = idx[i0:i1]; t = tiles[k]
        r0 = int(t//1000000)*TS; c0 = int(t%1000000)*TS
        win = rasterio.windows.Window(c0, r0, min(TS, W-c0), min(TS, H-r0))
        arr = ds.read(1, window=win)
        out[sel] = arr[rows[sel]-r0, cols[sel]-c0]
        k += 1
np.save(PART+".tmp.npy", out); os.replace(PART+".tmp.npy", PART)
open(PROG,"w").write(str(k))
if k >= len(starts):
    out[np.isfinite(out)&(out>1e30)] = 0.0; out[out<0] = 0.0
    np.save(f"{CK}/archdepth_{rp}.npy", out)
    print(rp, "DONE | flooded:", int((out>0).sum()), "| %.0fs"%(time.time()-t0), flush=True)
else:
    print(rp, f"progress {k}/{len(starts)} | %.0fs"%(time.time()-t0), flush=True)
