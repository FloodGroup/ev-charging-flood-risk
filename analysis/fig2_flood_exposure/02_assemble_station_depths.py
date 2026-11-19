#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 2 -- assemble the station-level depth table.

Joins the seven sampled depth vectors onto the station database and writes
the core intermediate product of the whole study:

  data/station_depths.parquet   2,291,903 stations x (lon, lat, country,
                                src, kind, d010...d500)

Figures 3, 4 and 5 and the SI all reuse this table, so run Figure 2 first
when rebuilding from source data. The script also prints the first real
exposure numbers as an immediate sanity check.
"""
import os, sys, gc
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG2_DATA")
RPS = ["010","020","050","075","100","200","500"]
ev = pd.read_parquet(P.get("DB_CHARGING"))[["lon","lat","country","src"]]
gs = pd.read_parquet(P.get("DB_FUEL"))[["lon","lat","country"]]
n_ev = len(ev)
lon = np.r_[ev.lon.values, gs.lon.values]; lat = np.r_[ev.lat.values, gs.lat.values]
country = pd.Categorical(np.r_[ev.country.values, gs.country.values])
src = pd.Categorical(np.r_[ev.src.values, np.repeat("GP", len(gs))])
kind = pd.Categorical(np.r_[np.repeat("ev", n_ev), np.repeat("gas", len(gs))])
del ev, gs; gc.collect()
cols = {"lon": lon, "lat": lat, "country": pa.array(country),
        "src": pa.array(src), "kind": pa.array(kind)}
for rp in RPS:
    cols[f"d{rp}"] = np.load(f"{DATA}/depth_RP{rp}.npy")
pq.write_table(pa.table(cols), f"{DATA}/station_depths.parquet")
print("SAVED station_depths.parquet")
ev_m = np.asarray(kind) == "ev"
print(f"{'RP':>5} {'EV%':>7} {'GAS%':>7} {'ratio':>6}")
for rp in RPS:
    d = cols[f"d{rp}"]
    e = 100*float((d[ev_m] > 0).mean()); g = 100*float((d[~ev_m] > 0).mean())
    print(f"RP{rp:>3} {e:7.2f} {g:7.2f} {e/g:6.3f}")
d = cols["d100"]
for th in [0.15, 0.5, 1.0]:
    e = 100*float((d[ev_m] >= th).mean()); g = 100*float((d[~ev_m] >= th).mean())
    print(f"RP100 depth>={th}m: EV {e:.2f}% GAS {g:.2f}% ratio {e/g:.3f}")
