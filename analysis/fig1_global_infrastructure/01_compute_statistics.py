#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1, step 1 -- global statistics of the two refuelling networks.

Reads the merged station database (Google Places worldwide + AMap for
mainland China; built by stage 0) and the IEA charging-point series, and
writes three summary products used by the drawing script and the manuscript:

  data/country_stats.csv              per-country station counts and ratios
  data/global_summary.json            global headline numbers
  data/crosscheck_stations_vs_iea.csv station counts audited against IEA
                                      2025 charging-point totals

Note: the 'status' literals below are values carried verbatim in the source
database. They translate as 'open for business' and 'temporarily closed'.
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

OUT = P.ensure_dir(P.get("FIG1_DATA"))
IEA = P.get("IEA_DIR")

ev = pd.read_parquet(P.get("DB_CHARGING"))
gs = pd.read_parquet(P.get("DB_FUEL"))

st_known = ev.status.notna()
g = dict(
    n_ev=int(len(ev)),
    n_ev_gp=int((ev.src == "GP").sum()),
    n_ev_amap=int((ev.src == "AMAP").sum()),
    n_gas=int(len(gs)),
    ev_status_known=int(st_known.sum()),
    ev_open=int((ev.status == "营业中").sum()),
    ev_temp_closed=int((ev.status == "暂时关闭").sum()),
    gas_open=int((gs.status == "营业中").sum()),
    gas_temp_closed=int((gs.status == "暂时关闭").sum()),
    pw_known=int(ev.maxkw.notna().sum()),
)
g["ev_closed_pct"] = round(100 * g["ev_temp_closed"] / g["ev_status_known"], 2)
g["gas_closed_pct"] = round(100 * g["gas_temp_closed"] / g["n_gas"], 2)
g["fast_share_global"] = round(100 * float((ev.maxkw > 22).sum()) / g["pw_known"], 1)
g["pw_coverage_pct"] = round(100 * g["pw_known"] / g["n_ev"], 1)
g["n_countries_ev"] = int(ev[ev.country != ""].country.nunique())

# ---- per-country table -----------------------------------------------------
ce = ev[ev.country != ""].groupby("country").agg(
    ev_n=("country", "size"),
    pw_known=("maxkw", lambda s: s.notna().sum()),
    fast_n=("maxkw", lambda s: (s > 22).sum())).reset_index()
cg = gs[gs.country != ""].groupby("country").agg(
    gas_n=("country", "size")).reset_index()
ct = ce.merge(cg, on="country", how="outer").fillna(0)
for c in ["ev_n", "pw_known", "fast_n", "gas_n"]:
    ct[c] = ct[c].astype(int)
ct["ev_per_gas"] = np.where(ct.gas_n > 0, ct.ev_n / ct.gas_n, np.nan)
ct["fast_share"] = np.where(ct.pw_known > 0, 100 * ct.fast_n / ct.pw_known, np.nan)
ct["pw_cov"] = np.where(ct.ev_n > 0, 100 * ct.pw_known / ct.ev_n, np.nan)
ct.sort_values("ev_n", ascending=False).to_csv(
    os.path.join(OUT, "country_stats.csv"), index=False)

# ---- cross-audit against IEA 2025 charging-point counts ---------------------
iea = pd.read_csv(os.path.join(IEA, "iea_ev_charging_points_2010_2025.csv"))
iea_c = iea[iea.year == 2025].groupby("region")["value"].sum()
rows = []
# South Korea is deliberately absent: Google Places coverage there is a
# small fraction of the real network (domestic providers dominate digital
# mapping), so a points-per-station audit would only measure source
# undercoverage; the manuscript discloses this in Methods and keeps Korea
# out of the six audited markets (SI Fig. 1a).
for our, theirs in [("China", "China"), ("United States of America", "USA"),
                    ("Germany", "Germany"), ("Netherlands", "Netherlands"),
                    ("France", "France"), ("Japan", "Japan")]:
    n_ = int(ct.loc[ct.country == our, "ev_n"].sum())
    p_ = float(iea_c.get(theirs, np.nan))
    if n_ and not np.isnan(p_):
        rows.append(dict(country=our, our_stations=n_, iea_points_2025=int(p_),
                         points_per_station=round(p_ / n_, 2)))
pd.DataFrame(rows).to_csv(
    os.path.join(OUT, "crosscheck_stations_vs_iea.csv"), index=False)
json.dump(g, open(os.path.join(OUT, "global_summary.json"), "w"), indent=1)
print(json.dumps(g, indent=1))
print(pd.DataFrame(rows).to_string())
