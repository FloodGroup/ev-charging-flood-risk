#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3, step 0 -- ingest the Open Charge Map registry archive.

Cleans the Kaggle-distributed OCM export (257,585 records), keeps dated,
geolocated stations registered 2010-2025 (dropping records whose status is
'Planned For Future Date' or 'Not Operational'; other statuses, including
unknown, are retained), maps country names onto Natural Earth conventions, and computes the per-country coverage-gate
table used to select the 15 well-covered markets:

  gates: n >= 1,200 dated stations; >= 6 years with meaningful additions;
         no single registration year holding > 55% of entries.

The recorded year is the database-registration year, a proxy for the
opening year; bulk imports mix station ages within a year, which biases
estimated trends towards zero, so flat results are conservative.

Outputs: data/archive_stations.parquet, data/coverage_table.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.ensure_dir(P.get("FIG3_DATA"))
SRC = P.get("OCM_ARCHIVE")
FIG2 = P.get("FIG2_DATA")

d = pd.read_csv(SRC, usecols=["Country","Latitude","Longitude","YearCreated",
                              "StatusType","Operator","MaxPowerKW"],
                low_memory=False, encoding="latin-1")
for c in ["Latitude","Longitude","YearCreated"]:
    d[c] = pd.to_numeric(d[c], errors="coerce")
n0 = len(d)
d = d[d.YearCreated.notna() & d.Latitude.between(-60, 80) & d.Longitude.between(-180, 180)]
d = d[~d.StatusType.astype(str).isin(["Planned For Future Date", "Not Operational"])]
d["year"] = d.YearCreated.astype(int)
d = d[(d.year >= 2010) & (d.year <= 2025)].copy()

# the Kaggle export writes "," as ";" inside country names -- map both forms
NE_MAP = {"United States": "United States of America", "Russian Federation": "Russia",
          "Czech Republic": "Czechia", "Korea, Republic of": "South Korea",
          "Korea; Republic Of": "South Korea", "Hong Kong": "Hong Kong S.A.R.",
          "Taiwan, Province of China": "Taiwan", "Taiwan; Province Of China": "Taiwan"}
d["country"] = d.Country.astype(str).map(lambda c: NE_MAP.get(c, c))
d = d.rename(columns={"Latitude": "lat", "Longitude": "lon"})
d[["lon","lat","year","country","Operator","MaxPowerKW"]].reset_index(drop=True) \
 .to_parquet(os.path.join(DATA, "archive_stations.parquet"))
print(f"usable {len(d):,} / {n0:,} raw")

# ---------- per-country QA ----------
ev_ct = pd.read_parquet(os.path.join(FIG2, "station_depths.parquet"), columns=["kind","country"])
merged_n = ev_ct[ev_ct.kind == "ev"].country.value_counts()
rows = []
for c, g in d.groupby("country"):
    s = g.year.value_counts(normalize=True)
    yrs_meaningful = int((g.year.value_counts() >= max(30, 0.01*len(g))).sum())
    rows.append(dict(country=c, n=len(g), max_year_share=round(float(s.max()), 3),
                     peak_year=int(s.idxmax()), spread_years=yrs_meaningful,
                     coverage_vs_merged=round(float(len(g)/merged_n.get(c, np.nan)), 3)
                     if merged_n.get(c, 0) else np.nan))
q = pd.DataFrame(rows).sort_values("n", ascending=False)
q["qa_pass"] = (q.n >= 1200) & (q.spread_years >= 6) & (q.max_year_share <= 0.55)
q.to_csv(os.path.join(DATA, "coverage_table.csv"), index=False)
print(q.head(20).to_string(index=False))
print("QA pass:", q[q.qa_pass].country.tolist())
