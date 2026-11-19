#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1, step 3 -- verification.

Recomputes every headline number that Figure 1 and the manuscript quote from
the station database, independently of step 1, and compares the results with
the values stored in data/global_summary.json and data/country_stats.csv.
Exits non-zero (and prints FAIL lines) on any mismatch.

Checks: global station counts (947,602 / 1,344,301), source split
(Google Places / AMap), the global charging-to-fuel ratio (0.70), the
IEA ten-year growth factor (x40) and the 2025 world total (7.0 M), and the
points-per-station cross-audit values quoted in the manuscript (six
audited markets; Korea is excluded -- see 01_compute_statistics.py).
"""
import os
import sys
import json
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG1_DATA")
failures = []

def check(name, got, want, tol=0.0):
    if isinstance(want, (int,)) and tol == 0:
        ok = int(got) == int(want)
    else:
        ok = abs(float(got) - float(want)) <= tol
    print(f"{'OK  ' if ok else 'FAIL'} {name}: {got} vs {want}")
    if not ok:
        failures.append(name)

ev = pd.read_parquet(P.get("DB_CHARGING"))
gs = pd.read_parquet(P.get("DB_FUEL"))
G = json.load(open(os.path.join(DATA, "global_summary.json")))

check("charging station count", len(ev), G["n_ev"])
check("fuel station count", len(gs), G["n_gas"])
check("Google Places share of charging DB", int((ev.src == "GP").sum()),
      G["n_ev_gp"])
check("AMap share of charging DB", int((ev.src == "AMAP").sum()),
      G["n_ev_amap"])
check("global charging-to-fuel ratio", round(len(ev) / len(gs), 2), 0.70,
      tol=0.005)

# IEA growth factor: world total 2025 / 2015 (charging points)
iea = pd.read_csv(os.path.join(P.get("IEA_DIR"), "iea_world_by_speed.csv"))
tot = iea.set_index("year")[["slow", "fast", "ultra"]].sum(axis=1)
growth = tot.loc[2025] / tot.loc[2015]
check("IEA ten-year growth factor (x40)", round(growth), 40, tol=1)
check("IEA 2025 world total (millions)", round(tot.loc[2025] / 1e6, 1), 7.0,
      tol=0.05)

# cross-audit values quoted in the manuscript (points per station)
cc = pd.read_csv(os.path.join(DATA, "crosscheck_stations_vs_iea.csv")
                 ).set_index("country")
for ctry, want in [("China", 16.4), ("France", 4.0),
                   ("United States of America", 3.8), ("Germany", 2.8),
                   ("Netherlands", 2.1), ("Japan", 1.2)]:
    got = cc.loc[ctry, "points_per_station"]
    check(f"points per station, {ctry}", round(float(got), 1), want, tol=0.05)

print("-" * 60)
if failures:
    print(f"RESULT: {len(failures)} check(s) FAILED")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
