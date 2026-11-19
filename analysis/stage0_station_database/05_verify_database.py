#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage 0, step 5 -- verification of the rebuilt station databases.

Compares the databases produced by steps 1-3 with the published counts:
947,602 charging stations (662,962 Google Places + 284,640 AMap) and
1,344,301 fuel stations; the ten largest charging markets; and the count
of stations without a country assignment (48 charging / 40 fuel).
Exits non-zero on any mismatch.

"""
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

failures = []

def check(name, got, want):
    ok = int(got) == int(want)
    print(f"{'OK  ' if ok else 'FAIL'} {name}: {got} vs {want}")
    if not ok:
        failures.append(name)

ev = pd.read_parquet(P.get("DB_CHARGING"))
gas = pd.read_parquet(P.get("DB_FUEL"))

check("charging stations", len(ev), 947602)
check("  of which Google Places", (ev.src == "GP").sum(), 662962)
check("  of which AMap (mainland China)", (ev.src == "AMAP").sum(), 284640)
check("fuel stations", len(gas), 1344301)
check("charging stations without country", (ev.country == "").sum(), 48)
check("fuel stations without country", (gas.country == "").sum(), 40)

# cross-check the largest markets against Figure 1's stored country table
ref_csv = os.path.join(P.get("FIG1_DATA"), "country_stats.csv")
if os.path.exists(ref_csv):
    ref = pd.read_csv(ref_csv).set_index("country")
    vc = ev.country.value_counts()
    for c in ref.sort_values("ev_n", ascending=False).head(10).index:
        check(f"stations in {c}", int(vc.get(c, 0)), int(ref.loc[c, "ev_n"]))
else:
    print("note: Figure 1 country_stats.csv not present yet -- top-market "
          "cross-check skipped")

print("-" * 60)
if failures:
    print(f"RESULT: {len(failures)} check(s) FAILED")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
