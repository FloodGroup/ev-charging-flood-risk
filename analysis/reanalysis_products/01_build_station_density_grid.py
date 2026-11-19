#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shareable re-analysis product: station density and flood exposure on a
0.25-degree grid.

The station-level databases cannot be redistributed (see
RESTRICTED_DATA_NOTICE.md). This script aggregates them to 0.25-degree cells,
which is the form released with the article:

  lon_w, lat_s                         south-west corner of the cell (degrees)
  charging_stations, fuel_stations     number of stations in the cell
  *_exposed_rp10 ... rp500             stations with modelled water depth > 0
                                       at each return period
  charging_expected_offline_rp100      expected chargers without power at the
                                       100-year flood (primary feeder)
  charging_expected_offline_backup_rp100   same, backup-feeder counterfactual

Inputs:  analysis/fig2_flood_exposure/data/station_depths.parquet
         analysis/fig4_grid_cascade/data/charger_cascade_rp100.parquet
Output:  analysis/reanalysis_products/data/station_density_exposure_0p25deg.csv

The column sums are checked against the totals of the article before the
file is written.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

CELL = 0.25
RPS = ["010", "020", "050", "075", "100", "200", "500"]
OUT = os.path.join(P.PACKAGE_ROOT, "analysis", "reanalysis_products", "data",
                   "station_density_exposure_0p25deg.csv")

depths = pd.read_parquet(os.path.join(P.get("FIG2_DATA"), "station_depths.parquet"))
cascade = pd.read_parquet(os.path.join(P.get("FIG4_DATA"), "charger_cascade_rp100.parquet"))

kinds = depths["kind"].astype(str)
labels = sorted(kinds.unique())
charging_label = [k for k in labels if k.lower() in ("ev", "charging", "charger")]
fuel_label = [k for k in labels if k.lower() in ("gas", "fuel", "gasoline")]
if len(charging_label) != 1 or len(fuel_label) != 1:
    raise SystemExit(f"unexpected station kinds: {labels}")

charging = depths[kinds == charging_label[0]].reset_index(drop=True)
fuel = depths[kinds == fuel_label[0]].reset_index(drop=True)
if len(charging) != len(cascade):
    raise SystemExit("charging stations and cascade table differ in length")
if not np.array_equal(cascade["direct"].values, (charging["d100"].values > 0)):
    raise SystemExit("cascade table is not aligned with the charging stations")


def cell(frame):
    ix = np.floor(frame["lon"].values / CELL).astype(np.int64)
    iy = np.floor(frame["lat"].values / CELL).astype(np.int64)
    return ix, iy


def aggregate(frame, prefix, extra=None):
    ix, iy = cell(frame)
    cols = {"ix": ix, "iy": iy, f"{prefix}_stations": 1}
    for rp in RPS:
        cols[f"{prefix}_exposed_rp{int(rp)}"] = (frame[f"d{rp}"].values > 0).astype(np.int64)
    for name, values in (extra or {}).items():
        cols[name] = values
    return pd.DataFrame(cols).groupby(["ix", "iy"], sort=True).sum()


ev = aggregate(charging, "charging", {
    "charging_expected_offline_rp100": cascade["p_off_primary"].values.astype(np.float64),
    "charging_expected_offline_backup_rp100": cascade["p_off_redundant_new"].values.astype(np.float64),
})
fu = aggregate(fuel, "fuel")
grid = ev.join(fu, how="outer").fillna(0).reset_index()
grid.insert(0, "lat_s", grid.pop("iy") * CELL)
grid.insert(0, "lon_w", grid.pop("ix") * CELL)
for c in grid.columns:
    if c.endswith("_stations") or "_exposed_" in c:
        grid[c] = grid[c].astype(np.int64)
for c in ("charging_expected_offline_rp100", "charging_expected_offline_backup_rp100"):
    grid[c] = grid[c].round(3)
grid = grid.sort_values(["lat_s", "lon_w"]).reset_index(drop=True)

checks = [
    ("charging stations", grid["charging_stations"].sum(), len(charging)),
    ("fuel stations", grid["fuel_stations"].sum(), len(fuel)),
    ("charging exposed RP100", grid["charging_exposed_rp100"].sum(), int((charging["d100"] > 0).sum())),
    ("fuel exposed RP100", grid["fuel_exposed_rp100"].sum(), int((fuel["d100"] > 0).sum())),
]
for name, got, want in checks:
    print(f"{'OK  ' if got == want else 'FAIL'} {name}: {got:,} vs {want:,}")
    if got != want:
        raise SystemExit("column sums do not reproduce the station totals")
exp = float(cascade["p_off_primary"].astype(np.float64).sum())
print(f"OK   expected offline RP100: {grid['charging_expected_offline_rp100'].sum():,.1f} vs {exp:,.1f}")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
grid.to_csv(P.fresh(OUT), index=False)
print(f"{len(grid):,} cells written to {OUT}")
print("GRID PRODUCT DONE")
