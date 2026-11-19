#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compile the station inventories. Assign country polygons by intersection and use the nearest polygon within 0.3 degrees for unmatched points. Retain Google charging records outside mainland China and AMap records assigned to mainland China. Standardise the Google fuel inventory separately."""
import os
import sys
import numpy as np
import pandas as pd
import pyogrio
import shapely
from shapely.strtree import STRtree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.ensure_dir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))

# Column names as exported by the collection tool (verbatim, in Chinese):
# station name, latitude, longitude, connector count, max power (kW),
# business status, place identifier.
GP_EV_COLS = {"lat": "纬度", "lon": "经度", "nconn": "充电接口数量",
              "maxkw": "最大功率(kW)", "status": "营业状态", "pid": "地点ID"}
GP_GAS_COLS = {"lat": "纬度", "lon": "经度", "status": "营业状态",
               "pid": "地点ID"}

def load_google(path, colmap):
    usecols = list(colmap.values())
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    df = df.rename(columns={v: k for k, v in colmap.items()})
    for c in ("lat", "lon"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    n0 = len(df)
    df = df.dropna(subset=["lat", "lon"]).drop_duplicates(subset="pid")
    print(f"{os.path.basename(path)}: {n0:,} raw -> {len(df):,} unique "
          "geolocated")
    return df.reset_index(drop=True)

def mainland_china_mask(lon, lat):
    """True for points assigned to mainland China.

    A point is assigned to the country polygon that contains it or, when it
    falls just outside every polygon (coastal and river-mouth sites), to the
    nearest polygon within 0.3 degrees. This is the same rule as step 3
    (03_assign_countries.py), so the replaced Google records and the retained
    AMap records are selected consistently with the final country labels.
    """
    ne = pyogrio.read_dataframe(P.get("NE10M"), columns=["NAME"])
    names = ne.NAME.values
    tree = STRtree(ne.geometry.values)
    pts = shapely.points(lon, lat)
    res = np.full(len(lon), "", dtype=object)
    CH = 200000
    for s in range(0, len(pts), CH):
        e = min(s + CH, len(pts))
        qi, ti = tree.query(pts[s:e], predicate="intersects")
        res[s + qi] = names[ti]
    un = np.where(res == "")[0]
    if len(un):
        ii, tt = tree.query_nearest(pts[un], max_distance=0.3, all_matches=False)
        res[un[ii]] = names[tt]
    return res == "China"

# ---- charging: Google worldwide minus mainland, plus AMap ------------------
gp = load_google(P.get("CHARGING_RAW_GOOGLE"), GP_EV_COLS)
gp["maxkw"] = pd.to_numeric(gp["maxkw"], errors="coerce")
gp["nconn"] = pd.to_numeric(gp["nconn"], errors="coerce")
gp["src"] = "GP"
in_chn = mainland_china_mask(gp.lon.values, gp.lat.values)
print(f"Google records inside mainland China (removed): {int(in_chn.sum()):,}")
gp = gp[~in_chn].reset_index(drop=True)

am = pd.read_parquet(os.path.join(DATA, "amap_stations_wgs84.parquet"))
am_in = mainland_china_mask(am.lon.values, am.lat.values)
print(f"AMap records outside mainland China (discarded): "
      f"{int((~am_in).sum()):,}")
am = am[am_in].reset_index(drop=True)

cols = ["lat", "lon", "status", "maxkw", "nconn", "src"]
merged = pd.concat([gp[cols], am[cols]], ignore_index=True)
merged.to_parquet(os.path.join(DATA, "charging_points_merged.parquet"),
                  index=False)
print(f"merged charging network: {len(merged):,} "
      f"(Google {int((merged.src=='GP').sum()):,} + "
      f"AMap {int((merged.src=='AMAP').sum()):,})")

# ---- fuel: single source, standardised --------------------------------------
gas = load_google(P.get("FUEL_RAW_GOOGLE"), GP_GAS_COLS)
gas["maxkw"] = np.nan
gas["nconn"] = np.nan
gas["src"] = "GP"
gas[cols].to_parquet(os.path.join(DATA, "fuel_points.parquet"), index=False)
print(f"fuel network: {len(gas):,}")
