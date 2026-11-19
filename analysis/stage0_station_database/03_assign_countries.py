#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage 0, step 3 -- assign a country and ISO3 code to every station.

Chunked STRtree point-in-polygon overlay on Natural Earth 10-m boundaries;
points that miss every polygon (coastal pixels, piers) are snapped to the
nearest boundary within 0.3 degrees. Writes the two station databases that
every figure pipeline starts from.

Outputs: shared/station_database/charging_stations_with_country.parquet
         shared/station_database/fuel_stations_with_country.parquet
"""
import os
import sys
import gc
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.strtree import STRtree
import shapely

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

world = gpd.read_file(P.get("NE10M"))
tree = STRtree(world.geometry.values)
names = world["NAME"].values
iso = world["ISO_A3_EH"].values

def assign(df):
    res_name = np.full(len(df), "", dtype=object)
    res_iso = np.full(len(df), "", dtype=object)
    CH = 120000
    for s in range(0, len(df), CH):
        e = min(s + CH, len(df))
        pts = shapely.points(df["lon"].values[s:e], df["lat"].values[s:e])
        qi, ti = tree.query(pts, predicate="intersects")
        res_name[s + qi] = names[ti]
        res_iso[s + qi] = iso[ti]
        del pts
        gc.collect()
        print("chunk", s, e, flush=True)
    un = np.where(res_name == "")[0]
    if len(un):
        pts = shapely.points(df["lon"].values[un], df["lat"].values[un])
        ii, tt = tree.query_nearest(pts, max_distance=0.3, all_matches=False)
        res_name[un[ii]] = names[tt]
        res_iso[un[ii]] = iso[tt]
    df = df.copy()
    df["country"] = res_name
    df["iso3"] = res_iso
    print("unmatched:", int((res_name == "").sum()), flush=True)
    return df

P.ensure_dir(P.get("DB_DIR"))
ev = pd.read_parquet(os.path.join(DATA, "charging_points_merged.parquet"))
ev = assign(ev)
ev.to_parquet(P.get("DB_CHARGING"), index=False)
print("saved charging database:", len(ev))

gas = pd.read_parquet(os.path.join(DATA, "fuel_points.parquet"))
gas = assign(gas)
gas.to_parquet(P.get("DB_FUEL"), index=False)
print("saved fuel database:", len(gas))
