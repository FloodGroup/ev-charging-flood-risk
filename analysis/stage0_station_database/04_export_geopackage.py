#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage 0, step 4 -- export the merged charging database as a GeoPackage
(openable directly in QGIS; EPSG:4326; fields country/iso3/src/status/
maxkw/nconn). Written in chunks to bound memory.

Output: shared/station_database/charging_stations_merged.gpkg
"""
import os
import sys
import gc
import pandas as pd
import geopandas as gpd
import pyogrio

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

out = P.get("DB_MERGED_GPKG")
if os.path.exists(out):
    os.remove(out)
df = pd.read_parquet(P.get("DB_CHARGING"))
n = len(df)
CH = 100000
for i, s in enumerate(range(0, n, CH)):
    e = min(s + CH, n)
    part = df.iloc[s:e]
    gdf = gpd.GeoDataFrame(
        part[["country", "iso3", "src", "status", "maxkw", "nconn"]].copy(),
        geometry=gpd.points_from_xy(part.lon, part.lat), crs="EPSG:4326")
    gdf["status"] = gdf["status"].astype(object).where(gdf["status"].notna(),
                                                       None)
    pyogrio.write_dataframe(gdf, out, layer="charging_stations",
                            driver="GPKG", append=(i > 0),
                            SPATIAL_INDEX="NO")
    del gdf, part
    gc.collect()
    print("chunk", s, e, flush=True)
print("WRITE DONE rows:", n, flush=True)
