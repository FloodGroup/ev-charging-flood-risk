#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage 0, step 1 -- AMap mainland-China chargers: GCJ-02 -> WGS-84.

Chinese map providers are required to serve coordinates in the GCJ-02
obfuscated datum. This step applies the standard published inverse
transformation (first-order inverse; accuracy 1-2 m, exact for practical
purposes) and maps the AMap fields onto the common station schema. The
mean correction is about 499 m -- without it, every subsequent flood-depth
sample for China would be systematically displaced.

Input:  source_data/stations_charging_amap/... (GCJ-02 CSV, n = 284,640
        unique station identifiers)
Output: data/amap_stations_wgs84.parquet
"""
import os
import sys
import math
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

OUT_DIR = P.ensure_dir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))

A = 6378245.0
EE = 0.00669342162296594323

def _transform_lat(x, y):
    r = (-100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y
         + 0.2 * np.sqrt(np.abs(x)))
    r += (20.0 * np.sin(6.0 * x * math.pi)
          + 20.0 * np.sin(2.0 * x * math.pi)) * 2.0 / 3.0
    r += (20.0 * np.sin(y * math.pi)
          + 40.0 * np.sin(y / 3.0 * math.pi)) * 2.0 / 3.0
    r += (160.0 * np.sin(y / 12.0 * math.pi)
          + 320.0 * np.sin(y * math.pi / 30.0)) * 2.0 / 3.0
    return r

def _transform_lon(x, y):
    r = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * np.sqrt(np.abs(x))
    r += (20.0 * np.sin(6.0 * x * math.pi)
          + 20.0 * np.sin(2.0 * x * math.pi)) * 2.0 / 3.0
    r += (20.0 * np.sin(x * math.pi)
          + 40.0 * np.sin(x / 3.0 * math.pi)) * 2.0 / 3.0
    r += (150.0 * np.sin(x / 12.0 * math.pi)
          + 300.0 * np.sin(x / 30.0 * math.pi)) * 2.0 / 3.0
    return r

def gcj02_to_wgs84(lon, lat):
    dlat = _transform_lat(lon - 105.0, lat - 35.0)
    dlon = _transform_lon(lon - 105.0, lat - 35.0)
    radlat = lat / 180.0 * math.pi
    magic = np.sin(radlat)
    magic = 1 - EE * magic * magic
    sqrtmagic = np.sqrt(magic)
    dlat = (dlat * 180.0) / ((A * (1 - EE)) / (magic * sqrtmagic) * math.pi)
    dlon = (dlon * 180.0) / (A / sqrtmagic * np.cos(radlat) * math.pi)
    return lon - dlon, lat - dlat   # first-order inverse (~1-2 m)

df = pd.read_csv(P.get("CHARGING_RAW_AMAP"))
print("rows:", len(df), "| unique id:", df.id.nunique())
df = df.dropna(subset=["longitude", "latitude"]).drop_duplicates(subset="id")
lon_w, lat_w = gcj02_to_wgs84(df.longitude.values, df.latitude.values)
out = pd.DataFrame(dict(lat=lat_w, lon=lon_w, status=pd.NA, maxkw=np.nan,
                        nconn=np.nan))
out["src"] = "AMAP"
out.to_parquet(os.path.join(OUT_DIR, "amap_stations_wgs84.parquet"),
               index=False)
print("converted & saved:", len(out), "| mean shift (m):",
      round(float(np.mean(np.hypot((df.longitude.values - lon_w) * 88000,
                                   (df.latitude.values - lat_w) * 111000))),
            1))
