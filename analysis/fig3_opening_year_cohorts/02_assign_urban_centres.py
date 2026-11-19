#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3, step 2 -- assign a GHS-UCDB urban-centre identifier to each
archive station (STRtree point-in-polygon; -1 = outside any urban centre).

Output: data/archive_metro.npy
"""
import os, sys
import numpy as np, pandas as pd, pyogrio, shapely
from shapely.strtree import STRtree
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG3_DATA")
UC = P.get("UCDB_GPKG")

st = pd.read_parquet(os.path.join(DATA, "archive_stations.parquet"))
uc = pyogrio.read_dataframe(UC, columns=["ID_HDC_G0", "UC_NM_MN"])
tree = STRtree(uc.geometry.values)
pts = shapely.points(st.lon.values, st.lat.values)
q = tree.query(pts, predicate="intersects")
metro = np.full(len(st), -1, np.int64)
metro[q[0]] = uc.ID_HDC_G0.values[q[1]].astype(np.int64)
np.save(os.path.join(DATA, "archive_metro.npy"), metro)
print("in-metro share %.1f%%" % (100*(metro >= 0).mean()), "| metros:", len(np.unique(metro[metro>=0])))
