#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, step 0 -- extract global substations and power plants.

Sources: the OpenStreetMap planet power layer (points + polygon centroids)
and the WRI Global Power Plant Database v1.3. Checkpointed by sub-step;
run with one argument per call, in this order:

  python 00_extract_power_nodes.py sub_pt     substation points
  python 00_extract_power_nodes.py sub_poly   substation polygon centroids
  python 00_extract_power_nodes.py plant      OSM power plants
  python 00_extract_power_nodes.py wri        WRI plant database
  python 00_extract_power_nodes.py merge      merge and classify

Voltage parsing: the OSM 'voltage' tag may hold several ';'-separated
values; the maximum numeric voltage is taken. Classes: transmission
>= 110 kV; primary 30-110 kV; distribution 1-30 kV; unknown (treated as
distribution-class in the fragility model, which errs towards failure).

Outputs: data/substations_global.parquet (852,587),
         data/plants_global.parquet (169,491)
"""
import os, sys, re, time, tempfile
import numpy as np, pandas as pd, pyogrio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

OUT = P.ensure_dir(P.get("FIG4_DATA"))
SRC = P.get("OSM_POWER")
TMP = os.path.join(tempfile.gettempdir(), "fig4_power_nodes")
os.makedirs(TMP, exist_ok=True)
_WRI = P.get("WRI_PLANTS")

def parse_volt(v):
    if v is None or (isinstance(v, float) and np.isnan(v)): return np.nan
    nums = re.findall(r"\d+", str(v))
    if not nums: return np.nan
    vals = [int(x) for x in nums if 100 <= int(x) <= 2_000_000]
    return max(vals) if vals else np.nan

def kind_of(volts):
    if np.isnan(volts): return "unknown"
    if volts >= 110_000: return "transmission"
    if volts >= 30_000: return "primary"
    if volts >= 1_000: return "distribution"
    return "lv"

def finish(df, name):
    df.to_parquet(f"{TMP}/{name}.parquet")
    print(name, len(df), "saved", flush=True)

step = sys.argv[1]
t0 = time.time()
if step == "sub_pt":
    d = pyogrio.read_dataframe(SRC, layer="power_points", where="power='substation'",
                               columns=["voltage", "osm_id"])
    df = pd.DataFrame(dict(lon=d.geometry.x, lat=d.geometry.y,
                           voltage=d["voltage"], osm_id=d["osm_id"], src="pt"))
    finish(df, "sub_pt")
elif step == "sub_poly":
    d = pyogrio.read_dataframe(SRC, layer="power_polygons", where="power='substation'",
                               columns=["voltage", "osm_id"])
    c = d.geometry.centroid
    df = pd.DataFrame(dict(lon=c.x, lat=c.y, voltage=d["voltage"],
                           osm_id=d["osm_id"], src="poly"))
    finish(df, "sub_poly")
elif step == "plant":
    rows = []
    for layer, tag in [("power_points", "pt"), ("power_polygons", "poly")]:
        d = pyogrio.read_dataframe(SRC, layer=layer, where="power='plant'",
                                   columns=["name", "osm_id"])
        g = d.geometry.centroid if tag == "poly" else d.geometry
        rows.append(pd.DataFrame(dict(lon=g.x, lat=g.y, name=d["name"],
                                      src="osm_" + tag)))
    df = pd.concat(rows, ignore_index=True)
    finish(df, "plant_osm")
elif step == "wri":
    w = pd.read_csv(_WRI,
                    usecols=["latitude", "longitude", "name", "capacity_mw", "country_long"],
                    low_memory=False)
    df = pd.DataFrame(dict(lon=w.longitude, lat=w.latitude, name=w.name,
                           cap_mw=w.capacity_mw, country=w.country_long, src="wri"))
    finish(df, "plant_wri")
elif step == "merge":
    sp = pd.read_parquet(f"{TMP}/sub_pt.parquet")
    so = pd.read_parquet(f"{TMP}/sub_poly.parquet")
    subs = pd.concat([sp, so], ignore_index=True)
    subs["volt_max"] = subs.voltage.map(parse_volt)
    subs["kind"] = subs.volt_max.map(kind_of)
    subs = subs.dropna(subset=["lon", "lat"]).reset_index(drop=True)
    subs["sub_id"] = np.arange(len(subs))
    subs[["sub_id", "lon", "lat", "volt_max", "kind", "osm_id", "src"]].to_parquet(
        os.path.join(OUT, "substations_global.parquet"))
    po = pd.read_parquet(f"{TMP}/plant_osm.parquet")
    pw = pd.read_parquet(f"{TMP}/plant_wri.parquet")
    pl = pd.concat([po, pw], ignore_index=True).dropna(subset=["lon", "lat"])
    pl["plant_id"] = np.arange(len(pl))
    pl.to_parquet(os.path.join(OUT, "plants_global.parquet"))
    print("substations:", len(subs), subs.kind.value_counts().to_dict())
    print("plants:", len(pl), pl.src.value_counts().to_dict())
print("%.1fs" % (time.time() - t0))
