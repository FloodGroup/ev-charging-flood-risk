#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, stage 5 -- export the complete global grid model as a
GIS-ready GeoPackage (openable directly in QGIS).

grid_model_global.gpkg layers
  power_grid_lines        3.68 M original network lines  (fid, source,
                          length_m, component)
  connector_substations   physical tie: substation -> attachment vertex
  connector_plants        physical tie: plant -> attachment vertex
  connector_chargers      physical tie: charger -> attachment vertex
  substations             852,587 points (class, site, country, depths,
                          P_fail at RP100/500)
  substation_sites        clustered physical sites (feeder identity)
  power_plants            169,491 points (+ nearest site over the network)
  chargers                947,602 points (feeder set, probabilities, mode)
  supply_links            logical charger -> feeder-site relations
                          (rank 1 = primary, rank 2 = redundant backup)

grid_model_light.gpkg: the same layers cropped to the four vignette
cities for local GIS inspection. Both files contain station locations.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

FIG4_DATA = P.get("FIG4_DATA")
FIG2_DATA = P.get("FIG2_DATA")
GRAPH = P.ensure_dir(os.path.join(FIG4_DATA, "checkpoints", "graph"))
OUT = FIG4_DATA
import json, time, gc
import numpy as np, pandas as pd, shapely, pyogrio
from scipy.stats import norm

GPKG_IN = P.get("GRIDFINDER")
GIS = P.ensure_dir(os.path.join(FIG4_DATA, "gis"))
OUTFILE = os.path.join(GIS, "grid_model_global.gpkg")
LIGHT   = os.path.join(GIS, "grid_model_light.gpkg")
THETA_DIST, THETA_TRANS, BETA = 0.5, 2.0, 0.6
t0 = time.time()

def pf(depth, majr):
    th = np.where(majr, THETA_TRANS, THETA_DIST)
    with np.errstate(divide="ignore"):
        z = np.log(np.maximum(depth, 1e-12) / th) / BETA
    return np.where(depth > 0, norm.cdf(z), 0.0)

import geopandas as gpd

# ---------------- load shared model pieces ----------------
zv = np.load(f"{GRAPH}/vertices.npz")
vlon, vlat, vcomp = zv["lon"], zv["lat"], zv["comp"]
feat = pd.read_parquet(f"{GRAPH}/features.parquet")
SS   = pd.read_parquet(f"{OUT}/sub_sites.parquet")
SITE = pd.read_parquet(f"{OUT}/site_table.parquet")
NETC = pd.read_parquet(f"{OUT}/charger_network_assignment.parquet")
NETP = pd.read_parquet(f"{OUT}/plant_network_assignment.parquet")

S = pd.read_parquet(f"{FIG4_DATA}/substation_depths.parquet")
SK = pd.read_parquet(f"{FIG4_DATA}/substations_global.parquet")
SC = pd.read_parquet(f"{FIG4_DATA}/substation_country.parquet")
S = S.merge(SK[["sub_id", "kind", "volt_max", "osm_id", "src"]], on="sub_id")
S = S.merge(SC, on="sub_id").merge(SS, on="sub_id")
za = np.load(f"{GRAPH}/attach_subs.npz");     sub_vid, sub_d = za["vid"], za["dist_m"]
za = np.load(f"{GRAPH}/attach_plants.npz");   pl_vid,  pl_d  = za["vid"], za["dist_m"]
za = np.load(f"{GRAPH}/attach_chargers.npz"); ch_vid,  ch_d  = za["vid"], za["dist_m"]

NSITE_ALL = int(SS.site_id.max()) + 1
site_major = np.zeros(NSITE_ALL, bool)
np.logical_or.at(site_major, S.site_id.values,
                 S.kind.isin(["transmission", "primary"]).values)
site_d100 = np.zeros(NSITE_ALL); site_d500 = np.zeros(NSITE_ALL)
np.maximum.at(site_d100, S.site_id.values, np.nan_to_num(S.d100.values.astype(float)))
np.maximum.at(site_d500, S.site_id.values, np.nan_to_num(S.d500.values.astype(float)))
site_p100 = pf(site_d100, site_major); site_p500 = pf(site_d500, site_major)

sd = pd.read_parquet(f"{FIG2_DATA}/station_depths.parquet")
ev = sd[sd.kind == "ev"].reset_index(drop=True)
del sd; gc.collect()

CAS = json.load(open(f"{OUT}/assignment_stats.json"))
R_SERV = CAS["r_serv_m"]; R_EUC = CAS["r_euclidean_fallback_m"]
EUC = pd.read_parquet(f"{FIG4_DATA}/charger_assignment.parquet")

rank2site = np.full(int(SITE.site_rank.max()) + 1, -1, np.int64)
rank2site[SITE.site_rank.values] = SITE.site_id.values
site_of_sub = np.zeros(int(SS.sub_id.max()) + 1, np.int64)
site_of_sub[SS.sub_id.values] = SS.site_id.values

d1 = NETC.d1_net_m.values; d2 = NETC.d2_net_m.values
net_ok = (NETC.s1_rank.values >= 0) & np.isfinite(d1) & (d1 <= R_SERV)
s2_ok  = net_ok & (NETC.s2_rank.values >= 0) & np.isfinite(d2) & (d2 <= R_SERV)
s1_site = np.where(net_ok, rank2site[np.maximum(NETC.s1_rank.values, 0)], -1)
s2_site = np.where(s2_ok, rank2site[np.maximum(NETC.s2_rank.values, 0)], -1)
fb = ~net_ok & (EUC.dist1_m.values <= R_EUC)
s1_site = np.where(net_ok, s1_site, np.where(fb, site_of_sub[EUC.sub1.values], -1))
mode = np.where(net_ok, "network",
                np.where(fb, "euclidean_fallback", "grid_independent"))

# site representative coordinates (rep member substation)
rep_lon = np.zeros(NSITE_ALL); rep_lat = np.zeros(NSITE_ALL)
rep_lon[site_of_sub[SK.sub_id.values]] = 0  # init
srt = SK.merge(SS, on="sub_id")
tmp = srt.assign(v=np.nan_to_num(srt.volt_max.values, nan=-1)) \
         .sort_values(["site_id", "v"], ascending=[True, False]) \
         .drop_duplicates("site_id")
rep_lon[tmp.site_id.values] = tmp.lon.values
rep_lat[tmp.site_id.values] = tmp.lat.values

def write(path, gdf, layer, append):
    pyogrio.write_dataframe(gdf, path, layer=layer, driver="GPKG",
                            append=append)
    print(f"  layer {layer}: {len(gdf):,} rows  {time.time()-t0:.0f}s",
          flush=True)

def lines_from_pairs(x0, y0, x1, y1):
    a = np.stack([np.c_[x0, y0], np.c_[x1, y1]], axis=1)
    return shapely.linestrings(a)

P.fresh(OUTFILE)
P.fresh(LIGHT)

# ---------------- 1: original power lines (chunked append) ----------------
info = pyogrio.read_info(GPKG_IN, layer="grid")
NF = int(info["features"])
CH = 400_000
for i0 in range(0, NF, CH):
    g = pyogrio.read_dataframe(GPKG_IN, layer="grid", skip_features=i0,
                               max_features=min(CH, NF - i0))
    sl = feat.iloc[i0:i0 + len(g)]
    g = gpd.GeoDataFrame(
        dict(line_id=sl.fid.values, source=g["source"].values,
             length_m=sl.length_m.values, component=sl.comp.values),
        geometry=g.geometry.values, crs="EPSG:4326")
    write(OUTFILE, g, "power_grid_lines", append=i0 > 0)
    del g; gc.collect()

# ---------------- 2-4: connectors ----------------
ok = sub_vid >= 0
g = gpd.GeoDataFrame(
    dict(sub_id=S.sub_id.values[ok], site_id=S.site_id.values[ok],
         attach_dist_m=sub_d[ok].astype(np.float32),
         component=vcomp[sub_vid[ok]]),
    geometry=lines_from_pairs(S.lon.values[ok], S.lat.values[ok],
                              vlon[sub_vid[ok]], vlat[sub_vid[ok]]),
    crs="EPSG:4326")
write(OUTFILE, g, "connector_substations", append=True)

PL = pd.read_parquet(f"{FIG4_DATA}/plants_global.parquet")   # do not shadow paths alias P
ok = pl_vid >= 0
g = gpd.GeoDataFrame(
    dict(plant_id=PL.plant_id.values[ok],
         attach_dist_m=pl_d[ok].astype(np.float32),
         component=vcomp[pl_vid[ok]]),
    geometry=lines_from_pairs(PL.lon.values[ok], PL.lat.values[ok],
                              vlon[pl_vid[ok]], vlat[pl_vid[ok]]),
    crs="EPSG:4326")
write(OUTFILE, g, "connector_plants", append=True)

ok = ch_vid >= 0
g = gpd.GeoDataFrame(
    dict(charger_idx=np.arange(len(ev), dtype=np.int64)[ok],
         attach_dist_m=ch_d[ok].astype(np.float32),
         component=vcomp[ch_vid[ok]]),
    geometry=lines_from_pairs(ev.lon.values[ok], ev.lat.values[ok],
                              vlon[ch_vid[ok]], vlat[ch_vid[ok]]),
    crs="EPSG:4326")
write(OUTFILE, g, "connector_chargers", append=True)
del g; gc.collect()

# ---------------- 5: substations ----------------
subs_pf = pf(S.d100.values.astype(float),
             S.kind.isin(["transmission", "primary"]).values)
g = gpd.GeoDataFrame(
    dict(sub_id=S.sub_id.values, site_id=S.site_id.values,
         kind=S.kind.values, volt_max=S.volt_max.values,
         osm_id=S.osm_id.values, country=S.country.values,
         continent=S.continent.values,
         attached=(sub_vid >= 0), attach_dist_m=sub_d.astype(np.float32),
         d010=S.d010.values, d020=S.d020.values, d050=S.d050.values,
         d075=S.d075.values, d100=S.d100.values, d200=S.d200.values,
         d500=S.d500.values,
         p_fail_rp100=subs_pf.astype(np.float32)),
    geometry=shapely.points(np.c_[S.lon.values, S.lat.values]), crs="EPSG:4326")
write(OUTFILE, g, "substations", append=True)

# ---------------- 6: sites ----------------
n_primary = np.zeros(NSITE_ALL, np.int64)
np.add.at(n_primary, np.maximum(s1_site, 0), (s1_site >= 0).astype(np.int64))
site_ids = np.arange(NSITE_ALL)
site_nmem = np.bincount(SS.site_id.values, minlength=NSITE_ALL)
att = np.zeros(NSITE_ALL, bool); att[SITE.site_id.values] = True
g = gpd.GeoDataFrame(
    dict(site_id=site_ids, n_members=site_nmem,
         major=site_major, attached=att,
         d100=site_d100.astype(np.float32), d500=site_d500.astype(np.float32),
         p_fail_rp100=site_p100.astype(np.float32),
         p_fail_rp500=site_p500.astype(np.float32),
         n_chargers_primary=n_primary),
    geometry=shapely.points(np.c_[rep_lon, rep_lat]), crs="EPSG:4326")
write(OUTFILE, g, "substation_sites", append=True)

# ---------------- 7: plants ----------------
pl_site = np.where((NETP.s1_rank.values >= 0),
                   rank2site[np.maximum(NETP.s1_rank.values, 0)], -1)
g = gpd.GeoDataFrame(
    dict(plant_id=PL.plant_id.values, name=PL.name.astype(str).values,
         src=PL.src.values, cap_mw=PL.cap_mw.values,
         attached=(pl_vid >= 0), attach_dist_m=pl_d.astype(np.float32),
         nearest_site=pl_site,
         net_dist_to_site_m=NETP.d1_net_m.values.astype(np.float32)),
    geometry=shapely.points(np.c_[PL.lon.values, PL.lat.values]), crs="EPSG:4326")
write(OUTFILE, g, "power_plants", append=True)

# ---------------- 8: chargers ----------------
p1 = np.where(s1_site >= 0, site_p100[np.maximum(s1_site, 0)], 0.0)
p2 = np.where(s2_site >= 0, site_p100[np.maximum(s2_site, 0)], 1.0)
direct = ev.d100.values > 0
p_grid = np.where(s1_site >= 0, p1, 0.0)
p_off = np.where(direct, 1.0, p_grid)
g = gpd.GeoDataFrame(
    dict(charger_idx=np.arange(len(ev), dtype=np.int64),
         country=ev.country.astype(str).values,
         d100=ev.d100.values,
         direct_rp100=direct,
         supply_mode=mode,
         s1_site=s1_site, d1_net_m=np.where(net_ok, d1, np.nan).astype(np.float32),
         s2_site=s2_site, d2_net_m=np.where(s2_ok, d2, np.nan).astype(np.float32),
         redundant=s2_ok,
         p_feeder_fail_rp100=p_grid.astype(np.float32),
         p_offline_rp100=p_off.astype(np.float32),
         p_offline_redundant_rp100=np.where(
             direct, 1.0, p_grid * np.where(s2_ok, p2, 1.0)).astype(np.float32)),
    geometry=shapely.points(np.c_[ev.lon.values, ev.lat.values]), crs="EPSG:4326")
write(OUTFILE, g, "chargers", append=True)

# ---------------- 9: supply links ----------------
m1 = s1_site >= 0
l1 = gpd.GeoDataFrame(
    dict(charger_idx=np.arange(len(ev), dtype=np.int64)[m1],
         site_id=s1_site[m1], rank=np.ones(int(m1.sum()), np.int16),
         mode=mode[m1],
         net_dist_m=np.where(net_ok, d1, np.nan)[m1].astype(np.float32)),
    geometry=lines_from_pairs(ev.lon.values[m1], ev.lat.values[m1],
                              rep_lon[s1_site[m1]], rep_lat[s1_site[m1]]),
    crs="EPSG:4326")
m2 = s2_site >= 0
l2 = gpd.GeoDataFrame(
    dict(charger_idx=np.arange(len(ev), dtype=np.int64)[m2],
         site_id=s2_site[m2], rank=np.full(int(m2.sum()), 2, np.int16),
         mode=np.full(int(m2.sum()), "network"),
         net_dist_m=d2[m2].astype(np.float32)),
    geometry=lines_from_pairs(ev.lon.values[m2], ev.lat.values[m2],
                              rep_lon[s2_site[m2]], rep_lat[s2_site[m2]]),
    crs="EPSG:4326")
links = pd.concat([l1, l2], ignore_index=True)
write(OUTFILE, links, "supply_links", append=True)
print(f"global model written: {os.path.getsize(OUTFILE)/1e9:.2f} GB", flush=True)

# ---------------- city-window extracts ----------------
WINDOWS = {}
for mid in [2288, 12703, 10912, 1048]:
    z = np.load(f"{FIG4_DATA}/vignette_v2_{mid}.npz")
    w, s_, e, n = [float(x) for x in z["bbox"]]
    WINDOWS[mid] = (w, s_, e, n)
print("light windows:", WINDOWS, flush=True)

if os.path.exists(LIGHT):
    os.remove(LIGHT)
first = True
for layer in ["power_grid_lines", "connector_substations", "connector_plants",
              "connector_chargers", "substations", "substation_sites",
              "power_plants", "chargers", "supply_links"]:
    parts = []
    for mid, (w, s_, e, n) in WINDOWS.items():
        gg = pyogrio.read_dataframe(OUTFILE, layer=layer, bbox=(w, s_, e, n))
        gg["vignette_metro"] = mid
        parts.append(gg)
    gg = pd.concat(parts, ignore_index=True)
    if len(gg):
        pyogrio.write_dataframe(gg, LIGHT, layer=layer, driver="GPKG",
                                append=not first)
        first = False
    print(f"  light {layer}: {len(gg):,}", flush=True)
print(f"light model: {os.path.getsize(LIGHT)/1e6:.1f} MB", flush=True)
print(f"DONE {time.time()-t0:.0f}s", flush=True)
