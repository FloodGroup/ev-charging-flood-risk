#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build local feeder maps for Düsseldorf, Tainan, Handan and Québec City. Each extract contains charger states, feeder sites, paths on the intact and residual networks, and the local flood mask. Keep the resulting station-level arrays private.

Usage: python 09_extract_city_vignettes.py <metro_id> <side_km> [full|standard]
Full mode reads the source layers. Standard mode reads locally cached geometry."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import norm
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

FIG4_DATA = P.get("FIG4_DATA")
FIG2_DATA = P.get("FIG2_DATA")
THETA_DIST, THETA_TRANS, BETA = 0.5, 2.0, 0.6
R_SERV, R_EUC = 10_000.0, 8_000.0
METRO = float(sys.argv[1])
MODE = sys.argv[3] if len(sys.argv) > 3 else "auto"

OUT_NPZ = os.path.join(FIG4_DATA, f"vignette_v2_{int(METRO)}.npz")
GEOMETRY_CACHE = os.path.join(FIG4_DATA, f"vignette_{int(METRO)}.npz")

# ---------- charger states under the network model (identical to stage 07) ----------
sd = pd.read_parquet(os.path.join(FIG2_DATA, "station_depths.parquet"),
                     columns=["lon", "lat", "kind", "d100"])
ev_mask = (sd.kind == "ev").values
metro = np.concatenate([np.load(os.path.join(FIG2_DATA, f"metro_p{i}.npy"))
                        for i in range(4)])[ev_mask]
ev = sd[ev_mask].reset_index(drop=True)

SS = pd.read_parquet(os.path.join(FIG4_DATA, "sub_sites.parquet"))
SITE = pd.read_parquet(os.path.join(FIG4_DATA, "site_table.parquet"))
S = pd.read_parquet(os.path.join(FIG4_DATA, "substation_depths.parquet")).merge(
    pd.read_parquet(os.path.join(FIG4_DATA, "substations_global.parquet"),
                    columns=["sub_id", "kind", "volt_max", "lon", "lat"]
                    ).rename(columns={"lon": "slon", "lat": "slat"}),
    on="sub_id").merge(SS, on="sub_id")
NSITE_ALL = int(SS.site_id.max()) + 1
site_major = np.zeros(NSITE_ALL, bool)
np.logical_or.at(site_major, S.site_id.values,
                 S.kind.isin(["transmission", "primary"]).values)
site_d100 = np.zeros(NSITE_ALL)
np.maximum.at(site_d100, S.site_id.values,
              np.nan_to_num(S.d100.values.astype(float)))
th = np.where(site_major, THETA_TRANS, THETA_DIST)
with np.errstate(divide="ignore"):
    site_p100 = np.where(site_d100 > 0,
                         norm.cdf(np.log(np.maximum(site_d100, 1e-12) / th)
                                  / BETA), 0.0)
tmp = S.assign(v=np.nan_to_num(S.volt_max.values, nan=-1)) \
       .sort_values(["site_id", "v"], ascending=[True, False]) \
       .drop_duplicates("site_id")
rep_lon = np.zeros(NSITE_ALL); rep_lat = np.zeros(NSITE_ALL)
rep_lon[tmp.site_id.values] = tmp.slon.values
rep_lat[tmp.site_id.values] = tmp.slat.values

NET = pd.read_parquet(os.path.join(FIG4_DATA,
                                   "charger_network_assignment.parquet"))
EUC = pd.read_parquet(os.path.join(FIG4_DATA, "charger_assignment.parquet"))
rank2site = np.full(int(SITE.site_rank.max()) + 1, -1, np.int64)
rank2site[SITE.site_rank.values] = SITE.site_id.values
site_of_sub = np.zeros(int(SS.sub_id.max()) + 1, np.int64)
site_of_sub[SS.sub_id.values] = SS.site_id.values
d1 = NET.d1_net_m.values; d2 = NET.d2_net_m.values
net_ok = (NET.s1_rank.values >= 0) & np.isfinite(d1) & (d1 <= R_SERV)
s2_ok = net_ok & (NET.s2_rank.values >= 0) & np.isfinite(d2) & (d2 <= R_SERV)
fb = ~net_ok & (EUC.dist1_m.values <= R_EUC)
s1_site = np.where(net_ok, rank2site[np.maximum(NET.s1_rank.values, 0)],
                   np.where(fb, site_of_sub[EUC.sub1.values], -1))
s2_site = np.where(s2_ok, rank2site[np.maximum(NET.s2_rank.values, 0)], -1)
p_grid = np.where(s1_site >= 0, site_p100[np.maximum(s1_site, 0)], 0.0)
p_grid2 = np.where(s2_ok, site_p100[np.maximum(s2_site, 0)], np.nan)
CAS = pd.read_parquet(os.path.join(FIG4_DATA,
                                   "charger_cascade_rp100.parquet"))
p_off_red_new = CAS.p_off_redundant_new.values
p_off_red_old = CAS.p_off_redundant_old.values

# ---------- this metro ----------
m = metro == METRO
lon, lat = ev.lon.values[m], ev.lat.values[m]
lat0 = np.deg2rad(lat.mean()); kx = np.cos(lat0) * 111.32; ky = 110.54

direct = ev.d100.values[m] > 0
pg = p_grid[m]
pg2 = p_grid2[m]
pred_new = p_off_red_new[m]
pred_old = p_off_red_old[m]
s1m = s1_site[m]
s2m = s2_site[m]
fbm = fb[m]
gonly = ~direct & (pg >= 0.5)
rescued = gonly & (pred_new < 0.5)
e_off = float(direct.sum() + pg[~direct].sum())
e_off_red = float(pred_new.sum())
lk = gonly & (s1m >= 0)
bk = rescued  # rescue links: charger -> backup site

# ---------- window geometry: cached (standard) or from sources (full) ----------
have_cache = os.path.exists(OUT_NPZ) or os.path.exists(GEOMETRY_CACHE)
if MODE == "auto":
    MODE = "standard" if have_cache else "full"

if MODE == "standard":
    src = OUT_NPZ if os.path.exists(OUT_NPZ) else GEOMETRY_CACHE
    with np.load(src) as zc:            # close before P.fresh() deletes/rewrites
        w, s, e, n = [float(v) for v in zc["bbox"]]
        osm_seg = zc["osm_seg"]; gf_seg = zc["gf_seg"]; arr = zc["flood"]
    side_km = (e - w) * kx
    geom_src = os.path.basename(src)
else:
    import pyogrio, shapely, rasterio
    OSM = P.get("OSM_POWER"); GF = P.get("GRIDFINDER")
    if len(sys.argv) > 2:                  # zoom on the grid-offline cluster
        side_km = float(sys.argv[2])
        src_pts = (gonly if gonly.sum() >= 5 else (gonly | direct))
        cx = float(np.median(lon[src_pts])); cy = float(np.median(lat[src_pts]))
    else:
        side_km = max((lon.max() - lon.min()) * kx,
                      (lat.max() - lat.min()) * ky) * 1.08 + 3.0
        cx, cy = (lon.min() + lon.max()) / 2, (lat.min() + lat.max()) / 2
    w, e = cx - side_km / kx / 2, cx + side_km / kx / 2
    s, n = cy - side_km / ky / 2, cy + side_km / ky / 2

    def segs(f, layer, cols, filt=None):
        g = pyogrio.read_dataframe(f, layer=layer, bbox=(w, s, e, n),
                                   columns=cols)
        if filt is not None:
            g = g[filt(g)]
        if not len(g):
            return np.zeros((0, 4), np.float32)
        cs = shapely.get_coordinates(g.geometry.values)
        nc = shapely.get_num_coordinates(g.geometry.values)
        idx = np.r_[0, np.cumsum(nc)]
        keep = np.ones(len(cs), bool); keep[idx[1:] - 1] = False
        a = np.where(keep)[0]
        out = np.c_[cs[a], cs[a + 1]].astype(np.float32)
        mid = (out[:, 0] + out[:, 2]) / 2; mya = (out[:, 1] + out[:, 3]) / 2
        return out[(mid >= w) & (mid <= e) & (mya >= s) & (mya <= n)]

    osm_seg = segs(OSM, "power_lines", ["power"])
    gf_seg = segs(GF, "grid", ["source"], lambda d: d["source"] == "gridfinder")
    with rasterio.open(P.jrc_raster("100")) as ds:
        win = rasterio.windows.from_bounds(w, s, e, n, ds.transform)
        arr = ds.read(1, window=win)
        arr = np.where(np.isfinite(arr) & (arr < 1e30) & (arr > 0), 1, 0
                       ).astype(np.uint8)
    geom_src = "sources"

# ---------- site table: window sites + every link-referenced site ----------
in_win = (rep_lon >= w) & (rep_lon <= e) & (rep_lat >= s) & (rep_lat <= n) & \
         (np.bincount(SS.site_id.values, minlength=NSITE_ALL) > 0)
ref = np.zeros(NSITE_ALL, bool)
ref[s1m[lk]] = True
ref[s2m[bk]] = True
keep = in_win | ref
site_ids = np.where(keep)[0]
pos = np.full(NSITE_ALL, -1, np.int64)          # site_id -> row in site arrays
pos[site_ids] = np.arange(len(site_ids))

n_dep = np.zeros(NSITE_ALL, np.int32)           # metro chargers per primary site
np.add.at(n_dep, np.maximum(s1m, 0), (s1m >= 0).astype(np.int32))
n_dep_dry = np.zeros(NSITE_ALL, np.int32)       # unflooded dependents
np.add.at(n_dep_dry, np.maximum(s1m, 0), ((s1m >= 0) & ~direct).astype(np.int32))
n_disp = np.zeros(NSITE_ALL, np.int32)          # displayed offline (P>=0.5)
np.add.at(n_disp, np.maximum(s1m, 0), lk.astype(np.int32))
n_resc = np.zeros(NSITE_ALL, np.int32)          # of those, backup survives
np.add.at(n_resc, np.maximum(s1m, 0), (lk & rescued).astype(np.int32))

PL = pd.read_parquet(os.path.join(FIG4_DATA, "plants_global.parquet"),
                     columns=["lon", "lat"])
plm = (PL.lon.values >= w) & (PL.lon.values <= e) & \
      (PL.lat.values >= s) & (PL.lat.values <= n)

# ---------- routed supply paths over the line graph ----------
# Trace, for every displayed offline charger, the shortest path along the
# gridfinder graph from its assigned feeder site (multi-source over the
# site's member-substation attachment vertices) to the charger's own
# attachment vertex, and aggregate a per-site FLOW (number of chargers
# carried) on every traversed segment. Uses the full graph checkpoints
# when present; otherwise the local city crop
# (checkpoints/graph/graph_crop_<metro>.npz, same content by construction).
GRAPH = os.path.join(FIG4_DATA, "checkpoints", "graph")
CROP = os.path.join(GRAPH, f"graph_crop_{int(METRO)}.npz")
pw, ph = (e - w) * 2.0, (n - s) * 2.0
if os.path.exists(os.path.join(GRAPH, "vertices.npz")):
    zv = np.load(os.path.join(GRAPH, "vertices.npz"))
    vlon_all, vlat_all = zv["lon"], zv["lat"]
    zs = np.load(os.path.join(GRAPH, "segments.npz"))
    ia_all, ib_all, ln_all = zs["ia"], zs["ib"], zs["length_m"]
    gm = (vlon_all >= w - pw) & (vlon_all <= e + pw) & \
         (vlat_all >= s - ph) & (vlat_all <= n + ph)
    vid = np.where(gm)[0].astype(np.int64)
    kseg = gm[ia_all] & gm[ib_all]
    g_lon, g_lat = vlon_all[vid], vlat_all[vid]
    g_ia, g_ib = ia_all[kseg].astype(np.int64), ib_all[kseg].astype(np.int64)
    g_len = ln_all[kseg].astype(np.float32)
    np.savez_compressed(P.fresh(CROP), vid=vid, lon=g_lon, lat=g_lat,
                        ia=g_ia, ib=g_ib, length_m=g_len)
    del zv, zs, vlon_all, vlat_all, ia_all, ib_all, ln_all
else:
    zc2 = np.load(CROP)
    vid, g_lon, g_lat = zc2["vid"], zc2["lon"], zc2["lat"]
    g_ia, g_ib, g_len = zc2["ia"], zc2["ib"], zc2["length_m"]

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as _dijkstra
NV = len(vid)
loc = {int(v): i for i, v in enumerate(vid)}
li = np.array([loc[int(a)] for a in g_ia], np.int64)
lj = np.array([loc[int(b)] for b in g_ib], np.int64)
A = csr_matrix((np.r_[g_len, g_len],
                (np.r_[li, lj], np.r_[lj, li])), shape=(NV, NV))

ATT_CH = np.load(os.path.join(GRAPH, "attach_chargers.npz"))["vid"]
ATT_SUB = np.load(os.path.join(GRAPH, "attach_subs.npz"))["vid"]
ev_idx = np.where(m)[0]                       # metro chargers -> ev rows
ch_vid_m = ATT_CH[ev_idx]                     # attachment vertex per charger
ch_loc = np.array([loc.get(int(v), -1) for v in ch_vid_m], np.int64)
con_ok = (ch_vid_m >= 0) & (ch_loc >= 0)
con_mask = lk & con_ok                        # displayed chargers with a spoke
con_seg = np.c_[lon[con_mask], lat[con_mask],
                g_lon[ch_loc[con_mask]], g_lat[ch_loc[con_mask]]
                ].astype(np.float32)
con_site = s1m[con_mask].astype(np.int64)


def _site_attach_locs(sgid):
    members = SS.sub_id.values[SS.site_id.values == sgid]
    mv = ATT_SUB[members]
    return np.unique([loc[int(v)] for v in mv
                      if int(v) >= 0 and int(v) in loc]).astype(np.int64)


def _blocked_graph(blocked_locs):
    """Drop all edges incident to blocked attachment vertices."""
    if not len(blocked_locs):
        return A
    blocked = np.asarray(blocked_locs, dtype=np.int64)
    coo = A.tocoo()
    keep = ~(np.isin(coo.row, blocked) | np.isin(coo.col, blocked))
    return csr_matrix((coo.data[keep], (coo.row[keep], coo.col[keep])),
                      shape=A.shape)


routed = np.zeros(len(lon), bool)
flow_keys = {}
sitecon = {}
d_routed, d_assign = [], []
for sgid in np.unique(s1m[lk]):
    src = _site_attach_locs(int(sgid))
    dep_ix = np.where(lk & (s1m == sgid) & con_ok)[0]
    if not len(src) or not len(dep_ix):
        continue
    dist, pred, _src = _dijkstra(A, directed=False, indices=src,
                                 return_predecessors=True, min_only=True)
    cnt = {}
    for ci in dep_ix:
        t = int(ch_loc[ci])
        if not np.isfinite(dist[t]):
            continue
        routed[ci] = True
        d_routed.append(float(dist[t]))
        d_assign.append(float(NET.d1_net_m.values[ev_idx[ci]]))
        u = t
        while pred[u] >= 0:
            v2 = int(pred[u])
            kkey = (v2, u) if v2 < u else (u, v2)
            cnt[kkey] = cnt.get(kkey, 0) + 1
            u = v2
    if cnt:
        ke = np.array(list(cnt.keys()), np.int64)
        flow_keys[int(sgid)] = (
            np.c_[g_lon[ke[:, 0]], g_lat[ke[:, 0]],
                  g_lon[ke[:, 1]], g_lat[ke[:, 1]]].astype(np.float32),
            np.array(list(cnt.values()), np.int32))
    # member-substation connectors of this site (the physical tie to the grid)
    members = SS.sub_id.values[SS.site_id.values == sgid]
    mv = ATT_SUB[members]
    okm = (mv >= 0) & np.isin(mv, vid)
    if okm.any():
        mloc = np.array([loc[int(v)] for v in mv[okm]], np.int64)
        SUBXY = S[["sub_id", "slon", "slat"]].drop_duplicates("sub_id") \
            .set_index("sub_id")
        mm = SUBXY.reindex(members[okm])
        sitecon[int(sgid)] = np.c_[mm.slon.values, mm.slat.values,
                                   g_lon[mloc], g_lat[mloc]].astype(np.float32)

if d_routed:
    dr, da = np.array(d_routed), np.array(d_assign)
    dd = np.abs(dr - da)
    print(f"  path check: {routed.sum()}/{int(lk.sum())} displayed chargers "
          f"routed | median |d_path - d1_net| = {np.median(dd):.0f} m "
          f"(connector lengths excluded from the traced path)", flush=True)

# ---------- routed rescue paths (charger -> backup along the lines) --------
# Prefer intact-network s2. If that hang is removed with the failed primary
# (shared attachment) or unreachable in the residual crop, fall back to the
# nearest other live site still reachable on the residual graph.
rescue_routed = np.zeros(len(lon), bool)
ch_bk_draw = np.full(len(lon), -1, np.int64)
rflow_keys = {}
rbk_con = {}
d_rescue = []
sub_ids_all = SS.sub_id.values.astype(np.int64)
vids_all = ATT_SUB[sub_ids_all]
in_crop = np.array([int(v) in loc for v in vids_all], bool)
cand_sites = np.unique(SS.site_id.values[(vids_all >= 0) & in_crop])
cand_sites = cand_sites[cand_sites >= 0]


def _add_bsitecon(bid):
    if bid in rbk_con:
        return
    b_members = SS.sub_id.values[SS.site_id.values == bid]
    bmv = ATT_SUB[b_members]
    okm = (bmv >= 0) & np.isin(bmv, vid)
    if not okm.any():
        return
    mloc = np.array([loc[int(v)] for v in bmv[okm]], np.int64)
    SUBXY = S[["sub_id", "slon", "slat"]].drop_duplicates("sub_id") \
        .set_index("sub_id")
    mm = SUBXY.reindex(b_members[okm])
    rbk_con[bid] = np.c_[mm.slon.values, mm.slat.values,
                         g_lon[mloc], g_lat[mloc]].astype(np.float32)


def _accumulate_path(pred, t, cnt):
    u = int(t)
    while pred[u] >= 0:
        v2 = int(pred[u])
        kkey = (v2, u) if v2 < u else (u, v2)
        cnt[kkey] = cnt.get(kkey, 0) + 1
        u = v2


for sgid in np.unique(s1m[bk]):
    sgid = int(sgid)
    if sgid < 0:
        continue
    blocked = _site_attach_locs(sgid)
    A_res = _blocked_graph(blocked)
    resc_ix = np.where(bk & (s1m == sgid) & con_ok)[0]
    if not len(resc_ix):
        continue
    cnt = {}

    # --- attempt fixed s2 ---
    for bid in np.unique(s2m[resc_ix].astype(np.int64)):
        bid = int(bid)
        dep_ix = resc_ix[s2m[resc_ix] == bid]
        if bid < 0 or bid == sgid:
            continue
        src = _site_attach_locs(bid)
        if len(blocked):
            src = src[~np.isin(src, blocked)]
        if not len(src):
            continue
        dist, pred, _ = _dijkstra(
            A_res, directed=False, indices=src,
            return_predecessors=True, min_only=True)
        for ci in dep_ix:
            t = int(ch_loc[ci])
            if t < 0 or not np.isfinite(dist[t]):
                continue
            rescue_routed[ci] = True
            ch_bk_draw[ci] = bid
            d_rescue.append(float(dist[t]))
            _accumulate_path(pred, t, cnt)
            _add_bsitecon(bid)

    # --- residual fallback for still-unrouted rescued chargers ---
    pending = resc_ix[~rescue_routed[resc_ix]]
    if len(pending):
        src_all, src_site = [], []
        for bid in cand_sites:
            bid = int(bid)
            if bid == sgid or site_p100[bid] >= 0.5:
                continue
            src = _site_attach_locs(bid)
            if len(blocked):
                src = src[~np.isin(src, blocked)]
            for u in src:
                src_all.append(int(u))
                src_site.append(bid)
        if src_all:
            src_all = np.asarray(src_all, np.int64)
            src_site = np.asarray(src_site, np.int64)
            # dedupe source nodes keeping first site label
            _, uniq_ix = np.unique(src_all, return_index=True)
            src_all, src_site = src_all[uniq_ix], src_site[uniq_ix]
            dist, pred, srcs = _dijkstra(
                A_res, directed=False, indices=src_all,
                return_predecessors=True, min_only=True)
            node_to_site = {int(u): int(b)
                            for u, b in zip(src_all, src_site)}
            for ci in pending:
                t = int(ch_loc[ci])
                if t < 0 or not np.isfinite(dist[t]):
                    continue
                snode = int(srcs[t])
                bid = node_to_site.get(snode, -1)
                if bid < 0 and 0 <= snode < len(src_all):
                    # some scipy builds return index into indices=
                    bid = int(src_site[snode])
                    snode = int(src_all[snode])
                if bid < 0:
                    continue
                rescue_routed[ci] = True
                ch_bk_draw[ci] = bid
                d_rescue.append(float(dist[t]))
                _accumulate_path(pred, t, cnt)
                _add_bsitecon(bid)

    if cnt:
        ke = np.array(list(cnt.keys()), np.int64)
        rflow_keys[sgid] = (
            np.c_[g_lon[ke[:, 0]], g_lat[ke[:, 0]],
                  g_lon[ke[:, 1]], g_lat[ke[:, 1]]].astype(np.float32),
            np.array(list(cnt.values()), np.int32))

if d_rescue:
    print(f"  rescue path check: {int(rescue_routed.sum())}/{int(bk.sum())} "
          f"rescued chargers routed on residual graph | median path "
          f"{np.median(d_rescue):.0f} m", flush=True)
elif int(bk.sum()):
    print(f"  rescue path check: 0/{int(bk.sum())} rescued chargers routed "
          f"(no residual path inside crop)", flush=True)

# Ensure residual-draw backup sites are present in the site table for panels.
extra_bk = np.unique(ch_bk_draw[ch_bk_draw >= 0])
missing = extra_bk[~np.isin(extra_bk, site_ids)]
if len(missing):
    site_ids = np.r_[site_ids, missing]

flow_npz = {}
for sgid, (segs_, n_) in flow_keys.items():
    flow_npz[f"flow_{sgid}_seg"] = segs_
    flow_npz[f"flow_{sgid}_n"] = n_
for sgid, segs_ in sitecon.items():
    flow_npz[f"sitecon_{sgid}"] = segs_
for sgid, (segs_, n_) in rflow_keys.items():
    flow_npz[f"rflow_{sgid}_seg"] = segs_
    flow_npz[f"rflow_{sgid}_n"] = n_
for bid, segs_ in rbk_con.items():
    flow_npz[f"bsitecon_{bid}"] = segs_

np.savez_compressed(
    P.fresh(OUT_NPZ),
    routed=routed, rescue_routed=rescue_routed,
    con_seg=con_seg, con_site=con_site,
    **flow_npz,
    metro=METRO, bbox=[w, s, e, n],
    # chargers
    ch_lon=lon, ch_lat=lat, direct=direct, gridonly=gonly,
    rescued=rescued, fallback=fbm,
    p_grid=pg.astype(np.float32), p_grid2=pg2.astype(np.float32),
    p_off_red_new=pred_new.astype(np.float32),
    p_off_red_old=pred_old.astype(np.float32),
    ch_s1=s1m.astype(np.int64), ch_s2=s2m.astype(np.int64),
    ch_bk_draw=ch_bk_draw.astype(np.int64),
    # supply links to the failed primary feeder site (grouped by site)
    link_lon0=lon[lk], link_lat0=lat[lk],
    link_lon1=rep_lon[s1m[lk]], link_lat1=rep_lat[s1m[lk]],
    link_site=s1m[lk].astype(np.int64),
    # rescue links retained for fallback / site-pair bookkeeping
    blink_lon0=lon[bk], blink_lat0=lat[bk],
    blink_lon1=rep_lon[np.maximum(
        np.where(ch_bk_draw[bk] >= 0, ch_bk_draw[bk], s2m[bk]), 0)],
    blink_lat1=rep_lat[np.maximum(
        np.where(ch_bk_draw[bk] >= 0, ch_bk_draw[bk], s2m[bk]), 0)],
    blink_site=np.where(ch_bk_draw[bk] >= 0, ch_bk_draw[bk], s2m[bk]
                        ).astype(np.int64),
    # feeder sites (window + link-referenced)
    site_id=site_ids.astype(np.int64),
    site_lon=rep_lon[site_ids], site_lat=rep_lat[site_ids],
    site_p=site_p100[site_ids].astype(np.float32),
    site_fail=(site_p100[site_ids] >= 0.5),
    site_major=site_major[site_ids],
    site_in_window=((rep_lon[site_ids] >= w) & (rep_lon[site_ids] <= e) &
                    (rep_lat[site_ids] >= s) & (rep_lat[site_ids] <= n)),
    site_n_dep=n_dep[site_ids], site_n_dep_dry=n_dep_dry[site_ids],
    site_n_disp=n_disp[site_ids], site_n_resc=n_resc[site_ids],
    # Attributes of sites within the map window.
    sub_lon=rep_lon[in_win], sub_lat=rep_lat[in_win],
    sub_fail=(site_p100[in_win] >= 0.5),
    sub_p=site_p100[in_win].astype(np.float32),
    # power plants inside the window (OSM + WRI merged set)
    plant_lon=PL.lon.values[plm], plant_lat=PL.lat.values[plm],
    # geometry
    flood=arr, osm_seg=osm_seg, gf_seg=gf_seg,
    stats=[float(direct.sum()), e_off, e_off_red])
nf = int((site_p100[site_ids] >= 0.5).sum())
print(f"metro {int(METRO)} [{MODE}, geometry: {geom_src}]: side {side_km:.0f} km"
      f" | ch {m.sum()} | sites {len(site_ids)} ({nf} failed)"
      f" | links {int(lk.sum())} (+{int(bk.sum())} rescue)"
      f" | plants {int(plm.sum())} | osm {len(osm_seg):,} gf {len(gf_seg):,} segs | flood {arr.shape}"
      f" | direct {int(direct.sum())} E[offline] {e_off:.1f}"
      f" A {e_off/max(direct.sum(),1):.2f}"
      f" -> route-aware {e_off_red/max(direct.sum(),1):.2f}"
      f" | rescued {int(rescued.sum())}/{int(gonly.sum())} displayed"
      f" | rescue_routed {int(rescue_routed.sum())}",
      flush=True)
