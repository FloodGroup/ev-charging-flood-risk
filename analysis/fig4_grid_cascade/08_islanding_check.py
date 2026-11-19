#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, stage 4b -- islanding robustness on the GLOBAL model (SI).

The headline cascade counts a charger as offline when its feeder site
fails. Islanding asks a stricter question: does the charger's grid
component still contain any live SOURCE (an attached power plant, or a
transmission-class site that has not failed)?  Failed sites isolate their
canonical geometric attachment vertex and all incident edges before connected
components are recomputed. Deterministic evaluation at RP100 uses
P_fail >= 0.5.

Output: data/islanding_global.json, data/islanding_country.csv
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
import numpy as np, pandas as pd
from scipy.stats import norm
from scipy.sparse.csgraph import connected_components
from _joint_failure_routing import failed_attachment_nodes, load_graph

THETA_DIST, THETA_TRANS, BETA = 0.5, 2.0, 0.6
R_SERV, R_EUC = 10_000.0, 8_000.0

za = np.load(f"{GRAPH}/attach_plants.npz");   pl_vid  = za["vid"]
za = np.load(f"{GRAPH}/attach_chargers.npz"); ch_vid  = za["vid"]

SS = pd.read_parquet(f"{OUT}/sub_sites.parquet")
S = pd.read_parquet(f"{FIG4_DATA}/substation_depths.parquet")
K = pd.read_parquet(f"{FIG4_DATA}/substations_global.parquet",
                    columns=["sub_id", "kind"])
S = S.merge(K, on="sub_id").merge(SS, on="sub_id")
NSITE_ALL = int(SS.site_id.max()) + 1
site_major = np.zeros(NSITE_ALL, bool)
np.logical_or.at(site_major, S.site_id.values,
                 S.kind.isin(["transmission", "primary"]).values)
site_d100 = np.zeros(NSITE_ALL)
np.maximum.at(site_d100, S.site_id.values, np.nan_to_num(S.d100.values.astype(float)))
th = np.where(site_major, THETA_TRANS, THETA_DIST)
with np.errstate(divide="ignore"):
    site_p100 = np.where(site_d100 > 0,
                         norm.cdf(np.log(np.maximum(site_d100, 1e-12) / th) / BETA), 0.0)
site_fail = site_p100 >= 0.5

# ---- isolate failed attachment vertices and recompute live components ----
A, GM = load_graph(GRAPH)
K_GRAPH = int(GM["K"])
kv = GM["kv"].astype(np.int64)
site_ids_rank = GM["att_sites"].astype(np.int64)
attach_k = GM["attach_k"].astype(np.int64)
attach_rank = GM["attach_rank"].astype(np.int64)
kidx = np.full(int(kv.max()) + 1, -1, np.int64)
kidx[kv] = np.arange(K_GRAPH)
fail_rank = site_fail[site_ids_rank]
failed_k = failed_attachment_nodes(GM, fail_rank)

Apost = A.copy()
for u in failed_k:
    Apost.data[Apost.indptr[u]:Apost.indptr[u + 1]] = 0.0
failed_node = np.zeros(Apost.shape[0], bool)
failed_node[failed_k] = True
Apost.data[failed_node[Apost.indices]] = 0.0
Apost.eliminate_zeros()
npre, comp_pre = connected_components(A, directed=False)
npost, comp_post = connected_components(Apost, directed=False)

source_pre = np.zeros(A.shape[0], bool)
source_post = np.zeros(A.shape[0], bool)
plant_ok = (pl_vid >= 0) & (pl_vid < len(kidx))
plant_k = kidx[pl_vid[plant_ok]]
plant_k = plant_k[plant_k >= 0]
source_pre[plant_k] = True
source_post[plant_k] = True
major_rank = site_major[site_ids_rank]
# Any live major site's attachment vertex counts as a source.
major_attach = np.unique(attach_k[major_rank[attach_rank]])
live_major_attach = np.unique(
    attach_k[major_rank[attach_rank] & ~fail_rank[attach_rank]])
source_pre[major_attach] = True
source_post[live_major_attach] = True
comp_live_pre = np.zeros(npre, bool)
comp_live_post = np.zeros(npost, bool)
comp_live_pre[comp_pre[source_pre]] = True
comp_live_post[comp_post[source_post]] = True

# ---- charger states ----
sd = pd.read_parquet(f"{FIG2_DATA}/station_depths.parquet")
evm = (sd.kind == "ev").values
ev = sd[evm].reset_index(drop=True)
NET = pd.read_parquet(f"{OUT}/charger_network_assignment.parquet")
EUC = pd.read_parquet(f"{FIG4_DATA}/charger_assignment.parquet")
SITE = pd.read_parquet(f"{OUT}/site_table.parquet")
rank2site = np.full(int(SITE.site_rank.max()) + 1, -1, np.int64)
rank2site[SITE.site_rank.values] = SITE.site_id.values
site_of_sub = np.zeros(int(SS.sub_id.max()) + 1, np.int64)
site_of_sub[SS.sub_id.values] = SS.site_id.values

d1 = NET.d1_net_m.values
net_ok = (NET.s1_rank.values >= 0) & np.isfinite(d1) & (d1 <= R_SERV)
fb = ~net_ok & (EUC.dist1_m.values <= R_EUC)
s1_site = np.where(net_ok, rank2site[np.maximum(NET.s1_rank.values, 0)],
                   np.where(fb, site_of_sub[EUC.sub1.values], -1))
direct = ev.d100.values > 0
feeder_fail = (s1_site >= 0) & site_fail[np.maximum(s1_site, 0)]

on_line = (ch_vid >= 0) & (ch_vid < len(kidx))
ch_k = np.where(on_line, kidx[np.maximum(ch_vid, 0)], -1)
on_line &= ch_k >= 0
was_live = on_line & comp_live_pre[comp_pre[np.maximum(ch_k, 0)]]
now_dead = was_live & ~comp_live_post[comp_post[np.maximum(ch_k, 0)]]
islanded = now_dead & ~feeder_fail & ~direct        # additional losses beyond feeder model

tot_service = float((direct | feeder_fail).sum())
tot_conn    = float((direct | feeder_fail | islanded).sum())
dsum = float(direct.sum())
g = dict(direct=int(dsum),
         service_total=tot_service, A_service=tot_service / dsum,
         conn_total=tot_conn, A_conn=tot_conn / dsum,
         islanded_add=int(islanded.sum()),
         delta_A=float((tot_conn - tot_service) / dsum),
         chargers_in_live_component_pct=float(100 * was_live.mean()),
         failed_site_buses=int(fail_rank.sum()),
         failed_attachment_vertices=int(len(failed_k)),
         components_pre=int(npre), components_post=int(npost),
         note=("deterministic at P_fail>=0.5, RP100; failed sites isolate "
               "their canonical attachment vertex and all incident edges"))
json.dump(g, open(P.fresh(f"{OUT}/islanding_global.json"), "w"), indent=1)
print(json.dumps(g, indent=1))

t = pd.DataFrame(dict(country=ev.country.values, direct=direct,
                      svc=direct | feeder_fail,
                      conn=direct | feeder_fail | islanded))
cy = t.groupby("country", observed=True).agg(
    n=("direct", "size"), direct=("direct", "sum"),
    svc=("svc", "sum"), conn=("conn", "sum"))
cy["A_service"] = cy.svc / cy.direct.replace(0, np.nan)
cy["A_conn"] = cy.conn / cy.direct.replace(0, np.nan)
cy["delta"] = cy.A_conn - cy.A_service
cy.to_csv(P.fresh(f"{OUT}/islanding_country.csv"))
top = cy[(cy.n >= 1000) & (cy.direct >= 100)].sort_values("delta", ascending=False)
print(top.head(15).round(3).to_string())
print("DONE")
