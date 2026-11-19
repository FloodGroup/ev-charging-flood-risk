#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify network attachment, feeder failure, the fixed-site analytical control and the residual-network counterfactual against the local analysis outputs."""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402


FIG4_DATA = P.get("FIG4_DATA")
FIG4_DIR = P.get("FIG4_DIR")
FIG2_DATA = P.get("FIG2_DATA")
GRAPH = os.path.join(FIG4_DATA, "checkpoints", "graph")
THETA_DIST, THETA_TRANS, BETA = 0.5, 2.0, 0.6
R_SERV, R_EUC = 10_000.0, 8_000.0

fails = []


def check(name, condition, detail=""):
    good = bool(condition)
    print(f"{'OK  ' if good else 'FAIL'} {name}"
          f"{': ' + str(detail) if detail != '' else ''}")
    if not good:
        fails.append(name)


def close(name, actual, expected, atol=1e-8):
    check(name, np.isclose(actual, expected, atol=atol, rtol=0),
          f"{actual} vs {expected}")


# ---------- topology-construction invariants (canonical site attachment) ----------
TA = json.load(open(os.path.join(FIG4_DATA, "site_topology_audit.json")))
check("canonical nearest-vertex attachment mode",
      TA.get("attachment_mode") == "canonical_nearest_vertex_hang")
check("one failure-cut hang per attached site",
      TA["n_failure_cut_attachments"] == TA["n_sites_attached"])
check("raw member links at least one per site",
      TA.get("n_member_attachments_raw", 0) >= TA["n_sites_attached"])
check("all geometric connectors within 10 km",
      TA["attachment_max_m"] <= R_SERV)
check("canonical attachments within 10 km",
      TA["canonical_attachment_max_m"] <= R_SERV)

CG = np.load(os.path.join(GRAPH, "contracted.npz"))
check("contracted source count matches attached sites",
      int(CG["NSRC"]) == TA["n_sites_attached"])
check("canonical buses survive contraction",
      np.all(CG["site_bus_k"] >= 0) and
      len(CG["site_bus_k"]) == int(CG["NSRC"]))
check("failure-cut arrays match one hang per site",
      len(CG["attach_k"]) == len(CG["attach_rank"]) ==
      TA["n_failure_cut_attachments"] == int(CG["NSRC"]))
check("attach_k equals site_bus_k (canonical-only policy)",
      np.array_equal(CG["attach_k"].astype(np.int64),
                     CG["site_bus_k"].astype(np.int64)))
check("attachment ranks are in range",
      np.all((CG["attach_rank"] >= 0) & (CG["attach_rank"] < int(CG["NSRC"]))))
check("contracted CSR is structurally valid",
      len(CG["indptr"]) == int(CG["K"]) + int(CG["NSRC"]) + 1 and
      int(CG["indptr"][-1]) == len(CG["indices"]) == len(CG["data"]))


# ---------- independent RP100 primary and fixed-site control recomputation ----------
S3 = json.load(open(os.path.join(FIG4_DATA, "cascade_stats.json")))
AS = json.load(open(os.path.join(FIG4_DATA, "assignment_stats.json")))
IS = json.load(open(os.path.join(FIG4_DATA, "islanding_global.json")))

sd = pd.read_parquet(os.path.join(FIG2_DATA, "station_depths.parquet"))
ev_mask = (sd.kind == "ev").values
metro = np.concatenate(
    [np.load(os.path.join(FIG2_DATA, f"metro_p{i}.npy")) for i in range(4)]
)[ev_mask]
ev = sd[ev_mask].reset_index(drop=True)

SS = pd.read_parquet(os.path.join(FIG4_DATA, "sub_sites.parquet"))
SITE = pd.read_parquet(os.path.join(FIG4_DATA, "site_table.parquet"))
SUB = pd.read_parquet(
    os.path.join(FIG4_DATA, "substation_depths.parquet")).merge(
        pd.read_parquet(
            os.path.join(FIG4_DATA, "substations_global.parquet"),
            columns=["sub_id", "kind"]), on="sub_id").merge(SS, on="sub_id")
NSITE_ALL = int(SS.site_id.max()) + 1
site_major = np.zeros(NSITE_ALL, bool)
np.logical_or.at(site_major, SUB.site_id.values,
                 SUB.kind.isin(["transmission", "primary"]).values)
site_d100 = np.zeros(NSITE_ALL)
np.maximum.at(site_d100, SUB.site_id.values,
              np.nan_to_num(SUB.d100.values.astype(float)))
theta = np.where(site_major, THETA_TRANS, THETA_DIST)
with np.errstate(divide="ignore"):
    pf = np.where(
        site_d100 > 0,
        norm.cdf(np.log(np.maximum(site_d100, 1e-12) / theta) / BETA),
        0.0)

NET = pd.read_parquet(
    os.path.join(FIG4_DATA, "charger_network_assignment.parquet"))
EUC = pd.read_parquet(os.path.join(FIG4_DATA, "charger_assignment.parquet"))
rank2site = np.full(int(SITE.site_rank.max()) + 1, -1, np.int64)
rank2site[SITE.site_rank.values] = SITE.site_id.values
site_of_sub = np.zeros(int(SS.sub_id.max()) + 1, np.int64)
site_of_sub[SS.sub_id.values] = SS.site_id.values
d1, d2 = NET.d1_net_m.values, NET.d2_net_m.values
net_ok = ((NET.s1_rank.values >= 0) & np.isfinite(d1) & (d1 <= R_SERV))
s2_geo = (net_ok & (NET.s2_rank.values >= 0) &
          np.isfinite(d2) & (d2 <= R_SERV))
KIND_CODE = {"transmission": 0, "primary": 1, "distribution": 2,
             "lv": 3, "unknown": 4}
subs_kind = pd.read_parquet(
    os.path.join(FIG4_DATA, "substations_global.parquet"),
    columns=["sub_id", "kind"])
kind_by_sub = dict(zip(subs_kind.sub_id.values, subs_kind.kind.values))
site_kind_rank = np.array(
    [KIND_CODE.get(kind_by_sub.get(int(sid), "unknown"), 4)
     for sid in SITE.rep_sub_id.values],
    dtype=np.int16)
s1_kind = np.where(
    net_ok, site_kind_rank[np.maximum(NET.s1_rank.values, 0)], -1)
s2_kind = np.where(
    s2_geo, site_kind_rank[np.maximum(NET.s2_rank.values, 0)], -1)
s2_ok = s2_geo & (s1_kind == s2_kind) & (s1_kind >= 0)
fallback = ~net_ok & (EUC.dist1_m.values <= R_EUC)
s1_site = np.where(
    net_ok, rank2site[np.maximum(NET.s1_rank.values, 0)],
    np.where(fallback, site_of_sub[EUC.sub1.values], -1))
s2_site = np.where(
    s2_ok, rank2site[np.maximum(NET.s2_rank.values, 0)], -1)
direct = ev.d100.values > 0
pg = np.where(s1_site >= 0, pf[np.maximum(s1_site, 0)], 0.0)
p_primary = np.where(direct, 1.0, pg)
p_old = np.where(
    direct, 1.0,
    pg * np.where(s2_ok, pf[np.maximum(s2_site, 0)], 1.0))

close("stored direct count", int(direct.sum()), S3["direct_rp100"], .5)
close("primary expected total", float(p_primary.sum()),
      S3["total_rp100"], 1e-5)
close("primary amplification", float(p_primary.sum() / direct.sum()),
      S3["A_rp100"], 1e-10)
close("fixed-site control", float(p_old.sum() / direct.sum()),
      S3["A_redundant_old_rp100"], 1e-10)
check("primary and backup sites are distinct",
      np.all(NET.s1_rank.values[s2_ok] != NET.s2_rank.values[s2_ok]))
close("network-assigned percentage", 100 * net_ok.mean(),
      AS["network_assigned_pct"], 1e-9)
close("redundant percentage", 100 * s2_ok.mean(),
      AS["redundant_pct"], 1e-9)


# ---------- alternate-site topology checks ----------
CAS = pd.read_parquet(
    os.path.join(FIG4_DATA, "charger_cascade_rp100.parquet"))
pnew = CAS.p_off_redundant_new.values.astype(float)
pold_mc = CAS.p_off_redundant_old_mc.values.astype(float)
check("alternate-search never exceeds primary-only loss",
      pnew.sum() <= p_primary.sum() + 1e-9)
check("at least one charger differs from the fixed-site control",
      np.any(np.abs(CAS.topology_delta.values) > 0))
close("stored route-aware amplification", float(pnew.sum() / direct.sum()),
      S3["A_redundant_rp100"], 1e-8)
close("stored topology delta",
      float((pnew - p_old).sum() / direct.sum()),
      S3["A_redundant_delta_rp100"], 1e-8)
close("independent-site control matches the fixed-site marginal",
      float(pold_mc.sum() / direct.sum()),
      S3["A_redundant_old_rp100"], 1e-10)
check("exact alternate-search evaluation (no MC sampling)",
      S3.get("exact_topology_backup", S3["mc_samples_rp100"] == 0) and
      S3["mc_se_rp100"] == 0.0)
check("dynamic alternate search enabled",
      S3.get("dynamic_alternate_search", False))
check("kind-matched backup enabled",
      S3.get("kind_matched_backup", False))

MC = json.load(open(os.path.join(FIG4_DATA, "cascade_mc_summary.json")))
check("all seven return periods evaluated", len(MC) == 7, list(MC))
for rp, row in sorted(MC.items()):
    check(f"{rp} exact evaluation flag", row.get("exact", True))
    check(f"{rp} zero sampling error", row["A_new_se"] == 0.0)
    check(f"{rp} dynamic alternate flag",
          row.get("dynamic_alternate_search", False))
    check(f"{rp} kind-matched backup flag",
          row.get("kind_matched_backup", False))


# ---------- optional path-certificate integrity ----------
cert_path = os.path.join(GRAPH, "backup_path_certificates.npz")
if os.path.exists(cert_path):
    CERT = np.load(cert_path)
    ps = CERT["path_sites"]
    check("certificate ranks are in range",
          np.all((ps == -1) | ((ps >= 0) & (ps < int(CG["NSRC"])))))
else:
    check("path certificates absent (optional)", True)


# ---------- deterministic islanding and vignette consistency ----------
check("islanding removed failed sites",
      IS["failed_site_buses"] == int(
          (pf[CG["att_sites"].astype(np.int64)] >= 0.5).sum()))
check("islanding isolated attachment vertices",
      IS.get("failed_attachment_vertices", 0) > 0 or
      IS["failed_site_buses"] == 0)
check("post-fault graph has at least as many components",
      IS["components_post"] >= IS["components_pre"])
check("islanding delta is nonnegative", IS["delta_A"] >= 0)

for mid, name in [(2288, "Dusseldorf"), (12703, "Tainan"),
                  (10912, "Handan"), (1048, "Quebec City")]:
    m = metro == mid
    z = np.load(os.path.join(FIG4_DATA, f"vignette_v2_{mid}.npz"))
    close(f"{name} direct", float(z["stats"][0]), int(direct[m].sum()), .5)
    close(f"{name} primary expected", float(z["stats"][1]),
          float(p_primary[m].sum()), 5e-2)
    close(f"{name} route-aware expected", float(z["stats"][2]),
          float(pnew[m].sum()), 5e-2)
    check(f"{name} rescued subset", np.all(
        ~z["rescued"].astype(bool) | z["gridonly"].astype(bool)))


# ---------- generated artifacts and synthetic regression ----------
fig_dirs = [FIG4_DIR] + ([os.path.join(P.FIGURES_ROOT, "Figure4")]
                         if P.FIGURES_ROOT else [])
rendered = [os.path.join(d, "Figure4.png") for d in fig_dirs
            if os.path.exists(os.path.join(d, "Figure4.png"))]
if rendered:
    check("rendered Figure 4 has content", os.path.getsize(rendered[0]) > 0, rendered[0])

test_script = os.path.join(os.path.dirname(__file__),
                           "13_test_site_failure_topology.py")
proc = subprocess.run([sys.executable, test_script],
                      capture_output=True, text=True)
check("synthetic radial/meshed topology regression",
      proc.returncode == 0, proc.stdout.strip() or proc.stderr.strip())

if fails:
    print(f"RESULT: {len(fails)} check(s) FAILED: {fails}")
    raise SystemExit(1)
print("ALL SITE-AWARE CHECKS DONE")
