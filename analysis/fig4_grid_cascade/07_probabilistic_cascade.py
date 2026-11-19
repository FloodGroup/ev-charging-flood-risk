#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, stage 4 -- fragility-based probabilistic cascade.

Substation-site fragility (lognormal CDF, literature-anchored):
    P_fail(d) = Phi( ln(d / theta) / beta ),  d > 0
    theta = 0.5 m (distribution/low-voltage/unknown)  |  2.0 m (transmission/primary)
    beta  = 0.6
Site depth = max over member substations (dual OSM representations);
site class = transmission/primary if any member is.

Charger states per return period:
    direct_i             deterministic (depth > 0, identical to Figure 2)
    p_grid_i             = P_fail(site1)                    [primary feeder]
    p_grid_redundant_i   = if a same-OSM-kind alternate site B remains
                           reachable within 10 km after isolating the
                           primary's attachment: P_fail(site1)*P_fail(B);
                           otherwise P_fail(site1).  Cross-kind sites are
                           rejected even if a network path exists.
    E[offline]           = sum_i [ direct_i ? 1 : p_grid_i ]

Feeder validity: post-fault network distance <= R_SERV (10 km). Chargers that cannot
reach a site over the lines fall back to the euclidean nearest substation
within 8 km (fallback radius, flagged); the rest are grid-independent.
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
from _joint_failure_routing import BoundedRouter, load_graph

RPS = ["010", "020", "050", "075", "100", "200", "500"]
R_SERV = 10_000.0
R_EUC  = 8_000.0
THETA_DIST, THETA_TRANS, BETA = 0.5, 2.0, 0.6
t0 = time.time()

# ---------------- charger table ----------------
sd = pd.read_parquet(f"{FIG2_DATA}/station_depths.parquet")
ev_mask = (sd.kind == "ev").values
metro = np.concatenate([np.load(f"{FIG2_DATA}/metro_p{i}.npy")
                        for i in range(4)])[ev_mask]
ev = sd[ev_mask].reset_index(drop=True)
NEV = len(ev)

# ---------------- substation sites ----------------
SS = pd.read_parquet(f"{OUT}/sub_sites.parquet")            # sub_id -> site_id
SITE = pd.read_parquet(f"{OUT}/site_table.parquet")          # attached sites (rank space)
S = pd.read_parquet(f"{FIG4_DATA}/substation_depths.parquet")
K = pd.read_parquet(f"{FIG4_DATA}/substations_global.parquet",
                    columns=["sub_id", "kind"])
S = S.merge(K, on="sub_id").merge(SS, on="sub_id")
NSITE_ALL = int(SS.site_id.max()) + 1
site_of_sub = np.zeros(int(SS.sub_id.max()) + 1, np.int64)
site_of_sub[SS.sub_id.values] = SS.site_id.values

sub_major = S.kind.isin(["transmission", "primary"]).values
site_major = np.zeros(NSITE_ALL, bool)
np.logical_or.at(site_major, S.site_id.values, sub_major)
site_depth = {}
for drp in RPS:
    d = np.zeros(NSITE_ALL, np.float64)
    np.maximum.at(d, S.site_id.values, np.nan_to_num(S[f"d{drp}"].values.astype(np.float64)))
    site_depth[drp] = d

# ---------------- feeder resolution ----------------
NET = pd.read_parquet(f"{OUT}/charger_network_assignment.parquet")
EUC = pd.read_parquet(f"{FIG4_DATA}/charger_assignment.parquet")
assert len(NET) == NEV == len(EUC)

rank2site = np.full(int(SITE.site_rank.max()) + 1, -1, np.int64)
rank2site[SITE.site_rank.values] = SITE.site_id.values

d1 = NET.d1_net_m.values; d2 = NET.d2_net_m.values
net_ok = (NET.s1_rank.values >= 0) & np.isfinite(d1) & (d1 <= R_SERV)
s2_geo = net_ok & (NET.s2_rank.values >= 0) & np.isfinite(d2) & (d2 <= R_SERV)

# Site OSM kind from the representative (max-voltage) member.  Alternate
# feeders must match the primary site's kind (voltage-class proxy).
KIND_CODE = {"transmission": 0, "primary": 1, "distribution": 2,
             "lv": 3, "unknown": 4}
subs_kind = pd.read_parquet(
    f"{FIG4_DATA}/substations_global.parquet", columns=["sub_id", "kind"])
kind_by_sub = dict(zip(subs_kind.sub_id.values, subs_kind.kind.values))
site_kind_rank = np.array(
    [KIND_CODE.get(kind_by_sub.get(int(sid), "unknown"), 4)
     for sid in SITE.rep_sub_id.values],
    dtype=np.int16)
assert len(site_kind_rank) == len(SITE)
s1_kind = np.where(
    net_ok, site_kind_rank[np.maximum(NET.s1_rank.values, 0)], -1)
s2_kind = np.where(
    s2_geo, site_kind_rank[np.maximum(NET.s2_rank.values, 0)], -1)
# Geographic second site counts as a fixed-product backup only if kinds match.
s2_ok = s2_geo & (s1_kind == s2_kind) & (s1_kind >= 0)
s1_site = np.where(net_ok, rank2site[np.maximum(NET.s1_rank.values, 0)], -1)
s2_site = np.where(s2_ok, rank2site[np.maximum(NET.s2_rank.values, 0)], -1)

fb = ~net_ok & (EUC.dist1_m.values <= R_EUC)
s1_site = np.where(net_ok, s1_site,
                   np.where(fb, site_of_sub[EUC.sub1.values], -1))
assigned = s1_site >= 0
s1i = np.maximum(s1_site, 0)
s2i = np.maximum(s2_site, 0)

astats = dict(
    n_chargers=int(NEV), n_subs=int(len(S)), n_sites=int(NSITE_ALL),
    n_sites_attached=int(len(SITE)),
    attached_to_line_pct=float(100 * (NET.attach_vid.values >= 0).mean()),
    attach_median_m=float(np.nanmedian(NET.attach_dist_m.values)),
    network_assigned_pct=float(100 * net_ok.mean()),
    euclidean_fallback_pct=float(100 * fb.mean()),
    grid_independent_pct=float(100 * (~assigned).mean()),
    redundant_pct=float(100 * s2_ok.mean()),
    d1_net_median_km=float(np.nanmedian(np.where(net_ok, d1, np.nan)) / 1000),
    d1_net_p95_km=float(np.nanpercentile(np.where(net_ok, d1, np.nan), 95) / 1000),
    r_serv_m=R_SERV, r_euclidean_fallback_m=R_EUC,
    theta_dist_m=THETA_DIST, theta_trans_m=THETA_TRANS, beta=BETA,
    site_radius_m=100.0)
json.dump(astats, open(P.fresh(f"{OUT}/assignment_stats.json"), "w"), indent=1)
print(json.dumps(astats, indent=1), flush=True)

# ---------------- fragility ----------------
def p_fail(drp, theta_scale=1.0, beta=BETA, threshold=False):
    d = site_depth[drp]
    if threshold:
        return np.where(site_major, d >= 2.0, d >= 0.5).astype(np.float64)
    th = np.where(site_major, THETA_TRANS, THETA_DIST) * theta_scale
    with np.errstate(divide="ignore"):
        z = np.log(np.maximum(d, 1e-12) / th) / beta
    return np.where(d > 0, norm.cdf(z), 0.0)

# ---------------- alternate-site routing after primary isolation ----------
att_ch = np.load(f"{GRAPH}/attach_chargers.npz")
graph_csr, graph_meta = load_graph(GRAPH)
K_GRAPH = int(graph_meta["K"])
NSRC_GRAPH = int(graph_meta["NSRC"])
kv = graph_meta["kv"].astype(np.int64)
kidx = np.full(int(kv.max()) + 1, -1, np.int64)
kidx[kv] = np.arange(K_GRAPH)
ch_vid = att_ch["vid"].astype(np.int64)
ok_attach = (ch_vid >= 0) & (ch_vid < len(kidx))
charger_node = np.full(NEV, -1, np.int64)
charger_node[ok_attach] = kidx[ch_vid[ok_attach]]
router = BoundedRouter(
    graph_csr, K_GRAPH,
    attach_k=graph_meta["attach_k"].astype(np.int64),
    attach_rank=graph_meta["attach_rank"].astype(np.int64))
assert NSRC_GRAPH == len(SITE)

topo_cache = {}
topo_audit_rows = []


def _any_alternate_after_primary_fail(charger_idx, failed_ranks):
    """Search a same-kind live site after isolating each failed primary."""
    out = np.zeros(len(charger_idx), bool)
    ranks = np.full(len(charger_idx), -1, np.int64)
    if not len(charger_idx):
        return out, ranks, dict(groups=0, visited=0, queries=0,
                                unique_starts=0)
    order = np.argsort(failed_ranks, kind="stable")
    fr = failed_ranks[order]
    starts = np.r_[0, np.where(np.diff(fr) != 0)[0] + 1, len(fr)]
    groups = visited = unique_starts = 0
    attach_m = NET.attach_dist_m.values.astype(np.float64)
    for a, b in zip(starts[:-1], starts[1:]):
        pos = order[a:b]
        fail_r = int(fr[a])
        fail_rank = np.zeros(NSRC_GRAPH, bool)
        fail_rank[fail_r] = True
        req_kind = int(site_kind_rank[fail_r])
        rr, rk, rst = router.route_any_alternate(
            charger_idx[pos], charger_node, attach_m, fail_rank,
            site_kind=site_kind_rank, required_kind=req_kind)
        out[pos] = rr
        ranks[pos] = rk
        groups += int(rst.get("groups", 0))
        visited += int(rst.get("visited", 0))
        unique_starts += int(rst.get("unique_starts", 0))
    return out, ranks, dict(groups=groups, visited=visited,
                            queries=int(len(charger_idx)),
                            unique_starts=unique_starts)


def independent_feeder_topology(drp, theta_scale=1.0, beta=BETA,
                                threshold=False):
    """Primary failure + residual same-kind alternate with P1×P2 product.

    When a same-kind alternate B remains reachable within 10 km after
    isolating the primary's attachment, offline probability is
    P_fail(primary)×P_fail(B).  Otherwise it equals P_fail(primary).
    The intact-network second-nearest same-kind site remains the analytic
    control (fixed P1×P2 without residual rerouting).
    """
    tmc = time.time()
    key = (str(drp), float(theta_scale), float(beta), bool(threshold))
    if key in topo_cache:
        return topo_cache[key]

    pf = p_fail(drp, theta_scale, beta, threshold)
    direct = ev[f"d{drp}"].values > 0
    dsum = float(direct.sum())
    p1 = np.where(assigned, pf[s1i], 0.0)
    p2 = np.where(s2_ok, pf[s2i], 1.0)
    p_old = np.where(direct, 1.0, p1 * p2)

    s1r = NET.s1_rank.values.astype(np.int64)
    ci = np.where((~direct) & net_ok & (s1r >= 0))[0]
    alt_ok = np.zeros(NEV, bool)
    alt_rank = np.full(NEV, -1, np.int64)
    reroute_queries = route_groups = route_visited = unique_starts = 0
    n_no_alternate = 0

    if len(ci):
        # Any other intact site within 10 km is necessary; kind match is
        # enforced inside the residual search (S2 may be wrong-kind while a
        # farther same-kind site still lies inside the radius).
        need = ci[s2_geo[ci]]
        no_s2 = ci[~s2_geo[ci]]
        alt_ok[no_s2] = False
        if len(need):
            found, found_rank, rst = _any_alternate_after_primary_fail(
                need, s1r[need])
            alt_ok[need] = found
            alt_rank[need] = found_rank
            reroute_queries += rst["queries"]
            route_groups += rst["groups"]
            route_visited += rst["visited"]
            unique_starts += rst["unique_starts"]
        n_no_alternate = int((~alt_ok[ci]).sum())

    # Product of independent feeder failures when a residual alternate exists;
    # otherwise the charger has only the primary (multiply by 1).
    p2_dyn = np.ones(NEV, np.float64)
    hit = alt_ok & (alt_rank >= 0)
    if hit.any():
        alt_site = rank2site[alt_rank[hit]]
        assert np.all(alt_site >= 0)
        p2_dyn[hit] = pf[alt_site]
    p_new = np.where(direct, 1.0, p1 * p2_dyn)
    p_new = np.clip(p_new, 0.0, 1.0)
    A_new = float(p_new.sum() / dsum)
    A_old = float(p_old.sum() / dsum)
    result = dict(
        marginal_new=p_new,
        marginal_old_mc=p_old,
        samples=0,
        seed=None,
        converged=True,
        exact=True,
        independent_feeder_failures=True,
        dynamic_alternate_search=True,
        kind_matched_backup=True,
        product_with_residual_alternate=True,
        fixed_backup=False,
        A_new=A_new,
        A_new_se=0.0,
        A_new_ci95=[A_new, A_new],
        A_old_mc=A_old,
        A_old_mc_se=0.0,
        reroute_queries=int(reroute_queries),
        route_groups=int(route_groups),
        route_visited_nodes=int(route_visited),
        unique_start_nodes=int(unique_starts),
        n_no_alternate=n_no_alternate,
        n_backup_topology_blocked=n_no_alternate,
        n_primary_topology_blocked=0,
        elapsed_seconds=float(time.time() - tmc),
        seconds_per_sample=None)
    topo_audit_rows.append(dict(
        rp=int(drp), theta_scale=float(theta_scale), beta=float(beta),
        threshold=bool(threshold), A_new=A_new, A_old=A_old,
        delta_A=float(A_new - A_old),
        reroute_queries=int(reroute_queries),
        n_no_alternate=n_no_alternate,
        elapsed_seconds=result["elapsed_seconds"]))
    print(f"exact kind-matched alt RP{drp}: A_backup={A_new:.5f} "
          f"(fixed-site control {A_old:.5f}, Δ={A_new - A_old:+.5f}) "
          f"searches={reroute_queries:,} no_alt={n_no_alternate:,}",
          flush=True)
    topo_cache[key] = result
    return result


# Entry point for single-period evaluation.
joint_backup_mc = independent_feeder_topology

if os.environ.get("EVFR_MC_PILOT") == "1":
    pilot = independent_feeder_topology("100")
    json.dump(
        {k: v for k, v in pilot.items()
         if k not in {"marginal_new", "marginal_old_mc"}},
        open(P.fresh(f"{OUT}/cascade_mc_pilot.json"), "w"), indent=1)
    print("INDEPENDENT-FEEDER TOPOLOGY PILOT DONE", flush=True)
    raise SystemExit(0)

# ---------------- cascade ----------------
def scenario(drp, theta_scale=1.0, beta=BETA, threshold=False, redundant=False):
    pf = p_fail(drp, theta_scale, beta, threshold)
    direct = ev[f"d{drp}"].values > 0
    pg = np.where(assigned, pf[s1i], 0.0)
    if redundant:
        mc = joint_backup_mc(drp, theta_scale, beta, threshold)
        p_off = mc["marginal_new"]
        pg = np.where(direct, 1.0, p_off)
    else:
        p_off = np.where(direct, 1.0, pg)
    return direct, pg, p_off

def summarize(direct, p_off):
    d = float(direct.sum()); tot = float(p_off.sum())
    return dict(direct=int(d), grid_add=float(tot - d), total=float(tot),
                A_total=float(tot / d))

def old_redundant_analytic(drp):
    """Analytical control with a fixed independent alternate feeder."""
    pf = p_fail(drp)
    direct = ev[f"d{drp}"].values > 0
    pg = np.where(
        assigned,
        pf[s1i] * np.where(s2_ok, pf[s2i], 1.0),
        0.0)
    return np.where(direct, 1.0, pg)

out = {}
for drp in RPS:
    dv, _, p = scenario(drp);                 out[f"fragility_RP{drp}"] = summarize(dv, p)
    dv, _, p = scenario(drp, redundant=True); out[f"fragility_RP{drp}_redundant"] = summarize(dv, p)
    dv, _, p = scenario(drp, threshold=True); out[f"threshold_RP{drp}"] = summarize(dv, p)
for tag, kw in [("beta04", dict(beta=0.4)), ("beta08", dict(beta=0.8)),
                ("soft", dict(theta_scale=0.75)), ("hard", dict(theta_scale=1.5))]:
    dv, _, p = scenario("100", **kw)
    out[f"sens_{tag}_RP100"] = summarize(dv, p)
out["meta"] = dict(
    **astats,
    fragility="lognormal CDF, theta=0.5/2.0 m, beta=0.6, site-level",
    redundant_model=(
        "after primary-site failure, isolate its canonical attachment "
        "vertex; if any other live site B of the SAME OSM kind remains "
        "reachable within 10 km, P_offline = P1×P_fail(B); otherwise "
        "P_offline = P1. The fixed-site control uses the intact second-nearest "
        "site only when kinds match"))
json.dump(out, open(P.fresh(f"{OUT}/cascade_global.json"), "w"), indent=1)

rows = []
for drp in RPS:
    r = dict(rp=int(drp))
    dv, _, p = scenario(drp)
    r["direct"] = int(dv.sum()); r["total"] = float(p.sum())
    r["A_fragility"] = float(p.sum() / dv.sum())
    _, _, p2 = scenario(drp, redundant=True)
    r["A_redundant"] = float(p2.sum() / dv.sum()); r["total_redundant"] = float(p2.sum())
    po = old_redundant_analytic(drp)
    r["A_redundant_old"] = float(po.sum() / dv.sum())
    mc = independent_feeder_topology(drp)
    r["A_redundant_old_mc"] = float(mc["A_old_mc"])
    r["mc_samples"] = int(mc["samples"])
    r["mc_se"] = float(mc["A_new_se"])
    for col, kw in [("A_threshold", dict(threshold=True)),
                    ("A_soft", dict(theta_scale=0.75)), ("A_hard", dict(theta_scale=1.5)),
                    ("A_beta04", dict(beta=0.4)), ("A_beta08", dict(beta=0.8))]:
        _, _, px = scenario(drp, **kw)
        r[col] = float(px.sum() / dv.sum())
    rows.append(r)
rp_df = pd.DataFrame(rows)
rp_df.to_csv(P.fresh(f"{OUT}/cascade_rp_curve.csv"), index=False)
pd.DataFrame(topo_audit_rows).to_csv(
    P.fresh(f"{OUT}/cascade_mc_convergence.csv"), index=False)
mc_summary = {
    f"RP{k[0]}": {kk: vv for kk, vv in val.items()
                  if kk not in {"marginal_new", "marginal_old_mc"}}
    for k, val in topo_cache.items()
    if k[1:] == (1.0, BETA, False)
}
json.dump(mc_summary, open(P.fresh(f"{OUT}/cascade_mc_summary.json"), "w"),
          indent=1)
print(rp_df.round(3).to_string(index=False), flush=True)

# ---------------- aggregations at RP100 ----------------
def agg(groups, direct, p_off, p_off_red):
    g = pd.DataFrame(dict(grp=groups, direct=direct.astype(float),
                          total=p_off, total_red=p_off_red))
    r = g.groupby("grp", observed=True).agg(
        n=("direct", "size"), direct=("direct", "sum"),
        total=("total", "sum"), total_redundant=("total_red", "sum"))
    r["direct"] = r.direct.astype(int)
    r["grid_add"] = r.total - r.direct
    r["A_total"] = r.total / r.direct.replace(0, np.nan)
    r["A_redundant"] = r.total_redundant / r.direct.replace(0, np.nan)
    r["direct_pct"] = 100 * r.direct / r.n
    return r

dv, _, p = scenario("100")
_, _, pr = scenario("100", redundant=True)
d5, _, p5 = scenario("500")
mc100 = independent_feeder_topology("100")
po100 = old_redundant_analytic("100")
pd.DataFrame(dict(
    direct=dv,
    p_off_primary=p.astype(np.float32),
    p_off_redundant_new=pr.astype(np.float32),
    p_off_redundant_old=po100.astype(np.float32),
    p_off_redundant_old_mc=mc100["marginal_old_mc"].astype(np.float32),
    topology_delta=(pr - po100).astype(np.float32),
    topology_delta_mc=(
        pr - mc100["marginal_old_mc"]).astype(np.float32))
).to_parquet(P.fresh(f"{OUT}/charger_cascade_rp100.parquet"))

cy = agg(ev.country.values, dv, p, pr)
cy5 = agg(ev.country.values, d5, p5, p5)[["total"]].rename(
    columns={"total": "total_rp500"})
sc = pd.read_parquet(f"{FIG4_DATA}/substation_country.parquet")
cont_map = sc[sc.country != ""].groupby("country").continent.agg(
    lambda x: x.mode().iloc[0])
cy["continent"] = [cont_map.get(c, "") for c in cy.index]
cy.join(cy5).to_csv(P.fresh(f"{OUT}/cascade_country_v2.csv"))
cont = pd.Series(ev.country.values).map(lambda c: cont_map.get(c, "")).values
agg(cont, dv, p, pr).to_csv(P.fresh(f"{OUT}/cascade_continent.csv"))
mm = np.where(metro >= 0, metro, -1)
mt = agg(mm, dv, p, pr)
mt = mt[mt.index != -1]
cm = pd.read_csv(f"{FIG2_DATA}/city_map_exposure.csv")[
    ["metro", "UC_NM_MN", "CTR_MN_NM"]]
mt.join(cm.set_index("metro"), how="left").to_csv(P.fresh(f"{OUT}/cascade_metro_named.csv"))

m = cy[(cy.n >= 1000) & (cy.direct >= 100)].dropna(subset=["A_total"])
corr = float(np.corrcoef(np.log(m.direct_pct), m.A_total)[0, 1])
r100 = rp_df.loc[rp_df.rp == 100].iloc[0]
stats = dict(
    A_rp100=float(r100.A_fragility),
    direct_rp100=int(r100.direct),
    total_rp100=float(r100.total),
    grid_add_rp100=float(r100.total - r100.direct),
    A_redundant_rp100=float(r100.A_redundant),
    A_redundant_old_rp100=float(r100.A_redundant_old),
    A_redundant_old_mc_rp100=float(r100.A_redundant_old_mc),
    A_redundant_delta_rp100=float(
        r100.A_redundant - r100.A_redundant_old),
    mc_samples_rp100=int(r100.mc_samples),
    mc_se_rp100=float(r100.mc_se),
    independent_feeder_failures=True,
    exact_topology_backup=True,
    dynamic_alternate_search=True,
    kind_matched_backup=True,
    fixed_backup=False,
    A_threshold_rp100=float(r100.A_threshold),
    A_band_rp100=[float(min(r100.A_hard, r100.A_beta04, r100.A_beta08)),
                  float(max(r100.A_soft, r100.A_beta04, r100.A_beta08))],
    mech_corr_logdirect_A=corr, n_countries_mech=int(len(m)),
    assigned_pct=float(100 * assigned.mean()),
    network_assigned_pct=astats["network_assigned_pct"],
    redundant_pct=astats["redundant_pct"])
json.dump(stats, open(P.fresh(f"{OUT}/cascade_stats.json"), "w"), indent=1)
print(json.dumps(stats, indent=1), flush=True)
print(f"DONE {time.time()-t0:.0f}s", flush=True)
