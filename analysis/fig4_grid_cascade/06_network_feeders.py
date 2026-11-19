#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4, stage 3 -- network feeder assignment over the grid graph.

Substation SITES: OSM maps many substations twice (point + polygon
centroid). Substations within 100 m are clustered into one physical site
(union-find); redundancy later requires two DISTINCT sites, so dual
representations can never masquerade as a backup feed.

Graph: contracted line graph (degree-2 chains folded, every attachment
vertex kept) + one virtual node per attached site, linked ONLY to that
site's canonical hang: the member with the smallest geometric attachment
distance to the line network (not every member hang, and not "min total
path among multiple entrances").

When a site later fails, that single canonical attachment vertex is
isolated together with all incident contracted edges (see BoundedRouter /
islanding).

Feeder search (exact, all C-speed scipy):
  pass 0      : multi-source Dijkstra over all sites -> D1, S1 labels
  bit passes  : sites split by bit b of their rank; the half that excludes
                S1(node) bounds the best OTHER site; min over bits = exact
                second-nearest DISTINCT site (any two ranks differ in a bit).

Output: FIG4_DATA/checkpoints/graph/contracted.npz, node_feeders.npz
        FIG4_DATA/sub_sites.parquet         sub_id -> site
        FIG4_DATA/site_table.parquet        site attributes
        FIG4_DATA/charger_network_assignment.parquet
        FIG4_DATA/plant_network_assignment.parquet
        FIG4_DATA/site_topology_audit.json
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

FIG4_DATA = P.get("FIG4_DATA")
FIG2_DATA = P.get("FIG2_DATA")
GRAPH = P.ensure_dir(os.path.join(FIG4_DATA, "checkpoints", "graph"))
OUT = FIG4_DATA
import time, gc, json
import numpy as np, pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

DJ_LIMIT = 60_000.0            # m network search horizon
SITE_RADIUS = 100.0            # m substation site clustering
ATTACH_MAX = 10_000.0          # m, same cap as step 05
RE = 6_371_000.0
t0 = time.time()

zv = np.load(f"{GRAPH}/vertices.npz")
deg = zv["deg"]; NV = len(deg)
zs = np.load(f"{GRAPH}/segments.npz")
ia, ib, slen = zs["ia"], zs["ib"], zs["length_m"].astype(np.float64)
original_line_length_m = float(slen.sum())

za = np.load(f"{GRAPH}/attach_subs.npz");     sub_vid, sub_d = za["vid"], za["dist_m"]
za = np.load(f"{GRAPH}/attach_plants.npz");   pl_vid,  pl_d  = za["vid"], za["dist_m"]
za = np.load(f"{GRAPH}/attach_chargers.npz"); ch_vid,  ch_d  = za["vid"], za["dist_m"]

# ---------------- substation sites (union-find over 100 m pairs) ----------
subs = pd.read_parquet(f"{FIG4_DATA}/substations_global.parquet",
                       columns=["sub_id", "lon", "lat", "volt_max", "kind"])
lo, la = np.deg2rad(subs.lon.values), np.deg2rad(subs.lat.values)
xyz = np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
pairs = cKDTree(xyz).query_pairs(2 * np.sin(SITE_RADIUS / RE / 2),
                                 output_type="ndarray")
parent = np.arange(len(subs))
def find(x):
    root = x
    while parent[root] != root:
        root = parent[root]
    while parent[x] != root:
        parent[x], x = root, parent[x]
    return root
for a, b in pairs:
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra
root = np.array([find(i) for i in range(len(subs))])
_, site_of_sub = np.unique(root, return_inverse=True)
NSITE_ALL = site_of_sub.max() + 1
print(f"substations {len(subs):,} -> sites {NSITE_ALL:,} "
      f"({time.time()-t0:.0f}s)", flush=True)
pd.DataFrame(dict(sub_id=subs.sub_id.values, site_id=site_of_sub)
             ).to_parquet(f"{OUT}/sub_sites.parquet")

# site attributes (class from best-documented member; representative = max volt)
sub_major = subs.kind.isin(["transmission", "primary"]).values
site_major = np.zeros(NSITE_ALL, bool)
np.logical_or.at(site_major, site_of_sub, sub_major)
volt = np.nan_to_num(subs.volt_max.values, nan=-1.0)
order = np.lexsort((-volt, site_of_sub))
first = np.r_[True, np.diff(site_of_sub[order]) != 0]
rep = np.zeros(NSITE_ALL, np.int64)
rep[site_of_sub[order][first]] = subs.sub_id.values[order][first]
site_n = np.bincount(site_of_sub, minlength=NSITE_ALL)

# ---------------- contraction ----------------
attach_mask = np.zeros(NV, bool)
for v in (sub_vid, pl_vid, ch_vid):
    attach_mask[v[v >= 0]] = True
keep = (deg != 2) | attach_mask
print(f"vertices {NV:,}  keep {keep.sum():,} ({100*keep.mean():.1f}%)  "
      f"{time.time()-t0:.0f}s", flush=True)

kk   = keep[ia] & keep[ib]
bnd  = keep[ia] ^ keep[ib]
inte = ~keep[ia] & ~keep[ib]

int_ids = np.where(~keep)[0]
ridx = np.full(NV, -1, np.int64)
ridx[int_ids] = np.arange(len(int_ids))
gi = coo_matrix((np.ones(int(inte.sum()), np.int8),
                 (ridx[ia[inte]], ridx[ib[inte]])),
                shape=(len(int_ids), len(int_ids)))
nch, chain = connected_components(gi, directed=False)
del gi; gc.collect()
body_len = np.zeros(max(nch, 1), np.float64)
np.add.at(body_len, chain[ridx[ia[inte]]], slen[inte])

b_ia, b_ib, b_w = ia[bnd], ib[bnd], slen[bnd]
b_keep  = np.where(keep[b_ia], b_ia, b_ib)
b_int   = np.where(keep[b_ia], b_ib, b_ia)
b_chain = chain[ridx[b_int]]
order = np.argsort(b_chain, kind="stable")
b_keep, b_chain, b_w = b_keep[order], b_chain[order], b_w[order]
starts = np.r_[0, np.where(np.diff(b_chain) != 0)[0] + 1, len(b_chain)]
cnt = np.diff(starts)
two = cnt == 2
i0 = starts[:-1][two]
e_u = b_keep[i0]; e_v = b_keep[i0 + 1]
e_w = b_w[i0] + b_w[i0 + 1] + body_len[b_chain[i0]]
print(f"chains {nch:,}  chain edges {int(two.sum()):,}  "
      f"degenerate {int((~two).sum()):,}  {time.time()-t0:.0f}s", flush=True)

kv = np.where(keep)[0]
K = len(kv)
kidx = np.full(NV, -1, np.int64)
kidx[kv] = np.arange(K)

eu = np.r_[kidx[ia[kk]], kidx[e_u]]
ev = np.r_[kidx[ib[kk]], kidx[e_v]]
ew = np.r_[slen[kk], e_w]
ok = eu != ev
eu, ev, ew = eu[ok], ev[ok], ew[ok]

# ---------------- virtual site nodes (canonical hang only) ----------------
# One physical site -> one grid hang: the member with smallest attach_dist.
# Virtual node connects ONLY to that canonical vertex (not every member).
sub_ok = (sub_vid >= 0) & np.isfinite(sub_d) & (sub_d <= ATTACH_MAX)
att_sites = np.unique(site_of_sub[sub_ok])
NSRC = len(att_sites)
site_rank = np.full(NSITE_ALL, -1, np.int64)
site_rank[att_sites] = np.arange(NSRC)
m = sub_ok.copy()

rows_ok = np.where(m)[0]
bo = np.lexsort((sub_d[rows_ok], site_of_sub[rows_ok]))
rows_sorted = rows_ok[bo]
bf = np.r_[True, np.diff(site_of_sub[rows_sorted]) != 0]
best_row = np.full(NSITE_ALL, -1, np.int64)
best_row[site_of_sub[rows_sorted[bf]]] = rows_sorted[bf]
site_bus_raw = sub_vid[best_row[att_sites]].astype(np.int64)
site_bus_k = kidx[site_bus_raw].astype(np.int64)
canon_dist = sub_d[best_row[att_sites]].astype(np.float32)
assert (site_bus_k >= 0).all()

# Audit-only: how often multi-member sites would have disagreed on hang point.
member_k = kidx[sub_vid[m]].astype(np.int64)
member_rank = site_rank[site_of_sub[m]].astype(np.int32)
member_dist = sub_d[m].astype(np.float32)
n_member_links = int(len(member_k))
n_sites_multi_hang = 0
if NSRC > 0:
    # per-site unique hang count among members
    order_m = np.argsort(member_rank, kind="stable")
    rk_s = member_rank[order_m]
    kk_s = member_k[order_m]
    starts = np.r_[0, np.where(np.diff(rk_s) != 0)[0] + 1, len(rk_s)]
    for i0, i1 in zip(starts[:-1], starts[1:]):
        if len(np.unique(kk_s[i0:i1])) > 1:
            n_sites_multi_hang += 1

virt_edges_u = K + np.arange(NSRC, dtype=np.int64)
virt_edges_v = site_bus_k.astype(np.int64)
virt_edges_w = canon_dist.astype(np.float64)
assert (virt_edges_v >= 0).all()
eu = np.r_[eu, virt_edges_u]
ev = np.r_[ev, virt_edges_v]
ew = np.r_[ew, virt_edges_w]
N = K + NSRC
virt = K + np.arange(NSRC)
print(f"graph: {N:,} nodes  {len(eu):,} edges  source sites {NSRC:,}  "
      f"(canonical hang only)  {time.time()-t0:.0f}s", flush=True)

# Failure cuts use the same single canonical hang per site.
attach_rank = np.arange(NSRC, dtype=np.int32)
attach_k = site_bus_k.astype(np.int32)
attach_raw = site_bus_raw.astype(np.int64)
attach_dist = canon_dist.astype(np.float32)

# Duplicate undirected pairs would be SUMMED by tocsr(); canonicalise each
# pair and deduplicate with MIN before constructing a symmetric CSR graph.
e0 = np.minimum(eu, ev).astype(np.int64)
e1 = np.maximum(eu, ev).astype(np.int64)
key = e0 * N + e1
order = np.argsort(key, kind="stable")
key, ew2 = key[order], ew[order]
first = np.r_[True, np.diff(key) != 0]
grp = np.cumsum(first) - 1
wmin = np.full(grp[-1] + 1, np.inf)
np.minimum.at(wmin, grp, ew2)
ku = (key[first] // N).astype(np.int64); kvv = (key[first] % N).astype(np.int64)
A = coo_matrix(
    (np.r_[wmin, wmin], (np.r_[ku, kvv], np.r_[kvv, ku])),
    shape=(N, N)).tocsr()
del key, ew2, grp, wmin, ku, kvv, eu, ev, ew, e0, e1
gc.collect()

np.savez(
    f"{GRAPH}/contracted.npz",
    kv=kv, K=K, NSRC=NSRC, att_sites=att_sites,
    site_bus_raw=site_bus_raw, site_bus_k=site_bus_k,
    attach_k=attach_k, attach_rank=attach_rank, attach_raw=attach_raw,
    indptr=A.indptr.astype(np.int64),
    indices=A.indices.astype(np.int32),
    data=A.data.astype(np.float32))

topology_audit = dict(
    attachment_mode="canonical_nearest_vertex_hang",
    virtual_site_policy=("one virtual node per site, linked only to the "
                         "member with smallest geometric attach_dist"),
    n_sites_all=int(NSITE_ALL),
    n_sites_attached=int(NSRC),
    n_member_attachments_raw=n_member_links,
    n_sites_with_multiple_member_hangs=int(n_sites_multi_hang),
    n_failure_cut_attachments=int(len(attach_k)),
    n_unique_attachment_vertices=int(len(np.unique(attach_k))),
    n_canonical_attachment_vertices=int(len(np.unique(site_bus_k))),
    attachment_median_m=float(np.median(attach_dist)),
    attachment_p95_m=float(np.percentile(attach_dist, 95)),
    attachment_max_m=float(np.max(attach_dist)),
    canonical_attachment_max_m=float(np.max(canon_dist)),
    n_original_vertices=int(NV),
    n_kept_vertices=int(K),
    n_original_edges=int(len(zs["ia"])),
    original_line_length_m=original_line_length_m,
    failure_isolation=("failed virtual site isolates its single canonical "
                       "attachment vertex and all incident contracted edges"))
json.dump(topology_audit, open(P.fresh(f"{OUT}/site_topology_audit.json"), "w"),
          indent=1)
print(f"canonical hangs: {len(attach_k):,} sites on "
      f"{len(np.unique(attach_k)):,} vertices  "
      f"(raw member links {n_member_links:,}; "
      f"multi-hang sites {n_sites_multi_hang:,})  "
      f"median {np.median(attach_dist):.0f} m  {time.time()-t0:.0f}s",
      flush=True)

def multi_dijkstra(source_nodes):
    d, _, srcs = dijkstra(A, directed=False, indices=source_nodes,
                          return_predecessors=True, min_only=True,
                          limit=DJ_LIMIT)
    return d, srcs

# ---------------- pass 0 ----------------
D1, S1n = multi_dijkstra(virt)
S1 = np.where(np.isfinite(D1), S1n - K, -1).astype(np.int64)   # rank 0..NSRC-1
print(f"pass0: reached {np.isfinite(D1[:K]).mean()*100:.1f}% of nodes  "
      f"{time.time()-t0:.0f}s", flush=True)

# ---------------- bit passes ----------------
D2 = np.full(N, np.inf)
S2 = np.full(N, -1, np.int64)
B2 = np.full(N, -1, np.int8)
nbits = int(np.ceil(np.log2(max(NSRC, 2))))
rk = np.arange(NSRC)
for b in range(nbits):
    m0 = ((rk >> b) & 1) == 0
    dA, sA = multi_dijkstra(virt[m0])
    dB, sB = multi_dijkstra(virt[~m0])
    s1_in_A = np.zeros(N, bool)
    valid = S1 >= 0
    s1_in_A[valid] = ((S1[valid] >> b) & 1) == 0
    candD = np.where(s1_in_A, dB, dA)
    candS = np.where(s1_in_A, sB, sA)
    upd = candD < D2
    D2[upd] = candD[upd]
    S2[upd] = candS[upd] - K
    B2[upd] = b
    print(f"bit {b+1}/{nbits}  {time.time()-t0:.0f}s", flush=True)
S2 = np.where(np.isfinite(D2), S2, -1)

np.savez(f"{GRAPH}/node_feeders.npz",
         D1=D1[:K].astype(np.float32), S1=S1[:K].astype(np.int32),
         D2=D2[:K].astype(np.float32), S2=S2[:K].astype(np.int32),
         B2=B2[:K])

# site table (rank space used by s1/s2 labels)
pd.DataFrame(dict(site_rank=np.arange(NSRC),
                  site_id=att_sites,
                  rep_sub_id=rep[att_sites],
                  n_members=site_n[att_sites],
                  major=site_major[att_sites],
                  bus_raw=site_bus_raw,
                  bus_k=site_bus_k,
                  attach_dist_m=canon_dist)
             ).to_parquet(f"{OUT}/site_table.parquet")

# ---------------- read out chargers & plants ----------------
def rank_arrays(vid, conn_d):
    okv = vid >= 0
    node = np.where(okv, kidx[np.maximum(vid, 0)], 0)
    d1 = np.where(okv, conn_d + D1[node], np.inf)
    s1 = np.where(okv & np.isfinite(d1), S1[node], -1)
    d2 = np.where(okv, conn_d + D2[node], np.inf)
    s2 = np.where(okv & np.isfinite(d2), S2[node], -1)
    return d1, s1, d2, s2

d1, s1, d2, s2 = rank_arrays(ch_vid, ch_d.astype(np.float64))
pd.DataFrame(dict(
    charger_idx=np.arange(len(ch_vid), dtype=np.int64),
    attach_vid=ch_vid.astype(np.int64),
    attach_dist_m=ch_d.astype(np.float32),
    s1_rank=s1.astype(np.int64),
    d1_net_m=np.where(np.isfinite(d1), d1, np.nan).astype(np.float32),
    s2_rank=s2.astype(np.int64),
    d2_net_m=np.where(np.isfinite(d2), d2, np.nan).astype(np.float32))
    ).to_parquet(f"{OUT}/charger_network_assignment.parquet")
fin = np.isfinite(d1)
print(f"chargers: line-attached {(ch_vid>=0).mean()*100:.2f}% | reach a site "
      f"{fin.mean()*100:.2f}% | median d1 "
      f"{np.nanmedian(np.where(fin, d1, np.nan))/1000:.2f} km | "
      f"2nd distinct site {np.isfinite(d2).mean()*100:.1f}%", flush=True)

d1, s1, d2, s2 = rank_arrays(pl_vid, pl_d.astype(np.float64))
pd.DataFrame(dict(
    plant_id=np.arange(len(pl_vid), dtype=np.int64),
    attach_vid=pl_vid.astype(np.int64),
    attach_dist_m=pl_d.astype(np.float32),
    s1_rank=s1.astype(np.int64),
    d1_net_m=np.where(np.isfinite(d1), d1, np.nan).astype(np.float32))
    ).to_parquet(f"{OUT}/plant_network_assignment.parquet")
print(f"DONE {time.time()-t0:.0f}s", flush=True)
