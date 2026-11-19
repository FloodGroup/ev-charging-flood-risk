#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Route chargers to reachable feeder sites on the residual network. A failed site isolates its attachment vertices and all incident edges. The production graph uses the attachment of the member nearest to the line."""
from __future__ import annotations

import gc
import heapq
import os
import time
from collections import defaultdict

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra


R_SERV = 10_000.0
MAX_PATH_SITES = 16


def load_graph(graph_dir):
    """Return the symmetric contracted CSR graph and mapping arrays."""
    z = np.load(os.path.join(graph_dir, "contracted.npz"))
    n = int(z["K"]) + int(z["NSRC"])
    graph = csr_matrix(
        (z["data"].astype(np.float64),
         z["indices"].astype(np.int32),
         z["indptr"].astype(np.int64)),
        shape=(n, n))
    return graph, {k: z[k] for k in z.files if k not in
                   {"data", "indices", "indptr"}}


def failed_attachment_nodes(meta, failed_sites):
    """Contracted k-indices isolated when the given site ranks fail."""
    failed_sites = np.asarray(failed_sites, dtype=bool)
    if "attach_k" in meta and "attach_rank" in meta:
        attach_k = np.asarray(meta["attach_k"], dtype=np.int64)
        attach_rank = np.asarray(meta["attach_rank"], dtype=np.int64)
        return np.unique(attach_k[failed_sites[attach_rank]])
    site_bus_k = np.asarray(meta["site_bus_k"], dtype=np.int64)
    return np.unique(site_bus_k[np.where(failed_sites)[0]])


def _attachment_site_lists(meta, K, NSRC):
    """Map each contracted vertex to the site ranks attached there."""
    lists = [[] for _ in range(K)]
    if "attach_k" in meta and "attach_rank" in meta:
        attach_k = np.asarray(meta["attach_k"], dtype=np.int64)
        attach_rank = np.asarray(meta["attach_rank"], dtype=np.int64)
        for k, r in zip(attach_k, attach_rank):
            if 0 <= k < K and 0 <= r < NSRC:
                lists[int(k)].append(int(r))
    else:
        site_bus_k = np.asarray(meta["site_bus_k"], dtype=np.int64)
        for r, k in enumerate(site_bus_k):
            if 0 <= k < K:
                lists[int(k)].append(int(r))
    # Deduplicate while preserving order.
    out = []
    for lst in lists:
        seen = set()
        uniq = []
        for r in lst:
            if r not in seen:
                seen.add(r)
                uniq.append(r)
        out.append(uniq)
    return out


def _multi_dijkstra(graph, sources, limit=60_000.0):
    return dijkstra(
        graph, directed=False, indices=sources,
        return_predecessors=True, min_only=True, limit=limit)


def _trace_path_sites(u, pred, K, sites_at, max_path_sites, path_row,
                      overflow_row, ensure_rank=None):
    nsite = 0
    steps = 0
    seen_ranks = set()
    while 0 <= u < K:
        for sr in sites_at[u]:
            if sr in seen_ranks:
                continue
            seen_ranks.add(sr)
            if nsite < max_path_sites:
                path_row[nsite] = sr
            else:
                overflow_row[0] = True
            nsite += 1
        p = int(pred[u])
        if p < 0:
            break
        u = p
        steps += 1
        if steps > 100_000:
            raise RuntimeError("predecessor cycle while tracing path")
    # Rare predecessor ties can miss the labelled destination attachment;
    # the assigned feeder is always on the route by construction.
    if ensure_rank is not None and int(ensure_rank) >= 0:
        sr = int(ensure_rank)
        if sr not in seen_ranks:
            if nsite < max_path_sites:
                path_row[nsite] = sr
            else:
                overflow_row[0] = True
            nsite += 1
    return nsite


def build_backup_certificates(graph_dir, charger_attach_vid,
                              out_path, max_path_sites=MAX_PATH_SITES):
    """Build compact site-node certificates for every fixed backup route.

    The bit that won the exact second-nearest search is stored by step 06.
    We rerun one bit partition at a time, retain its predecessor forest only
    while tracing the chargers won by that bit, then release it.  This keeps
    memory bounded while preserving the exact intact-network route.
    """
    t0 = time.time()
    graph, meta = load_graph(graph_dir)
    K, NSRC = int(meta["K"]), int(meta["NSRC"])
    kv = meta["kv"].astype(np.int64)
    virt = K + np.arange(NSRC, dtype=np.int64)
    sites_at = _attachment_site_lists(meta, K, NSRC)

    nf = np.load(os.path.join(graph_dir, "node_feeders.npz"))
    S1 = nf["S1"].astype(np.int64)
    S2 = nf["S2"].astype(np.int64)
    B2 = nf["B2"].astype(np.int16)
    kidx = np.full(int(kv.max()) + 1, -1, np.int64)
    kidx[kv] = np.arange(K)
    ok_attach = (charger_attach_vid >= 0) & (charger_attach_vid < len(kidx))
    ch_node = np.full(len(charger_attach_vid), -1, np.int64)
    ch_node[ok_attach] = kidx[charger_attach_vid[ok_attach]]

    path_sites = np.full(
        (len(charger_attach_vid), max_path_sites), -1, np.int32)
    path_count = np.zeros(len(charger_attach_vid), np.uint8)
    overflow = np.zeros(len(charger_attach_vid), bool)

    nbits = int(np.ceil(np.log2(max(NSRC, 2))))
    ranks = np.arange(NSRC)
    for b in range(nbits):
        m0 = ((ranks >> b) & 1) == 0
        ci = np.where(
            (ch_node >= 0) & (B2[ch_node.clip(min=0)] == b) &
            (S1[ch_node.clip(min=0)] >= 0) &
            (S2[ch_node.clip(min=0)] >= 0))[0]
        if not len(ci):
            continue
        n0 = ch_node[ci]
        s1_in_0 = ((S1[n0] >> b) & 1) == 0
        for use_zero_sources, choose in ((True, ~s1_in_0),
                                         (False, s1_in_0)):
            cc = ci[choose]
            if not len(cc):
                continue
            sources = virt[m0] if use_zero_sources else virt[~m0]
            _, pred, _ = _multi_dijkstra(graph, sources)
            for charger_i in cc:
                ov = [False]
                n0 = int(ch_node[charger_i])
                nsite = _trace_path_sites(
                    n0, pred, K, sites_at,
                    max_path_sites, path_sites[charger_i], ov,
                    ensure_rank=int(S2[n0]))
                overflow[charger_i] = ov[0]
                path_count[charger_i] = min(nsite, 255)
            del pred
            gc.collect()
        print(f"certificate bit {b+1}/{nbits}: {len(ci):,} chargers  "
              f"{time.time()-t0:.0f}s", flush=True)

    np.savez_compressed(
        out_path,
        path_sites=path_sites,
        path_count=path_count,
        overflow=overflow,
        charger_node=ch_node,
        max_path_sites=np.int32(max_path_sites))
    print(f"backup certificates: {(path_count>0).sum():,} paths, "
          f"overflow {overflow.sum():,}, {time.time()-t0:.0f}s", flush=True)
    return path_sites, path_count, overflow, ch_node


def build_primary_certificates(graph_dir, charger_attach_vid,
                               out_path, max_path_sites=MAX_PATH_SITES):
    """Build compact site-node certificates for intact primary routes."""
    t0 = time.time()
    graph, meta = load_graph(graph_dir)
    K, NSRC = int(meta["K"]), int(meta["NSRC"])
    kv = meta["kv"].astype(np.int64)
    virt = K + np.arange(NSRC, dtype=np.int64)
    sites_at = _attachment_site_lists(meta, K, NSRC)
    nf = np.load(os.path.join(graph_dir, "node_feeders.npz"))
    S1 = nf["S1"].astype(np.int64)

    kidx = np.full(int(kv.max()) + 1, -1, np.int64)
    kidx[kv] = np.arange(K)
    ok_attach = (charger_attach_vid >= 0) & (charger_attach_vid < len(kidx))
    ch_node = np.full(len(charger_attach_vid), -1, np.int64)
    ch_node[ok_attach] = kidx[charger_attach_vid[ok_attach]]
    _, pred, _ = _multi_dijkstra(graph, virt)

    path_sites = np.full(
        (len(charger_attach_vid), max_path_sites), -1, np.int32)
    path_count = np.zeros(len(charger_attach_vid), np.uint8)
    overflow = np.zeros(len(charger_attach_vid), bool)
    ci = np.where((ch_node >= 0) & (S1[ch_node.clip(min=0)] >= 0))[0]
    for charger_i in ci:
        ov = [False]
        n0 = int(ch_node[charger_i])
        nsite = _trace_path_sites(
            n0, pred, K, sites_at,
            max_path_sites, path_sites[charger_i], ov,
            ensure_rank=int(S1[n0]))
        overflow[charger_i] = ov[0]
        path_count[charger_i] = min(nsite, 255)

    np.savez_compressed(
        out_path, path_sites=path_sites, path_count=path_count,
        overflow=overflow, charger_node=ch_node,
        max_path_sites=np.int32(max_path_sites))
    print(f"primary certificates: {(path_count>0).sum():,} paths, "
          f"transit {(path_count>1).sum():,}, overflow {overflow.sum():,}, "
          f"{time.time()-t0:.0f}s", flush=True)
    return path_sites, path_count, overflow, ch_node


class BoundedRouter:
    """Exact fixed-target rerouting with failed attachment-vertex cuts."""

    def __init__(self, graph, K, site_bus_k=None, attach_k=None,
                 attach_rank=None):
        self.indptr = graph.indptr
        self.indices = graph.indices
        self.weights = graph.data
        self.K = int(K)
        self.n = graph.shape[0]
        if attach_k is not None and attach_rank is not None:
            self.attach_k = np.asarray(attach_k, dtype=np.int64)
            self.attach_rank = np.asarray(attach_rank, dtype=np.int64)
            self.site_bus_k = None
        else:
            self.attach_k = None
            self.attach_rank = None
            self.site_bus_k = np.asarray(site_bus_k, dtype=np.int64)

    def _bus_failed(self, failed_sites):
        bus_failed = np.zeros(self.K, bool)
        failed_sites = np.asarray(failed_sites, dtype=bool)
        if self.attach_k is not None:
            bus_failed[self.attach_k[failed_sites[self.attach_rank]]] = True
        else:
            bus_failed[self.site_bus_k[np.where(failed_sites)[0]]] = True
        return bus_failed

    def route_groups(self, charger_indices, charger_nodes, target_ranks,
                     charger_attach_m, failed_sites, limit=R_SERV):
        """Return availability for affected charger-to-fixed-site queries.

        Queries are grouped by fixed target.  One bounded Dijkstra settles
        all charger nodes requested for that target, and stops once all are
        found or the 10 km frontier is exhausted.  Failed attachment
        vertices cannot be entered or left, so every incident edge is cut.
        """
        out = np.zeros(len(charger_indices), bool)
        if not len(charger_indices):
            return out, dict(groups=0, visited=0)
        bus_failed = self._bus_failed(failed_sites)

        order = np.argsort(target_ranks, kind="stable")
        rr = target_ranks[order]
        starts = np.r_[0, np.where(np.diff(rr) != 0)[0] + 1, len(rr)]
        total_visited = 0
        groups = 0
        for a, b in zip(starts[:-1], starts[1:]):
            pos = order[a:b]
            rank = int(rr[a])
            if failed_sites[rank]:
                continue
            source = self.K + rank
            node_to_pos = defaultdict(list)
            max_graph_limit = 0.0
            for q in pos:
                node = int(charger_nodes[charger_indices[q]])
                if node < 0 or bus_failed[node]:
                    continue
                node_to_pos[node].append(int(q))
                max_graph_limit = max(
                    max_graph_limit,
                    float(limit - charger_attach_m[charger_indices[q]]))
            if not node_to_pos or max_graph_limit < 0:
                continue
            groups += 1
            dist = {source: 0.0}
            heap = [(0.0, source)]
            remaining = set(node_to_pos)
            while heap and remaining:
                du, u = heapq.heappop(heap)
                if du != dist.get(u) or du > max_graph_limit:
                    continue
                total_visited += 1
                if u in remaining:
                    for q in node_to_pos[u]:
                        ci = charger_indices[q]
                        if du + charger_attach_m[ci] <= limit:
                            out[q] = True
                    remaining.remove(u)
                if u < self.K and bus_failed[u]:
                    continue
                for p in range(self.indptr[u], self.indptr[u + 1]):
                    v = int(self.indices[p])
                    if v >= self.K and v != source:
                        continue
                    if v < self.K and bus_failed[v]:
                        continue
                    nd = du + float(self.weights[p])
                    if nd <= max_graph_limit and nd < dist.get(v, np.inf):
                        dist[v] = nd
                        heapq.heappush(heap, (nd, v))
        return out, dict(groups=groups, visited=total_visited)

    def nearest_live_site(self, source_node, bus_failed, failed_sites,
                          limit, site_kind=None, required_kind=None):
        """Nearest residual live site: (distance, site_rank) or (inf, -1).

        If required_kind is set, only virtual sites whose site_kind matches
        are accepted as service targets (voltage / OSM-kind compatibility).
        """
        source_node = int(source_node)
        if source_node < 0 or source_node >= self.K or bus_failed[source_node]:
            return np.inf, -1
        failed_sites = np.asarray(failed_sites, dtype=bool)
        dist = {source_node: 0.0}
        heap = [(0.0, source_node)]
        while heap:
            du, u = heapq.heappop(heap)
            if du != dist.get(u) or du > limit:
                continue
            if u >= self.K:
                r = u - self.K
                if failed_sites[r]:
                    continue
                if (required_kind is not None and site_kind is not None and
                        int(site_kind[r]) != int(required_kind)):
                    continue
                return du, int(r)
            if bus_failed[u]:
                continue
            for p in range(self.indptr[u], self.indptr[u + 1]):
                v = int(self.indices[p])
                if v < self.K and bus_failed[v]:
                    continue
                if v >= self.K and failed_sites[v - self.K]:
                    continue
                nd = du + float(self.weights[p])
                if nd <= limit and nd < dist.get(v, np.inf):
                    dist[v] = nd
                    heapq.heappush(heap, (nd, v))
        return np.inf, -1

    def nearest_live_site_dist(self, source_node, bus_failed, failed_sites,
                               limit, site_kind=None, required_kind=None):
        """Shortest residual distance from a physical node to any live site."""
        d, _ = self.nearest_live_site(
            source_node, bus_failed, failed_sites, limit,
            site_kind=site_kind, required_kind=required_kind)
        return d

    def route_any_alternate(self, charger_indices, charger_nodes,
                            charger_attach_m, failed_sites, limit=R_SERV,
                            site_kind=None, required_kind=None):
        """Reachability and site_rank of a compatible live residual backup.

        Used after the assigned primary site is isolated.  Optional
        required_kind enforces OSM-kind / voltage-class matching with the
        primary feeder site.  Returns (ok, ranks, stats); ranks[i] is the
        site_rank of the nearest accepted alternate, or -1 if none.
        """
        out = np.zeros(len(charger_indices), bool)
        ranks = np.full(len(charger_indices), -1, np.int64)
        if not len(charger_indices):
            return out, ranks, dict(groups=0, visited=0, unique_starts=0)
        bus_failed = self._bus_failed(failed_sites)
        node_budget = {}
        node_members = defaultdict(list)
        for q, ci in enumerate(charger_indices):
            node = int(charger_nodes[ci])
            am = float(charger_attach_m[ci])
            if node < 0 or node >= self.K or bus_failed[node]:
                continue
            budget = float(limit - am)
            if budget < 0:
                continue
            node_members[node].append((q, am))
            prev = node_budget.get(node)
            if prev is None or budget > prev:
                node_budget[node] = budget
        visited = 0
        for node, budget in node_budget.items():
            d_site, rank = self.nearest_live_site(
                node, bus_failed, failed_sites, budget,
                site_kind=site_kind, required_kind=required_kind)
            visited += 1
            if not np.isfinite(d_site) or rank < 0:
                continue
            for q, am in node_members[node]:
                if d_site + am <= limit:
                    out[q] = True
                    ranks[q] = rank
        return out, ranks, dict(groups=1, visited=visited,
                                unique_starts=len(node_budget))
