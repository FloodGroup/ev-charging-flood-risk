#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Synthetic regression tests for attachment-vertex cuts and fixed rerouting."""
import numpy as np
from scipy.sparse import csr_matrix

from _joint_failure_routing import BoundedRouter


def graph(edges, k=3, nsrc=2):
    u = np.array([e[0] for e in edges], np.int64)
    v = np.array([e[1] for e in edges], np.int64)
    w = np.array([e[2] for e in edges], np.float64)
    return csr_matrix(
        (np.r_[w, w], (np.r_[u, v], np.r_[v, u])),
        shape=(k + nsrc, k + nsrc))


def query(a, failed, target=1, attach_k=None, attach_rank=None):
    if attach_k is None:
        attach_k = np.array([1, 2])
        attach_rank = np.array([0, 1])
    router = BoundedRouter(
        a, K=3, attach_k=attach_k, attach_rank=attach_rank)
    ok, stats = router.route_groups(
        charger_indices=np.array([0]),
        charger_nodes=np.array([0]),
        target_ranks=np.array([target]),
        charger_attach_m=np.array([0.0]),
        failed_sites=np.asarray(failed, bool))
    return bool(ok[0]), stats


# Physical 0--1--2; source virtual nodes 3->1 (primary), 4->2 (backup).
# Failing the primary isolates attachment vertex 1 and cuts its incident edges.
radial = graph([(0, 1, 1000), (1, 2, 1000), (3, 1, 1), (4, 2, 1)])
ok, _ = query(radial, [True, False])
assert not ok, "backup path must not traverse the failed primary attachment"

# A loop around failed attachment 1 gives an exact route to the same backup.
meshed = graph([
    (0, 1, 1000), (1, 2, 1000), (0, 2, 2500),
    (3, 1, 1), (4, 2, 1)])
ok, _ = query(meshed, [True, False])
assert ok, "meshed network must reroute around the failed attachment"

# A third site is not silently substituted when the fixed backup fails.
ok, _ = query(meshed, [False, True])
assert not ok, "failed fixed backup must remain unavailable"

# Service-radius enforcement applies to the post-failure route.
long_loop = graph([
    (0, 1, 1000), (1, 2, 1000), (0, 2, 10_500),
    (3, 1, 1), (4, 2, 1)])
ok, _ = query(long_loop, [True, False])
assert not ok, "post-failure route beyond 10 km must be unavailable"

# Multi-member attachment (router API): failing site 0 isolates BOTH hangs.
# The grid model stores one canonical hang per site only.
multi = graph([
    (0, 1, 1000), (1, 2, 1000),
    (3, 0, 1), (3, 1, 1), (4, 2, 1)], k=3, nsrc=2)
ok, _ = query(
    multi, [True, False],
    attach_k=np.array([0, 1, 2]),
    attach_rank=np.array([0, 0, 1]))
assert not ok, "all member attachment vertices of a failed site must be cut"

# Synthetic regression: primary isolation then any-site residual search.
# Radial 0-1-2 with primary on 1: isolating attachment 1 also blocks the
# only corridor to the other site, so no alternate remains.
radial = graph([(0, 1, 1000), (1, 2, 1000), (3, 1, 1), (4, 2, 1)])
router = BoundedRouter(
    radial, K=3,
    attach_k=np.array([1, 2]), attach_rank=np.array([0, 1]))
ok, *_ = router.route_any_alternate(
    charger_indices=np.array([0]),
    charger_nodes=np.array([0]),
    charger_attach_m=np.array([0.0]),
    failed_sites=np.array([True, False]))
assert not ok[0], "radial cut at the primary attachment must leave no alternate"

# Meshed bypass around the failed attachment still reaches another site.
meshed2 = graph([
    (0, 1, 1000), (1, 2, 1000), (0, 2, 2500),
    (3, 1, 1), (4, 2, 1)])
router2 = BoundedRouter(
    meshed2, K=3,
    attach_k=np.array([1, 2]), attach_rank=np.array([0, 1]))
ok, *_ = router2.route_any_alternate(
    charger_indices=np.array([0]),
    charger_nodes=np.array([0]),
    charger_attach_m=np.array([0.0]),
    failed_sites=np.array([True, False]))
assert ok[0], "meshed residual graph must reach another live site"

# Kind-matched alternate: only a same-kind live site counts.
meshed2 = graph([
    (0, 1, 1000), (1, 2, 1000), (0, 2, 2500),
    (3, 1, 1), (4, 2, 1)])
router2 = BoundedRouter(
    meshed2, K=3,
    attach_k=np.array([1, 2]), attach_rank=np.array([0, 1]))
site_kind = np.array([0, 1], np.int16)  # different kinds
ok, *_ = router2.route_any_alternate(
    charger_indices=np.array([0]),
    charger_nodes=np.array([0]),
    charger_attach_m=np.array([0.0]),
    failed_sites=np.array([True, False]),
    site_kind=site_kind, required_kind=0)
assert not ok[0], "different-kind residual site must be rejected"

site_kind_same = np.array([0, 0], np.int16)
ok, *_ = router2.route_any_alternate(
    charger_indices=np.array([0]),
    charger_nodes=np.array([0]),
    charger_attach_m=np.array([0.0]),
    failed_sites=np.array([True, False]),
    site_kind=site_kind_same, required_kind=0)
assert ok[0], "same-kind residual site must be accepted"

print("ALL SITE-FAILURE TOPOLOGY TESTS PASSED")
