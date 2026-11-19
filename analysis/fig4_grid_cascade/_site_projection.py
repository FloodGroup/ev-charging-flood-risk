#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact bounded spatial search for site-to-line projection.

Every line segment is indexed in every 0.1-degree cell touched by its
longitude/latitude bounding box.  A site searches expanding cell rings and
stops only when the geodesic lower bound to the unsearched exterior exceeds
its current best distance.  Bounding-box replication is conservative:
segments may be visited unnecessarily, but a potentially closer segment
cannot be omitted.
"""
from __future__ import annotations

import gc
import time

import numpy as np


RE = 6_371_000.0


def build_segment_grid(vlon, vlat, ia, ib, cell_deg=0.1):
    """Build a conservative segment-bounding-box grid index."""
    t0 = time.time()
    nx_world = int(round(360.0 / cell_deg))
    ny_world = int(round(180.0 / cell_deg))
    a = vlon[ia].astype(np.float64)
    b = vlon[ib].astype(np.float64)
    d = b - a
    b = np.where(d > 180.0, b - 360.0,
                 np.where(d < -180.0, b + 360.0, b))
    ix0 = np.floor((np.minimum(a, b) + 180.0) / cell_deg).astype(np.int32)
    ix1 = np.floor((np.maximum(a, b) + 180.0) / cell_deg).astype(np.int32)
    iy0 = np.floor(
        (np.minimum(vlat[ia], vlat[ib]) + 90.0) / cell_deg
    ).astype(np.int32)
    iy1 = np.floor(
        (np.maximum(vlat[ia], vlat[ib]) + 90.0) / cell_deg
    ).astype(np.int32)
    iy0 = np.clip(iy0, 0, ny_world - 1)
    iy1 = np.clip(iy1, 0, ny_world - 1)
    sx = (ix1 - ix0 + 1).astype(np.int32)
    sy = (iy1 - iy0 + 1).astype(np.int32)
    counts = (sx.astype(np.int64) * sy).astype(np.int32)
    total = int(counts.sum(dtype=np.int64))
    if total >= np.iinfo(np.int32).max:
        raise MemoryError("segment-grid assignments exceed int32 range")

    edge = np.repeat(np.arange(len(ia), dtype=np.int32), counts)
    prefix = (np.cumsum(counts, dtype=np.int64) - counts).astype(np.int32)
    local = np.arange(total, dtype=np.int32) - np.repeat(prefix, counts)
    ex = sx[edge]
    cx = (ix0[edge] + local % ex) % nx_world
    cy = iy0[edge] + local // ex
    cell = (cy * nx_world + cx).astype(np.int32)
    order = np.argsort(cell, kind="stable")
    cell, edge = cell[order], edge[order]
    unique, starts, ncell = np.unique(
        cell, return_index=True, return_counts=True)
    start_lookup = np.full(nx_world * ny_world, -1, np.int32)
    count_lookup = np.zeros(nx_world * ny_world, np.int32)
    start_lookup[unique] = starts.astype(np.int32)
    count_lookup[unique] = ncell.astype(np.int32)
    audit = dict(
        cell_deg=float(cell_deg),
        nx_world=int(nx_world),
        ny_world=int(ny_world),
        assignments=int(total),
        assignments_per_edge=float(total / len(ia)),
        max_cells_per_edge=int(counts.max()),
        build_seconds=float(time.time() - t0))
    del a, b, d, ix0, ix1, iy0, iy1, sx, sy, counts
    del prefix, local, ex, cx, cy, cell, order, unique, starts, ncell
    gc.collect()
    return edge, start_lookup, count_lookup, audit


def _ring_cells(ix, iy, ring, nx_world, ny_world):
    if ring == 0:
        return [iy * nx_world + ix]
    cells = []
    y0, y1 = iy - ring, iy + ring
    for x in range(ix - ring, ix + ring + 1):
        xx = x % nx_world
        if 0 <= y0 < ny_world:
            cells.append(y0 * nx_world + xx)
        if y1 != y0 and 0 <= y1 < ny_world:
            cells.append(y1 * nx_world + xx)
    for y in range(iy - ring + 1, iy + ring):
        if not 0 <= y < ny_world:
            continue
        cells.append(y * nx_world + (ix - ring) % nx_world)
        cells.append(y * nx_world + (ix + ring) % nx_world)
    return cells


def _outside_lower_bound(lon, lat, ix, iy, ring, cell_deg):
    """Tangent-plane lower bound from point to searched-cell exterior."""
    left = -180.0 + (ix - ring) * cell_deg
    right = -180.0 + (ix + ring + 1) * cell_deg
    bottom = -90.0 + (iy - ring) * cell_deg
    top = -90.0 + (iy + ring + 1) * cell_deg
    coslat = max(float(np.cos(np.deg2rad(lat))), 1e-6)
    vals = [
        max(lon - left, 0.0) * coslat,
        max(right - lon, 0.0) * coslat,
    ]
    if bottom > -90.0:
        vals.append(max(lat - bottom, 0.0))
    if top < 90.0:
        vals.append(max(top - lat, 0.0))
    return min(vals) * RE * np.pi / 180.0


def _haversine(lon0, lat0, lon1, lat1):
    lo0, la0 = np.deg2rad(lon0), np.deg2rad(lat0)
    lo1, la1 = np.deg2rad(lon1), np.deg2rad(lat1)
    dlo = (lo1 - lo0 + np.pi) % (2 * np.pi) - np.pi
    dla = la1 - la0
    h = np.sin(dla / 2) ** 2 + np.cos(la0) * np.cos(la1) * \
        np.sin(dlo / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


def refine_site_projections(vlon, vlat, ia, ib, site_lon, site_lat,
                            best_edge, best_t, best_dist,
                            best_plon, best_plat, cell_deg=0.1):
    """Refine seed projections to the globally nearest indexed segment."""
    edge_index, starts, counts, audit = build_segment_grid(
        vlon, vlat, ia, ib, cell_deg)
    nx_world = int(round(360.0 / cell_deg))
    ny_world = int(round(180.0 / cell_deg))
    seed_edge = best_edge.copy()
    seed_dist = best_dist.copy()
    max_ring = 0
    tested = 0
    rad_m = RE * np.pi / 180.0
    t0 = time.time()

    for i in range(len(site_lon)):
        ix = int(np.floor((site_lon[i] + 180.0) / cell_deg)) % nx_world
        iy = int(np.floor((site_lat[i] + 90.0) / cell_deg))
        iy = min(max(iy, 0), ny_world - 1)
        ring = 0
        while True:
            chunks = []
            for cell in _ring_cells(ix, iy, ring, nx_world, ny_world):
                st = int(starts[cell])
                if st >= 0:
                    chunks.append(edge_index[st:st + int(counts[cell])])
            if chunks:
                ee = chunks[0] if len(chunks) == 1 else np.concatenate(chunks)
                tested += len(ee)
                cosla = max(float(np.cos(np.deg2rad(site_lat[i]))), 1e-6)
                dlo_a = (vlon[ia[ee]] - site_lon[i] + 180.0) % 360.0 - 180.0
                dlo_b = (vlon[ib[ee]] - site_lon[i] + 180.0) % 360.0 - 180.0
                ax = dlo_a * cosla * rad_m
                ay = (vlat[ia[ee]] - site_lat[i]) * rad_m
                bx = dlo_b * cosla * rad_m
                by = (vlat[ib[ee]] - site_lat[i]) * rad_m
                dx, dy = bx - ax, by - ay
                den = dx * dx + dy * dy
                tt = np.where(
                    den > 0, -(ax * dx + ay * dy) / den, 0.0)
                tt = np.clip(tt, 0.0, 1.0)
                px, py = ax + tt * dx, ay + tt * dy
                dd = np.hypot(px, py)
                j = int(np.argmin(dd))
                if float(dd[j]) < best_dist[i]:
                    best_dist[i] = float(dd[j])
                    best_edge[i] = int(ee[j])
                    best_t[i] = float(tt[j])
                    best_plon[i] = (
                        (site_lon[i] + px[j] / (cosla * rad_m) + 180.0)
                        % 360.0) - 180.0
                    best_plat[i] = site_lat[i] + py[j] / rad_m
            bound = _outside_lower_bound(
                site_lon[i], site_lat[i], ix, iy, ring, cell_deg)
            if bound > best_dist[i]:
                break
            ring += 1
            max_ring = max(max_ring, ring)
            if ring > 64:
                raise RuntimeError("segment-grid nearest search exceeded 64 rings")
        if (i + 1) % 100_000 == 0:
            print(f"  exact projection {i+1:,}/{len(site_lon):,} "
                  f"{time.time()-t0:.0f}s", flush=True)

    # Connector output follows the model's great-circle distance definition.
    best_dist[:] = _haversine(
        site_lon, site_lat, best_plon, best_plat)
    changed = best_edge != seed_edge
    audit.update(
        search_seconds=float(time.time() - t0),
        candidate_projections_tested=int(tested),
        candidate_tests_per_site=float(tested / len(site_lon)),
        max_search_ring=int(max_ring),
        sites_changed_from_nearest_vertex_incident_seed=int(changed.sum()),
        sites_changed_pct=float(100 * changed.mean()),
        median_improvement_m=float(
            np.median(np.maximum(seed_dist - best_dist, 0.0))),
        p95_improvement_m=float(
            np.percentile(np.maximum(seed_dist - best_dist, 0.0), 95)))
    return best_edge, best_t, best_dist, best_plon, best_plat, audit
