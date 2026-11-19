#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resolve analysis inputs and outputs relative to the repository. EVFR_SOURCE_ROOT selects the input tree. Each named path also accepts an EVFR_<KEY> override."""

import os

# ---------------------------------------------------------------------------
# Package root: this file lives at <package>/shared/paths.py
# ---------------------------------------------------------------------------
PACKAGE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _root(*parts):
    return os.path.join(PACKAGE_ROOT, *parts)

# ---------------------------------------------------------------------------
# Input data root. Default: <package>/source_data. When the package sits in an
# archive that keeps the inputs in a sibling folder (<archive>/02_input_data
# next to <archive>/03_analysis), that folder is used instead. The environment
# variable EVFR_SOURCE_ROOT overrides both.
# ---------------------------------------------------------------------------
def _first_dir(candidates, default):
    for c in candidates:
        if os.path.isdir(c):
            return c
    return default

SOURCE_ROOT = os.environ.get("EVFR_SOURCE_ROOT") or _first_dir(
    (_root("source_data"),
     os.path.join(os.path.dirname(PACKAGE_ROOT), "02_input_data")),
    _root("source_data"))

# Folder that holds the rendered figures (one sub-folder per figure), if any.
FIGURES_ROOT = os.environ.get("EVFR_FIGURES_ROOT") or _first_dir(
    (_root("figures"),
     os.path.join(os.path.dirname(PACKAGE_ROOT), "05_figures"),
     os.path.join(os.path.dirname(PACKAGE_ROOT), "04_figures_and_source_data")),
    None)

def _src(*parts):
    return os.path.join(SOURCE_ROOT, *parts)

# ---------------------------------------------------------------------------
# Canonical locations (inputs under SOURCE_ROOT, everything else under the
# package root, unless overridden)
# ---------------------------------------------------------------------------
_DEFAULTS = {
    # -- source data: station point databases (restricted licence, see notice)
    "CHARGING_RAW_GOOGLE": _src("stations_charging_google",
                                 "google_places_charging_stations.csv"),
    "CHARGING_RAW_AMAP":   _src("stations_charging_amap",
                                 "amap_charging_stations_china_gcj02.csv"),
    "FUEL_RAW_GOOGLE":     _src("stations_fuel_google",
                                 "google_places_fuel_stations.csv"),
    "OCM_ARCHIVE":         _src("stations_registry_ocm",
                                 "ocm_charging_stations_2010_2026.csv"),

    # -- source data: flood hazard (JRC global river flood maps, 90 m)
    "JRC_DIR":             _src("flood_hazard_jrc"),

    # -- source data: power system
    "OSM_POWER":           _src("power_grid_osm",
                                 "osm_power_layer.gpkg"),
    "GRIDFINDER":          _src("power_grid_gridfinder",
                                 "gridfinder_mv_lines.gpkg"),
    "WRI_PLANTS":          _src("power_plants_wri",
                                 "wri_global_power_plant_database.csv"),

    # -- source data: socioeconomic and protection layers
    "WB_INCOME":           _src("socioeconomic_worldbank",
                                 "worldbank_income_groups_fy2026.json"),
    "OWID_POPULATION":     _src("socioeconomic_worldbank",
                                 "owid_population_2023.csv"),
    "KUMMU_ADM0":          _src("socioeconomic_gdp_kummu",
                                 "kummu_adm0_gdp_per_capita.csv"),
    "KUMMU_ADM2_TIF":      _src("socioeconomic_gdp_kummu",
                                 "kummu_adm2_gdp_per_capita_1990_2024.tif"),
    "FLOPROS_SHP_DIR":     _src("flood_protection_flopros",
                                 "FLOPROS_shp_V1"),
    "FLOPROS_TIF":         _src("flood_protection_flopros",
                                 "flopros_merl_riv_2p5arcmin.tif"),

    # -- source data: boundaries, urban centres, statistics
    "NE10M":               _src("boundaries_natural_earth",
                                 "ne_10m_admin_0_countries.geojson"),
    "NE110M":              _src("boundaries_natural_earth",
                                 "ne_110m_countries_simplified.geojson"),
    "UCDB_GPKG":           _src("urban_centres_ghs_ucdb",
                                 "GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg"),
    "IEA_DIR":             _src("statistics_iea"),

    # -- source data: map backgrounds (SRTM tiles and OSM extracts, cached)
    "MAP_BG_CACHE":        _src("map_backgrounds_srtm_osm"),

    # -- shared station database (built by stage 0, consumed by every figure)
    "DB_DIR":              _root("shared", "station_database"),
    "DB_CHARGING":         _root("shared", "station_database",
                                 "charging_stations_with_country.parquet"),
    "DB_FUEL":             _root("shared", "station_database",
                                 "fuel_stations_with_country.parquet"),
    "DB_CHARGING_GOOGLE_ONLY": _root("shared", "station_database",
                                 "charging_stations_google_only.parquet"),
    "DB_MERGED_GPKG":      _root("shared", "station_database",
                                 "charging_stations_merged.gpkg"),

    # -- per-module output folders
    "FIG1_DIR":            _root("analysis", "fig1_global_infrastructure"),
    "FIG1_DATA":           _root("analysis", "fig1_global_infrastructure", "data"),
    "FIG2_DIR":            _root("analysis", "fig2_flood_exposure"),
    "FIG2_DATA":           _root("analysis", "fig2_flood_exposure", "data"),
    "FIG3_DIR":            _root("analysis", "fig3_opening_year_cohorts"),
    "FIG3_DATA":           _root("analysis", "fig3_opening_year_cohorts", "data"),
    "FIG4_DIR":            _root("analysis", "fig4_grid_cascade"),
    "FIG4_DATA":           _root("analysis", "fig4_grid_cascade", "data"),
    "FIG5_DIR":            _root("analysis", "fig5_inequality_protection"),
    "FIG5_DATA":           _root("analysis", "fig5_inequality_protection", "data"),
    "SI_TABLES_DIR":       _root("analysis", "si_tables"),
    "SI_TABLES_DATA":      _root("analysis", "si_tables", "data"),
}

def get(key):
    """Return the resolved location for *key* (environment override first)."""
    env = os.environ.get("EVFR_" + key)
    if env:
        return env
    try:
        return _DEFAULTS[key]
    except KeyError:
        raise KeyError(f"paths.py: unknown location key '{key}'")

def jrc_raster(rp):
    """Water-depth GeoTIFF for return period *rp* ('010' ... '500')."""
    return os.path.join(get("JRC_DIR"), f"RP{rp}.tif")

def ensure_dir(path):
    """Create *path* (a directory) if missing; return it unchanged."""
    os.makedirs(path, exist_ok=True)
    return path

def ensure_parent(path):
    """Create the parent directory of *path* if missing; return *path*."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path

def fresh(path):
    """Return *path* after removing any existing file there.

    Outputs are always written to a fresh file: overwriting in place can
    fail on Windows when a previous copy carries a read-only attribute or
    a restrictive ACL (for example after a file transfer), while deleting
    and recreating always works. A missing file or an undeletable parent
    is reported by the subsequent write, not here."""
    try:
        if os.path.isfile(path):
            os.remove(path)
    except OSError:
        pass
    return path

def describe():
    """Return '(key, resolved path, exists?)' rows for every location."""
    rows = []
    for k in sorted(_DEFAULTS):
        p = get(k)
        rows.append((k, p, os.path.exists(p)))
    return rows

if __name__ == "__main__":
    print(f"Package root: {PACKAGE_ROOT}")
    print(f"Input data  : {SOURCE_ROOT}")
    print(f"Figures     : {FIGURES_ROOT}\n")
    for k, p, ok in describe():
        print(f"  [{'ok' if ok else 'MISSING'}] {k:26s} {p}")
