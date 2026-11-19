"""Export the public figure tables from local analysis outputs.

Use --local-maps to copy the six city extracts needed for Figures 3 and 4
into their respective data/private directories.
"""

import argparse
import csv
import json
from pathlib import Path
import shutil

from matplotlib.figure import Figure
import numpy as np
import pandas as pd

from shared import paths, plot_style


# Only these files are copied intact. Other public products are built below.
PUBLIC_TABLES = (
    ("Figure1", "FIG1_DATA", ("global_summary.json", "country_stats.csv")),
    ("Figure1", "IEA_DIR", ("iea_world_by_speed.csv",)),
    ("Figure2", "FIG2_DATA", (
        "fig02_stats.json", "fig02_stats_v2.json", "rp_curve.csv", "depth_curve.csv",
        "metro_ratios.csv", "city_map_exposure.csv", "city_map_exposure_gas.csv",
        "country_scatter.csv", "exposed_stock.csv",
    )),
    ("Figure3", "FIG3_DATA", (
        "cohort_stats.json", "cohort_curves.csv", "country_slopes.csv", "metro_slopes.csv",
    )),
    ("Figure4", "FIG4_DATA", ("cascade_country_v2.csv", "cascade_metro_named.csv")),
    ("Figure5", "FIG5_DATA", (
        "income_group_exposure.csv", "concentration_curve.csv", "concentration_indices.json",
        "protection_by_country.csv", "protection_bins_income.csv", "protection_by_income.csv",
        "quintile_within_country.csv", "quintile_pooled.csv", "inheritance_country.csv",
        "fig05_stats.json",
    )),
    ("Supplementary_Figure1", "FIG1_DATA", ("crosscheck_stations_vs_iea.csv",)),
    ("Supplementary_Figure1", "FIG3_DATA", ("coverage_table.csv",)),
    ("Supplementary_Figure2", "FIG4_DATA", ("islanding_global.json", "islanding_country.csv")),
    ("Supplementary_Figure3", "FIG5_DATA", ("protection_sensitivity.json",)),
)

PUBLIC_NAMES = {
    "fig02_stats_v2.json": "exposure_decomposition.json",
    "cascade_country_v2.csv": "cascade_by_country.csv",
}
CASCADE_FIELDS = (
    "A_rp100", "direct_rp100", "total_rp100", "grid_add_rp100", "A_redundant_rp100",
    "A_threshold_rp100", "A_band_rp100", "mech_corr_logdirect_A", "n_countries_mech",
)
CURVE_FIELDS = (
    "rp", "direct", "total", "A_fragility", "A_redundant", "total_redundant",
    "A_threshold", "A_soft", "A_hard", "A_beta04", "A_beta08",
)


def copy_public_tables(figures):
    figures.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(paths.get("NE110M"), figures / "world_countries_110m.geojson")
    for figure, source_key, names in PUBLIC_TABLES:
        source = Path(paths.get(source_key))
        destination = figures / figure / "data"
        destination.mkdir(parents=True, exist_ok=True)
        for name in names:
            shutil.copyfile(source / name, destination / PUBLIC_NAMES.get(name, name))


def export_station_density(destination):
    """Figure 1 uses occupied hexagons in Robinson coordinates."""
    extent = (plot_style.XL, plot_style.XR, plot_style.YB, plot_style.YT)
    for key, network in (("DB_CHARGING", "charging"), ("DB_FUEL", "fuel")):
        stations = pd.read_parquet(paths.get(key), columns=["lon", "lat"])
        x, y = plot_style.to_rob(stations.lon.to_numpy(), stations.lat.to_numpy())
        axes = Figure().subplots()
        cells = axes.hexbin(x, y, gridsize=440, mincnt=1, extent=extent)
        centres = cells.get_offsets()
        counts = pd.DataFrame({
            "x": np.round(centres[:, 0], 1),
            "y": np.round(centres[:, 1], 1),
            "stations": cells.get_array().astype(int),
        })
        counts.to_csv(destination / f"hexbin_{network}_stations.csv", index=False)
        vertices = cells.get_paths()[0].vertices[:6]
        polygon = [[round(float(x), 2), round(float(y), 2)] for x, y in vertices]
        with (destination / f"hexagon_{network}.json").open("w", encoding="utf-8") as target:
            json.dump(polygon, target)


def export_cascade_tables(destination):
    """Select the global statistics, scenario columns and public city centroids for Figure 4."""
    source = Path(paths.get("FIG4_DATA"))
    with (source / "cascade_stats.json").open(encoding="utf-8") as stream:
        stats = json.load(stream)
    selected = {key: stats[key] for key in CASCADE_FIELDS}
    with (destination / "cascade_stats.json").open("w", encoding="utf-8") as target:
        json.dump(selected, target, indent=2)

    curve = "cascade_rp_curve.csv"
    with (source / curve).open(encoding="utf-8", newline="") as stream:
        with (destination / curve).open("w", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=CURVE_FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(csv.DictReader(stream))

    columns = ["metro", "GCPNT_LAT", "GCPNT_LON"]
    cities = pd.read_csv(Path(paths.get("FIG2_DATA")) / "city_map_exposure.csv", usecols=columns)
    cities.to_csv(destination / "urban_centre_coordinates.csv", columns=columns, index=False)


def distance_curve(distances_m):
    """Sample up to 4,000 ranks. All chargers contribute to the denominator."""
    distances_km = np.sort(distances_m[np.isfinite(distances_m)]) / 1000.0
    cumulative_share = np.arange(1, len(distances_km) + 1) / len(distances_m)
    positions = np.linspace(0, len(distances_km) - 1, min(4000, len(distances_km))).astype(int)
    return pd.DataFrame({
        "distance_km": distances_km[positions],
        "cumulative_share": cumulative_share[positions],
    })


def export_feeder_distances(destination):
    """Supplementary Figure 1 reports two distance distributions and their summaries."""
    assignments = pd.read_parquet(
        Path(paths.get("FIG4_DATA")) / "charger_network_assignment.parquet",
        columns=["attach_dist_m", "d1_net_m"],
    )
    attachment = assignments["attach_dist_m"].to_numpy()
    network = assignments["d1_net_m"].to_numpy()
    distance_curve(attachment).to_csv(destination / "cdf_charger_to_power_line.csv", index=False)
    distance_curve(network).to_csv(destination / "cdf_charger_to_feeder_site.csv", index=False)

    within_limit = network <= 10_000
    summary = {
        "n_chargers": len(assignments),
        "median_attach_km": float(np.nanmedian(attachment)) / 1000.0,
        "median_network_km_within_10km": float(np.nanmedian(np.where(within_limit, network, np.nan))) / 1000.0,
        "within_10km_pct": float(within_limit.mean() * 100),
    }
    with (destination / "coupling_summary.json").open("w", encoding="utf-8") as target:
        json.dump(summary, target, indent=1)


def copy_local_maps(figures):
    for figure, source_key, source_prefix, target_prefix, cities in (
        ("Figure3", "FIG3_DATA", "citymap", "citymap", (2046, 1677)),
        ("Figure4", "FIG4_DATA", "vignette_v2", "feeder_map", (2288, 12703, 10912, 1048)),
    ):
        destination = figures / figure / "data" / "private"
        destination.mkdir(parents=True, exist_ok=True)
        source = Path(paths.get(source_key))
        for city in cities:
            shutil.copyfile(source / f"{source_prefix}_{city}.npz", destination / f"{target_prefix}_{city}.npz")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figures_directory", type=Path)
    parser.add_argument("--local-maps", action="store_true", help="also copy licensed city extracts to data/private")
    args = parser.parse_args()
    figures = args.figures_directory
    copy_public_tables(figures)
    export_station_density(figures / "Figure1" / "data")
    export_cascade_tables(figures / "Figure4" / "data")
    export_feeder_distances(figures / "Supplementary_Figure1" / "data")
    if args.local_maps:
        copy_local_maps(figures)
    print("figure source data exported to", figures)


if __name__ == "__main__":
    main()
