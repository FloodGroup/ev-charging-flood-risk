# Data access

The analysis reads Google Places and AMap inventories supplied by the researcher under the respective provider terms. Keep names, addresses, identifiers, coordinates, station attributes, per-site flood depths, charger–feeder assignments and route geometries in local storage. The opening-year registry and its city extracts also require independent access and a licence appropriate to the intended use.

| Public figure material | Local inputs and products |
|:--|:--|
| Counts at hexagon centres, annual point stocks and national totals | Station inventories and station-database exports |
| City and country statistics, public urban-centre centroids, grouped opening cohorts | Per-station registrations, flood depths, raster samples and city-map extracts |
| Expected outages by area and return period, cumulative distance distributions | Substation assignments, individual charger distances and network checkpoints |
| Income, protection and wealth-group summaries | Geocoded station covariates, high-resolution source layers and licensed city map outputs |

## Using locally held data

1. Configure the named paths in [Data sources](DATA_SOURCES.md). With the required inputs, use `python run_all.py full` to build the station database, analyses and local maps.
2. If analysis outputs already exist, run `python export_figure_source_data.py figures` to prepare the aggregate plotting tables. To draw city panels locally, use `python export_figure_source_data.py figures --local-maps` and `python analysis/map_backgrounds/01_build_map_backgrounds.py figures` before running the plot scripts.
3. Licensed city extracts and terrain or street layers are read from `figures/Figure3/data/private/` and `figures/Figure4/data/private/`. The figure guides include the complete figures as raster previews from the manuscript. The private directories and locally generated map outputs are excluded from Git. Aggregate country and urban-centre tables remain in the figure’s `data/` folder.

The 0.25° density and exposure table is a separate reanalysis product. Its generation code is provided here. Each row reports counts or expected outages for a cell. The coordinates in this product mark grid cells.
