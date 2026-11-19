# Grid dependence amplifies flood impacts on the global electric vehicle charging network

Analysis code and figure data for the article by Zhi-Yong Long, Rui-Yi Yang, Long Jiang, Lu Chang and Huan-Feng Duan.

The study compares **947,602 public charging stations** with **1,344,301 fuel stations** across seven river-flood return periods. Under the 100-year event, 143,264 charging stations flood directly. Accounting for failures at their feeder substations raises expected outages to 170,142. The backup-feeder counterfactual reduces the amplification factor from 1.19 to 1.13.

[Redraw the composite figures](#redraw-the-composite-figures) · [Full reproduction](#reproduce-the-full-analysis-from-local-data) · [Figure guides](#figure-guides) · [Input sources](DATA_SOURCES.md) · [Data access](RESTRICTED_DATA_NOTICE.md) · [Citation](CITATION.cff)

## Getting started

From the repository root, create a Python 3.12 environment and install the pinned dependencies.

```bash
python -m venv .venv
```

Activate it with `.venv\Scripts\Activate.ps1` in PowerShell or `source .venv/bin/activate` in a POSIX shell.

```bash
python -m pip install -r requirements.txt
```

You can check the installed dependencies with `python check_environment.py`.

### Redraw the composite figures

To redraw the composites using the figure data already in this repository, run this command from the repository root.

```bash
python run_all.py figures
```

This mode reads the prepared figure tables directly. Each plotting script writes PNG, PDF, SVG and PDF-compatible Illustrator files beside the script. The runner records completed and skipped figures in `run_logs/`.

> [!NOTE]
> The public tables are sufficient for Figures 1, 2 and 5 and Supplementary Figures 1–3. Figures 3 and 4 also require local station-level map extracts. The runner skips those two composites until the required files are present. Supplementary Figure 4 is supplied as a static schematic.

### Reproduce the full analysis from local data

With the licensed station inventories and the other inputs listed in [Data sources](DATA_SOURCES.md) available locally, run the full pipeline from the repository root.

```bash
python run_all.py full
```

The runner executes the analyses and verification steps, reusing available checkpoints. It then prepares the city-map backgrounds, automatically calls `export_figure_source_data.py` to prepare the plotting tables and local city extracts, and draws the composites.

## Figure guides

| Figure | Scientific question | Inputs for the composite |
|:--|:--|:--|
| [1 · Infrastructure](figures/Figure1/README.md) | Where are charging and fuel stations concentrated? | Public tables |
| [2 · Direct flood exposure](figures/Figure2/README.md) | How much of the exposure gap persists within the same city? | Public tables |
| [3 · Opening cohorts](figures/Figure3/README.md) | Does floodplain share change with registration year? | Local city maps required for e–f |
| [4 · Grid-mediated outages](figures/Figure4/README.md) | How many dry charging sites lose power through their feeder? | Local feeder maps required for a–d |
| [5 · Protection and inequality](figures/Figure5/README.md) | Where does exposure remain after accounting for protection? | Public tables |
| [S1 · Input audits](figures/Supplementary_Figure1/README.md) | How well do the station, cohort and feeder records support the analysis? | Public tables |
| [S2 · Cascade robustness](figures/Supplementary_Figure2/README.md) | How sensitive are losses to fragility, backup supply and islanding? | Public tables |
| [S3 · Protection sensitivity](figures/Supplementary_Figure3/README.md) | How does effective exposure vary with event magnitude? | Public tables |
| [S4 · Study workflow](figures/Supplementary_Figure4/README.md) | How do the inventories, hazard model and protection layers fit together? | Static workflow schematic |

![Study workflow, from the two station inventories through flood and grid modelling to the distribution of effective exposure](figures/Supplementary_Figure4/Supplementary_Figure4.png)

## Run the analysis

The [path settings](shared/paths.py) locate data on the researcher’s computer. The figure-data exporter uses the [Robinson projection and map extent](shared/plot_style.py) to aggregate station coordinates for Figure 1. The [composite plotting functions](figures/plot_style.py) style and export the published figures.

[Input sources](DATA_SOURCES.md) lists the inputs, local paths and source acknowledgements. Place individually obtained datasets under `source_data/`, or set `EVFR_SOURCE_ROOT` to a separate input directory. Inspect individual overrides and resolved paths with `python shared/paths.py`.

| Additional analysis modes | Inputs and outputs |
|:--|:--|
| `python run_all.py standard` | Recompute statistics and run checks using local station databases, intermediate tables and network checkpoints. |
| `python run_all.py standard fig3` | Run the cohort module against prepared local inputs. |
| `python run_all.py full --fresh` | Move existing checkpoints aside, then rerun the full analysis. |

### From station records to figure data

1. [Station database](analysis/stage0_station_database/) converts AMap coordinates and combines the charging inventories by geographical coverage. The fuel network provides the reference for [direct flood exposure](analysis/fig2_flood_exposure/).
2. The [cohort analysis](analysis/fig3_opening_year_cohorts/), [grid simulation](analysis/fig4_grid_cascade/) and [protection adjustment](analysis/fig5_inequality_protection/) link local inputs to flood depths and city assignments. [Global infrastructure statistics](analysis/fig1_global_infrastructure/) provide network totals and the IEA cross-check.
3. [Supplementary tables](analysis/si_tables/) assemble Tables 2–8. The [grid aggregation script](analysis/reanalysis_products/01_build_station_density_grid.py) creates the 0.25° reanalysis table from local records. The [figure-data exporter](export_figure_source_data.py) copies or aggregates the fields needed for plotting.
4. [Independent tile resampling](analysis/fig2_flood_exposure/resample_public_jrc_tiles.py) checks the mosaicked hazard values against JRC tiles using a locally held station table. The [topology test](analysis/fig4_grid_cascade/13_test_site_failure_topology.py) runs on synthetic networks with `python analysis/fig4_grid_cascade/13_test_site_failure_topology.py`.

## Licence and citation

Code is distributed under the [MIT licence](LICENSE). The authors license their own contributions to the aggregate figure products under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Third-party material remains subject to the terms listed in [Data sources and attribution](DATA_SOURCES.md). Data derived from OpenStreetMap carry the ODbL attribution. Cite the article when using this repository. [Citation metadata](CITATION.cff)
