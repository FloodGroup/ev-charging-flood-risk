#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary tables, step 1 -- build the source tables for Supplementary
Tables 2-8 of the manuscript from the analysis-module statistics.

Supplementary Table 1 (the dataset inventory) is documentation rather than
a computation; it is maintained in DATA_SOURCES.md and reproduced verbatim
in the manuscript. Every other supplementary table is derived here, one
CSV per table, from the same statistics files that the figures draw from:

  ST2  country-level exposure of both networks     <- fig2 country_exposure
  ST3  opening-year cohort trends (15 markets)     <- fig3 country_slopes
  ST4  grid-cascade amplification by country       <- fig4 cascade_country_v2
  ST5  exposure/protection by income group         <- fig5 income tables
  ST6  pre-transition fuel networks (45 markets)   <- fig5 inheritance_country
  ST7  the 30 largest charging urban centres       <- fig2 city_map tables
  ST8  top-20 urban-centre amplification           <- fig4 cascade_metro_named

Values are written at full precision; the manuscript tables round for
display (percentages to 1 dp, A to 2 dp). 02_verify_si_tables.py asserts
the rounded values against the numbers printed in the supplementary
information document.

Outputs: data/si_table{2..8}_*.csv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

OUT = P.ensure_dir(P.get("SI_TABLES_DATA"))
FIG2 = P.get("FIG2_DATA")
FIG3 = P.get("FIG3_DATA")
FIG4 = P.get("FIG4_DATA")
FIG5 = P.get("FIG5_DATA")

def save(df, name):
    p = P.fresh(os.path.join(OUT, name))
    df.to_csv(p, index=False)
    print(f"{name}: {len(df)} rows")

# ---------------------------------------------------------------------------
# ST2 -- country-level flood exposure of the two networks (top 30 by chargers)
# ---------------------------------------------------------------------------
ce = pd.read_csv(os.path.join(FIG2, "country_exposure.csv"))
cs = pd.read_csv(os.path.join(FIG2, "country_scatter.csv"))
st2 = ce.merge(cs[["country", "continent"]], on="country", how="left")
st2 = st2.sort_values("n_ev", ascending=False).head(30)
st2 = st2[["country", "continent", "n_ev", "n_gas", "ev_pct", "gas_pct",
           "ratio"]]
save(st2, "si_table2_country_exposure_top30.csv")

# ---------------------------------------------------------------------------
# ST3 -- opening-year cohort trends in the 15 well-covered markets
# ---------------------------------------------------------------------------
ss = pd.read_csv(os.path.join(FIG3, "country_slopes.csv"))
st3 = ss.sort_values("slope_pp_yr", ascending=False)[
    ["country", "n", "slope_pp_yr", "lo", "hi",
     "fe_slope_pp_yr", "fe_lo", "fe_hi", "newest_pct", "stock_gas_pct"]]
save(st3, "si_table3_cohort_trends.csv")

# ---------------------------------------------------------------------------
# ST4 -- grid-cascade amplification by country (>=1,000 chargers and >=100
#        directly flooded), ordered by amplification
# ---------------------------------------------------------------------------
cy = pd.read_csv(os.path.join(FIG4, "cascade_country_v2.csv"))
st4 = cy[(cy.n >= 1000) & (cy.direct >= 100)].copy()
st4["A"] = st4.total / st4.direct
st4["A_backup"] = st4.total_redundant / st4.direct
st4 = st4.sort_values("A", ascending=False)[
    ["grp", "continent", "n", "direct", "grid_add", "total",
     "A", "A_backup", "direct_pct"]].rename(columns={"grp": "country"})
save(st4, "si_table4_cascade_by_country.csv")

# ---------------------------------------------------------------------------
# ST5 -- exposure, protection and effective risk by World Bank income group
# ---------------------------------------------------------------------------
ig = pd.read_csv(os.path.join(FIG5, "income_group_exposure.csv"))
pb = pd.read_csv(os.path.join(FIG5, "protection_by_income.csv"))
ig = ig[ig.inc_grp != "ALL"].merge(
    pb[["inc_grp", "pct_protected", "prot_median_at_exposed"]],
    on="inc_grp", how="left")
st5 = ig[["inc_grp", "label", "ev_n", "ev_pct", "gas_pct", "ratio",
          "ev_eff", "share_ev_eff_stock", "pct_protected",
          "prot_median_at_exposed"]]
save(st5, "si_table5_income_groups.csv")

# ---------------------------------------------------------------------------
# ST6 -- fuel-network exposure and protection in the pre-transition markets
#        (fewer than 0.1 chargers per fuel station; >=1,000 fuel stations)
# ---------------------------------------------------------------------------
ic = pd.read_csv(os.path.join(FIG5, "inheritance_country.csv"))
st6 = ic[(ic.ev_per_gas < 0.1) & (ic.n_gas >= 1000)].copy()
st6["pop_m"] = st6.pop_nat / 1e6
st6 = st6.sort_values("gas_pct", ascending=False)[
    ["country", "iso3", "inc_grp", "n_gas", "gas_pct", "gas_eff_pct",
     "prot_median_gas_exposed", "pop_m"]]
save(st6, "si_table6_pretransition_fuel.csv")

# ---------------------------------------------------------------------------
# ST7 -- flood exposure of the 30 largest charging urban centres
# ---------------------------------------------------------------------------
cm_ev = pd.read_csv(os.path.join(FIG2, "city_map_exposure.csv"))
st7 = cm_ev.sort_values("n_ev", ascending=False).head(30)
st7 = st7[["UC_NM_MN", "CTR_MN_NM", "n_ev", "ev_pct", "n_gas", "gas_pct"]]
st7.columns = ["urban_centre", "country", "n_ev", "ev_pct", "n_gas",
               "gas_pct"]
save(st7, "si_table7_city_exposure_top30.csv")

# ---------------------------------------------------------------------------
# ST8 -- the 20 urban centres with the highest grid amplification
#        (>=50 directly flooded chargers)
# ---------------------------------------------------------------------------
mt = pd.read_csv(os.path.join(FIG4, "cascade_metro_named.csv"))
st8 = mt[(mt.direct >= 50) & mt.A_total.notna()].copy()
st8["A"] = st8.total / st8.direct
st8 = st8.sort_values("A", ascending=False).head(20)[
    ["UC_NM_MN", "CTR_MN_NM", "n", "direct", "grid_add", "total", "A"]]
st8.columns = ["urban_centre", "country", "n", "direct", "grid_add",
               "total", "A"]
save(st8, "si_table8_metro_amplification_top20.csv")

print("SI TABLES BUILT")
