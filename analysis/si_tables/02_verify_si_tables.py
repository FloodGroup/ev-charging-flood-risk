#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary tables, step 2 -- verification against the manuscript.

Re-reads the CSVs written by step 1 and asserts, at the manuscript's
display rounding, the values printed in Supplementary Tables 2-8 of the
supplementary information document (spot rows, orderings, totals and every
headline the main text quotes from a table). Exits non-zero on any
mismatch, so a rebuilt package cannot silently diverge from the submitted
tables.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("SI_TABLES_DATA")
_fails = []

def ok(nm, a, b, tol=5e-3):
    good = abs(float(a) - float(b)) <= tol
    print(f"{'OK  ' if good else 'FAIL'} {nm}: {a} vs {b}")
    if not good:
        _fails.append(nm)

def oks(nm, a, b):
    good = str(a) == str(b)
    print(f"{'OK  ' if good else 'FAIL'} {nm}: {a} vs {b}")
    if not good:
        _fails.append(nm)

# ---- ST2: top-30 country exposure ----------------------------------------
t = pd.read_csv(os.path.join(DATA, "si_table2_country_exposure_top30.csv"))
ok("ST2 rows", len(t), 30, .5)
oks("ST2 first country", t.country.iloc[0], "China")
r = t[t.country == "China"].iloc[0]
ok("ST2 China chargers", r.n_ev, 284640, .5)
ok("ST2 China fuel", r.n_gas, 196756, .5)
ok("ST2 China ev%", round(r.ev_pct, 1), 25.8, .05)
ok("ST2 China ratio", round(r.ratio, 2), 0.96, .005)
r = t[t.country == "United Kingdom"].iloc[0]
ok("ST2 UK ratio x1.98", round(r.ratio, 2), 1.98, .005)
r = t[t.country == "Vietnam"].iloc[0]
ok("ST2 Vietnam ev% (48.0)", round(r.ev_pct, 1), 48.0, .05)

# ---- ST3: cohort trends ----------------------------------------------------
t = pd.read_csv(os.path.join(DATA, "si_table3_cohort_trends.csv"))
ok("ST3 rows (15 markets)", len(t), 15, .5)
oks("ST3 first country", t.country.iloc[0], "United Kingdom")
r = t.iloc[0]
ok("ST3 UK n", r.n, 27157, .5)
ok("ST3 UK slope +0.73", round(r.slope_pp_yr, 2), 0.73, .005)
ok("ST3 UK CI low", round(r.lo, 2), 0.63, .005)
ok("ST3 UK CI high", round(r.hi, 2), 0.83, .005)
ok("ST3 UK FE slope +0.48", round(r.fe_slope_pp_yr, 2), 0.48, .005)
ok("ST3 median slope -0.09", round(t.slope_pp_yr.median(), 2), -0.09, .005)
au = t[t.country == "Australia"].iloc[0]
ok("ST3 Australia has no FE column",
   0 if pd.isna(au.fe_slope_pp_yr) else 1, 0, .5)

# ---- ST4: cascade by country ----------------------------------------------
t = pd.read_csv(os.path.join(DATA, "si_table4_cascade_by_country.csv"))
ok("ST4 rows (35 countries)", len(t), 35, .5)
oks("ST4 first country", t.country.iloc[0], "Switzerland")
r = t.iloc[0]
ok("ST4 Switzerland A x1.77", round(r.A, 2), 1.77, .005)
ok("ST4 Switzerland backup x1.41", round(r.A_backup, 2), 1.41, .005)
r = t[t.country == "China"].iloc[0]
ok("ST4 China grid-added 8,725", round(r.grid_add), 8725, .5)
ok("ST4 China total 82,231", round(r.total), 82231, .5)
ok("ST4 grid-added column sums to <= global 26,878",
   t.grid_add.sum() <= 26878.5, True, 0)

# ---- ST5: income groups ----------------------------------------------------
t = pd.read_csv(os.path.join(DATA, "si_table5_income_groups.csv")
                ).set_index("inc_grp")
ok("ST5 rows", len(t), 4, .5)
ok("ST5 HIC chargers", t.loc["HIC", "ev_n"], 585013, .5)
ok("ST5 HIC ratio x1.39", round(t.loc["HIC", "ratio"], 2), 1.39, .005)
ok("ST5 UMC effective 70,617", t.loc["UMC", "ev_eff"], 70617, .5)
ok("ST5 UMC share 88.4",
   round(t.loc["UMC", "share_ev_eff_stock"], 1), 88.4, .05)
ok("ST5 LIC effective count is 2", t.loc["LIC", "ev_eff"], 2, .5)
ok("ST5 HIC protected 85.3",
   round(t.loc["HIC", "pct_protected"], 1), 85.3, .05)
ok("ST5 UMC median protection RP20",
   t.loc["UMC", "prot_median_at_exposed"], 20, .5)

# ---- ST6: pre-transition markets -------------------------------------------
t = pd.read_csv(os.path.join(DATA, "si_table6_pretransition_fuel.csv"))
ok("ST6 rows (45 countries)", len(t), 45, .5)
oks("ST6 first country", t.country.iloc[0], "Bangladesh")
r = t.iloc[0]
ok("ST6 Bangladesh exposed 47.1", round(r.gas_pct, 1), 47.1, .05)
ok("ST6 Bangladesh effective 46.8", round(r.gas_eff_pct, 1), 46.8, .05)
r = t[t.country == "Indonesia"].iloc[0]
ok("ST6 Indonesia effective 0.0", round(r.gas_eff_pct, 1), 0.0, .05)
ok("ST6 population 3.0 bn", round(t.pop_m.sum() / 1000, 1), 3.0, .05)
w = (t.gas_pct * t.n_gas).sum() / t.n_gas.sum()
ok("ST6 weighted exposure 9.9%", round(w, 1), 9.9, .06)

# ---- ST7: largest charging centres -----------------------------------------
t = pd.read_csv(os.path.join(DATA, "si_table7_city_exposure_top30.csv"))
ok("ST7 rows", len(t), 30, .5)
oks("ST7 first centre", t.urban_centre.iloc[0], "London")
r = t.iloc[0]
ok("ST7 London chargers", r.n_ev, 20907, .5)
ok("ST7 London ev% 11.8", round(r.ev_pct, 1), 11.8, .05)
ok("ST7 London fuel 866", r.n_gas, 866, .5)
r = t[t.urban_centre == "Shanghai"].iloc[0]
ok("ST7 Shanghai ev% 77.8", round(r.ev_pct, 1), 77.8, .05)
r = t[t.urban_centre == "Dordrecht"].iloc[0]
ok("ST7 Dordrecht ev% 75.9", round(r.ev_pct, 1), 75.9, .05)

# ---- ST8: top-20 amplification ---------------------------------------------
t = pd.read_csv(os.path.join(DATA,
                             "si_table8_metro_amplification_top20.csv"))
ok("ST8 rows", len(t), 20, .5)
oks("ST8 first centre", t.urban_centre.iloc[0], "Dusseldorf")
r = t.iloc[0]
ok("ST8 Dusseldorf A x3.03", round(r.A, 2), 3.03, .005)
ok("ST8 Dusseldorf direct 72", r.direct, 72, .5)
ok("ST8 Dusseldorf total 218", round(r.total), 218, .5)
oks("ST8 last centre (Berlin)", t.urban_centre.iloc[-1], "Berlin")
ok("ST8 Berlin A x1.46", round(t.A.iloc[-1], 2), 1.46, .005)

print("-" * 60)
if _fails:
    print(f"RESULT: {len(_fails)} check(s) FAILED")
    sys.exit(1)
print("ALL CHECKS PASSED")
