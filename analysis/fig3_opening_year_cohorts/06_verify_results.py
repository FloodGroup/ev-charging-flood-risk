#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3, step 6 -- verification of every number shown on the figure.

Independently recomputes the cohort statistics from the archive tables and
compares them with data/cohort_stats.json and the stored slope tables:
country count (15), station count (218,560), the UK slope (+0.73 pp/yr),
the median national slope (-0.09), significance counts (national,
fixed-effects and city level), the newest-cohort benchmark comparison
(7 of 15), the early/late shares of the mapped cities Paris and Malaga
(6.3 -> 13.3% and 46.2 -> 18.5%) and of Vancouver, which is quoted in the
text (33.7 -> 18.8%), the city trend annotations drawn on panels c/e/f,
consistency of the city-map extracts, and the 80-city panel (median
+0.00). Prints ALL CHECKS DONE only if every check passed; exits non-zero
otherwise.
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG3_DATA")
_fails = []
def ok(nm, a, b, tol=5e-3):
    good = abs(a-b) <= tol
    print(f"{'OK ' if good else 'FAIL'} {nm}: {a} vs {b}")
    if not good: _fails.append(nm)

S4 = json.load(open(os.path.join(DATA, "cohort_stats.json")))
ss = pd.read_csv(f"{DATA}/country_slopes.csv")
st = pd.read_parquet(f"{DATA}/archive_depths.parquet")
meta = pd.read_parquet(f"{DATA}/archive_stations.parquet", columns=["year","country"])
if len(st) != len(meta):
    raise SystemExit(f"archive tables out of sync: depths {len(st):,} rows vs "
                     f"stations {len(meta):,} rows -- re-run stages 0-2 together")
st = pd.concat([st.reset_index(drop=True), meta.reset_index(drop=True)], axis=1)
st["metro"] = np.load(f"{DATA}/archive_metro.npy")
st["exp"] = (st.d100 > 0).astype(float)
qa = pd.read_csv(f"{DATA}/coverage_table.csv")
st = st[st.country.isin(qa[qa.qa_pass].country)]

ok("n_countries 15", S4["n_countries"], 15, .5)
ok("n_stations", S4["n_stations"], len(st), .5)
# UK slope recompute
g = st[st.country=="United Kingdom"]
x = g.year.values.astype(float); y = g.exp.values
b = ((x-x.mean())*(y-y.mean())).sum()/((x-x.mean())**2).sum()
ok("UK slope +0.733", 100*b, float(ss.loc[ss.country=="United Kingdom","slope_pp_yr"].iloc[0]), 1e-6)
ok("median slope -0.09", S4["slope_median"], float(ss.slope_pp_yr.median()), 1e-9)
ok("pos_sig 1", S4["pos_sig"], int((ss.lo>0).sum()), .5)
ok("neg_sig 2", S4["neg_sig"], int((ss.hi<0).sum()), .5)
ok("newest_below_gas 7", S4["newest_below_gas"], int((ss.newest_pct<ss.stock_gas_pct).sum()), .5)
# city examples
for mid, nm, e_, l_, mapped in [(2046,"Paris",6.3,13.3,True),(1677,"Malaga",46.2,18.5,True),
                                (187,"Vancouver",33.7,18.8,False)]:
    g = st[st.metro==mid]
    early = g[g.year<=2018].exp.mean()*100; late = g[g.year>=2022].exp.mean()*100
    ok(f"{nm} early {e_}%", early, e_, .06); ok(f"{nm} late {l_}%", late, l_, .06)
    if mapped:
        z = np.load(f"{DATA}/citymap_{mid}.npz")["stats"]
        ok(f"{nm} npz early", float(z[1]), early, 1e-6); ok(f"{nm} npz late", float(z[3]), late, 1e-6)
        ok(f"{nm} npz cohort sizes", float(z[0] + z[2]), float((g.year<=2018).sum() + (g.year>=2022).sum()), .5)
g = st[st.metro==1677]
ok("Malaga cohort sizes 13 and 265", float((g.year<=2018).sum()*1000 + (g.year>=2022).sum()), 13265.0, .5)
mm = pd.read_csv(f"{DATA}/metro_slopes.csv")
ok("n_metros 80", S4["n_metros"], len(mm), .5)
ok("metro median 0.00", S4["metro_slope_median"], float(mm.slope_pp_yr.median()), 1e-9)
# significance counts quoted on panels b/c and the drawn city annotations
ok("fe_pos_sig", S4["fe_pos_sig"], int((ss.fe_lo > 0).sum()), .5)
ok("fe_neg_sig", S4["fe_neg_sig"], int((ss.fe_hi < 0).sum()), .5)
ok("metro_pos_sig", S4["metro_pos_sig"], int((mm.lo > 0).sum()), .5)
ok("metro_neg_sig", S4["metro_neg_sig"], int((mm.hi < 0).sum()), .5)
ok("Paris trend annotation (+0.60)",
   float(mm.loc[mm.metro == 2046, "slope_pp_yr"].iloc[0]), 0.60, .005)
ok("Malaga trend annotation (-4.14)",
   float(mm.loc[mm.metro == 1677, "slope_pp_yr"].iloc[0]), -4.14, .005)
ok("Vancouver trend quoted in the text (-2.07)",
   float(mm.loc[mm.metro == 187, "slope_pp_yr"].iloc[0]), -2.07, .005)
_rob = f"{DATA}/operator_robustness_centres.csv"
if os.path.exists(_rob):
    rb = pd.read_csv(_rob).set_index("unit")
    ok("Malaga trend since 2019 (-3.64)", float(rb.loc["Málaga", "trend_since_2019"]), -3.64, .005)
    ok("robust city trends (3 of 6)", float((rb.significance_changes == 0).sum()), 3.0, .5)
if _fails:
    print(f"RESULT: {len(_fails)} check(s) FAILED")
    sys.exit(1)
print("ALL CHECKS DONE")
