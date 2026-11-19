#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3, step 3 -- cohort exposure statistics for the QA-passing markets.

Exposure is identical to Figure 2 (90-m cell depth > 0 m at RP100). For
each of the 15 well-covered markets the script computes:

  * five opening-cohort exposure shares with Wilson 95% intervals
    (<=2015, 2016-18, 2019-21, 2022-23, 2024-25);
  * the national trend as a station-level linear probability model
      x_i = alpha + beta * y_i + e_i
    where x_i indicates floodplain location and y_i the registration year,
    with heteroskedasticity-robust (HC1) standard errors; beta is the
    drift in percentage points per year;
  * the within-city fixed-effects variant (Frisch-Waugh demeaning by
    urban centre, over centres with >= 50 dated stations);
  * city-level trends for the 80 centres with >= 200 dated stations and
    >= 5 cohort years.

The fuel benchmark is each country's full current fuel-network exposure
from the Figure 2 pipeline.

Outputs: data/cohort_curves.csv, country_slopes.csv, metro_slopes.csv,
         cohort_stats.json
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG3_DATA")
FIG2 = P.get("FIG2_DATA")

st = pd.read_parquet(f"{DATA}/archive_depths.parquet")
meta = pd.read_parquet(f"{DATA}/archive_stations.parquet", columns=["year","country"])
if len(st) != len(meta):
    raise SystemExit(f"archive tables out of sync: depths {len(st):,} rows vs "
                     f"stations {len(meta):,} rows -- re-run stages 0-2 together")
st = pd.concat([st.reset_index(drop=True), meta.reset_index(drop=True)], axis=1)
st["metro"] = np.load(f"{DATA}/archive_metro.npy")
st["exp"] = (st.d100 > 0).astype(float)
qa = pd.read_csv(f"{DATA}/coverage_table.csv")
COUNTRIES = qa[qa.qa_pass].country.tolist()
st = st[st.country.isin(COUNTRIES)].copy()

COHORTS = [(2010, 2015, "≤2015"), (2016, 2018, "2016–18"), (2019, 2021, "2019–21"),
           (2022, 2023, "2022–23"), (2024, 2025, "2024–25")]

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n
    den = 1 + z*z/n
    c = (p + z*z/(2*n)) / den
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return (100*(c-h), 100*(c+h))

def ols_slope(y, x, w=None):
    """LPM slope + HC1 robust SE (bivariate)."""
    x = x.astype(float); y = y.astype(float)
    xm = x - x.mean(); ym = y - y.mean()
    sxx = (xm*xm).sum()
    if sxx == 0: return np.nan, np.nan
    b = (xm*ym).sum() / sxx
    e = ym - b*xm
    n = len(y); k = 2
    se = np.sqrt((n/(n-k)) * (xm*xm*e*e).sum() / sxx**2)
    return b, se

rows_c, rows_s = [], []
ref = pd.read_csv(f"{FIG2}/country_exposure.csv").set_index("country")
for c, g in st.groupby("country"):
    for y0, y1, lab in COHORTS:
        gg = g[(g.year >= y0) & (g.year <= y1)]
        k = int(gg.exp.sum()); n = len(gg)
        lo, hi = wilson(k, n)
        rows_c.append(dict(country=c, cohort=lab, mid=(y0+y1)/2, n=n,
                           exp_pct=100*k/n if n else np.nan, lo=lo, hi=hi))
    b, se = ols_slope(g.exp.values, g.year.values)
    # within-metro FE (metros with >=50 archive stations in this country)
    gm = g[g.metro >= 0].copy()
    cnt = gm.metro.value_counts()
    gm = gm[gm.metro.isin(cnt[cnt >= 50].index)]
    if len(gm) >= 300 and gm.metro.nunique() >= 3:
        ye = gm.exp - gm.groupby("metro").exp.transform("mean")
        xe = gm.year - gm.groupby("metro").year.transform("mean")
        bf, sef = ols_slope(ye.values + 0.0, xe.values + 0.0)
    else:
        bf, sef = np.nan, np.nan
    rows_s.append(dict(country=c, n=len(g), slope_pp_yr=100*b, se=100*se,
                       lo=100*(b-1.96*se), hi=100*(b+1.96*se),
                       fe_slope_pp_yr=100*bf if bf==bf else np.nan,
                       fe_lo=100*(bf-1.96*sef) if bf==bf else np.nan,
                       fe_hi=100*(bf+1.96*sef) if bf==bf else np.nan,
                       n_fe=len(gm), n_metros_fe=int(gm.metro.nunique()),
                       stock_ev_pct=float(ref.ev_pct.get(c, np.nan)),
                       stock_gas_pct=float(ref.gas_pct.get(c, np.nan))))
cc = pd.DataFrame(rows_c); ss = pd.DataFrame(rows_s).sort_values("n", ascending=False)
newest = cc[cc.cohort == "2024–25"].set_index("country").exp_pct
ss["newest_pct"] = ss.country.map(newest)
cc.to_csv(f"{DATA}/cohort_curves.csv", index=False)
ss.to_csv(f"{DATA}/country_slopes.csv", index=False)

# ---------- metro-level slopes ----------
gm = st[st.metro >= 0]
rows_m = []
for (c, m), g in gm.groupby(["country", "metro"]):
    if len(g) >= 200 and g.year.nunique() >= 5:
        b, se = ols_slope(g.exp.values, g.year.values)
        rows_m.append(dict(country=c, metro=m, n=len(g), exp_pct=100*g.exp.mean(),
                           slope_pp_yr=100*b, lo=100*(b-1.96*se), hi=100*(b+1.96*se)))
mm = pd.DataFrame(rows_m)
names = pd.read_csv(f"{FIG2}/city_map_exposure.csv")[["metro","UC_NM_MN","CTR_MN_NM","GCPNT_LAT","GCPNT_LON"]]
mm = mm.merge(names, on="metro", how="left").sort_values("slope_pp_yr")
mm.to_csv(f"{DATA}/metro_slopes.csv", index=False)

below_gas = int((ss.newest_pct < ss.stock_gas_pct).sum())
pos_sig = int(((ss.lo > 0)).sum()); neg_sig = int(((ss.hi < 0)).sum())
fe_pos = int((ss.fe_lo > 0).sum()); fe_neg = int((ss.fe_hi < 0).sum())
out = dict(n_countries=len(ss), n_stations=int(ss.n.sum()),
           slope_median=float(ss.slope_pp_yr.median()),
           pos_sig=pos_sig, neg_sig=neg_sig, fe_pos_sig=fe_pos, fe_neg_sig=fe_neg,
           fe_slope_median=float(ss.fe_slope_pp_yr.median()),
           n_metros=len(mm), metro_pos_sig=int((mm.lo > 0).sum()),
           metro_neg_sig=int((mm.hi < 0).sum()),
           metro_slope_median=float(mm.slope_pp_yr.median()) if len(mm) else np.nan,
           newest_below_gas=below_gas)
json.dump(out, open(os.path.join(DATA, "cohort_stats.json"), "w"), indent=1)
print(ss[["country","n","slope_pp_yr","lo","hi","fe_slope_pp_yr","fe_lo","fe_hi","stock_ev_pct","stock_gas_pct"]].round(3).to_string(index=False))
print(json.dumps(out, indent=1))
