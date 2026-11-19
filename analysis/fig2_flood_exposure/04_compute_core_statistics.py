#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 4 -- core exposure statistics from the sampled depths.

Computes, entirely from data/station_depths.parquet:
  * global exposure by return period (depth > 0 m) and by depth threshold
    at RP100, with a log-delta-method 95% CI for the charging:fuel ratio
      R = p1/p2,  SE(ln R) = sqrt((1-p1)/(n1 p1) + (1-p2)/(n2 p2));
  * the per-country table at RP100 (countries with >= 1,000 stations of
    each network);
  * within-urban-centre paired ratios (centres with >= 50 stations of each
    network and at least one exposed station of each).

Outputs: data/fig02_stats.json, country_exposure.csv, metro_ratios.csv,
         rp_curve.csv, depth_curve.csv
"""
import os, sys, json
import numpy as np, pandas as pd, pyogrio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG2_DATA")
UC = P.get("UCDB_GPKG")
RPS = ["010","020","050","075","100","200","500"]

sd = pd.read_parquet(f"{DATA}/station_depths.parquet")
sd["metro"] = np.concatenate([np.load(f"{DATA}/metro_p{i}.npy") for i in range(4)])
ev = sd.kind == "ev"

def ratio_ci(p1, n1, p2, n2):
    import math
    lr = math.log(p1 / p2)
    se = math.sqrt((1-p1)/(n1*p1) + (1-p2)/(n2*p2))
    return math.exp(lr - 1.96*se), math.exp(lr + 1.96*se)

n1, n2 = int(ev.sum()), int((~ev).sum())
rows = []
for rp in RPS:
    p1 = float((sd.loc[ev, f"d{rp}"] > 0).mean()); p2 = float((sd.loc[~ev, f"d{rp}"] > 0).mean())
    lo, hi = ratio_ci(p1, n1, p2, n2)
    rows.append(dict(rp=int(rp), ev_pct=100*p1, gas_pct=100*p2, ratio=p1/p2, lo=lo, hi=hi))
rp_curve = pd.DataFrame(rows); rp_curve.to_csv(f"{DATA}/rp_curve.csv", index=False)

drows = []
for th in [0.0, 0.15, 0.3, 0.5, 0.75, 1.0, 1.5]:
    cmp_ = (sd["d100"] > 0) if th == 0 else (sd["d100"] >= th)
    p1 = float(cmp_[ev].mean()); p2 = float(cmp_[~ev].mean())
    lo, hi = ratio_ci(p1, n1, p2, n2)
    drows.append(dict(thresh=th, ev_pct=100*p1, gas_pct=100*p2, ratio=p1/p2, lo=lo, hi=hi))
depth_curve = pd.DataFrame(drows); depth_curve.to_csv(f"{DATA}/depth_curve.csv", index=False)

exp100 = sd["d100"] > 0
ct = sd.groupby(["country", "kind"], observed=True)["d100"].agg(
    n="size", n_exp=lambda s: int((s > 0).sum())).reset_index()
piv = ct.pivot(index="country", columns="kind", values=["n", "n_exp"])
piv.columns = ["n_ev", "n_gas", "exp_ev", "exp_gas"]
piv = piv.fillna(0).astype(int)
piv = piv[(piv.n_ev >= 1000) & (piv.n_gas >= 1000)].copy()
piv["ev_pct"] = 100 * piv.exp_ev / piv.n_ev
piv["gas_pct"] = 100 * piv.exp_gas / piv.n_gas
piv["ratio"] = np.where(piv.gas_pct > 0, piv.ev_pct / piv.gas_pct, np.nan)
piv.sort_values("n_ev", ascending=False).to_csv(f"{DATA}/country_exposure.csv")

m = sd[sd.metro >= 0]
mg = m.groupby(["metro", "kind"], observed=True)["d100"].agg(
    n="size", n_exp=lambda s: int((s > 0).sum())).reset_index()
mp = mg.pivot(index="metro", columns="kind", values=["n", "n_exp"])
mp.columns = ["n_ev", "n_gas", "exp_ev", "exp_gas"]
mp = mp.fillna(0).astype(int)
mp = mp[(mp.n_ev >= 50) & (mp.n_gas >= 50) & ((mp.exp_ev + mp.exp_gas) > 0)].copy()
mp["ev_pct"] = 100 * mp.exp_ev / mp.n_ev
mp["gas_pct"] = 100 * mp.exp_gas / mp.n_gas
valid = (mp.exp_gas > 0) & (mp.exp_ev > 0)
mp["ratio"] = np.where(mp.gas_pct > 0, mp.ev_pct / mp.gas_pct, np.nan)
names = pyogrio.read_dataframe(UC, columns=["ID_HDC_G0","UC_NM_MN","CTR_MN_NM"], read_geometry=False)
mp = mp.merge(names.rename(columns={"ID_HDC_G0": "metro"}).set_index("metro"),
              left_index=True, right_index=True, how="left")
mp.to_csv(f"{DATA}/metro_ratios.csv")
mr = mp.loc[valid, "ratio"].values
share_gt1 = float((mr > 1).mean()) * 100
med = float(np.median(mr))

amap = sd[(sd.kind == "ev") & (sd.src == "AMAP")]
gp_ev = sd[(sd.kind == "ev") & (sd.src == "GP")]
stats = dict(
    n_ev=n1, n_gas=n2,
    ev_exp100=int((sd.loc[ev, "d100"] > 0).sum()),
    gas_exp100=int((sd.loc[~ev, "d100"] > 0).sum()),
    rp100_ev=rp_curve.loc[rp_curve.rp == 100, "ev_pct"].iloc[0],
    rp100_gas=rp_curve.loc[rp_curve.rp == 100, "gas_pct"].iloc[0],
    rp100_ratio=rp_curve.loc[rp_curve.rp == 100, "ratio"].iloc[0],
    rp100_ci=[rp_curve.loc[rp_curve.rp == 100, "lo"].iloc[0],
              rp_curve.loc[rp_curve.rp == 100, "hi"].iloc[0]],
    d05_ev=depth_curve.loc[depth_curve.thresh == 0.5, "ev_pct"].iloc[0],
    d05_gas=depth_curve.loc[depth_curve.thresh == 0.5, "gas_pct"].iloc[0],
    d05_ratio=depth_curve.loc[depth_curve.thresh == 0.5, "ratio"].iloc[0],
    n_metros=int(valid.sum()), metro_share_gt1=share_gt1, metro_median=med,
    amap_exp100=100*float((amap.d100 > 0).mean()),
    gpev_exp100=100*float((gp_ev.d100 > 0).mean()),
    n_outside_bounds=int(sd.d100.isna().sum()),
)
json.dump(stats, open(f"{DATA}/fig02_stats.json", "w"), indent=1)
print(json.dumps(stats, indent=1))
print("\ntop-12 countries by ev_n:")
print(piv.sort_values("n_ev", ascending=False).head(12)[["n_ev","n_gas","ev_pct","gas_pct","ratio"]].round(2).to_string())
