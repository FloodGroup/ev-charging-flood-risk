#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 5 -- city map tables, continental aggregates and the
Kitagawa decomposition.

The Kitagawa (1955) identity splits the global exposure gap over country
strata c into a composition term (where each network is located) and a
siting term (how exposed each network is within a country):

  gap = sum_c (w1c - w2c) * (p1c + p2c)/2    composition
      + sum_c (w1c + w2c)/2 * (p1c - p2c)    within-country siting

where w is a country's share of the respective network (1 = charging,
2 = fuel) and p its exposed share. Strata cover countries hosting both
networks (100% of chargers, 99.8% of fuel stations).

Outputs: data/city_map_exposure.csv (centres with >= 20 chargers),
         city_map_ratio.csv, continent_exposure.csv, fig02_stats_v2.json
"""
import os, sys, json
import numpy as np, pandas as pd, pyogrio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG2_DATA")
UC = P.get("UCDB_GPKG")
sd = pd.read_parquet(f"{DATA}/station_depths.parquet", columns=["lon","lat","country","kind","d100"])
sd["metro"] = np.concatenate([np.load(f"{DATA}/metro_p{i}.npy") for i in range(4)])
sd["exp"] = sd.d100 > 0
ev = sd.kind == "ev"

# ---------- city tables (with UCDB centroids) ----------
uc = pyogrio.read_dataframe(UC, columns=["ID_HDC_G0","UC_NM_MN","CTR_MN_NM","GCPNT_LAT","GCPNT_LON"],
                            read_geometry=False).rename(columns={"ID_HDC_G0":"metro"})
m = sd[sd.metro >= 0]
g = m.groupby(["metro","kind"], observed=True)["exp"].agg(n="size", ne="sum").reset_index()
p = g.pivot(index="metro", columns="kind", values=["n","ne"])
p.columns = ["n_ev","n_gas","e_ev","e_gas"]; p = p.fillna(0).astype(int)
p["ev_pct"] = np.where(p.n_ev>0, 100*p.e_ev/p.n_ev, np.nan)
p["gas_pct"] = np.where(p.n_gas>0, 100*p.e_gas/p.n_gas, np.nan)
p = p.merge(uc.set_index("metro"), left_index=True, right_index=True, how="left")
a_map = p[p.n_ev >= 20].copy()
a_map.to_csv(f"{DATA}/city_map_exposure.csv")
c_map = p[(p.n_ev>=50)&(p.n_gas>=50)&(p.e_ev>0)&(p.e_gas>0)].copy()
c_map["ratio"] = c_map.ev_pct/c_map.gas_pct
c_map.to_csv(f"{DATA}/city_map_ratio.csv")
print("city tables:", len(a_map), "exposure dots;", len(c_map), "ratio dots")

# ---------- continents via Natural Earth attribute ----------
ne = pyogrio.read_dataframe(P.get("NE10M"),
                            columns=["NAME","CONTINENT"], read_geometry=False)
cont = dict(zip(ne.NAME, ne.CONTINENT))
sd["continent"] = sd.country.map(cont).fillna("")
ct = sd[sd.continent != ""].groupby(["continent","kind"], observed=True)["exp"].agg(n="size", ne="sum").reset_index()
cp = ct.pivot(index="continent", columns="kind", values=["n","ne"])
cp.columns = ["n_ev","n_gas","e_ev","e_gas"]; cp = cp.fillna(0).astype(int)
cp["ev_pct"] = 100*cp.e_ev/cp.n_ev; cp["gas_pct"] = 100*cp.e_gas/cp.n_gas
cp["ratio"] = cp.ev_pct/cp.gas_pct
cp.to_csv(f"{DATA}/continent_exposure.csv")
print(cp.round(2).to_string())

# ---------- Kitagawa decomposition (country strata) ----------
cc = sd.groupby(["country","kind"], observed=True)["exp"].agg(n="size", ne="sum").reset_index()
w = cc.pivot(index="country", columns="kind", values=["n","ne"])
w.columns = ["n_ev","n_gas","e_ev","e_gas"]; w = w.fillna(0)
both = w[(w.n_ev>0)&(w.n_gas>0)].copy()
only = w[~((w.n_ev>0)&(w.n_gas>0))]
cover_ev = both.n_ev.sum()/w.n_ev.sum(); cover_gas = both.n_gas.sum()/w.n_gas.sum()
w1 = both.n_ev/both.n_ev.sum(); w2 = both.n_gas/both.n_gas.sum()
p1 = both.e_ev/both.n_ev; p2 = both.e_gas/both.n_gas
E1 = float((w1*p1).sum())*100; E2 = float((w2*p2).sum())*100
comp = float(((w1-w2)*(p1+p2)/2).sum())*100
sit  = float((((w1+w2)/2)*(p1-p2)).sum())*100
print(f"Kitagawa: EV {E1:.2f} vs GAS {E2:.2f} | gap {E1-E2:.2f}pp = composition {comp:.2f}pp + siting {sit:.2f}pp")
print(f"coverage: EV {100*cover_ev:.2f}% GAS {100*cover_gas:.2f}% of stations in both-type countries")

out = dict(kita_ev=E1, kita_gas=E2, kita_gap=E1-E2, kita_comp=comp, kita_sit=sit,
           kita_comp_share=100*comp/(E1-E2), cover_ev=100*cover_ev, cover_gas=100*cover_gas,
           n_cities_a=len(a_map), n_cities_c=len(c_map))
cont_order = ["Asia","Europe","North America","South America","Africa","Oceania"]
out["continents"] = {k: dict(ev=float(cp.loc[k,"ev_pct"]), gas=float(cp.loc[k,"gas_pct"]),
                             n_ev=int(cp.loc[k,"n_ev"]), n_gas=int(cp.loc[k,"n_gas"]))
                     for k in cont_order if k in cp.index}
json.dump(out, open(f"{DATA}/fig02_stats_v2.json","w"), indent=1)
print(json.dumps(out, indent=1)[:600])
