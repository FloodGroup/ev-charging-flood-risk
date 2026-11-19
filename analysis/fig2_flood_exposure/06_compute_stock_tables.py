#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 6 -- fuel-network city table, exposed-stock table and the
country scatter table.

Produces the remaining panel tables: the fuel-station city map (centres
with >= 20 fuel stations, mirroring the charger map), the exposed-stock
ranking (every country's count of chargers inside the 100-year floodplain,
for panel g) and the country scatter table (country_exposure.csv joined
with Natural Earth continent attributes, for panel e).

Outputs: data/city_map_exposure_gas.csv, exposed_stock.csv,
         country_scatter.csv
"""
import os, sys
import numpy as np, pandas as pd, pyogrio
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG2_DATA")
UC = P.get("UCDB_GPKG")
sd = pd.read_parquet(f"{DATA}/station_depths.parquet", columns=["country", "kind", "d100"])
sd["metro"] = np.concatenate([np.load(f"{DATA}/metro_p{i}.npy") for i in range(4)])
sd["exp"] = sd.d100 > 0

# ---------- gasoline city table (mirror of city_map_exposure.csv which is n_ev>=20) ----------
uc = pyogrio.read_dataframe(UC, columns=["ID_HDC_G0", "UC_NM_MN", "CTR_MN_NM",
                                         "GCPNT_LAT", "GCPNT_LON"],
                            read_geometry=False).rename(columns={"ID_HDC_G0": "metro"})
m = sd[sd.metro >= 0]
g = m.groupby(["metro", "kind"], observed=True)["exp"].agg(n="size", ne="sum").reset_index()
p = g.pivot(index="metro", columns="kind", values=["n", "ne"])
p.columns = ["n_ev", "n_gas", "e_ev", "e_gas"]
p = p.fillna(0).astype(int)
p["ev_pct"] = np.where(p.n_ev > 0, 100 * p.e_ev / p.n_ev, np.nan)
p["gas_pct"] = np.where(p.n_gas > 0, 100 * p.e_gas / p.n_gas, np.nan)
p = p.merge(uc.set_index("metro"), left_index=True, right_index=True, how="left")
gmap = p[p.n_gas >= 20].copy()
gmap.to_csv(f"{DATA}/city_map_exposure_gas.csv")
print("gas city dots:", len(gmap), "| ev city dots (existing):", (p.n_ev >= 20).sum())

# in-metro coverage shares for caption
for k, col in [("ev", "n_ev"), ("gas", "n_gas")]:
    tot = (sd.kind == k).sum()
    print(f"{k}: {100 * p[col].sum() / tot:.1f}% of stations inside UCDB metros (n={tot:,})")

# ---------- exposed stock by country (for panel g) ----------
stk = (sd[(sd.kind == "ev") & sd.exp].groupby("country", observed=True).size()
         .sort_values(ascending=False).rename("exp_ev").reset_index())
total_exp = int(sd[sd.kind == "ev"].exp.sum())
stk["share"] = 100 * stk.exp_ev / total_exp
stk.to_csv(f"{DATA}/exposed_stock.csv", index=False)
print(f"exposed stock: total {total_exp:,} | China+NL share "
      f"{stk[stk.country.isin(['China','Netherlands'])].share.sum():.1f}%")

# ---------- country scatter table (add continent) ----------
ne = pyogrio.read_dataframe(P.get("NE10M"),
                            columns=["NAME", "CONTINENT"], read_geometry=False)
cont = dict(zip(ne.NAME, ne.CONTINENT))
cy = pd.read_csv(f"{DATA}/country_exposure.csv")
cy["continent"] = cy.country.map(cont)
cy.to_csv(f"{DATA}/country_scatter.csv", index=False)
print(len(cy), "countries |", cy.continent.isna().sum(), "unmapped")
print(cy.nlargest(8, "n_ev")[["country", "continent", "n_ev", "ev_pct", "gas_pct", "ratio"]]
        .round(1).to_string(index=False))
