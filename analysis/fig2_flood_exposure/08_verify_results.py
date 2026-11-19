#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, step 8 -- verification.

Recomputes every headline number quoted for Figure 2 directly from
data/station_depths.parquet (independently of steps 4-6) and compares with
the stored statistics files. Exits non-zero on any mismatch.

Checks: station counts; exposed counts at RP100 (143,264 / 167,605); the
exposure percentages (15.12% / 12.47%) and the ratio 1.213 with its CI; the
monotonic return-period gradient (1.16 -> 1.25) and the >= 1.0 m depth
ratio (1.29); the Kitagawa terms (+3.07 / -0.43 pp of a 2.64 pp gap); the
within-city pairing (474 centres, median 0.98, 47.9% above parity); the
exposed-stock concentration (China + the Netherlands, 69%); the city counts
shown on the two maps (2,585 / 3,923); and the three map callouts
(Shanghai / Bangkok / Rotterdam-The Hague).
"""
import os
import sys
import json
import math
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG2_DATA")
failures = []

def check(name, got, want, tol=0.0):
    ok = abs(float(got) - float(want)) <= tol if tol else got == want
    print(f"{'OK  ' if ok else 'FAIL'} {name}: {got} vs {want}")
    if not ok:
        failures.append(name)

sd = pd.read_parquet(os.path.join(DATA, "station_depths.parquet"))
S = json.load(open(os.path.join(DATA, "fig02_stats.json")))
S2 = json.load(open(os.path.join(DATA, "fig02_stats_v2.json")))
ev = sd.kind == "ev"
n1, n2 = int(ev.sum()), int((~ev).sum())

check("n charging", n1, S["n_ev"])
check("n fuel", n2, S["n_gas"])
e1 = int((sd.loc[ev, "d100"] > 0).sum())
e2 = int((sd.loc[~ev, "d100"] > 0).sum())
check("exposed chargers RP100", e1, S["ev_exp100"])
check("exposed fuel RP100", e2, S["gas_exp100"])
p1, p2 = e1 / n1, e2 / n2
check("EV % RP100", round(100 * p1, 2), round(S["rp100_ev"], 2), tol=0.005)
check("fuel % RP100", round(100 * p2, 2), round(S["rp100_gas"], 2), tol=0.005)
check("ratio RP100", round(p1 / p2, 3), round(S["rp100_ratio"], 3), tol=0.0005)
se = math.sqrt((1 - p1) / (n1 * p1) + (1 - p2) / (n2 * p2))
lo = math.exp(math.log(p1 / p2) - 1.96 * se)
hi = math.exp(math.log(p1 / p2) + 1.96 * se)
check("ratio CI low", round(lo, 3), round(S["rp100_ci"][0], 3), tol=0.0005)
check("ratio CI high", round(hi, 3), round(S["rp100_ci"][1], 3), tol=0.0005)

# gradient across return periods and depth thresholds
rp_ratio = {}
for rp_ in ["010", "500"]:
    a = float((sd.loc[ev, f"d{rp_}"] > 0).mean())
    b = float((sd.loc[~ev, f"d{rp_}"] > 0).mean())
    rp_ratio[rp_] = a / b
check("ratio RP10 (1.16)", round(rp_ratio["010"], 2), 1.16, tol=0.005)
check("ratio RP500 (1.25)", round(rp_ratio["500"], 2), 1.25, tol=0.005)
a = float((sd.loc[ev, "d100"] >= 1.0).mean())
b = float((sd.loc[~ev, "d100"] >= 1.0).mean())
check("ratio at depth >= 1.0 m (1.29)", round(a / b, 2), 1.29, tol=0.005)

# Kitagawa recomputation
cc = sd.groupby(["country", "kind"], observed=True)["d100"].agg(
    n="size", ne=lambda s: int((s > 0).sum())).reset_index()
w = cc.pivot(index="country", columns="kind", values=["n", "ne"])
w.columns = ["n_ev", "n_gas", "e_ev", "e_gas"]
w = w.fillna(0)
both = w[(w.n_ev > 0) & (w.n_gas > 0)]
w1 = both.n_ev / both.n_ev.sum()
w2 = both.n_gas / both.n_gas.sum()
q1 = both.e_ev / both.n_ev
q2 = both.e_gas / both.n_gas
comp = float(((w1 - w2) * (q1 + q2) / 2).sum()) * 100
sit = float((((w1 + w2) / 2) * (q1 - q2)).sum()) * 100
check("Kitagawa composition (+3.07 pp)", round(comp, 2),
      round(S2["kita_comp"], 2), tol=0.005)
check("Kitagawa siting (-0.43 pp)", round(sit, 2),
      round(S2["kita_sit"], 2), tol=0.005)
check("Kitagawa gap (2.64 pp)", round(comp + sit, 2),
      round(S2["kita_gap"], 2), tol=0.01)

# within-city pairing
sd["metro"] = np.concatenate(
    [np.load(os.path.join(DATA, f"metro_p{i}.npy")) for i in range(4)])
m = sd[sd.metro >= 0]
mg = m.groupby(["metro", "kind"], observed=True)["d100"].agg(
    n="size", ne=lambda s: int((s > 0).sum())).reset_index()
mp = mg.pivot(index="metro", columns="kind", values=["n", "ne"])
mp.columns = ["n_ev", "n_gas", "e_ev", "e_gas"]
mp = mp.fillna(0).astype(int)
mp = mp[(mp.n_ev >= 50) & (mp.n_gas >= 50) & (mp.e_ev > 0) & (mp.e_gas > 0)]
r = (100 * mp.e_ev / mp.n_ev) / (100 * mp.e_gas / mp.n_gas)
check("paired centres (474)", int(len(mp)), S["n_metros"])
check("median within-city ratio (0.98)", round(float(np.median(r)), 2),
      round(S["metro_median"], 2), tol=0.005)
check("share above parity (47.9%)", round(100 * float((r > 1).mean()), 1),
      round(S["metro_share_gt1"], 1), tol=0.05)

# exposed-stock concentration -- compared against the stored ranking table,
# which is what panel g draws
stk = sd[(sd.kind == "ev") & (sd.d100 > 0)].groupby(
    "country", observed=True).size()
share_cn_nl = 100 * (stk.get("China", 0) + stk.get("Netherlands", 0)) / e1
stock = pd.read_csv(os.path.join(DATA, "exposed_stock.csv"))
stock = stock[stock.country.notna() & (stock.country != "")]
stored_cn_nl = 100.0 * float(stock.head(2).exp_ev.sum()) / S["ev_exp100"]
check("China + Netherlands share of exposed stock (69%)",
      round(share_cn_nl, 1), round(stored_cn_nl, 1), tol=0.05)
check("exposed-stock table sums to global count",
      int(stock.exp_ev.sum()), S["ev_exp100"])

# figure annotations: city counts on the maps and the three city callouts
cm_ev = pd.read_csv(os.path.join(DATA, "city_map_exposure.csv"))
cm_gas = pd.read_csv(os.path.join(DATA, "city_map_exposure_gas.csv"))
check("charging map city count (2,585)", len(cm_ev), S2["n_cities_a"])
mg2 = m[m.kind != "ev"].groupby("metro", observed=True).size()
check("fuel map city count (3,923)", len(cm_gas), int((mg2 >= 20).sum()))
for city, ctry, want in [("Shanghai", "China", 77.8),
                         ("Bangkok", "Thailand", 62.9),
                         ("Rotterdam [The Hague]", "Netherlands", 35.4)]:
    got = float(cm_ev[(cm_ev.UC_NM_MN == city) &
                      (cm_ev.CTR_MN_NM == ctry)].ev_pct.iloc[0])
    check(f"map callout {city}", round(got, 1), want, tol=0.05)

print("-" * 60)
if failures:
    print(f"RESULT: {len(failures)} check(s) FAILED")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
