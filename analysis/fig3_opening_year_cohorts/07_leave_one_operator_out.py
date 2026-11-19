"""
Figure 3, step 7 -- leave-one-operator-out robustness of cohort trends.

For every country that passes the coverage gates and for every urban centre
with a significant city-level trend, the linear cohort trend (percentage
points per year, HC1 errors) is re-estimated after removing each operator
with at least 20 dated stations in turn. The output lists the full-sample
trend, the largest change produced by a single operator and the number of
operators whose removal changes the significance class of the trend.
For urban centres the table also gives the size and floodplain share of the
early (up to 2018) and late (2022-25) cohorts and the trend re-estimated from
stations registered since 2019 only, which shows how much a trend depends on
a small early cohort.

Outputs: data/operator_robustness_countries.csv
         data/operator_robustness_centres.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared"))
import paths as P  # noqa: E402

DATA = P.get("FIG3_DATA")
GENERIC = ("NA", "(Unknown Operator)", "(Business Owner at Location)")

st = pd.read_parquet(f"{DATA}/archive_depths.parquet")
meta = pd.read_parquet(f"{DATA}/archive_stations.parquet", columns=["year", "country", "Operator"])
st = pd.concat([st.reset_index(drop=True), meta.reset_index(drop=True)], axis=1)
st["metro"] = np.load(f"{DATA}/archive_metro.npy")
st["exp"] = (st.d100 > 0).astype(float)
st["op"] = st.Operator.fillna("NA")


def slope(g):
    x = g.year.values.astype(float)
    y = 100 * g.exp.values
    X = np.c_[np.ones(len(x)), x - x.mean()]
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ b
    n = len(x)
    XtXi = np.linalg.inv(X.T @ X)
    V = XtXi @ (X.T * (e ** 2)) @ X @ XtXi * n / (n - 2)
    se = np.sqrt(V[1, 1])
    return b[1], b[1] - 1.96 * se, b[1] + 1.96 * se


def sig(s):
    return 1 if s[1] > 0 else (-1 if s[2] < 0 else 0)


def loo(g, label):
    full = slope(g)
    vc = g.op.value_counts()
    ops = [o for o in vc.index if vc[o] >= 20 and o not in GENERIC]
    res = [(o, int(vc[o])) + slope(g[g.op != o]) for o in ops]
    worst = max(res, key=lambda r: abs(r[2] - full[0])) if res else None
    flips = [r for r in res if sig(r[2:]) != sig(full)]
    return dict(unit=label, n=len(g), trend=full[0], lo=full[1], hi=full[2], operators_tested=len(ops),
                most_influential_operator=worst[0] if worst else "", operator_stations=worst[1] if worst else 0,
                trend_without=worst[2] if worst else np.nan, lo_without=worst[3] if worst else np.nan,
                hi_without=worst[4] if worst else np.nan, significance_changes=len(flips))


ss = pd.read_csv(f"{DATA}/country_slopes.csv")
pd.DataFrame([loo(st[st.country == c], c) for c in ss.country]).round(3).to_csv(
    f"{DATA}/operator_robustness_countries.csv", index=False)
mm = pd.read_csv(f"{DATA}/metro_slopes.csv")
rows = []
for r in mm.itertuples():
    g = st[st.metro == r.metro]
    if sig(slope(g)) != 0:
        row = loo(g, r.UC_NM_MN)
        early, late, recent = g[g.year <= 2018], g[g.year >= 2022], g[g.year >= 2019]
        rec = slope(recent)
        row.update(n_early=len(early), share_early=100 * early.exp.mean(), n_late=len(late), share_late=100 * late.exp.mean(),
                   n_since_2019=len(recent), trend_since_2019=rec[0], lo_since_2019=rec[1], hi_since_2019=rec[2])
        rows.append(row)
pd.DataFrame(rows).round(3).to_csv(f"{DATA}/operator_robustness_centres.csv", index=False)
print("leave-one-operator-out tables written")
