#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 5, step 2 -- core inequality and protection statistics.

Definitions (identical benchmark to Figure 2):
  exposed:    RP100 cell depth > 0 m
  effective:  exposed AND FLOPROS MerL_Riv < 100, i.e. the 100-year event
              exceeds the local protection standard

The script begins with a hard QA gate: it must reproduce the Figure 2
headline (charging 15.12%, fuel 12.47%, ratio 1.213) from the covariate
table, otherwise it stops -- this guards against any upstream drift.

Outputs: income_group_exposure.csv, concentration_curve.csv,
concentration_indices.json, protection_by_income/country/bins_income.csv,
protection_sensitivity.json, quintile_within_country/pooled/slopes.csv,
inheritance_country/bins.csv, fig05_stats.json, qa_fig02_check.json
"""
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P

DATA = P.get("FIG5_DATA")
INC_ORDER = ['HIC', 'UMC', 'LMC', 'LIC']
INC_LABEL = {'HIC': 'High income', 'UMC': 'Upper-middle', 'LMC': 'Lower-middle', 'LIC': 'Low income'}

st = pd.read_parquet(os.path.join(DATA, 'fig05_station_covariates.parquet'))
st['exp'] = (st['d100'] > 0)
st['eff'] = st['exp'] & (st['prot_rp'] < 100)
ev = st[st['kind'] == 'ev']
gas = st[st['kind'] == 'gas']

# ---------- 0. QA gate: reproduce the Figure 2 headline ----------
qa = {
    'ev_n': int(len(ev)), 'gas_n': int(len(gas)),
    'ev_exp': int(ev['exp'].sum()), 'gas_exp': int(gas['exp'].sum()),
    'ev_pct': round(100 * ev['exp'].mean(), 2), 'gas_pct': round(100 * gas['exp'].mean(), 2),
    'ratio': round(float(ev['exp'].mean() / gas['exp'].mean()), 3),
}
# Match the Figure 2 counts for the study inventory. A local sensitivity run
# with coarsened coordinates can opt into the percentage tolerances below.
if os.environ.get('EVFR_COARSE_OK') == '1':
    assert qa['ev_n'] == 947602 and qa['gas_n'] == 1344301, qa
    assert abs(qa['ev_pct'] - 15.12) < 0.1 and abs(qa['ratio'] - 1.213) < 0.01, qa
else:
    assert qa['ev_n'] == 947602 and qa['gas_n'] == 1344301, qa
    assert qa['ev_exp'] == 143264 and abs(qa['ratio'] - 1.213) < 0.001, qa
json.dump(qa, open(os.path.join(DATA, 'qa_fig02_check.json'), 'w'), indent=1)
print('QA vs FIG02: PASS', qa)

def ci_cols(x1, n1, x2, n2):
    """Wilson 95% CIs for each share + log-delta 95% CI for the ratio."""
    from math import sqrt, log, exp as mexp
    z = 1.96
    def wilson(x, n):
        if n == 0: return (np.nan, np.nan)
        p = x / n; d = 1 + z*z/n
        c = (p + z*z/(2*n)) / d
        h = z*sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
        return (100*(c-h), 100*(c+h))
    lo1, hi1 = wilson(x1, n1); lo2, hi2 = wilson(x2, n2)
    out = {'ev_pct_lo': lo1, 'ev_pct_hi': hi1, 'gas_pct_lo': lo2, 'gas_pct_hi': hi2}
    if min(x1, x2) > 0:
        se = sqrt(1/x1 - 1/n1 + 1/x2 - 1/n2)
        r = (x1/n1)/(x2/n2)
        out['ratio_lo'] = r*mexp(-z*se); out['ratio_hi'] = r*mexp(z*se)
    else:
        out['ratio_lo'] = np.nan; out['ratio_hi'] = np.nan
    return out

# ---------- 1. exposure by income group ----------
rows = []
glob_pop = st.drop_duplicates('iso3')['pop_nat'].sum()
for g in INC_ORDER + ['ALL']:
    sub = st if g == 'ALL' else st[st['inc_grp'] == g]
    e, f = sub[sub['kind'] == 'ev'], sub[sub['kind'] == 'gas']
    pop = sub.drop_duplicates('iso3')['pop_nat'].sum()
    rows.append({
        'inc_grp': g, 'label': INC_LABEL.get(g, 'All'),
        'ev_n': len(e), 'gas_n': len(f),
        'ev_exp': int(e['exp'].sum()), 'gas_exp': int(f['exp'].sum()),
        'ev_pct': 100 * e['exp'].mean(), 'gas_pct': 100 * f['exp'].mean(),
        'ratio': float(e['exp'].mean() / f['exp'].mean()),
        'ev_eff': int(e['eff'].sum()), 'gas_eff': int(f['eff'].sum()),
        'ev_eff_pct': 100 * e['eff'].mean(),
        'share_ev_stock': 100 * len(e) / len(ev),
        'share_ev_exp_stock': 100 * e['exp'].sum() / ev['exp'].sum(),
        'share_ev_eff_stock': 100 * e['eff'].sum() / ev['eff'].sum(),
        'pop_m': pop / 1e6, 'share_pop': 100 * pop / glob_pop,
        'exp_per_million': e['exp'].sum() / (pop / 1e6) if pop > 0 else np.nan,
        **ci_cols(int(e['exp'].sum()), len(e), int(f['exp'].sum()), len(f)),
    })
pd.DataFrame(rows).to_csv(os.path.join(DATA, 'income_group_exposure.csv'), index=False)

# ---------- 2. concentration curves (countries by GDP per capita) ----------
cn = st.groupby(['iso3', 'kind']).agg(n=('exp', 'size'), exp=('exp', 'sum'), eff=('eff', 'sum')).reset_index()
cn = cn.pivot(index='iso3', columns='kind', values=['n', 'exp', 'eff']).fillna(0)
cn.columns = [f'{a}_{b}' for a, b in cn.columns]
meta = st.drop_duplicates('iso3').set_index('iso3')[['gdppc_nat', 'pop_nat', 'inc_grp', 'country']]
cn = cn.join(meta).dropna(subset=['gdppc_nat', 'pop_nat']).sort_values('gdppc_nat')
for col in ['pop_nat', 'n_ev', 'exp_ev', 'eff_ev', 'n_gas', 'exp_gas']:
    cn[f'cum_{col}'] = cn[col].cumsum() / cn[col].sum()
cn.reset_index().to_csv(os.path.join(DATA, 'concentration_curve.csv'), index=False)

def conc_index(xw, yw):
    """Kakwani concentration index: countries in ascending GDP per capita;
    x = cumulative population share, y = cumulative burden share. A curve
    below the diagonal (burden concentrated at the rich end) gives C > 0
    (pro-rich); above it, C < 0 (pro-poor)."""
    x = np.concatenate([[0], np.cumsum(xw) / np.sum(xw)])
    y = np.concatenate([[0], np.cumsum(yw) / np.sum(yw)])
    return float(1 - 2 * np.trapezoid(y, x))

ci = {}
for name, col in [('ev_stock', 'n_ev'), ('ev_exposed', 'exp_ev'), ('ev_effective', 'eff_ev'),
                  ('gas_stock', 'n_gas'), ('gas_exposed', 'exp_gas'), ('gas_effective', 'eff_gas')]:
    ci[name] = round(conc_index(cn['pop_nat'].values, cn[col].values), 3)
json.dump(ci, open(os.path.join(DATA, 'concentration_indices.json'), 'w'), indent=1)

# ---------- 3. protection adjustment ----------
pb = []
for g in INC_ORDER:
    e = ev[ev['inc_grp'] == g]
    ex = e[e['exp']]
    pb.append({
        'inc_grp': g, 'nominal': int(e['exp'].sum()), 'effective': int(e['eff'].sum()),
        'pct_protected': 100 * (1 - e['eff'].sum() / max(e['exp'].sum(), 1)),
        'prot_median_at_exposed': float(ex['prot_rp'].median()) if len(ex) else np.nan,
        'share_nominal': 100 * e['exp'].sum() / ev['exp'].sum(),
        'share_effective': 100 * e['eff'].sum() / ev['eff'].sum(),
        'gas_nominal': int(gas[gas['inc_grp'] == g]['exp'].sum()),
        'gas_effective': int(gas[gas['inc_grp'] == g]['eff'].sum()),
    })
pd.DataFrame(pb).to_csv(os.path.join(DATA, 'protection_by_income.csv'), index=False)

top = ev[ev['exp']].groupby('iso3').size().nlargest(15).index
pc = []
for c in top:
    e = ev[ev['iso3'] == c]
    ex = e[e['exp']]
    pc.append({'iso3': c, 'country': e['country'].iat[0], 'inc_grp': e['inc_grp'].iat[0],
               'nominal': int(e['exp'].sum()), 'effective': int(e['eff'].sum()),
               'prot_median_at_exposed': float(ex['prot_rp'].median()),
               'ev_n': len(e), 'ev_pct': 100 * e['exp'].mean(), 'eff_pct': 100 * e['eff'].mean()})
pd.DataFrame(pc).to_csv(os.path.join(DATA, 'protection_by_country.csv'), index=False)

bins = [0, 10, 50, 100, 1000, 1e9]
blab = ['<10', '10-50', '50-100', '100-1000', '>=1000']
pbn = []
for g in INC_ORDER:
    ex = ev[(ev['inc_grp'] == g) & ev['exp']]
    cut = pd.cut(ex['prot_rp'], bins=bins, labels=blab, right=False)
    v = cut.value_counts(normalize=True).reindex(blab) * 100
    pbn.append({'inc_grp': g, 'n_exposed': len(ex), **{f'p_{b}': v[b] for b in blab}})
pd.DataFrame(pbn).to_csv(os.path.join(DATA, 'protection_bins_income.csv'), index=False)

# sensitivity: event return periods x protection; variant treating prot=0
# as missing (replaced by the national non-zero median)
sens = {}
for rp, dcol in [(50, 'd050'), (100, 'd100'), (200, 'd200'), (500, 'd500')]:
    isev = st['kind'] == 'ev'
    nom = int(((st[dcol] > 0) & isev).sum())
    effn = int(((st[dcol] > 0) & (st['prot_rp'] < rp) & isev).sum())
    sens[f'RP{rp}'] = {'nominal': nom, 'effective': effn,
                       'pct_protected': round(100 * (1 - effn / nom), 1)}
med0 = st[st['prot_rp'] > 0].groupby('iso3')['prot_rp'].median()
alt = st['prot_rp'].where(st['prot_rp'] > 0, st['iso3'].map(med0))
fl_alt = st['exp'] & (alt < 100) & (st['kind'] == 'ev')
sens['zero_as_missing'] = {'effective': int(fl_alt.sum())}
json.dump(sens, open(os.path.join(DATA, 'protection_sensitivity.json'), 'w'), indent=1)

# ---------- 4. within-country quintiles (admin-2 GDP per capita) ----------
qrows, srows = [], []
big = [c for c, n in ev.groupby('iso3').size().items() if n >= 2000]
for c in big:
    e = ev[ev['iso3'] == c].copy()
    if e['gdppc_adm2'].nunique() < 8:
        continue
    r = e['gdppc_adm2'].rank(method='first')
    e['q'] = pd.qcut(r, 5, labels=[1, 2, 3, 4, 5]).astype(int)
    per = e.groupby('q')['exp'].agg(['mean', 'size'])
    for q, row in per.iterrows():
        qrows.append({'iso3': c, 'q': q, 'ev_pct': 100 * row['mean'], 'n': int(row['size'])})
    x = per.index.values.astype(float)
    y = 100 * per['mean'].values
    slope = np.polyfit(x, y, 1)[0]
    srows.append({'iso3': c, 'country': e['country'].iat[0], 'inc_grp': e['inc_grp'].iat[0],
                  'n_ev': len(e), 'slope_pp_per_q': slope,
                  'q5_minus_q1': y[-1] - y[0], 'q1_pct': y[0], 'q5_pct': y[-1]})
pd.DataFrame(qrows).to_csv(os.path.join(DATA, 'quintile_within_country.csv'), index=False)
pd.DataFrame(srows).to_csv(os.path.join(DATA, 'quintile_slopes.csv'), index=False)

# pooled summary: quintiles assigned within each country, then pooled
# globally; one series each for charging and fuel
pool = []
for kind in ['ev', 'gas']:
    k = st[st['kind'] == kind]
    parts = []
    for c, sub in k.groupby('iso3'):
        if len(sub) < 2000 or sub['gdppc_adm2'].nunique() < 8:
            continue
        r = sub['gdppc_adm2'].rank(method='first')
        q = pd.qcut(r, 5, labels=[1, 2, 3, 4, 5]).astype(int)
        parts.append(pd.DataFrame({'q': q, 'exp': sub['exp']}))
    allq = pd.concat(parts)
    for q, subq in allq.groupby('q'):
        pool.append({'kind': kind, 'q': int(q), 'pct': 100 * subq['exp'].mean(), 'n': len(subq)})
pd.DataFrame(pool).to_csv(os.path.join(DATA, 'quintile_pooled.csv'), index=False)

# ---------- 5. transition inheritance ----------
ct = st.groupby(['iso3', 'kind']).agg(n=('exp', 'size'), exp=('exp', 'sum'), eff=('eff', 'sum')).reset_index()
ct = ct.pivot(index='iso3', columns='kind', values=['n', 'exp', 'eff']).fillna(0)
ct.columns = [f'{a}_{b}' for a, b in ct.columns]
ct = ct.join(meta)
gex = gas[gas['exp']].groupby('iso3')['prot_rp'].median().rename('prot_median_gas_exposed')
ct = ct.join(gex)
ct['ev_per_gas'] = ct['n_ev'] / ct['n_gas'].replace(0, np.nan)
ct['gas_pct'] = 100 * ct['exp_gas'] / ct['n_gas'].replace(0, np.nan)
ct['gas_eff_pct'] = 100 * ct['eff_gas'] / ct['n_gas'].replace(0, np.nan)
ct['ev_pct'] = 100 * ct['exp_ev'] / ct['n_ev'].replace(0, np.nan)
ct = ct[ct['n_gas'] >= 1000]
# Natural Earth assigns ISO3 '-99' to several disputed/miscellaneous
# territories; pooled they would form one mixed pseudo-country row in the
# country table (labelled by whichever territory sorts first, i.e.
# 'Kosovo'). The published COUNTRY table therefore excludes the bucket; its
# stations remain in every aggregate (global, income-group, maturity-bin and
# Spearman) statistic, including those for territories without an ISO3 code.
ct[ct.index != '-99'].reset_index().to_csv(
    os.path.join(DATA, 'inheritance_country.csv'), index=False)

# maturity quintiles (weighted by fuel-station count)
v = ct.dropna(subset=['ev_per_gas']).sort_values('ev_per_gas').copy()
v['cw'] = v['n_gas'].cumsum() / v['n_gas'].sum()
v['mq'] = np.minimum((v['cw'] * 5).astype(int) + 1, 5)
ib = v.groupby('mq').apply(lambda d: pd.Series({
    'gas_pct_w': 100 * d['exp_gas'].sum() / d['n_gas'].sum(),
    'gas_eff_pct_w': 100 * d['eff_gas'].sum() / d['n_gas'].sum(),
    'eff_share_of_exposed': 100 * d['eff_gas'].sum() / max(d['exp_gas'].sum(), 1),
    'ev_per_gas_max': d['ev_per_gas'].max(), 'n_countries': len(d),
    'prot_median': d['prot_median_gas_exposed'].median(),
    'pop_m': d['pop_nat'].sum() / 1e6}), include_groups=False)
ib.reset_index().to_csv(os.path.join(DATA, 'inheritance_bins.csv'), index=False)

zone = v[v['ev_per_gas'] < 0.1]
zstats = {'n_countries': int(len(zone)), 'pop_b': round(float(zone['pop_nat'].sum()) / 1e9, 2),
          'gas_pct_w': round(100 * float(zone['exp_gas'].sum() / zone['n_gas'].sum()), 1),
          'gas_eff_pct_w': round(100 * float(zone['eff_gas'].sum() / zone['n_gas'].sum()), 1),
          'prot_median': round(float(zone['prot_median_gas_exposed'].median()), 1)}

from scipy.stats import spearmanr
rho, p = spearmanr(v['ev_per_gas'], v['gas_pct'])

# ---------- 6. headline summary ----------
inc_df = pd.DataFrame(rows).set_index('inc_grp')
pbdf = pd.DataFrame(pb).set_index('inc_grp')
chn = ev[ev['iso3'] == 'CHN']
nld = ev[ev['iso3'] == 'NLD']
stats = {
    'qa': qa,
    'income_ratio_gradient': {g: round(float(inc_df.loc[g, 'ratio']), 2) for g in INC_ORDER},
    'ev_pct_by_income': {g: round(float(inc_df.loc[g, 'ev_pct']), 1) for g in INC_ORDER},
    'share_exposed_HIC_UMC': round(float(inc_df.loc[['HIC', 'UMC'], 'share_ev_exp_stock'].sum()), 1),
    'share_pop_HIC_UMC': round(float(inc_df.loc[['HIC', 'UMC'], 'share_pop'].sum()), 1),
    'conc_indices': ci,
    'protection': {
        'global_nominal': int(ev['exp'].sum()), 'global_effective': int(ev['eff'].sum()),
        'global_pct_protected': round(100 * (1 - ev['eff'].sum() / ev['exp'].sum()), 1),
        'share_effective_by_income': {g: round(float(pbdf.loc[g, 'share_effective']), 1) for g in INC_ORDER},
        'share_nominal_by_income': {g: round(float(pbdf.loc[g, 'share_nominal']), 1) for g in INC_ORDER},
        'CHN': {'nominal': int(chn['exp'].sum()), 'effective': int(chn['eff'].sum())},
        'NLD': {'nominal': int(nld['exp'].sum()), 'effective': int(nld['eff'].sum())},
    },
    'inheritance': {'spearman_maturity_vs_gasexp': round(float(rho), 3), 'p': float(p),
                    'bins': ib.round(2).to_dict('index')},
    'pretransition_zone': zstats,
    'quintile_pooled': {f"{r['kind']}_q{r['q']}": round(float(r['pct']), 2) for r in pool},
}
json.dump(stats, open(os.path.join(DATA, 'fig05_stats.json'), 'w'), ensure_ascii=False, indent=1)
print('STATS WRITTEN')
print('conc:', ci)
print('inherit:', ib.round(2).to_dict('index'))
