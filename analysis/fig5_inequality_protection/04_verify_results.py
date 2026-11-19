#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute exposure, income shares, protection and concentration statistics from the local station covariates. EVFR_FIG5_BASELINE optionally selects a separate reference directory for comparing outputs."""
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P

DATA = P.get("FIG5_DATA")
BASE = os.environ.get("EVFR_FIG5_BASELINE", "")
# optional: point EVFR_FIG5_BASELINE at a reference snapshot of this module's
# data folder to report drift; when unset the comparison is skipped.

st = pd.read_parquet(os.path.join(DATA, 'fig05_station_covariates.parquet'))
st['exp'] = st['d100'] > 0
st['eff'] = st['exp'] & (st['prot_rp'] < 100)
ev, gas = st[st['kind'] == 'ev'], st[st['kind'] == 'gas']
S = json.load(open(os.path.join(DATA, 'fig05_stats.json')))
ci = json.load(open(os.path.join(DATA, 'concentration_indices.json')))

fails = []


def chk(name, got, want, tol=0.0):
    ok = (abs(got - want) <= tol) if tol else (got == want)
    print(('PASS' if ok else 'FAIL'), name, '| got', got, '| want', want)
    if not ok:
        fails.append(name)


# ---------- A: must match manuscript / physical exposure pipeline ----------
chk('ev_n', len(ev), 947602)
chk('gas_n', len(gas), 1344301)
chk('ev_exp', int(ev['exp'].sum()), 143264)
chk('gas_exp', int(gas['exp'].sum()), 167605)
chk('ratio', round(float(ev['exp'].mean() / gas['exp'].mean()), 3), 1.213)

for g, want in [('HIC', 1.39), ('UMC', 1.48), ('LMC', 0.69), ('LIC', 0.25)]:
    e, f = ev[ev['inc_grp'] == g], gas[gas['inc_grp'] == g]
    chk(f'ratio_{g}', round(float(e['exp'].mean() / f['exp'].mean()), 2), want)

hu_exp = ev[ev['inc_grp'].isin(['HIC', 'UMC'])]['exp'].sum() / ev['exp'].sum() * 100
chk('share_exposed_HIC_UMC(97.3)', round(float(hu_exp), 1), 97.3)

chk('global_effective(79908)', int(ev['eff'].sum()), 79908)
chk('pct_protected(44.2)', round(100 * (1 - ev['eff'].sum() / ev['exp'].sum()), 1), 44.2)
chk('CHN_nominal', int(ev[(ev['iso3'] == 'CHN')]['exp'].sum()), 73506)
chk('CHN_effective', int(ev[(ev['iso3'] == 'CHN')]['eff'].sum()), 65876)
chk('NLD_effective(0)', int(ev[(ev['iso3'] == 'NLD')]['eff'].sum()), 0)
effsum = ev['eff'].sum()
for g, want in [('HIC', 11.2), ('UMC', 88.4), ('LMC', 0.4)]:
    chk(f'eff_share_{g}', round(100 * float(ev[ev['inc_grp'] == g]['eff'].sum() / effsum), 1), want)

# Aggregate inputs consumed by the Figure 5 composite.
for name in (
    'income_group_exposure.csv', 'concentration_curve.csv', 'concentration_indices.json',
    'protection_by_country.csv', 'protection_bins_income.csv', 'quintile_within_country.csv',
    'quintile_pooled.csv', 'inheritance_country.csv', 'fig05_stats.json',
):
    path = os.path.join(DATA, name)
    chk(f'plot input {name}', os.path.isfile(path) and os.path.getsize(path) > 0, True)

# ---------- B: internal consistency for GDP/pop-dependent quantities ----------
cn = st.groupby(['iso3', 'kind']).agg(n=('exp', 'size'), exp=('exp', 'sum'), eff=('eff', 'sum')).reset_index()
cn = cn.pivot(index='iso3', columns='kind', values=['n', 'exp', 'eff']).fillna(0)
cn.columns = [f'{a}_{b}' for a, b in cn.columns]
meta = st.drop_duplicates('iso3').set_index('iso3')[['gdppc_nat', 'pop_nat']]
cn = cn.join(meta).dropna(subset=['gdppc_nat', 'pop_nat']).sort_values('gdppc_nat')


def conc(yw):
    x = np.concatenate([[0], np.cumsum(cn['pop_nat'].values) / cn['pop_nat'].sum()])
    y = np.concatenate([[0], np.cumsum(yw) / np.sum(yw)])
    return round(float(1 - 2 * np.trapezoid(y, x)), 3)


for name, col in [('ev_stock', 'n_ev'), ('ev_exposed', 'exp_ev'), ('ev_effective', 'eff_ev'),
                  ('gas_stock', 'n_gas'), ('gas_exposed', 'exp_gas'), ('gas_effective', 'eff_gas')]:
    chk(f'C_{name}_vs_this_run', conc(cn[col].values), ci[name], 0.001)

pop = st.drop_duplicates('iso3')
hu_pop = pop[pop['inc_grp'].isin(['HIC', 'UMC'])]['pop_nat'].sum() / pop['pop_nat'].sum() * 100
chk('share_pop_HIC_UMC_vs_stats', round(float(hu_pop), 1),
    round(float(S['share_pop_HIC_UMC']), 1), 0.05)

qp = pd.read_csv(os.path.join(DATA, 'quintile_pooled.csv')).set_index(['kind', 'q'])
for kind, q in [('ev', 1), ('ev', 5), ('gas', 1), ('gas', 5)]:
    parts = []
    k = st[st['kind'] == kind]
    for c, sub in k.groupby('iso3'):
        if len(sub) < 2000 or sub['gdppc_adm2'].nunique() < 8:
            continue
        r = sub['gdppc_adm2'].rank(method='first')
        qq = pd.qcut(r, 5, labels=[1, 2, 3, 4, 5]).astype(int)
        parts.append(pd.DataFrame({'q': qq, 'exp': sub['exp']}))
    allq = pd.concat(parts)
    got = round(100 * float(allq[allq['q'] == q]['exp'].mean()), 2)
    chk(f'quintile_{kind}_q{q}_vs_csv', got, round(float(qp.loc[(kind, q), 'pct']), 2), 0.01)

ctv = st.groupby(['iso3', 'kind']).agg(n=('exp', 'size'), exp=('exp', 'sum'), eff=('eff', 'sum')).reset_index()
ctv = ctv.pivot(index='iso3', columns='kind', values=['n', 'exp', 'eff']).fillna(0)
ctv.columns = [f'{a}_{b}' for a, b in ctv.columns]
ctv = ctv.join(meta)
ctv = ctv[ctv['n_gas'] >= 1000]
ctv['epg'] = ctv['n_ev'] / ctv['n_gas']
zone = ctv[ctv['epg'] < 0.1]
chk('zone_n_countries_vs_stats', int(len(zone)), int(S['pretransition_zone']['n_countries']))
chk('zone_pop_b_vs_stats', round(float(zone['pop_nat'].sum()) / 1e9, 2),
    float(S['pretransition_zone']['pop_b']), 0.01)

# ---------- optional comparison vs a reference snapshot ----------
cmp = {'baseline': BASE, 'deltas': {}}
if BASE and os.path.exists(os.path.join(BASE, 'fig05_stats.json')):
    Bold = json.load(open(os.path.join(BASE, 'fig05_stats.json'), encoding='utf-8'))
    cold = json.load(open(os.path.join(BASE, 'concentration_indices.json'), encoding='utf-8'))
    for k in ['ev_stock', 'ev_exposed', 'ev_effective', 'gas_stock', 'gas_exposed', 'gas_effective']:
        cmp['deltas'][f'C_{k}'] = {
            'baseline': cold[k], 'new': ci[k],
            'diff': round(ci[k] - cold[k], 3),
        }
    cmp['deltas']['share_pop_HIC_UMC'] = {
        'baseline': Bold['share_pop_HIC_UMC'], 'new': S['share_pop_HIC_UMC'],
        'diff': round(S['share_pop_HIC_UMC'] - Bold['share_pop_HIC_UMC'], 3),
    }
    cmp['deltas']['share_exposed_HIC_UMC'] = {
        'baseline': Bold['share_exposed_HIC_UMC'], 'new': S['share_exposed_HIC_UMC'],
        'diff': round(S['share_exposed_HIC_UMC'] - Bold['share_exposed_HIC_UMC'], 3),
    }
    cmp['deltas']['pretransition_zone'] = {
        'baseline': Bold['pretransition_zone'], 'new': S['pretransition_zone'],
    }
    cmp['deltas']['quintile_pooled'] = {
        'baseline': Bold.get('quintile_pooled'), 'new': S.get('quintile_pooled'),
    }
    cmp['deltas']['protection_unchanged'] = {
        'baseline': Bold['protection']['global_effective'],
        'new': S['protection']['global_effective'],
        'same': Bold['protection']['global_effective'] == S['protection']['global_effective'],
    }
    json.dump(cmp, open(os.path.join(DATA, 'compare_vs_fig5.json'), 'w'),
              ensure_ascii=False, indent=1)
    print('WROTE compare_vs_fig5.json')
    print('C deltas:', cmp['deltas'])
else:
    print('SKIP baseline compare (no', BASE, ')')

print()
print('ALL PASS' if not fails else f'FAILURES: {fails}')
sys.exit(1 if fails else 0)
