#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Attach national income and population, subnational GDP and flood-protection standards to each station. Use the 2025 World Bank WDI estimates, with 2024 as fallback. Kummu supplies missing national GDP and admin-2 GDP. OWID supplies the population of Taiwan. FLOPROS supplies the local river-protection standard."""
import os, sys, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

T0 = time.time()
DATA = P.ensure_dir(P.get("FIG5_DATA"))
CKPT = P.ensure_dir(os.path.join(DATA, "checkpoints"))

CK_BASE = os.path.join(CKPT, "ck_base.parquet")
FLO_TIF = P.get("FLOPROS_TIF")
FLO_SHP = os.path.join(P.get("FLOPROS_SHP_DIR"), "FLOPROS_shp_V1.shp")
OUT = os.path.join(DATA, "fig05_station_covariates.parquet")
LOG = os.path.join(DATA, "build_log.json")

# Local World Bank WDI exports
_NEW = P.SOURCE_ROOT
WB_GDP_CSV = os.environ.get(
    "EVFR_WB_GDP_CSV",
    os.path.join(_NEW, "API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_33034",
                 "API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_33034.csv"),
)
WB_POP_CSV = os.environ.get(
    "EVFR_WB_POP_CSV",
    os.path.join(_NEW, "API_SP.POP.TOTL_DS2_en_csv_v2_33112",
                 "API_SP.POP.TOTL_DS2_en_csv_v2_33112.csv"),
)
WB_GDP_META = os.path.join(os.path.dirname(WB_GDP_CSV),
                           "Metadata_Country_" + os.path.basename(WB_GDP_CSV))

INCOME_PATCH = {'TWN': 'HIC'}
NAME_PATCH = {'Kosovo': 'UMC', 'Somaliland': 'LIC', 'Jersey': 'HIC', 'Guernsey': 'HIC',
              'Åland': 'HIC', 'St. Pierre and Miquelon': 'HIC'}
RES = 1 / 24.0


def _read_wdi_wide(csv_path):
    """World Bank DataBank CSV: 4 preamble rows, then Country Code + year columns."""
    df = pd.read_csv(csv_path, skiprows=4)
    if "Unnamed: 70" in df.columns:
        df = df.drop(columns=["Unnamed: 70"])
    return df


def _country_codes_from_meta(meta_path):
    """Prefer rows with a non-empty Region (economies); drop aggregates."""
    if not os.path.exists(meta_path):
        return None
    meta = pd.read_csv(meta_path)
    m = meta[meta["Region"].notna() & (meta["Region"].astype(str).str.len() > 0)]
    return set(m["Country Code"].astype(str))


def load_wb_gdp_series():
    """Return dict iso3 -> value and per-iso3 source label (wb_2025 / wb_2024)."""
    df = _read_wdi_wide(WB_GDP_CSV)
    codes = _country_codes_from_meta(WB_GDP_META)
    if codes is not None:
        df = df[df["Country Code"].isin(codes)].copy()
    latest = df.get("2025", pd.Series(np.nan, index=df.index))
    val = latest.combine_first(df["2024"])
    src = np.where(latest.notna(), "wb_2025",
                   np.where(df["2024"].notna(), "wb_2024", "missing"))
    out_v = dict(zip(df["Country Code"].astype(str), val.astype(float)))
    out_s = dict(zip(df["Country Code"].astype(str), src))
    # drop NaN values from value dict
    out_v = {k: v for k, v in out_v.items() if pd.notna(v)}
    return out_v, out_s


def load_wb_pop_series():
    df = _read_wdi_wide(WB_POP_CSV)
    codes = _country_codes_from_meta(WB_GDP_META)  # same country list
    if codes is not None:
        df = df[df["Country Code"].isin(codes)].copy()
    latest = df.get("2025", pd.Series(np.nan, index=df.index))
    val = latest.combine_first(df["2024"])
    src = np.where(latest.notna(), "wb_2025",
                   np.where(df["2024"].notna(), "wb_2024", "missing"))
    out_v = {k: v for k, v in zip(df["Country Code"].astype(str), val.astype(float))
             if pd.notna(v)}
    out_s = dict(zip(df["Country Code"].astype(str), src))
    return out_v, out_s


def stage_A():
    import pyarrow.parquet as pq
    st = pq.read_table(os.path.join(P.get('FIG2_DATA'), 'station_depths.parquet')).to_pandas()
    m = []
    for key in ['DB_CHARGING', 'DB_FUEL']:
        t = pq.read_table(P.get(key), columns=['country', 'iso3']).to_pandas()
        m.append(t.drop_duplicates())
    iso = pd.concat(m).drop_duplicates()
    iso = iso.groupby('country')['iso3'].agg(lambda s: s.mode().iat[0]).to_dict()
    st['iso3'] = st['country'].map(iso)

    wb = json.load(open(P.get('WB_INCOME')))[1]
    inc = {c['id']: c['incomeLevel']['id'] for c in wb if c['incomeLevel']['id'] != 'NA'}
    inc.update(INCOME_PATCH)
    st['inc_grp'] = st['iso3'].map(inc)
    fix = st['inc_grp'].isna() & st['country'].isin(NAME_PATCH)
    st.loc[fix, 'inc_grp'] = st.loc[fix, 'country'].map(NAME_PATCH)

    # --- national GDP: World Bank WDI (2025 then 2024), Kummu fallback ---
    gdp_wb, gdp_src = load_wb_gdp_series()
    st['gdppc_nat'] = st['iso3'].map(gdp_wb)
    st['gdppc_nat_source'] = st['iso3'].map(gdp_src)
    n_wb = int(st['gdppc_nat'].notna().sum())
    a0 = pd.read_csv(P.get('KUMMU_ADM0'))
    kummu = a0.set_index('iso3')['2024'].to_dict()
    need = st['gdppc_nat'].isna() & st['iso3'].notna()
    st.loc[need, 'gdppc_nat'] = st.loc[need, 'iso3'].map(kummu)
    filled_k = int(need.sum() - st.loc[need, 'gdppc_nat'].isna().sum())
    st.loc[need & st['gdppc_nat'].notna(), 'gdppc_nat_source'] = 'fallback_kummu_2024'
    st['gdppc_nat_source'] = st['gdppc_nat_source'].fillna('missing')

    # --- population: WB SP.POP.TOTL (2025 then 2024); OWID 2023 only for TWN ---
    pop_wb, pop_src = load_wb_pop_series()
    st['pop_nat'] = st['iso3'].map(pop_wb)
    st['pop_nat_source'] = st['iso3'].map(pop_src)
    owid = pd.read_csv(P.get('OWID_POPULATION'))
    owid = owid[(owid['Year'] == 2023) & owid['Code'].notna()
                & ~owid['Code'].fillna('').str.startswith('OWID')]
    owid_map = owid.set_index('Code')['Population'].to_dict()
    # Only Taiwan is patched from OWID; other WB-missing codes stay missing.
    twn = (st['iso3'] == 'TWN') & st['pop_nat'].isna()
    st.loc[twn, 'pop_nat'] = st.loc[twn, 'iso3'].map(owid_map)
    st.loc[twn & st['pop_nat'].notna(), 'pop_nat_source'] = 'fallback_owid_2023_TWN_only'
    st['pop_nat_source'] = st['pop_nat_source'].fillna('missing')
    n_twn = int(twn.sum())

    st['col5'] = np.clip(((st['lon'].values + 180) / (1 / 12.)).astype(np.int32), 0, 4319)
    st['row5'] = np.clip(((90 - st['lat'].values) / (1 / 12.)).astype(np.int32), 0, 2159)
    st['col2'] = np.clip(((st['lon'].values + 180) / RES).astype(np.int32), 0, 8639)
    st['row2'] = np.clip(((90 - st['lat'].values) / RES).astype(np.int32), 0, 4319)

    src_counts = st['gdppc_nat_source'].value_counts(dropna=False).to_dict()
    pop_counts = st['pop_nat_source'].value_counts(dropna=False).to_dict()
    log = {
        'variant': 'fig5',
        'wb_gdp_csv': WB_GDP_CSV,
        'wb_pop_csv': WB_POP_CSV,
        'wb_gdp_indicator': 'NY.GDP.PCAP.PP.CD',
        'wb_pop_indicator': 'SP.POP.TOTL',
        'n_total': int(len(st)),
        'iso3_unmatched': int(st['iso3'].isna().sum()),
        'income_unmatched_n': int(st['inc_grp'].isna().sum()),
        'income_unmatched_top': st.loc[st['inc_grp'].isna(), 'country'].value_counts().head(12).to_dict(),
        'gdppc_nat_unmatched_n': int(st['gdppc_nat'].isna().sum()),
        'pop_unmatched_n': int(st['pop_nat'].isna().sum()),
        'gdppc_nat_source_counts': {str(k): int(v) for k, v in src_counts.items()},
        'pop_nat_source_counts': {str(k): int(v) for k, v in pop_counts.items()},
        'gdppc_nat_matched_wb': n_wb,
        'gdppc_nat_fallback_kummu_n_stations': filled_k,
        'pop_nat_fallback_owid_TWN_n_stations': n_twn,
        'note': ('GDP: WB 2025 then 2024, else Kummu 2024 for missing ISO3. '
                 'Population: SP.POP.TOTL 2025 then 2024; OWID 2023 only for TWN.'),
    }
    json.dump(log, open(LOG, 'w'), ensure_ascii=False, indent=1)
    tmp = CK_BASE + '.tmp'
    st.to_parquet(tmp, index=False)
    os.replace(tmp, CK_BASE)
    print('stage A done', round(time.time() - T0, 1), 's |', log)


def stage_B():
    import geopandas as gpd, rasterio
    from rasterio import features
    from rasterio.transform import from_origin
    g = gpd.read_file(FLO_SHP, engine='pyogrio', columns=['MerL_Riv', 'geometry'])
    tr = from_origin(-180, 90, RES, RES)
    arr = features.rasterize(((geom, v) for geom, v in zip(g.geometry, g['MerL_Riv'])),
                             out_shape=(4320, 8640), transform=tr, fill=np.nan,
                             all_touched=True, dtype='float32')
    interior = features.rasterize(((geom, v) for geom, v in zip(g.geometry, g['MerL_Riv'])),
                                  out_shape=(4320, 8640), transform=tr, fill=np.nan,
                                  all_touched=False, dtype='float32')
    arr = np.where(np.isnan(interior), arr, interior)
    prof = dict(driver='GTiff', height=4320, width=8640, count=1, dtype='float32',
                crs='EPSG:4326', transform=tr, compress='deflate', nodata=np.nan)
    tmp = FLO_TIF + '.tmp.tif'
    with rasterio.open(tmp, 'w', **prof) as dst:
        dst.write(arr, 1)
    os.replace(tmp, FLO_TIF)
    print('stage B done', round(time.time() - T0, 1), 's | valid px:', int(np.isfinite(arr).sum()))


def stage_C():
    import rasterio
    st = pd.read_parquet(CK_BASE)
    with rasterio.open(P.get('KUMMU_ADM2_TIF')) as r:
        band = r.read(35)
    st['gdppc_adm2'] = band[st['row5'].values, st['col5'].values].astype('float64')
    del band
    with rasterio.open(FLO_TIF) as r:
        flo = r.read(1)
    st['prot_rp'] = flo[st['row2'].values, st['col2'].values].astype('float64')
    del flo
    n_g = int(st['gdppc_adm2'].isna().sum())
    st['gdppc_adm2'] = st['gdppc_adm2'].fillna(st['gdppc_nat'])
    n_p = int(st['prot_rp'].isna().sum())
    med = st.groupby('iso3')['prot_rp'].transform('median')
    st['prot_rp_filled'] = st['prot_rp'].isna()
    st['prot_rp'] = st['prot_rp'].fillna(med)
    log = json.load(open(LOG))
    log.update({'gdppc_adm2_filled_from_nat': n_g, 'prot_rp_filled_from_country_median': n_p,
                'prot_rp_still_nan': int(st['prot_rp'].isna().sum())})
    json.dump(log, open(LOG, 'w'), ensure_ascii=False, indent=1)
    st = st.drop(columns=['col5', 'row5', 'col2', 'row2'])
    tmp = OUT + '.tmp'
    st.to_parquet(tmp, index=False)
    os.replace(tmp, OUT)
    print('stage C done', round(time.time() - T0, 1), 's |', log)


if __name__ == '__main__':
    if not os.path.exists(CK_BASE):
        stage_A()
    elif not os.path.exists(FLO_TIF):
        stage_B()
    elif not os.path.exists(OUT):
        stage_C()
    else:
        print('DONE')
