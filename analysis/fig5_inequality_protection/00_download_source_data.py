#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Retrieve the World Bank income classification, Kummu GDP, FLOPROS and OWID population files into the local input tree. Existing files are reused. Supply the national WDI GDP and population exports separately. Input locations and source licences are listed in DATA_SOURCES.md."""
import os
import sys
import zipfile

import requests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "shared"))
import paths as P  # noqa: E402

ZEN = "https://zenodo.org/api/records/18429133/files"
JOBS = [
    (P.get("WB_INCOME"),
     "https://api.worldbank.org/v2/country?format=json&per_page=400"),
    (P.get("KUMMU_ADM2_TIF"),
     f"{ZEN}/rast_adm2_gdp_perCapita_1990_2024.tif/content"),
    (P.get("KUMMU_ADM0"),
     f"{ZEN}/tabulated_adm0_gdp_perCapita.csv/content"),
    (os.path.join(os.path.dirname(P.get("FLOPROS_TIF")),
                  "flopros_supplement.zip"),
     "https://nhess.copernicus.org/articles/16/1049/2016/"
     "nhess-16-1049-2016-supplement.zip"),
    (P.get("OWID_POPULATION"),
     "https://ourworldindata.org/grapher/population.csv"),
]

if __name__ == "__main__":
    for out, url in JOBS:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            print("skip (exists):", os.path.basename(out))
            continue
        print("fetch:", os.path.basename(out), flush=True)
        with requests.get(url, stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(out, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
    # unpack the FLOPROS shapefile if it is not present yet
    shp = os.path.join(P.get("FLOPROS_SHP_DIR"), "FLOPROS_shp_V1.shp")
    if not os.path.exists(shp):
        zpath = os.path.join(os.path.dirname(P.get("FLOPROS_TIF")),
                             "flopros_supplement.zip")
        with zipfile.ZipFile(zpath) as z:
            for nm in z.namelist():
                if "FLOPROS_shp_V1/" in nm and not nm.endswith("/"):
                    tgt = os.path.join(P.get("FLOPROS_SHP_DIR"),
                                       os.path.basename(nm))
                    os.makedirs(os.path.dirname(tgt), exist_ok=True)
                    with open(tgt, "wb") as f:
                        f.write(z.read(nm))
        print("FLOPROS shapefile unpacked")
    print("DONE")
