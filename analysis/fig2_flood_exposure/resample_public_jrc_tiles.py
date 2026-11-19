"""
Independent check of the flood-depth sampling.

The main pipeline samples global mosaics of the JRC river flood hazard maps.
This script repeats the sampling from the public 10-degree tiles of the same
dataset and compares the exposed shares of both networks at one return period.

From the repository root:
  python analysis/fig2_flood_exposure/resample_public_jrc_tiles.py <tile_directory> <return_period>

<tile_directory> holds the tiles ID*_RP<return_period>_depth.tif downloaded from
https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard/RP<return_period>/
The station table is read from the Figure 2 data folder (station_depths.parquet).
"""
import glob, os, sys
import numpy as np, pandas as pd, rasterio
from rasterio.transform import rowcol
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared"))
import paths as P  # noqa: E402

tiles, rp = sys.argv[1], int(sys.argv[2])
a = pd.read_parquet(os.path.join(P.get("FIG2_DATA"), "station_depths.parquet"))
lon, lat = a.lon.values, a.lat.values
new = np.zeros(len(a), dtype=np.float32)
covered = np.zeros(len(a), dtype=bool)
for f in sorted(glob.glob(os.path.join(tiles, "*.tif"))):
    with rasterio.open(f) as ds:
        b = ds.bounds
        m = (lon >= b.left) & (lon < b.right) & (lat > b.bottom) & (lat <= b.top)
        if not m.any():
            continue
        r, c = rowcol(ds.transform, lon[m], lat[m])
        r = np.clip(np.asarray(r), 0, ds.height - 1)
        c = np.clip(np.asarray(c), 0, ds.width - 1)
        v = ds.read(1)[r, c].astype(np.float32)
        v[(v == ds.nodata) | ~np.isfinite(v) | (v >= 1e30) | (v < 0)] = 0.0
        new[m] = v
        covered[m] = True
old = a[f"d{rp:03d}"].values > 0
ev = (a.kind == "ev").values
print(f"return period {rp}: stations outside tile coverage {int((~covered).sum()):,}")
for name, m in (("charging", ev), ("fuel", ~ev)):
    print(f"  {name}: pipeline {100 * old[m].mean():.3f}%  re-sampled {100 * (new[m] > 0).mean():.3f}%  "
          f"flag agreement {100 * (old[m] == (new[m] > 0)).mean():.3f}%")
print(f"  ratio: pipeline {old[ev].mean() / old[~ev].mean():.4f}  re-sampled "
      f"{(new[ev] > 0).mean() / (new[~ev] > 0).mean():.4f}")
