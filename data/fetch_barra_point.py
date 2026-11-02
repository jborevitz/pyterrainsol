"""Pull a BARRA-C2 time series for one point from NCI's THREDDS server.

Reproduces data/barra_avonside_20y_20min.npz.  Nothing here is specific to
Avonside beyond the default coordinates.

    python fetch_barra_point.py --lat -36.435477 --lon 148.725680 \
        --start 2006 --end 2025 --out barra_point.npz

Two practical notes, both learned the hard way:

* netCDF4/HDF5 is NOT thread-safe over OPeNDAP -- concurrent reads collide on
  the DAP client's temp files.  Parallelise with separate *processes*, and note
  that `concurrent.futures.ProcessPoolExecutor` needs POSIX semaphores that some
  sandboxes deny; running N copies of this script with --shard is the portable
  way.
* BARRA's stored time coordinate carries +/- 1 s rounding, so raw timestamps look
  like 02:39:59.  Snap to the sampling interval before asserting regularity.
"""

import argparse, os, time
import numpy as np

TMPL = ("https://thredds.nci.org.au/thredds/dodsC/ob53/output/reanalysis/AUST-04/"
        "BOM/ERA5/historical/hres/BARRA-C2/v1/{freq}/{v}/latest/"
        "{v}_AUST-04_ERA5_historical_hres_BOM_BARRA-C2_v1_{freq}_{ym}-{ym}.nc")
STEP = {"20min": 1200, "1hr": 3600}


def fetch_month(var, ym, lat, lon, freq="20min", cache="parts", retries=3):
    """Nearest-gridpoint series for one variable-month.  Cached on disk."""
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, f"{var}_{ym}.npz")
    if os.path.exists(path):
        return path
    import xarray as xr
    for attempt in range(retries):
        try:
            ds = xr.open_dataset(TMPL.format(v=var, ym=ym, freq=freq), engine="netcdf4")
            da = ds.sel(lat=lat, lon=lon, method="nearest")[var].load()
            ds.close()
            np.savez(path, t=da.time.values.astype("datetime64[s]").astype("int64"),
                     x=da.values.astype("float32"),
                     glat=float(da.lat), glon=float(da.lon))
            return path
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(5 * (attempt + 1))


def consolidate(varnames, yms, out, freq="20min", cache="parts"):
    cols, tref = {}, None
    for v in varnames:
        ts = [np.load(os.path.join(cache, f"{v}_{ym}.npz")) for ym in yms]
        t = np.concatenate([z["t"] for z in ts])
        step = STEP[freq]
        t = ((t + step // 2) // step) * step        # undo the +/- 1 s rounding
        if tref is None:
            tref = t
            if not np.all(np.diff(t) == step):
                raise ValueError("series has gaps or duplicates after snapping")
        elif not np.array_equal(t, tref):
            raise ValueError(f"{v} time axis differs from {varnames[0]}")
        cols[v] = np.concatenate([z["x"] for z in ts])
    np.savez_compressed(out, time=tref, **cols)
    return tref, cols


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lat", type=float, default=-36.435477)
    p.add_argument("--lon", type=float, default=148.725680)
    p.add_argument("--start", type=int, default=2006)
    p.add_argument("--end", type=int, default=2025)
    p.add_argument("--freq", default="20min", choices=sorted(STEP))
    p.add_argument("--vars", nargs="+", default=["uas", "vas", "ua50m", "va50m"])
    p.add_argument("--out", default="barra_point.npz")
    p.add_argument("--shard", type=int, nargs=2, metavar=("I", "N"),
                   help="fetch only every Nth job (offset I); run N copies, then "
                        "re-run once without --shard to consolidate")
    a = p.parse_args()

    yms = [f"{y}{m:02d}" for y in range(a.start, a.end + 1) for m in range(1, 13)]
    jobs = [(v, ym) for v in a.vars for ym in yms]
    if a.shard:
        i, n = a.shard
        jobs = [j for k, j in enumerate(jobs) if k % n == i]
    for k, (v, ym) in enumerate(jobs, 1):
        fetch_month(v, ym, a.lat, a.lon, a.freq)
        if k % 40 == 0:
            print(f"{k}/{len(jobs)}", flush=True)
    if a.shard:
        print("shard done; re-run without --shard to consolidate")
        return
    t, cols = consolidate(a.vars, yms, a.out, a.freq)
    print(f"{a.out}: {len(t):,} steps, {len(cols)} variables")


if __name__ == "__main__":
    main()
