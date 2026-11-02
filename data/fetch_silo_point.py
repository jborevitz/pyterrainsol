"""Pull SILO daily climate for one point from the DataDrill service.

Reproduces data/silo_avonside_20y_daily.npz.

    python fetch_silo_point.py --username YOUR_REGISTERED_EMAIL \
        --lat -36.435477 --lon 148.725680 --start 20060101 --finish 20251231 \
        --out silo_point.npz

SILO is free but requires a registered email as the `username` parameter; an
anonymous request returns a "missing essential parameters" page rather than an
error code, so check the first line of the response before parsing. Register at
https://www.longpaddock.qld.gov.au/silo/ and supply your own address -- never
reuse one found in example code.

Licence and citation: https://www.longpaddock.qld.gov.au/silo/
Data is supplied by the Queensland Department of Environment, Science and
Innovation. This script retrieves it; it redistributes nothing.
"""

import argparse
import io
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

BASE = "https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php"
ALL_CODES = "RXNJVDESCLFTAPWMHG"          # every available variable in one call
VARIABLES = ("daily_rain", "max_temp", "min_temp", "radiation", "vp", "vp_deficit",
             "evap_pan", "evap_syn", "evap_comb", "evap_morton_lake", "et_short_crop",
             "et_tall_crop", "et_morton_actual", "et_morton_potential", "et_morton_wet",
             "mslp", "rh_tmax", "rh_tmin")


def fetch(username, lat, lon, start, finish, codes=ALL_CODES, timeout=180):
    """Return SILO's CSV response as a DataFrame. Raises if the service complains."""
    q = dict(format="csv", username=username, password="apirequest",
             lat=lat, lon=lon, start=start, finish=finish, comment=codes)
    with urllib.request.urlopen(BASE + "?" + urllib.parse.urlencode(q), timeout=timeout) as r:
        txt = r.read().decode()
    if not txt.lstrip().lower().startswith("latitude"):
        raise RuntimeError("SILO did not return a data table: " + txt.strip()[:200])
    df = pd.read_csv(io.StringIO(txt))
    df.columns = [c.strip() for c in df.columns]
    return df


def to_npz(df, path, variables=VARIABLES):
    """Write the numeric columns as float32 plus an int32 day index."""
    date = pd.to_datetime(df["YYYY-MM-DD"])
    out = {v: pd.to_numeric(df[v], errors="coerce").to_numpy(np.float32)
           for v in variables if v in df.columns}
    out["date"] = date.values.astype("datetime64[D]").astype("int32")
    out["lat"] = np.float64(df["latitude"].iloc[0])
    out["lon"] = np.float64(df["longitude"].iloc[0])
    meta = str(df["metadata"].iloc[0]) if "metadata" in df.columns else ""
    if "elevation=" in meta:
        out["elevation_m"] = np.float64(meta.split("elevation=")[1].split("m")[0])
    np.savez_compressed(path, **out)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--username", required=True, help="your SILO registration email")
    p.add_argument("--lat", type=float, default=-36.435477)
    p.add_argument("--lon", type=float, default=148.725680)
    p.add_argument("--start", default="20060101")
    p.add_argument("--finish", default="20251231")
    p.add_argument("--out", default="silo_point.npz")
    a = p.parse_args()
    df = fetch(a.username, a.lat, a.lon, a.start, a.finish)
    out = to_npz(df, a.out)
    print(f"{a.out}: {len(df)} days, {len(out) - 3} variables")


if __name__ == "__main__":
    main()
