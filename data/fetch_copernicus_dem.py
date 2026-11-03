"""Rebuild snowy_monaro_dem30_mga55.tif from the Copernicus GLO-30 DEM.

Nine 1-degree tiles are mosaicked, reprojected to MGA zone 55 at 30 m, and
cropped to the largest window that contains no reprojection nodata.  The result
is the regional terrain used for upwind fetch and exposure work around
845 Avonside, where the airborne LiDAR covers only 4 x 2 km.

Copernicus DEM is free to use with attribution to ESA / Airbus Defence and
Space.  The committed GeoTIFF is a convenience extract, not an endorsed
product; regenerate it here if you need a different window.

    python fetch_copernicus_dem.py [--out PATH] [--keep-tiles]
"""
from __future__ import annotations
import argparse, pathlib, urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.transform import Affine, array_bounds
from rasterio.warp import calculate_default_transform, reproject, Resampling

HOST = "https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com"
TILES = [(la, lo) for la in ("S36", "S37", "S38") for lo in ("E147", "E148", "E149")]
WINDOW = (147.80, -36.98, 149.30, -35.82)      # west, south, east, north
RES = 30.0
CRS = "EPSG:28355"


def download(dest: pathlib.Path) -> list[pathlib.Path]:
    dest.mkdir(parents=True, exist_ok=True)

    def one(job):
        la, lo = job
        name = f"Copernicus_DSM_COG_10_{la}_00_{lo}_00_DEM"
        path = dest / f"{name}.tif"
        if not (path.exists() and path.stat().st_size > 10_000_000):
            urllib.request.urlretrieve(f"{HOST}/{name}/{name}.tif", path)
        return path

    with ThreadPoolExecutor(4) as pool:
        return list(pool.map(one, TILES))


def mosaic_reproject(paths):
    srcs = [rasterio.open(p) for p in paths]
    west, south, east, north = WINDOW
    arr, src_tr = merge(srcs, bounds=(west, south, east, north))
    src_crs, src_nodata = srcs[0].crs, srcs[0].nodata
    bounds = array_bounds(arr.shape[1], arr.shape[2], src_tr)
    dst_tr, w, h = calculate_default_transform(
        src_crs, CRS, arr.shape[2], arr.shape[1], *bounds, resolution=RES)
    dem = np.full((h, w), np.nan, np.float32)
    reproject(arr[0], dem, src_transform=src_tr, src_crs=src_crs,
              dst_transform=dst_tr, dst_crs=CRS, resampling=Resampling.bilinear,
              src_nodata=src_nodata, dst_nodata=np.nan)
    for s in srcs:
        s.close()
    return dem, dst_tr


def trim(dem, transform, tol=0.9999, step=8):
    """Shrink to the largest centred window essentially free of nodata.

    Reprojecting a lon/lat mosaic into a projected CRS rotates it, leaving
    empty corners.  Trim whichever edge currently carries the most nodata
    until the window is clean.
    """
    valid = np.isfinite(dem)
    r0, r1, c0, c1 = 0, dem.shape[0], 0, dem.shape[1]
    while valid[r0:r1, c0:c1].mean() < tol:
        edges = {"t": valid[r0, c0:c1].mean(), "b": valid[r1 - 1, c0:c1].mean(),
                 "l": valid[r0:r1, c0].mean(), "r": valid[r0:r1, c1 - 1].mean()}
        worst = min(edges, key=edges.get)
        r0, r1, c0, c1 = (r0 + step, r1, c0, c1) if worst == "t" else \
                         (r0, r1 - step, c0, c1) if worst == "b" else \
                         (r0, r1, c0 + step, c1) if worst == "l" else \
                         (r0, r1, c0, c1 - step)
    out_tr = Affine(transform.a, transform.b, transform.c + c0 * transform.a,
                    transform.d, transform.e, transform.f + r0 * transform.e)
    return dem[r0:r1, c0:c1], out_tr


def write(dem, transform, path):
    """Int16 decimetres: 0.1 m precision is far below the product's own error."""
    bad = ~np.isfinite(dem)
    dm = np.where(bad, -32768, np.rint(np.nan_to_num(dem) * 10.0)).astype(np.int16)
    with rasterio.open(path, "w", driver="GTiff", height=dm.shape[0], width=dm.shape[1],
                       count=1, dtype="int16", crs=CRS, transform=transform,
                       nodata=-32768, compress="deflate", predictor=2, zlevel=9,
                       tiled=True, blockxsize=512, blockysize=512) as ds:
        ds.write(dm, 1)
        ds.scales = (0.1,)
        ds.update_tags(1, UNITS="m AHD, stored in decimetres (apply scale 0.1)")
        ds.update_tags(SOURCE="Copernicus GLO-30 DEM, nine 1-degree tiles",
                       LICENCE="Copernicus DEM — free to use with attribution to ESA / Airbus DS")
    return dm.shape, int(bad.sum())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="snowy_monaro_dem30_mga55.tif")
    ap.add_argument("--tiles", default="_tiles")
    ap.add_argument("--keep-tiles", action="store_true")
    args = ap.parse_args()

    tiles = download(pathlib.Path(args.tiles))
    dem, tr = mosaic_reproject(tiles)
    dem, tr = trim(dem, tr)
    shape, nodata = write(dem, tr, args.out)
    print(f"{args.out}: {shape[0]} x {shape[1]} at {RES:.0f} m, {nodata} nodata cells")
    if not args.keep_tiles:
        for p in tiles:
            p.unlink()
        pathlib.Path(args.tiles).rmdir()
