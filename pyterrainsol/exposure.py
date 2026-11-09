"""Upwind terrain: profiles, horizon angles and the Sx exposure index.

Sx (Winstral) asks a single question of a point: looking upwind, what is the
largest angle above horizontal subtended by the terrain?  Positive means
sheltered, negative means the ground falls away and the point is exposed.  It
is the same geometry as a solar horizon angle, pointed at the wind instead of
at the sun, so :func:`horizon_angle` here and ``pylidarterrain.horizon_angles``
compute the same quantity for different purposes.

Two things the index does not do, worth knowing before trusting it:

* **It saturates close in.**  The angle is dominated by whatever is nearest,
  because the vertical offset is divided by distance.  At 845 Avonside the
  terrain within 700 m subtends 12.4 deg and a 2200 m range 41 km away adds
  1.6 deg, so a long search distance costs time and changes almost nothing.
  Run a short fetch on fine terrain for shelter; use long fetch only to ask
  what the air has crossed, not whether the point is sheltered.
* **It is a straight line.**  That is the right geometry for flow forced over
  a barrier, and only an approximation for flow steered along a valley, where
  the whole point is that the air bends.

Distances and elevations must be in the same units, so use a projected CRS.
"""

import numpy as np

from scipy.ndimage import map_coordinates


def _sample(dem, transform, x, y, order=1):
    """Bilinear sample of ``dem`` at projected coordinates.

    The half-cell shift is not cosmetic: an affine transform maps a coordinate
    to the cell's *corner* index, while array interpolation treats index 0 as
    the *centre* of cell 0.  Omit it and every profile is displaced half a cell
    along both axes.
    """
    inv = ~transform
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    cols = inv.a * x + inv.b * y + inv.c - 0.5
    rows = inv.d * x + inv.e * y + inv.f - 0.5
    return map_coordinates(dem, [rows, cols], order=order, mode="nearest")


def upwind_profile(dem, transform, x, y, bearing, max_distance=5000.0, step=None):
    """Terrain profile along a ray from ``(x, y)`` toward ``bearing``.

    ``bearing`` is the compass direction the ray points, degrees clockwise from
    north -- for an upwind profile pass the direction the wind comes FROM, which
    is the meteorological convention already used in :mod:`pyterrainsol.wind`.

    Returns ``(distance, elevation)`` with distance starting at zero at the
    point itself.
    """
    if step is None:
        step = abs(transform.a)
    d = np.arange(0.0, max_distance + step, step)
    a = np.radians(bearing)
    return d, _sample(dem, transform, x + d * np.sin(a), y + d * np.cos(a))


def horizon_angle(dem, transform, x, y, bearing, max_distance=5000.0, step=None,
                  z_offset=0.0, return_distance=False):
    """Largest elevation angle of the terrain along one upwind ray, in degrees.

    Positive is terrain above the point (shelter), negative is terrain falling
    away (exposure).  ``z_offset`` raises the observer -- use it to ask about
    the air a few metres above the surface rather than the ground itself.

    With ``return_distance``, also returns how far away the controlling terrain
    is, which is the honest way to report whether a long search mattered.
    """
    d, z = upwind_profile(dem, transform, x, y, bearing, max_distance, step)
    z0 = z[0] + z_offset
    with np.errstate(invalid="ignore", divide="ignore"):
        ang = np.degrees(np.arctan2(z[1:] - z0, d[1:]))
    i = int(np.nanargmax(ang))
    return (float(ang[i]), float(d[1:][i])) if return_distance else float(ang[i])


def sx(dem, transform, x, y, wind_from, half_window=15.0, n_rays=5,
       max_distance=5000.0, step=None, z_offset=0.0):
    """Winstral Sx: horizon angle averaged over an azimuth window.

    A single ray is noisy because one spur can control it, so Sx is averaged
    over ``+/- half_window`` degrees about ``wind_from``.  ``n_rays`` rays are
    spread evenly across that window, including both edges.
    """
    offsets = np.linspace(-half_window, half_window, int(n_rays))
    angles = [horizon_angle(dem, transform, x, y, wind_from + o,
                            max_distance=max_distance, step=step, z_offset=z_offset)
              for o in offsets]
    return float(np.mean(angles))


def sx_rose(dem, transform, x, y, bearings, weights=None, **kwargs):
    """Sx for a set of upwind bearings, and optionally their weighted mean.

    Pass ``weights`` -- typically the fraction of strong-wind hours arriving
    from each bearing -- to collapse the rose to one number per point.  Weights
    are normalised; without them the plain mean is returned alongside.
    """
    bearings = np.asarray(bearings, float)
    vals = np.array([sx(dem, transform, x, y, b, **kwargs) for b in bearings])
    if weights is None:
        return vals, float(vals.mean())
    w = np.asarray(weights, float)
    if w.shape != bearings.shape:
        raise ValueError("weights must match bearings")
    return vals, float(np.sum(vals * w) / np.sum(w))


def sx_raster(dem, res, bearings, weights=None, max_distance=500.0, step=None,
              sub_rays=(-7.5, 0.0, 7.5)):
    """Sx for every cell of a DEM, weighted across a wind rose.

    The point functions above are fine for a handful of sites; a 3 m farm grid
    needs the whole field at once, so this walks each ray outward one step at a
    time and keeps a running maximum over the entire array.  Cost is linear in
    ``len(bearings) * len(sub_rays) * max_distance / step``.

    ``dem`` must be north-up with rows increasing southward -- the usual raster
    convention -- and must contain no NaNs, since one propagates into every
    cell downwind of it.  ``weights`` is the share of strong wind arriving from
    each bearing; omit it for a plain mean.  Each bearing is split across
    ``sub_rays`` offsets, which stops a single spur controlling a whole sector.

    Returns degrees, same shape as ``dem``: positive is sheltered, negative is
    ground that falls away upwind.
    """
    dem = np.asarray(dem, float)
    if not np.all(np.isfinite(dem)):
        raise ValueError("dem contains non-finite values; fill them first")
    bearings = np.atleast_1d(np.asarray(bearings, float))
    w = np.ones_like(bearings) if weights is None else np.asarray(weights, float)
    if w.shape != bearings.shape:
        raise ValueError("weights must match bearings")
    if step is None:
        step = 2.0 * res
    # Start at the ADJACENT cell, not at `step`.  Sx is dominated by whatever is
    # nearest, so skipping the first cell to save one sample changed the 3 m
    # Avonside field by up to 21 deg on gully walls -- a far bigger error than
    # coarsening `step` anywhere further out.
    dist = np.arange(res, max_distance + 0.5 * step, step)
    rows, cols = np.indices(dem.shape, dtype=np.float32)

    out = np.zeros(dem.shape, np.float32)
    for bearing, weight in zip(bearings, w):
        for sub in sub_rays:
            a = np.radians(bearing + sub)
            dr = -np.cos(a) / res          # rows increase southward
            dc = np.sin(a) / res
            run = np.full(dem.shape, -90.0, np.float32)
            for d in dist:
                up = map_coordinates(dem, [rows + dr * d, cols + dc * d],
                                     order=1, mode="nearest")
                np.maximum(run, np.degrees(np.arctan2(up - dem, d)), out=run)
            out += (weight / len(sub_rays)) * run
    return out / w.sum()


# --- a physical terrain-to-wind-speed transfer -----------------------------
#
# Everything above measures SHELTER: Sx is the angle to the highest upwind
# obstruction, so its zero means "nothing upwind" and it has no mechanism for
# speed-up above the free stream.  Converting it to a wind multiplier requires
# an arbitrary re-centring, and the resulting values above 1.0 assert an
# acceleration the index never measured.
#
# The functions below take the other route, the one with a physical anchor:
# local slope in the wind direction and local curvature, after Liston & Elder
# (2006, J. Hydrometeorology 7, 217-234) -- the MicroMet terrain weighting.
# Convex ground and windward faces accelerate the flow, concavities and lee
# faces decelerate it, and the scheme is two-sided, so it nearly conserves the
# mean over any area containing both.


def slope_aspect(dem, res):
    """Slope (radians) and downslope aspect (degrees from north) of `dem`.

    Aspect follows the compass convention used throughout the package: the
    direction the ground FACES, so a west-facing slope reads about 270.
    """
    z = np.asarray(dem, float)
    gy, gx = np.gradient(z, res, res)
    dzdx, dzdy = gx, -gy                       # to east-positive, north-positive
    slope = np.arctan(np.hypot(dzdx, dzdy))
    # The gradient points uphill, so downhill is (-dzdx, -dzdy).  Azimuth is
    # clockwise from north, i.e. atan2(east component, north component).
    # Swapping these two arguments rotates aspect by exactly 90 degrees, which
    # silently ZEROES any slope-in-wind-direction term computed from it: the
    # result stays symmetric about a crest and so looks plausible rather than
    # broken.  Pinned by test_slope_aspect_compass_convention.
    aspect = (np.degrees(np.arctan2(-dzdx, -dzdy)) + 360.0) % 360.0
    return slope, aspect


def curvature(dem, res, length_scale=150.0):
    """MicroMet topographic curvature at `length_scale` metres.

    The mean, over four direction pairs (W-E, S-N and both diagonals), of the
    cell's elevation less the mean of its two neighbours at that separation,
    normalised by twice the separation.  Convex ground (ridges, crests, spur
    noses) is POSITIVE; concave ground (gullies, hollows) is negative.

    `length_scale` should be roughly half the wavelength of the terrain
    features that matter.  At 845 Avonside the result is insensitive to it:
    the west-minus-east aspect contrast holds between 0.14 and 0.19 for any
    value from 30 to 600 m, so the choice sets detail, not the answer.

    Neighbours are fetched with edge replication, NOT ``np.roll``: at a 150 m
    offset on a 3 m grid the stencil reaches 50 cells, and a wrapped neighbour
    imports terrain from the opposite side of the window.  At 845 Avonside the
    farm mask touches the bottom edge of the DEM window with no margin at all,
    where wrapping corrupted 13% of farm cells by up to 0.294 against a true
    field range of 0.122 -- a solid block of extreme values in one corner.
    """
    z = np.asarray(dem, float)
    k = max(1, int(round(length_scale / res)))
    zp_ = np.pad(z, k, mode="edge")
    ny, nx = z.shape
    out = np.zeros_like(z)
    for dr, dc in ((0, k), (k, 0), (k, k), (k, -k)):
        a = zp_[k + dr:k + dr + ny, k + dc:k + dc + nx]
        b = zp_[k - dr:k - dr + ny, k - dc:k - dc + nx]
        out += (z - 0.5 * (a + b)) / (2.0 * np.hypot(dr, dc) * res)
    return out / 4.0


def micromet_weight(dem, res, bearings, weights=None, length_scale=150.0,
                    gamma_slope=0.58, gamma_curvature=0.42, scale_mask=None):
    """Wind speed multiplier from terrain, after Liston & Elder (2006).

    ``W = 1 + gamma_slope * Omega_s + gamma_curvature * Omega_c``, where
    Omega_s is the slope in the wind direction and Omega_c the curvature, each
    scaled to +/- 0.5 by dividing by twice its maximum absolute value.  The
    published weights 0.58 and 0.42 sum to one, so W is bounded by 0.5 and 1.5.

    `bearings` are the directions the wind comes FROM, in degrees; `weights`
    their relative frequency (equal if omitted), which is how a wind rose is
    folded in.  `scale_mask` restricts the normalisation to a region of
    interest -- important when `dem` is a padded window, since otherwise a few
    extreme cells outside the area you care about set the scaling for all of it.

    Returns an array shaped like `dem`.  Note that the scaling is
    DOMAIN-DEPENDENT by construction: the same terrain normalised over a
    different extent gives different numbers, so W is a relative redistribution
    of a reference wind, not an absolute speed.
    """
    z = np.asarray(dem, float)
    b = np.atleast_1d(np.asarray(bearings, float))
    w = np.ones_like(b) if weights is None else np.atleast_1d(np.asarray(weights, float))
    if w.shape != b.shape:
        raise ValueError("bearings and weights must have the same shape")

    slope, aspect = slope_aspect(z, res)
    curv = curvature(z, res, length_scale)
    m = np.ones(z.shape, bool) if scale_mask is None else np.asarray(scale_mask, bool)
    smax = np.abs(slope[m]).max()
    cmax = np.abs(curv[m]).max()
    if smax == 0 or cmax == 0:
        raise ValueError("flat DEM: terrain weighting is undefined")
    omega_c = curv / (2.0 * cmax)

    out = np.zeros_like(z)
    for bearing, weight in zip(b, w):
        # positive where the ground rises INTO the wind
        omega_s = slope * np.cos(np.deg2rad(bearing) - np.deg2rad(aspect)) / (2.0 * smax)
        out += weight * (1.0 + gamma_slope * omega_s + gamma_curvature * omega_c)
    return out / w.sum()
