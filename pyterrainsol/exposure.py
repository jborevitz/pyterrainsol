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
