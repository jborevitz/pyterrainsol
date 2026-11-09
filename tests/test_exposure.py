import numpy as np
import pyterrainsol.exposure as micromet
import pytest
from rasterio.transform import Affine

from pyterrainsol import upwind_profile, horizon_angle, sx, sx_rose, sx_raster

RES = 30.0
# 100 x 100 grid, 30 m cells, origin at the top-left of a projected CRS
TR = Affine(RES, 0.0, 500_000.0, 0.0, -RES, 6_000_000.0)
# centre of cell (50, 50): an affine maps this to index 50.0 exactly, so rays
# cast from here land on cell centres and the expected values stay exact
X0, Y0 = 500_000.0 + 50.5 * RES, 6_000_000.0 - 50.5 * RES


def ramp(slope_per_m, bearing_up):
    """Plane tilted so that elevation rises at a known rate toward `bearing_up`."""
    cols, rows = np.meshgrid(np.arange(100), np.arange(100))
    x = 500_000.0 + (cols + 0.5) * RES
    y = 6_000_000.0 - (rows + 0.5) * RES
    a = np.radians(bearing_up)
    return 1000.0 + slope_per_m * ((x - X0) * np.sin(a) + (y - Y0) * np.cos(a))


def test_profile_reads_a_known_gradient():
    dem = ramp(0.10, 270.0)                      # rises 10% toward the west
    d, z = upwind_profile(dem, TR, X0, Y0, 270.0, max_distance=600.0)
    assert d[0] == 0.0 and z[0] == pytest.approx(1000.0, abs=1e-6)
    assert np.allclose(np.diff(z) / np.diff(d), 0.10, atol=1e-6)


def test_horizon_angle_matches_the_plane_it_was_built_from():
    dem = ramp(np.tan(np.radians(12.0)), 300.0)
    assert horizon_angle(dem, TR, X0, Y0, 300.0, max_distance=600.0) == pytest.approx(12.0, abs=0.01)


def test_reported_distance_points_at_the_controlling_ridge():
    """On a plain with one ridge, the angle and its range are both known."""
    dem = np.full((100, 100), 1000.0)
    ridge_col = 50 + int(round(300.0 / RES))           # 300 m east of the point
    dem[:, ridge_col] = 1100.0
    ang, dist = horizon_angle(dem, TR, X0, Y0, 90.0, max_distance=600.0,
                              return_distance=True)
    assert dist == pytest.approx(300.0, abs=RES)
    assert ang == pytest.approx(np.degrees(np.arctan2(100.0, 300.0)), abs=0.01)


def test_downhill_fetch_is_negative():
    """Looking up a ray that falls away must report exposure, not shelter."""
    dem = ramp(np.tan(np.radians(12.0)), 300.0)
    assert horizon_angle(dem, TR, X0, Y0, 120.0, max_distance=600.0) < -11.0


def test_raising_the_observer_reduces_shelter():
    dem = ramp(np.tan(np.radians(12.0)), 300.0)
    low = horizon_angle(dem, TR, X0, Y0, 300.0, max_distance=600.0)
    high = horizon_angle(dem, TR, X0, Y0, 300.0, max_distance=600.0, z_offset=50.0)
    assert high < low


def test_sx_window_averages_rays_and_is_symmetric_on_a_plane():
    dem = ramp(np.tan(np.radians(12.0)), 300.0)
    one = horizon_angle(dem, TR, X0, Y0, 300.0, max_distance=600.0)
    win = sx(dem, TR, X0, Y0, 300.0, half_window=15.0, n_rays=5, max_distance=600.0)
    assert win < one                      # off-axis rays climb more slowly
    assert win == pytest.approx(
        sx(dem, TR, X0, Y0, 300.0, half_window=15.0, n_rays=5, max_distance=600.0))


def test_rose_weighting_follows_the_weights():
    dem = ramp(np.tan(np.radians(12.0)), 300.0)
    bear = np.array([300.0, 120.0])
    vals, plain = sx_rose(dem, TR, X0, Y0, bear, max_distance=600.0)
    _, upwind_heavy = sx_rose(dem, TR, X0, Y0, bear, weights=[0.99, 0.01],
                              max_distance=600.0)
    assert plain == pytest.approx(vals.mean())
    assert upwind_heavy == pytest.approx(vals[0], abs=0.3)


def test_sx_raster_agrees_with_the_point_function():
    dem = ramp(np.tan(np.radians(12.0)), 300.0)
    field = sx_raster(dem, RES, [300.0], max_distance=600.0, step=30.0)
    point = sx(dem, TR, X0, Y0, 300.0, half_window=7.5, n_rays=3,
               max_distance=600.0, step=30.0)
    assert field[50, 50] == pytest.approx(point, abs=0.05)


def test_sx_raster_rejects_holes():
    dem = ramp(0.1, 270.0); dem[10, 10] = np.nan
    with pytest.raises(ValueError):
        sx_raster(dem, RES, [270.0], max_distance=120.0)


def test_sx_raster_weighting_collapses_to_the_dominant_bearing():
    dem = ramp(np.tan(np.radians(12.0)), 300.0)
    both = sx_raster(dem, RES, [300.0, 120.0], weights=[0.99, 0.01],
                     max_distance=600.0, step=30.0)
    one = sx_raster(dem, RES, [300.0], max_distance=600.0, step=30.0)
    assert both[50, 50] == pytest.approx(one[50, 50], abs=0.4)


def test_sx_raster_samples_the_adjacent_cell():
    """A one-cell step adjacent to the point must be seen, not stepped over."""
    dem = np.full((40, 40), 100.0)
    dem[:, 19] = 110.0                       # wall one cell west of column 20
    field = sx_raster(dem, RES, [270.0], max_distance=300.0, step=4 * RES,
                      sub_rays=(0.0,))       # one ray, so the value is exact
    assert field[20, 20] == pytest.approx(np.degrees(np.arctan2(10.0, RES)), abs=0.01)


# --- MicroMet terrain weighting (Liston & Elder 2006) ----------------------

def _ridge(n=241, res=10.0, height=120.0, half_width=300.0):
    """A smooth north-south ridge, so east and west faces are symmetric."""
    x = (np.arange(n) - n // 2) * res
    prof = height / (1.0 + (x / half_width) ** 2)
    return np.tile(prof, (n, 1))


def test_slope_aspect_compass_convention():
    # Ground falling to the west must read aspect ~270 (the way it faces).
    n, res = 41, 10.0
    z = np.tile(np.arange(n) * res * 0.2, (n, 1))      # rises to the east
    slope, aspect = micromet.slope_aspect(z, res)
    assert np.degrees(slope[20, 20]) == pytest.approx(np.degrees(np.arctan(0.2)), abs=0.5)
    assert aspect[20, 20] == pytest.approx(270.0, abs=1.0)


def test_curvature_sign_convex_positive():
    z = _ridge()
    c = micromet.curvature(z, 10.0, length_scale=100.0)
    mid = z.shape[1] // 2
    assert c[120, mid] > 0                      # crest is convex
    assert c[120, mid - 40] < 0                 # the flank-to-plain hollow is concave


def test_curvature_does_not_wrap_at_the_array_edge():
    """A wrapped neighbour imports terrain from the far side of the window.

    The Avonside farm mask touches the bottom edge of its DEM window, so this
    is not hypothetical: wrapping corrupted 13% of farm cells by up to 0.294
    against a true field range of 0.122.  A ramp has zero curvature throughout
    its interior.  Edge replication still leaves a small residual in the
    boundary band -- the replicated neighbour flattens, which is the honest
    answer where there is no data beyond the edge -- but wrapping imports
    terrain from the far side and the residual is four times larger.
    """
    n, res, k = 81, 10.0, 20
    ramp = np.tile(np.arange(n)[:, None] * 3.0, (1, n))      # rises to the south
    c = micromet.curvature(ramp, res, length_scale=k * res)
    assert np.abs(c[k + 1:-k - 1]).max() == pytest.approx(0.0, abs=1e-12)
    assert np.abs(c[0]).max() < 0.05

    def _wrapped(z):
        out = np.zeros_like(z)
        for dr, dc in ((0, k), (k, 0), (k, k), (k, -k)):
            zp = np.roll(np.roll(z, -dr, 0), -dc, 1)
            zm = np.roll(np.roll(z, dr, 0), dc, 1)
            out += (z - 0.5 * (zp + zm)) / (2.0 * np.hypot(dr, dc) * res)
        return out / 4.0

    assert np.abs(_wrapped(ramp)[0]).max() > 3 * np.abs(c[0]).max()


def test_micromet_windward_exceeds_lee_on_a_symmetric_ridge():
    z = _ridge()
    res, mid = 10.0, z.shape[1] // 2
    W = micromet.micromet_weight(z, res, [270.0], length_scale=100.0)   # wind from the west
    west_face, east_face = W[120, mid - 15], W[120, mid + 15]
    assert west_face > 1.0 > east_face
    # and the asymmetry must reverse with the wind
    We = micromet.micromet_weight(z, res, [90.0], length_scale=100.0)
    assert We[120, mid + 15] > 1.0 > We[120, mid - 15]
    assert We[120, mid + 15] == pytest.approx(west_face, rel=1e-6)


def test_micromet_stays_within_the_published_bounds():
    # gamma_slope + gamma_curvature = 1 and both Omegas are scaled to +/-0.5,
    # so the multiplier cannot leave [0.5, 1.5] whatever the terrain.
    z = _ridge(height=400.0, half_width=80.0)
    W = micromet.micromet_weight(z, 10.0, [270.0, 315.0], [0.6, 0.4])
    assert W.min() >= 0.5 and W.max() <= 1.5


def test_micromet_normalises_the_rose_weights():
    # A real wind rose rarely sums to 1 (ours sums to 0.977).  The result must
    # not depend on that, or every multiplier is silently scaled down.
    z = _ridge()
    raw = [0.493, 0.395, 0.038, 0.030, 0.021]
    brg = [292.5, 315.0, 270.0, 337.5, 225.0]
    a = micromet.micromet_weight(z, 10.0, brg, raw)
    b = micromet.micromet_weight(z, 10.0, brg, np.array(raw) / sum(raw))
    assert np.allclose(a, b, atol=1e-12)
    assert abs(np.median(a) - 1.0) < 0.2        # centred near 1, not near 0.977


def test_micromet_scale_mask_changes_the_normalisation():
    # The scaling is domain-dependent by construction; the mask is how you
    # keep a padded window's extreme cells from setting it for the farm.
    # Mask a flank well away from the crest, so the masked region's maximum
    # slope and curvature are genuinely smaller than the whole window's.
    z = _ridge()
    m = np.zeros(z.shape, bool); m[100:140, 10:40] = True
    full = micromet.micromet_weight(z, 10.0, [270.0], length_scale=100.0)
    masked = micromet.micromet_weight(z, 10.0, [270.0], length_scale=100.0, scale_mask=m)
    assert not np.allclose(full[m], masked[m])
    # the masked normalisation must spread that region over a wider range
    assert np.ptp(masked[m]) > np.ptp(full[m])


def test_micromet_rejects_a_flat_dem():
    with pytest.raises(ValueError):
        micromet.micromet_weight(np.zeros((30, 30)), 10.0, [270.0])


def test_micromet_rejects_mismatched_weights():
    with pytest.raises(ValueError):
        micromet.micromet_weight(_ridge(), 10.0, [270.0, 315.0], [1.0])
