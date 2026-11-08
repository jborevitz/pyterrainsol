import numpy as np
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
