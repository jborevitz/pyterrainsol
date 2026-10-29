import numpy as np, pandas as pd, pathlib, pytest
from pyterrainsol import par_fraction, par_from_components, daily_light_integral
from pyterrainsol.geometry import solar_position

HERE = pathlib.Path(__file__).parent


def test_par_fraction_bounded_by_the_two_components():
    b = np.array([10.0, 1.0]); d = np.array([1.0, 10.0])
    f = par_fraction(b, d)
    assert (f >= 0.43 - 1e-9).all() and (f <= 0.60 + 1e-9).all()
    assert f[1] > f[0]        # diffuse-dominated ground is richer in PAR


def test_published_par_fraction_range():
    g = pd.read_csv(HERE / "farm_grid3m_sample.csv")
    assert g.par_fraction.min() >= 0.450 and g.par_fraction.max() <= 0.481


def test_dli_conversion():
    assert daily_light_integral(1.0) == pytest.approx(4.57)


def test_southern_hemisphere_sun_is_north_at_noon():
    alt, az = solar_position(172, 0.0, -36.4355)      # winter solstice, noon
    assert alt > 0
    deg = np.degrees(az) % 360
    assert deg < 15 or deg > 345, f"midday sun should be due north, got {deg:.0f} deg"
