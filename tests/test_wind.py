"""Pin the wind conventions and the published Avonside climatology."""

import numpy as np
import pytest

from pyterrainsol.wind import (
    components, directional_constancy, exceedance, fit_harmonic, harmonic_design,
    harmonic_summary, roughness_length, scale_speed, speed_direction, weibull_fit,
)


def test_a_westerly_reads_270():
    # u > 0 is wind blowing toward the east, i.e. coming FROM the west.
    spd, drc = speed_direction(5.0, 0.0)
    assert spd == pytest.approx(5.0) and drc == pytest.approx(270.0)


def test_each_cardinal_direction():
    for u, v, expect in [(0, -1, 0), (-1, 0, 90), (0, 1, 180), (1, 0, 270)]:
        assert speed_direction(u, v)[1] == pytest.approx(expect)


def test_components_round_trip():
    spd = np.array([1.0, 4.5, 12.0]); drc = np.array([15.0, 200.0, 337.5])
    assert speed_direction(*components(spd, drc))[0] == pytest.approx(spd)
    assert speed_direction(*components(spd, drc))[1] == pytest.approx(drc)


def test_constancy_bounds():
    steady_u, steady_v = components(np.full(100, 6.0), np.full(100, 290.0))
    assert directional_constancy(steady_u, steady_v) == pytest.approx(1.0, abs=1e-9)
    box = np.linspace(0, 360, 361)[:-1]
    assert directional_constancy(*components(np.full(360, 6.0), box)) < 1e-9


def test_roughness_reproduces_the_published_avonside_value():
    # 20 years of BARRA at the farm cell: median u50/u10 = 1.515.
    assert roughness_length(1.515) == pytest.approx(0.441, abs=0.002)
    assert roughness_length(1.498) == pytest.approx(0.395, abs=0.002)


def test_roughness_rises_with_the_ratio_and_rejects_bad_input():
    assert roughness_length(1.60) > roughness_length(1.45)
    for bad in (1.0, 0.9):
        with pytest.raises(ValueError):
            roughness_length(bad)


def test_scale_speed_is_consistent_with_the_inversion():
    # Recover the ratio the roughness came from.
    z0 = roughness_length(1.515)
    assert scale_speed(1.0, 10.0, 50.0, z0) == pytest.approx(1.515, abs=1e-6)


def test_design_matrix_width_is_the_tensor_product():
    d = np.arange(0, 365, 7.0); h = np.linspace(0, 24, len(d))
    assert harmonic_design(d, h, 1, 1).shape[1] == 9
    assert harmonic_design(d, h, 2, 2).shape[1] == 25
    assert harmonic_design(d, h, 2, 3).shape[1] == 35


def test_fit_recovers_a_known_signal():
    rng = np.random.default_rng(0)
    d = rng.uniform(0, 365, 20000); h = rng.uniform(0, 24, 20000)
    truth = (10.0 + 2.0 * np.cos(2 * np.pi * (d - 275) / 365.25)
             + 3.0 * np.cos(2 * np.pi * (h - 15) / 24))
    coef, r2 = fit_harmonic(truth, d, h)
    s = harmonic_summary(coef)
    assert r2 > 0.999
    assert s["mean"] == pytest.approx(10.0, abs=0.05)
    assert s["annual_amplitude"] == pytest.approx(2.0, abs=0.05)
    assert s["annual_peak_day"] == pytest.approx(275, abs=2)
    assert s["diurnal_amplitude"] == pytest.approx(3.0, abs=0.05)
    assert s["diurnal_peak_hour"] == pytest.approx(15.0, abs=0.1)


def test_summary_rejects_the_wrong_order():
    with pytest.raises(ValueError):
        harmonic_summary(np.zeros(25))


def test_weibull_shape_is_scale_free():
    # The result that settles whether a low morning k is a units artefact.
    rng = np.random.default_rng(1)
    x = rng.weibull(1.6, 40000) * 9.0
    k1, c1 = weibull_fit(x)
    k2, c2 = weibull_fit(x * 2.0)
    # Exact in theory; scipy's optimiser converges from a scale-dependent
    # starting point, so agreement is ~0.3% rather than to machine precision.
    assert k2 == pytest.approx(k1, rel=0.01)
    assert c2 == pytest.approx(2 * c1, rel=0.01)
    assert k1 == pytest.approx(1.6, abs=0.05)


def test_exceedance_matches_the_published_contrast():
    # Spring mornings are windier in the tail despite a LOWER scale than
    # summer afternoons, because their shape parameter is lower.
    spring = exceedance(1.41, 12.0, 30.0)
    summer = exceedance(2.61, 15.7, 30.0)
    assert spring == pytest.approx(0.0262, abs=0.002)
    assert summer == pytest.approx(0.0044, abs=0.002)
    assert spring > 5 * summer


def test_exceedance_is_monotone():
    assert exceedance(2.0, 12.0, 20.0) > exceedance(2.0, 12.0, 40.0)
    assert exceedance(2.0, 18.0, 30.0) > exceedance(2.0, 12.0, 30.0)
