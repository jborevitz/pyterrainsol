"""Pin the published Avonside demand numbers against the formulae.

Priestley-Taylor and the thermodynamics it rests on.  The Penman-Monteith
tests live in the companion package pyterrainwind, alongside the function.
"""

import numpy as np
import pytest

from pyterrainsol.demand import (
    clearness_index, energy_partition, net_longwave, net_radiation,
    priestley_taylor, psychrometric_constant, saturation_vapour_pressure,
    svp_slope,
)

FARM_MEAN_ELEVATION = 1040.0


FARM_MEAN_ELEVATION = 1040.0


def test_saturation_vapour_pressure_at_zero():
    assert saturation_vapour_pressure(0.0) == pytest.approx(0.6108, abs=1e-4)


def test_psychrometric_constant_falls_with_altitude():
    sea = psychrometric_constant(0.0)
    farm = psychrometric_constant(FARM_MEAN_ELEVATION)
    assert sea == pytest.approx(0.0674, abs=5e-4)
    assert farm == pytest.approx(0.0596, abs=5e-4)
    assert farm < sea


def test_partition_term_published_monthly_values():
    # WorldClim 2.1 monthly means at -36.4355, 148.7257
    july = energy_partition(3.5, FARM_MEAN_ELEVATION)
    january = energy_partition(17.5, FARM_MEAN_ELEVATION)
    assert july == pytest.approx(0.483, abs=0.002)
    assert january == pytest.approx(0.680, abs=0.002)
    assert january / july == pytest.approx(1.41, abs=0.02)


def test_altitude_raises_the_evaporating_fraction():
    # The same air temperature evaporates more efficiently up the hill.
    high = energy_partition(10.0, FARM_MEAN_ELEVATION)
    low = energy_partition(10.0, 0.0)
    assert high - low == pytest.approx(0.027, abs=0.005)


def test_svp_slope_rises_steeply_with_temperature():
    # Clausius-Clapeyron: roughly a doubling from 3.5 C to 17.5 C.
    assert svp_slope(17.5) / svp_slope(3.5) == pytest.approx(2.27, abs=0.05)


def test_longwave_loss_is_near_constant_through_the_year():
    # Published range at Avonside: 4.37 to 4.81 MJ m-2 d-1.
    loss = net_longwave(t_min_c=[-1.8, 9.9], t_max_c=[8.9, 25.1], clearness=[0.780, 0.796])
    assert loss.min() > 4.0 and loss.max() < 5.2


def test_winter_net_radiation_is_nearly_all_eaten_by_longwave():
    # June: 7.2 MJ m-2 d-1 arrives, 0.87 survives as net radiation.
    loss = net_longwave(-0.6, 10.0, 0.763)
    rn = net_radiation(7.2, loss)
    assert rn == pytest.approx(0.87, abs=0.15)
    assert rn < 1.0


def test_priestley_taylor_reproduces_published_january_demand():
    et = priestley_taylor(15.19, 17.5, FARM_MEAN_ELEVATION)
    assert et == pytest.approx(5.28, abs=0.05)


def test_priestley_taylor_reproduces_published_june_demand():
    et = priestley_taylor(0.87, 4.7, FARM_MEAN_ELEVATION)
    assert et == pytest.approx(0.23, abs=0.03)


def test_negative_net_radiation_gives_zero_not_condensation():
    assert priestley_taylor(-3.0, 5.0, FARM_MEAN_ELEVATION) == 0.0


def test_net_radiation_is_not_clipped():
    # A deeply shaded midwinter cell really does run a radiation deficit.
    assert net_radiation(2.0, 4.66) < 0.0


def test_clearness_index_matches_the_published_annual_mean():
    assert clearness_index(19.95, 26.26) == pytest.approx(0.760, abs=0.01)


def test_demand_is_monotonic_in_both_drivers():
    base = priestley_taylor(10.0, 10.0, FARM_MEAN_ELEVATION)
    assert priestley_taylor(12.0, 10.0, FARM_MEAN_ELEVATION) > base
    assert priestley_taylor(10.0, 14.0, FARM_MEAN_ELEVATION) > base
