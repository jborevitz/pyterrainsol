"""Pin the published Avonside demand numbers against the formulae."""

import numpy as np
import pytest

from pyterrainsol.demand import (
    clearness_index, energy_partition, net_longwave, net_radiation,
    penman_monteith, vapour_pressure_deficit, wind_at_2m,
    priestley_taylor, psychrometric_constant, saturation_vapour_pressure,
    svp_slope,
)

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


# --- Penman-Monteith, the wind-aware companion to Priestley-Taylor ----------

def test_wind_at_2m_matches_fao_table():
    # FAO-56 Table 4: the 10 m conversion factor is 0.748.
    assert wind_at_2m(1.0, 10.0) == pytest.approx(0.748, abs=0.001)
    # 2 m is the reference height, so the factor is 1 to within rounding
    # (4.87/ln(67.8*2-5.42) = 1.00022, not exactly unity).
    assert wind_at_2m(3.3, 2.0) == pytest.approx(3.3, rel=1e-3)


def test_vpd_uses_mean_of_saturation_not_saturation_of_mean():
    # The SVP curve is convex, so es(Tmin),es(Tmax) averaged always exceeds
    # es of the mean temperature.  A build that uses the latter under-reads.
    tmin, tmax = 2.0, 24.0
    vpd = vapour_pressure_deficit(tmin, tmax, 0.0)
    naive = saturation_vapour_pressure(0.5 * (tmin + tmax))
    assert vpd > naive
    # Never negative, even when the air is supersaturated in the record.
    assert vapour_pressure_deficit(5.0, 6.0, 3.0) == 0.0


def test_penman_monteith_assembles_the_fao56_equation():
    """Check the ET0 assembly against the equation written out by hand.

    Delta and gamma are typed as literals here rather than taken from the
    module, so this exercises the assembly of FAO-56 eq. 6 independently of
    svp_slope and psychrometric_constant (which have their own tests).
    """
    T, Z, u2, Rn, G, vpd = 30.2, 2.0, 2.0, 14.33, 0.14, 1.15
    delta, gamma = 0.2458, 0.06735          # at T = 30.2 C and 2 m altitude
    expected = (0.408 * delta * (Rn - G)
                + gamma * (900.0 / (T + 273.0)) * u2 * vpd) \
        / (delta + gamma * (1.0 + 0.34 * u2))
    assert expected == pytest.approx(5.246, abs=0.002)   # pins the arithmetic
    assert penman_monteith(Rn, T, Z, u2, vpd, soil_heat=G) \
        == pytest.approx(expected, rel=1e-3)


def test_penman_monteith_split_sums_to_total():
    rad, aero = penman_monteith(12.0, 18.0, 861.0, 2.4, 0.9, split=True)
    total = penman_monteith(12.0, 18.0, 861.0, 2.4, 0.9)
    assert rad + aero == pytest.approx(total, rel=1e-12)
    assert aero > 0.0


def test_aerodynamic_term_vanishes_in_still_air_but_radiation_does_not():
    # With u2 = 0 Penman-Monteith keeps only the radiation term, which is the
    # structural difference from Priestley-Taylor: PT has no u at all.
    rad, aero = penman_monteith(12.0, 18.0, 861.0, 0.0, 0.9, split=True)
    assert aero == 0.0
    assert rad > 0.0


def test_wind_raises_demand_monotonically():
    u = np.array([0.0, 1.0, 2.0, 4.0, 8.0])
    et = penman_monteith(12.0, 18.0, 861.0, u, 0.9)
    assert np.all(np.diff(et) > 0)


def test_priestley_taylor_ignores_wind_by_construction():
    # Guards the "two options" contract: PT must stay a wind-free baseline, so
    # anyone adding a u argument to it has made an error.
    import inspect
    assert "wind" not in inspect.signature(priestley_taylor).parameters


def test_silo_vapour_pressure_unit_trap():
    """SILO publishes vp in hPa; feeding it straight in zeroes the aero term.

    This is a real failure mode: es - ea goes negative, the clip sends VPD to
    zero, and ET0 silently loses its entire aerodynamic term rather than
    raising.  The test documents the symptom so it is recognised next time.
    """
    tmin, tmax, vp_hpa = 4.8, 18.4, 9.93
    wrong = vapour_pressure_deficit(tmin, tmax, vp_hpa)
    right = vapour_pressure_deficit(tmin, tmax, vp_hpa / 10.0)
    assert wrong == 0.0
    assert right == pytest.approx(0.495, abs=0.01)
