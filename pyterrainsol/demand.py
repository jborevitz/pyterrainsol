"""Atmospheric evaporative demand from terrain-resolved shortwave radiation.

The spatial field comes from :mod:`pyterrainsol.geometry` -- slope, aspect and
horizon shading at the resolution of the DEM.  Everything in this module is the
*temporal* half: how a given amount of energy is split between heating the air
and evaporating water, which depends on air temperature and barometric pressure
and varies over the year rather than across a paddock.

Formulae follow FAO Irrigation and Drainage Paper 56 (Allen et al. 1998) except
the final step, which is Priestley and Taylor (1972) rather than
Penman-Monteith.  Priestley-Taylor has no aerodynamic term: it assumes the
surface is well watered and that advection contributes a fixed fraction,
absorbed into ``alpha``.  It is therefore the right instrument when the spatial
signal you resolve is radiation, and the wrong one when wind dominates --
notably in winter, when net radiation is near zero but the air is still moving.
See the package README for the size of that bias.
"""

import numpy as np

LATENT_HEAT = 2.45          # MJ kg-1, latent heat of vaporisation near 20 C
STEFAN = 4.903e-9           # MJ K-4 m-2 d-1
ALPHA_PT = 1.26             # Priestley-Taylor coefficient, well-watered surface
ALBEDO_GRASS = 0.23         # FAO-56 reference short crop


def saturation_vapour_pressure(t_c):
    """Saturation vapour pressure (kPa) at air temperature `t_c` in degrees C."""
    t = np.asarray(t_c, float)
    return 0.6108 * np.exp(17.27 * t / (t + 237.3))


def svp_slope(t_c):
    """Delta, the slope d(es)/dT in kPa per degree C (FAO-56 eq. 13)."""
    t = np.asarray(t_c, float)
    return 4098.0 * saturation_vapour_pressure(t) / (t + 237.3) ** 2


def psychrometric_constant(elevation_m):
    """Gamma (kPa per degree C) from barometric pressure at `elevation_m`.

    Altitude matters: at 1040 m the pressure is about 89.6 kPa, so gamma is
    0.0596 rather than the sea-level 0.0674, which raises the evaporating
    fraction of any given energy input by roughly three percentage points.
    """
    z = np.asarray(elevation_m, float)
    pressure = 101.3 * ((293.0 - 0.0065 * z) / 293.0) ** 5.26
    return 0.000665 * pressure


def energy_partition(t_c, elevation_m):
    """Delta / (Delta + gamma): the fraction of available energy that evaporates.

    Rises steeply with temperature because saturation vapour pressure follows
    Clausius-Clapeyron.  Cold air has little room to accept water, so energy
    delivered to a cold surface mostly warms the air instead.
    """
    d = svp_slope(t_c)
    return d / (d + psychrometric_constant(elevation_m))


def net_longwave(t_min_c, t_max_c, clearness, vapour_pressure=None):
    """Outgoing net longwave radiation, MJ m-2 d-1 (FAO-56 eq. 39).

    `clearness` is Rs/Rso for an *unobstructed* site.  Pass a scalar or a daily
    series, not a per-cell field: terrain shading is not cloud, and feeding a
    shaded cell's own low radiation in here would wrongly treat a hillside as
    an overcast sky.

    `vapour_pressure` (kPa) defaults to the saturation value at `t_min_c`, the
    FAO-56 fallback when humidity is unmeasured.
    """
    t_min = np.asarray(t_min_c, float)
    t_max = np.asarray(t_max_c, float)
    ea = saturation_vapour_pressure(t_min) if vapour_pressure is None else np.asarray(vapour_pressure, float)
    temperature_term = ((t_max + 273.16) ** 4 + (t_min + 273.16) ** 4) / 2.0
    return (STEFAN * temperature_term * (0.34 - 0.14 * np.sqrt(ea))
            * (1.35 * np.asarray(clearness, float) - 0.35))


def net_radiation(shortwave, longwave_loss, albedo=ALBEDO_GRASS):
    """Net radiation, MJ m-2 d-1, from incoming shortwave and the longwave loss.

    `shortwave` is the per-cell actual-sky value; `longwave_loss` comes from
    :func:`net_longwave` and is spatially uniform for the reason given there.
    The result is deliberately not clipped -- midwinter net radiation on a
    shaded slope is genuinely negative, and hiding that would misrepresent it.
    """
    return (1.0 - albedo) * np.asarray(shortwave, float) - np.asarray(longwave_loss, float)


def priestley_taylor(net_rad, t_c, elevation_m, alpha=ALPHA_PT, soil_heat=0.0):
    """Potential evapotranspiration, mm per day.

    Demand for a well-watered surface -- an upper bound on what the atmosphere
    can take, not what it does take.  Negative values are clipped to zero.
    """
    available = np.asarray(net_rad, float) - np.asarray(soil_heat, float)
    et = alpha * energy_partition(t_c, elevation_m) * available / LATENT_HEAT
    return np.maximum(et, 0.0)


def clearness_index(observed_shortwave, clear_sky_shortwave):
    """Rs/Rso from an observed radiation climatology and the clear-sky model.

    Both arguments must describe a flat unobstructed surface at the same place
    and averaging period.  Values outside roughly 0.4-0.85 mean the two are
    not comparable rather than that the sky is unusual.
    """
    return np.asarray(observed_shortwave, float) / np.asarray(clear_sky_shortwave, float)
