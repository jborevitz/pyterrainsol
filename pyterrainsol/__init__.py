"""Terrain-corrected solar radiation: beam, diffuse, PAR and daily light integral.

Takes slope, aspect, horizon angles and sky view factor (e.g. from
`pylidarterrain`) and returns daily shortwave on the inclined, shaded surface,
then the evaporative demand that follows from it by Priestley-Taylor.

Why this rather than a temperature-based PET proxy: Hargreaves and friends run
on extraterrestrial radiation and temperature range, so they are structurally
blind to aspect and return identical values for a north and a south slope at
the same elevation.  If aspect is supposed to drive evaporative demand, the
radiation has to be computed on the actual surface.

Wind is deliberately absent.  The aerodynamic term, the Sx and MicroMet
exposure fields and the wind climatology live in the companion package
`pyterrainwind`, which depends on this one.  A solar-only user never installs
them.
"""
from .geometry import solar_position, daily_radiation
from .spectral import par_fraction, par_from_components, daily_light_integral
from .demand import (saturation_vapour_pressure, svp_slope, psychrometric_constant,
                     energy_partition, net_longwave, net_radiation, priestley_taylor,
                     clearness_index)

__all__ = ["solar_position","daily_radiation","par_fraction",
           "par_from_components","daily_light_integral",
           "saturation_vapour_pressure","svp_slope","psychrometric_constant",
           "energy_partition","net_longwave","net_radiation","priestley_taylor",
           "clearness_index"]
__version__ = "0.2.0"
