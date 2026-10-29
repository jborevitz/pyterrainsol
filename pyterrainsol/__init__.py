"""Terrain-corrected solar radiation: beam, diffuse, PAR and daily light integral.

Takes slope, aspect, horizon angles and sky view factor (e.g. from
`pylidarterrain`) and returns daily shortwave on the inclined, shaded surface.

Why this rather than a temperature-based PET proxy: Hargreaves and friends run
on extraterrestrial radiation and temperature range, so they are structurally
blind to aspect and return identical values for a north and a south slope at
the same elevation.  If aspect is supposed to drive evaporative demand, the
radiation has to be computed on the actual surface.
"""
from .geometry import solar_position, daily_radiation
from .spectral import par_fraction, par_from_components, daily_light_integral

__all__ = ["solar_position","daily_radiation","par_fraction",
           "par_from_components","daily_light_integral"]
__version__ = "0.1.0"
