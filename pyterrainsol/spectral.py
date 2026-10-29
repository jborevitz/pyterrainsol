"""Splitting shortwave into the photosynthetic and evaporative fractions.

A single fixed PAR fraction makes the photosynthesis map and the
evapotranspiration map the same picture in different units.  They separate
only because diffuse light is spectrally richer in PAR than direct beam, so
the split must be applied to the two components separately.  The resulting
spread is real but modest -- on steep terrain the PAR share of total shortwave
varies by only a few percent between sunlit and shaded ground.
"""
from __future__ import annotations
import numpy as np

PAR_FRACTION_BEAM = 0.43      # energy fraction, 400-700 nm, direct beam
PAR_FRACTION_DIFFUSE = 0.60   # diffuse skylight is bluer, so richer in PAR
UMOL_PER_JOULE_PAR = 4.57     # micromol photons per joule of PAR


def par_from_components(beam, diffuse,
                        f_beam=PAR_FRACTION_BEAM, f_diffuse=PAR_FRACTION_DIFFUSE):
    """PAR energy (same units as inputs) from beam and diffuse shortwave."""
    return f_beam * np.asarray(beam) + f_diffuse * np.asarray(diffuse)


def par_fraction(beam, diffuse, **kw):
    """PAR as a fraction of total shortwave, per cell."""
    b = np.asarray(beam); d = np.asarray(diffuse)
    return par_from_components(b, d, **kw) / (b + d)


def daily_light_integral(par_mj):
    """mol photons m-2 d-1 from PAR in MJ m-2 d-1."""
    return np.asarray(par_mj) * UMOL_PER_JOULE_PAR
