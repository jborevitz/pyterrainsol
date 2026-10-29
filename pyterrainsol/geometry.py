"""Solar position and clear-sky radiation on an inclined, shaded surface."""
from __future__ import annotations
import numpy as np

SOLAR_CONSTANT = 1367.0      # W m-2
DEFAULT_TAU = 0.75           # broadband clear-sky atmospheric transmittance
DIFFUSE_FRACTION = 0.30      # of the beam depleted by the atmosphere


def solar_position(doy: int, hour: float, lat_deg: float):
    """Return (altitude, azimuth) in radians; azimuth 0 = north, clockwise.

    `hour` is solar time, -12..12.  Altitude <= 0 means the sun is down.
    """
    lat = np.radians(lat_deg)
    dec = np.radians(23.45) * np.sin(2 * np.pi * (284 + doy) / 365.0)
    H = np.radians(15.0 * hour)
    sin_alt = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(H)
    alt = np.arcsin(np.clip(sin_alt, -1, 1))
    cos_az = (np.sin(dec) * np.cos(lat) - np.cos(dec) * np.sin(lat) * np.cos(H)) / np.cos(alt)
    az = np.arccos(np.clip(cos_az, -1, 1))
    if hour > 0:
        az = 2 * np.pi - az          # afternoon sun swings west
    return alt, az


def daily_radiation(slope_deg, aspect_deg, horizon, sky_view, lat_deg,
                    days=None, dt_hours=0.5, tau=DEFAULT_TAU):
    """Daily beam and diffuse shortwave on the surface, MJ m-2 d-1.

    Parameters are 1-D arrays over cells; `horizon` is (n_azimuth, n_cells)
    in degrees.  Returns (beam, diffuse) each shaped (n_days, n_cells).
    """
    slope = np.radians(np.asarray(slope_deg, float))
    aspect = np.radians(np.asarray(aspect_deg, float))
    hor = np.asarray(horizon, float)
    svf = np.asarray(sky_view, float)
    n_az = hor.shape[0]
    if days is None:
        days = np.arange(1, 366)
    days = np.asarray(days)
    lat = np.radians(lat_deg)
    beam = np.zeros((len(days), slope.size), dtype=np.float32)
    diff = np.zeros((len(days), slope.size), dtype=np.float32)
    hours = np.arange(-12, 12, dt_hours) + dt_hours / 2
    for i, n in enumerate(days):
        dec = np.radians(23.45) * np.sin(2 * np.pi * (284 + n) / 365.0)
        ecc = 1 + 0.033 * np.cos(2 * np.pi * n / 365.0)
        for h in hours:
            H = np.radians(15.0 * h)
            sa = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(H)
            if sa <= 0.01:
                continue
            alt = np.arcsin(sa)
            cos_az = (np.sin(dec) * np.cos(lat) - np.cos(dec) * np.sin(lat) * np.cos(H)) / np.cos(alt)
            az = np.arccos(np.clip(cos_az, -1, 1))
            if h > 0:
                az = 2 * np.pi - az
            airmass = 1.0 / max(sa, 0.05)
            bn = SOLAR_CONSTANT * ecc * tau ** airmass          # direct normal
            cosi = np.clip(np.cos(slope) * sa +
                           np.sin(slope) * np.cos(alt) * np.cos(az - aspect), 0, None)
            k = int(round(np.degrees(az) / (360.0 / n_az))) % n_az
            lit = (np.degrees(alt) > hor[k]).astype(np.float32)
            beam[i] += (bn * cosi * lit) * dt_hours * 3600 / 1e6
            dh = SOLAR_CONSTANT * ecc * sa * DIFFUSE_FRACTION * (1 - tau ** airmass)
            diff[i] += (dh * svf) * dt_hours * 3600 / 1e6
    return beam, diff
