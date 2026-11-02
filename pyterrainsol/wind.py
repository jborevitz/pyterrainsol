"""Wind: components, circular statistics, the log profile, and climatology.

Wind reaches evapotranspiration only through the aerodynamic term of
Penman-Monteith, which scales as ``u x VPD``.  Priestley-Taylor has no such
term, so nothing here changes :mod:`pyterrainsol.demand` on its own -- these
are the pieces a Penman-Monteith step would need, plus the climatology tools
used to characterise a site's wind before trusting any of it.

Three conventions, each of which has cost someone a day:

* **Direction is the direction wind comes FROM**, degrees clockwise from north.
  A westerly (blowing toward the east, ``u > 0``) reads 270.
* **Never average directions.**  Average the components, or use
  :func:`directional_constancy` to find out whether a mean direction means
  anything at all.
* **Reanalysis timestamps are usually UTC.**  Shift to local time before any
  diurnal analysis or the afternoon maximum lands at dawn.
"""

import numpy as np

KMH_PER_MS = 3.6
VON_KARMAN = 0.41


def speed_direction(u, v):
    """Speed and meteorological direction from eastward/northward components.

    Returns ``(speed, direction)`` in the units of the input and in degrees.
    Direction is where the wind comes FROM: 0 is northerly, 90 easterly,
    180 southerly, 270 westerly.
    """
    u = np.asarray(u, float); v = np.asarray(v, float)
    return np.hypot(u, v), (270.0 - np.degrees(np.arctan2(v, u))) % 360.0


def components(speed, direction_from_deg):
    """Inverse of :func:`speed_direction`."""
    speed = np.asarray(speed, float)
    theta = np.radians(270.0 - np.asarray(direction_from_deg, float))
    return speed * np.cos(theta), speed * np.sin(theta)


def directional_constancy(u, v):
    """Vector-mean speed divided by scalar-mean speed, in [0, 1].

    1 means the wind is always from the same quarter; 0 means it averages to
    nothing and a "prevailing direction" is meaningless.  Use it to decide
    whether a mean direction is worth reporting, and as the opacity channel
    when plotting a direction climatology -- fading what is directionless is
    more honest than fading what is merely slow.
    """
    u = np.asarray(u, float); v = np.asarray(v, float)
    scalar = np.hypot(u, v).mean()
    if scalar == 0:
        return 0.0
    return float(np.hypot(u.mean(), v.mean()) / scalar)


def roughness_length(ratio, z_lower=10.0, z_upper=50.0):
    """Roughness length z0 from the speed ratio between two heights.

    Inverts the neutral logarithmic profile, ``u(z) = (u*/k) ln(z/z0)``:

        z0 = exp[ (ln z_upper - ratio * ln z_lower) / (1 - ratio) ]

    `ratio` is ``u(z_upper) / u(z_lower)`` and must exceed 1.

    Fit this on near-neutral conditions only -- windy, well-mixed, daytime.
    A stable nocturnal boundary layer steepens the profile and inverts to a
    spuriously large z0; the ground does not get rougher after dark.

    If the two speeds come from a reanalysis rather than instruments, what you
    recover is the *model's* effective roughness, which is a property of its
    land-surface configuration and not a measurement of the landscape.
    """
    r = np.asarray(ratio, float)
    if np.any(r <= 1.0):
        raise ValueError("ratio must exceed 1; the upper level must be windier")
    return np.exp((np.log(z_upper) - r * np.log(z_lower)) / (1.0 - r))


def scale_speed(speed, z_from, z_to, z0, displacement=0.0):
    """Move a wind speed between heights along a neutral log profile."""
    z0 = np.asarray(z0, float)
    num = np.log((z_to - displacement) / z0)
    den = np.log((z_from - displacement) / z0)
    return np.asarray(speed, float) * num / den


def harmonic_design(day_of_year, hour, n_seasonal=1, n_diurnal=1):
    """Tensor-product harmonic design matrix: seasonal x diurnal, with interaction.

    Columns are every seasonal term times every diurnal term, so the width is
    ``(1 + 2*n_seasonal) * (1 + 2*n_diurnal)``.  The interaction columns are
    what let the diurnal amplitude vary through the year.

    `day_of_year` may be fractional.  Orders above 2 x 2 buy almost nothing at
    a synoptically dominated site and start fitting one year's weather.
    """
    d = np.asarray(day_of_year, float); h = np.asarray(hour, float)
    seas = [np.ones_like(d)]
    for j in range(1, n_seasonal + 1):
        seas += [np.cos(2 * np.pi * j * d / 365.25), np.sin(2 * np.pi * j * d / 365.25)]
    diur = [np.ones_like(h)]
    for k in range(1, n_diurnal + 1):
        diur += [np.cos(2 * np.pi * k * h / 24.0), np.sin(2 * np.pi * k * h / 24.0)]
    return np.column_stack([a * b for a in seas for b in diur])


def fit_harmonic(values, day_of_year, hour, n_seasonal=1, n_diurnal=1):
    """Least-squares fit of :func:`harmonic_design`. Returns ``(coef, r_squared)``.

    Fit the *components* separately for direction, never the angle itself.
    Expect a low R-squared wherever fronts rather than the sun set the wind.
    """
    X = harmonic_design(day_of_year, hour, n_seasonal, n_diurnal)
    y = np.asarray(values, float)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    denom = np.sum((y - y.mean()) ** 2)
    return coef, float(1.0 - np.sum(resid ** 2) / denom) if denom else np.nan


def harmonic_summary(coef):
    """Interpretable quantities from a 1x1 (nine-term) harmonic fit.

    Returns mean, annual amplitude and peak day, diurnal amplitude and peak
    hour, and the diurnal amplitude and peak hour evaluated at midsummer and
    midwinter, which is where the interaction terms show themselves.
    """
    coef = np.asarray(coef, float)
    if coef.size != 9:
        raise ValueError(f"expected a nine-term 1x1 fit, got {coef.size} coefficients")
    # Column order follows harmonic_design's tensor product, seasonal-major:
    #   [1, cos h, sin h, cos d, cos d cos h, cos d sin h, sin d, sin d cos h, sin d sin h]
    m = coef[0]
    dc, ds = coef[1], coef[2]            # diurnal main effect
    ac, as_ = coef[3], coef[6]           # annual main effect
    i1, i2 = coef[4], coef[5]            # cos d x diurnal
    i3, i4 = coef[7], coef[8]            # sin d x diurnal
    amp = lambda c, s: float(np.hypot(c, s))
    day = lambda c, s: float((np.degrees(np.arctan2(s, c)) / 360 * 365.25) % 365.25)
    hr = lambda c, s: float((np.degrees(np.arctan2(s, c)) / 360 * 24) % 24)
    out = {"mean": float(m), "annual_amplitude": amp(ac, as_), "annual_peak_day": day(ac, as_),
           "diurnal_amplitude": amp(dc, ds), "diurnal_peak_hour": hr(dc, ds)}
    for label, doy in (("midsummer", 0.0), ("midwinter", 181.0)):
        C, S = np.cos(2 * np.pi * doy / 365.25), np.sin(2 * np.pi * doy / 365.25)
        c_, s_ = dc + i1 * C + i3 * S, ds + i2 * C + i4 * S
        out[f"diurnal_amplitude_{label}"] = amp(c_, s_)
        out[f"diurnal_peak_hour_{label}"] = hr(c_, s_)
    return out


def weibull_fit(speed, floor=0.1):
    """Maximum-likelihood Weibull shape ``k`` and scale ``c``, location fixed at 0.

    `k` is **dimensionless and scale-free**: multiply every speed by two and
    `k` does not move, only `c` does.  So a low `k` is never an artefact of a
    low background -- it means the distribution is relatively more spread,
    which at this kind of site comes from a mixture of calm and synoptic
    regimes rather than from a long tail alone.

    Values at or below `floor` are dropped; an exact zero has no likelihood.
    """
    from scipy.stats import weibull_min
    x = np.asarray(speed, float)
    x = x[np.isfinite(x) & (x > floor)]
    if x.size < 50:
        raise ValueError(f"need at least 50 usable speeds, got {x.size}")
    k, _loc, c = weibull_min.fit(x, floc=0)
    return float(k), float(c)


def exceedance(k, c, threshold):
    """Fraction of time above `threshold` for a Weibull(k, c): exp(-(x/c)^k).

    Sensitive to `k`, which is why ranking sites or seasons by mean speed can
    order their extremes backwards.
    """
    return float(np.exp(-((float(threshold) / float(c)) ** float(k))))
