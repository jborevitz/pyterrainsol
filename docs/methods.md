# pyterrainsol — methods

Conventions, traps and limitations. Findings from the block this was developed on
are in the [README](../README.md).

## Why not a temperature-based proxy

Hargreaves and similar run on extraterrestrial radiation and temperature range, so
they are structurally blind to aspect: a north slope and a south slope at the same
elevation return identical values. If aspect is meant to drive evaporative demand,
the radiation has to be computed on the actual inclined, shaded surface.

## PAR is barely a separate map

A single fixed PAR fraction makes the photosynthesis and evaporation maps identical
in different units. They separate only because diffuse light is spectrally richer in
PAR (~0.60) than direct beam (~0.43), so the split is applied per component. The
resulting spread is real but small. Say so rather than overselling it.

## Radiation is a map; temperature is a clock

The Priestley-Taylor partition term Δ/(Δ+γ) is the fraction of available energy that
leaves the surface as water rather than as hot air. Because saturation vapour
pressure follows Clausius-Clapeyron, it climbs steeply with temperature.

That term varies **in time and barely in space**: on any given day every cell in a
small catchment shares an air temperature to within what the elevation range buys at
the lapse rate. Radiation is the opposite. This is why a coarse climate product
costs little — it supplies only the part that does not vary across a farm — and why
the DEM has to supply all the rest.

Altitude runs the wrong way for water. Lower pressure means a smaller psychrometric
constant, so a given energy input evaporates a slightly larger fraction. At 1040 m
γ is 0.0596 against the sea-level 0.0674, about three percentage points on the
partition term year-round.

Net longwave is close to a **constant tax** rather than a seasonal one: the warmer
summer surface emits more, but the drier winter air absorbs less, and the two
largely cancel, while incoming shortwave collapses in winter. Most of the gap
between the seasonal swing in demand and the swing in light comes from this term,
not from the partition term.

Terrain shading is **not cloud**. Feeding a shaded cell's own low radiation into the
longwave cloudiness factor treats a hillside as an overcast sky; pass a spatially
uniform clearness instead.

## Wind, and why it is not in the demand model

Priestley-Taylor has **no aerodynamic term**; α = 1.26 absorbs a typical advection
contribution, so bolting wind onto it would double-count. Wind enters properly only
through Penman-Monteith, whose aerodynamic term scales as **u × VPD** — a wind speed
without humidity multiplies an unknown and contributes nothing.

`pyterrainsol.wind` therefore provides the pieces rather than the model.

Three conventions, each of which has cost someone a day:

- **Direction is where wind comes FROM**, degrees clockwise from north. A westerly
  (blowing toward the east, `u > 0`) reads 270.
- **Never average directions.** Average the components. Use `directional_constancy`
  — vector-mean over scalar-mean — to find out whether a mean direction means
  anything before quoting one.
- **Reanalysis timestamps are usually UTC.** Shift to local time before any diurnal
  analysis, or the afternoon maximum lands at dawn.

### Fitting a wind climatology, and how to be wrong at it

- **Do not calibrate on one year.** At a synoptically dominated site the best year in
  a twenty-year record can have half again the explanatory power of the pooled fit,
  and the largest annual amplitude. Fit the pool; use single years to measure spread.
- **Keep the model small.** Annual and diurnal sine/cosine pairs plus their
  interaction — nine terms — is usually the ceiling. The interaction is what lets the
  diurnal amplitude vary through the year and is often worth as much as the annual
  cycle itself. Higher orders fit weather.
- **Expect a low R².** Where fronts rather than the sun set the wind, a climatology
  explains a small fraction of variance by construction. That is a property of the
  place, not a failure of the fit.
- **Report amplitudes and phases, not coefficients**, and check which is which. A
  tensor product can be laid out seasonal-major or diurnal-major, and reading the
  wrong pair silently swaps the annual and daily cycles.
- **Weibull `k` is dimensionless and scale-free.** Doubling every speed moves `c` and
  leaves `k` alone, so a low `k` is never an artefact of a low background — it means
  a genuinely wider relative spread, typically a mixture of calm and synoptic
  regimes. Rank sites or seasons by mean speed and you can get their extremes
  backwards, since exceedance is `exp(-(x/c)^k)` and depends sharply on `k`.
- **Fit roughness on near-neutral steps only.** A stable nocturnal boundary layer
  steepens the two-height profile and inverts to a spuriously large z0. And if the
  two speeds come from a reanalysis you recover the *model's* effective roughness,
  a property of its land-surface configuration rather than a measurement of the
  landscape.

### Plotting a direction climatology

Use a **cyclic** colormap so 359° and 1° are the same colour and north carries no
seam, and a compass ring rather than a linear bar — a bar still implies the ends are
far apart. Note that `twilight_shifted` is near-white at its midpoint, so masked
cells must take a tone that appears nowhere in the ramp or they read as southerlies.
Scale opacity by directional constancy rather than by speed: fading what is
directionless is more honest than fading what is merely slow, and a hard speed cut
erases real light breezes.

## What this does not do

- **No wind in the demand model**, for the reason above. The midwinter floor is a
  lower bound: net radiation approaches zero while the air keeps moving.
- **No weather.** A clear-sky model scaled by a clearness index has no individual
  cloudy days — these are climatological days, not a particular year.
- **No water limit.** This is demand, not actual evapotranspiration: an upper bound
  on what the atmosphere can take from a well-watered surface. What it does take is
  a soil-water question and belongs in a different package.
- **Humidity is inferred** where not supplied, as saturation vapour pressure at the
  daily minimum temperature — the FAO-56 fallback.
