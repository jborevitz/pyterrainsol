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

## Wind is deliberately absent

Priestley-Taylor has **no aerodynamic term**; α = 1.26 absorbs a typical advection
contribution, so bolting a wind speed onto it would double-count. Wind enters
properly only through Penman-Monteith, whose aerodynamic term scales as
**u × VPD** — a wind speed without humidity multiplies an unknown and contributes
nothing.

That model, the terrain exposure fields that feed it and the wind climatology that
forces it are in [`pyterrainwind`](https://github.com/jborevitz/pyterrainwind),
which depends on this package for net radiation. The split is along the dependency:
wind needs the solar chain, the solar chain never needs wind. A solar-only user
installs this package alone and is not asked to take a position on wind.

## What this does not do

- **No wind in the demand model**, for the reason above. The midwinter floor is a
  lower bound: net radiation approaches zero while the air keeps moving. Install
  `pyterrainwind` if that bound matters to you — in this landscape it is worth up
  to 61% of midwinter demand.
- **No weather.** A clear-sky model scaled by a clearness index has no individual
  cloudy days — these are climatological days, not a particular year.
- **No water limit.** This is demand, not actual evapotranspiration: an upper bound
  on what the atmosphere can take from a well-watered surface. What it does take is
  a soil-water question and belongs in a different package.
- **Humidity is inferred** where not supplied, as saturation vapour pressure at the
  daily minimum temperature — the FAO-56 fallback.
