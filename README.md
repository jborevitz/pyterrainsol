# pyterrainsol

Terrain-corrected solar radiation: beam, diffuse, PAR and daily light integral.

```python
from pyterrainsol import (daily_radiation, par_from_components, daily_light_integral,
                          net_longwave, net_radiation, priestley_taylor)

beam, diffuse = daily_radiation(slope_deg, aspect_deg, horizon, sky_view,
                                lat_deg=-36.4355)          # (365, n_cells) MJ/m2/d
par = par_from_components(beam, diffuse)
dli = daily_light_integral(par)                            # mol photons /m2/d

loss = net_longwave(t_min_c, t_max_c, clearness)           # (365,)  -- not per cell
et   = priestley_taylor(net_radiation(clearness * (beam + diffuse), loss),
                        t_mean_c[:, None], elevation_m)    # (365, n_cells) mm/d
```

Takes slope, aspect, horizon angles and sky view factor — e.g. from
[`pylidarterrain`](https://github.com/jborevitz/pylidarterrain) — and returns daily
shortwave on the inclined, shaded surface.

## Why not a temperature-based proxy

Hargreaves and similar run on extraterrestrial radiation and temperature range, so
they are structurally blind to aspect: a north slope and a south slope at the same
elevation return identical values. If aspect is meant to drive evaporative demand,
the radiation has to be computed on the actual inclined, shaded surface.

## PAR and shortwave are nearly the same map

A single fixed PAR fraction makes the photosynthesis and evaporation maps identical
in different units. They separate only because diffuse light is spectrally richer in
PAR (~0.60) than direct beam (~0.43), so the split is applied per component. The
resulting spread is real but modest. Say so rather than overselling it.

## Worked example — 47 ha in the Snowy Monaro, NSW

Computed on a **3 m** LiDAR DEM (51,467 cells) with horizon shading at 32 azimuths.

![Annual and winter insolation at 3 m](docs/insolation_annual_winter_3m.png)

*Annual total shortwave (3477–8890 MJ m⁻² yr⁻¹, a 2.6-fold span) and the winter
solstice daily light integral (4.4–31.8 mol photons m⁻² d⁻¹, 7.3-fold). On the summer
solstice the same ground spans only 1.48-fold, 47.5–70.1 mol m⁻² d⁻¹. Terrain sorts
the energy load overwhelmingly in winter — so summer differences between paddocks come
from water, not from aspect.*

![PAR and its share of shortwave](docs/par_and_fraction.png)

*PAR as a fraction of total shortwave spans only 0.450–0.481 across the block.
Shaded ground earns about 7% more photosynthetic light per unit of evaporative load,
but the two maps remain roughly 93% the same picture.*

## From light to evaporative demand

Radiation is a map; temperature is a clock. The Priestley-Taylor partition term
Δ/(Δ+γ) is the fraction of available energy that leaves the surface as water
rather than as hot air, and because saturation vapour pressure follows
Clausius-Clapeyron it climbs steeply with temperature. At Avonside it runs from
**0.483 in July to 0.680 in January** — a 1.41-fold swing. On any given day every
cell on the block shares that value to within the 2 °C the 309 m elevation range
buys, so it varies in time and barely in space, which is the mirror image of
radiation. That is why a 10 km climate product costs nothing here: it only ever
supplies the part that does not vary across the farm.

Altitude helps the wrong way. At 1040 m the pressure is 89.6 kPa, so γ is 0.0596
rather than the sea-level 0.0674, raising the evaporating fraction of any given
energy input by about three percentage points year-round.

### Two years in the life of one hillside

![A year of light](docs/avonside_light_year.gif)

*Incoming shortwave, every day of the year, 3 m grid. Clear-sky geometry scaled
by a monthly clearness index of 0.74–0.80 (annual mean 0.760) derived from
WorldClim's observed radiation climatology. The property boundary is drawn over
the full catchment rather than clipped to it. **The animated GIF steps every 6th day**,
which keeps it under 6 MB but makes the dates jump about a week at a time; the daily version is
[avonside_light_year.mp4](https://raw.githubusercontent.com/jborevitz/pyterrainsol/main/docs/avonside_light_year.mp4) — GitHub cannot play
video inline, so that link opens the file directly and will either play in the browser
or download, depending on the browser.*

![A year of evaporative demand](docs/avonside_et_year.gif)

*The same days as potential evapotranspiration in mm, on the same ramp as the annual
map below — **blue is where water stays, pale is where it leaves**. January bleaches
almost white as the whole block gives up water; July holds deep blue, with the
north-facing slopes showing as the pale patches that are still losing some. Note the **square-root colour
scale**: midwinter demand is 0.1–1 mm/day against a summer 5–6, and on a linear
scale the whole of winter is a blank screen. As above the GIF steps every 6th day; the daily
version is [avonside_et_year.mp4](https://raw.githubusercontent.com/jborevitz/pyterrainsol/main/docs/avonside_et_year.mp4).*

![Annual potential evapotranspiration at 3 m](docs/annual_et_3m.png)

*The year summed into one map, on the same ramp as the movie. Potential ET spans **318-1103 mm/yr** across 46.3 ha,
median 900, interquartile range 758-1012. The on-farm histogram is **bimodal**, not a
smooth gradient: aspect sorts the block into two populations rather than a continuum.
On slopes steeper than 20 degrees, north-facing ground demands **1071 mm/yr against
599 for south-facing - a 472 mm gap**, which at this site is most of a year's rainfall.
Two paddocks the same size, the same soil and the same sky are being asked for very
different amounts of water.*

The pair is worth watching together because they are not the same movie. Farm-median
light swings **3.7-fold** between summer and winter; demand swings far harder, because
the longwave loss of 4.4–4.8 MJ m⁻² d⁻¹ is a near-constant subtraction. In June only
**0.87 of 7.2 MJ m⁻² d⁻¹** survives as net radiation. Multiply what is left by a
partition term that has itself fallen by a third and midwinter demand nearly vanishes.

Over the block, potential ET runs **318–1103 mm/yr, median 900**; a flat unobstructed
cell at 1040 m gives 953 mm/yr. Watch July in the second movie: the north-facing
slopes are the only ground losing appreciable water, which is the same effect that
gives those faces an early spring.

## What this does not do

- **No wind.** Priestley-Taylor has no aerodynamic term; α = 1.26 absorbs a typical
  advection contribution. That is defensible when net radiation dominates and weakest
  exactly where this model puts its minimum — midwinter, when Rn is near zero but the
  air still moves. Treat the winter floor as a lower bound. Wind enters properly only
  through Penman-Monteith, which needs humidity and wind speed this package does not
  take.
- **No weather.** A clear-sky model scaled by a monthly clearness index has no
  individual cloudy days, so these are climatological days, not a particular year.
- **No water limit.** This is demand, not actual evapotranspiration — an upper bound
  on what the atmosphere can take from a well-watered surface. What it does take is a
  soil water question and belongs in a different package.
- **Humidity is inferred**, as saturation vapour pressure at the daily minimum
  temperature, the FAO-56 fallback when nothing is measured.

Climate inputs are WorldClim 2.1 monthly normals at 5 arc-minutes — roughly 10 km, a
single cell over this farm, and a 1970–2000 baseline.

## What the wind actually does here

The model above has no wind term, so it is worth being concrete about the thing
it is missing. `data/` holds twenty years of BARRA-C2 wind at the nearest grid
cell, 10 m and 50 m, at 20-minute steps.

![A year of wind, speed and direction](docs/wind_year_20min.png)

*2025 at 20-minute resolution: 365 days across, 72 steps of the day up, in local
AEST. Speed above; below, the direction the wind comes from on a **cyclic**
colormap so that 359° and 1° are the same colour and north carries no seam, with
a compass ring rather than a bar for the same reason. Direction opacity scales
with speed — 0.18 at 4 km/h rising to opaque at 20 km/h — because direction is
meaningless when nothing is moving, but a hard cut erases real light breezes.*

Two regimes are visible and they matter in opposite directions for water. The
dark columns through spring are a dry, fast north-westerly; of the steps above
36 km/h, 89% come from WNW through NNW. The pale pre-noon wedges from December
to March are a slow, humid easterly that rotates to westerly by afternoon — a
daily rotation absent in June and July. The first raises evaporative demand
above what this package reports, the second lowers it. A model with no
aerodynamic term sees neither, so its error changes sign with the season.

### Fitting it, and what survives twenty years

![Twenty years of diurnal cycle and gust shape](docs/wind_20y_structure.png)

A tensor-product harmonic model — annual and diurnal sine and cosine pairs plus
their interaction, nine terms — explains **10.0% of the variance in speed over
2006–2025**. Fronts, not the sun, set the wind here: the fitted climatology has
a standard deviation of 2.01 km/h against a residual of 6.02, so synoptic
weather is three times the amplitude of the entire annual plus daily cycle.
(Total spread is 6.34 km/h, and 2.01² + 6.02² recovers it.)

Fitting single years first was a mistake worth recording. **2025 alone gave
15.3%, and it ranks first of twenty on both goodness of fit and annual
amplitude** — the most seasonally organised year in the record. Anything
calibrated on it would have been calibrated on an outlier.

What is reproducible, and what is not:

| quantity | pooled | year-to-year sd | verdict |
|---|---|---|---|
| mean speed | 10.56 km/h | 0.36 | solid |
| diurnal amplitude | 2.33 km/h | 0.19 | solid |
| diurnal peak hour | 15.7 h AEST | 0.26 h | solid |
| diurnal amplitude, midsummer | 3.86 km/h | 0.40 | solid |
| annual amplitude | 1.10 km/h | **0.40** | weak |
| annual peak day | day 275 | **28 days** | do not interpret |

The diurnal cycle is the predictable part — the left panel is twenty years of
fitted daily cycles, and they lie almost on top of each other within a season
while the seasonal amplitude contrast is large and consistent. The annual cycle
in mean speed is not: its amplitude is the same size as its own year-to-year
spread, and the fitted peak day ranges from day 216 to day 341 across the
record. Pooling *lowers* the annual amplitude relative to the mean of per-year
fits, which only happens when the phase wanders.

### Gust shape is the robust part, and the mean misranks it

Weibull fits by season and six-hour block hold up far better than the harmonics.
Shape `k` per bin varies year to year by only 0.07–0.24, and the structure in the
right panel is the same on twenty years as on one: **mornings carry the lowest
`k` in all four seasons** (1.36–1.63), meaning the most episodic, fattest-tailed
wind. Afternoons are steadier, because convective mixing puts a floor under them.

This is why wind has to enter as a distribution. Summer afternoons have a
*higher* scale than spring mornings — 15.7 against 12.0 km/h — and yet **six
times fewer hours above 30 km/h**, 0.44% against 2.62%, because their `k` is
2.61 against 1.41. Rank the two by mean wind speed and you get the damage risk
backwards.

One honest correction from going to twenty years: 2025's tails were about twice
the climatology. Spring mornings spent 5.0% of hours above 30 km/h in 2025
against 2.62% over the full record. The *shape* was right; the *rate* was not.

Parameters are in `data/wind_weibull_20y.csv` and `data/wind_harmonic_by_year.csv`.
