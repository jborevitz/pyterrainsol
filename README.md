# pyterrainsol

**How much light reaches each patch of ground, and how much water the air pulls
back out of it.**

A plant does not experience elevation. It experiences the light that falls on the
particular slope it is rooted in, and the drying power of the air around it. On flat
ground those are much the same everywhere, and one weather station describes a whole
paddock. On a hillside they are not. Two shrubs fifty metres apart on opposite sides
of a spur can be living in different climates, and the difference is large enough to
decide which one survives a dry summer.

This package takes a bare-earth surface — slope, aspect, how much of the sky each
cell can see, and what the surrounding ridges block — and returns, for every cell and
every day of the year:

- **incoming sunlight on the real surface**, inclined and shaded, rather than on an
  imaginary flat plane;
- **photosynthetically active radiation and the daily light integral**, the quantity
  a growth model or a grower actually wants;
- **evaporative demand** by Priestley-Taylor — how much water the atmosphere can take
  from a well-watered surface, which is the ceiling on transpiration.

Wind is deliberately absent, because it is a separate decision rather than a missing
feature: aspect-driven light is predictable geometry and matters at any site with
relief, while wind needs a terrain-to-speed transfer function that is scheme-
dependent and only partly constrained. If that matters to you, add
[`pyterrainwind`](https://github.com/jborevitz/pyterrainwind), which builds on this
package. If it does not, you never have to think about it.

## Where this sits

Six small packages, one chain. Each does one thing and each is useful on its own.

| package | reads | gives you |
|---|---|---|
| [pylidarterrain](https://github.com/jborevitz/pylidarterrain) | classified airborne LiDAR | the ground surface a plant actually stands on: slope, aspect, how much sky it can see |
| **pyterrainsol** (this one) | that surface + daily weather | **light and evaporative demand, cell by cell, day by day** |
| [pyterrainwind](https://github.com/jborevitz/pyterrainwind) | that surface + reanalysis wind | wind exposure, and the extra demand it carries |
| [pychannelform](https://github.com/jborevitz/pychannelform) | that surface | where the water goes: catchments, channels, long profiles |
| [pylandstore](https://github.com/jborevitz/pylandstore) | channels + demand | where to hold water, and what it costs to move it uphill |
| [pylandgraze](https://github.com/jborevitz/pylandgraze) | all of the above + a herd | where the animals should be, hour by hour and season by season |

## Quickstart

```python
from pyterrainsol import (daily_radiation, par_from_components, daily_light_integral,
                          net_longwave, net_radiation, priestley_taylor)

beam, diffuse = daily_radiation(slope_deg, aspect_deg, horizon, sky_view,
                                lat_deg=-36.4355)          # (365, n_cells) MJ/m2/d
dli = daily_light_integral(par_from_components(beam, diffuse))   # mol photons/m2/d

loss = net_longwave(t_min_c, t_max_c, clearness)           # (365,) -- not per cell
et   = priestley_taylor(net_radiation(clearness * (beam + diffuse), loss),
                        t_mean_c[:, None], elevation_m)    # (365, n_cells) mm/d
```

Slope, aspect, horizon angles and sky view come from
[`pylidarterrain`](https://github.com/jborevitz/pylidarterrain).

## What it found at 845 Avonside

46.3 ha in the Snowy Monaro, NSW, 888-1197 m, on a 3 m LiDAR surface with horizon
shading computed at 32 compass directions. Rainfall averages **537 mm/yr** over
2006-2025 against **1043 mm/yr** of reference evapotranspiration — the air asks for
roughly twice what the sky delivers, every year. Read every number below against that
deficit.

### Terrain sorts the light in winter and barely at all in summer

![Annual and winter insolation at 3 m](docs/insolation_annual_winter_3m.png)

*Annual total sunlight spans 3477-8890 MJ m⁻² yr⁻¹ across the one block, a 2.6-fold
range. On the winter solstice the daily light integral spans 4.4-31.8 mol photons
m⁻² d⁻¹ — **7.3-fold**. On the summer solstice the same ground spans only 1.48-fold,
47.5-70.1 mol m⁻² d⁻¹.*

This is the single most useful thing the radiation model says about managing a
hillside. **In midsummer every part of the block is within about 50% of the same
light, so summer differences between paddocks are differences in water, not aspect.**
In winter the same ground spans seven-fold, which is what sets where frost lingers,
when growth restarts in spring, and which faces can carry anything winter-active at
all.

![PAR and its share of shortwave](docs/par_and_fraction.png)

*The photosynthetic share of sunlight spans only 0.450-0.481. Shaded ground earns
about **7% more photosynthetic light per unit of evaporative load**, because diffuse
sky light is spectrally richer in PAR than direct beam — a real effect, and a small
one. The two maps remain roughly 93% the same picture.*

### Demand swings much harder than light

The fraction of available energy that leaves as water vapour rather than as warm air
runs **0.483 in July to 0.680 in January**, a 1.41-fold swing. Light swings 3.7-fold.
Demand swings harder than either, and the reason is a term that is easy to forget:
the ground radiates heat away to the sky at 4.4-4.8 MJ m⁻² d⁻¹ almost regardless of
season. **In June only 0.87 of 7.2 MJ m⁻² d⁻¹ survives that tax as net radiation.**
Midwinter is close to radiatively shut, which is why what little demand remains then
is driven by moving air rather than by sun — the subject of the companion package.

![A year of light](docs/avonside_light_year.gif)

*Incoming sunlight, every day, on the 3 m grid. The boundary is drawn over the whole
catchment rather than clipped to it. This GIF steps every 6th day to stay under 6 MB;
the daily version is
[avonside_light_year.mp4](https://raw.githubusercontent.com/jborevitz/pyterrainsol/main/docs/avonside_light_year.mp4).*

### The block divides into four parts, and 365 mm separates them

![Annual potential evapotranspiration at 3 m](docs/annual_et_3m.png)

*Demand spans **318-1103 mm/yr** within the one property. The histogram is
**bimodal** rather than a smooth gradient: aspect sorts the block into two
populations. On slopes steeper than 20°, north-facing ground demands **1071 mm/yr
against 599 for south-facing — a 472 mm gap, 88% of a year's rainfall, inside one
fence line.***

Recomputing on daily SILO weather, with per-cell temperature lapsed from the grid
point, shifts the range to **310-1075 mm/yr, median 862**. By the four terrain parts
(defined in [`pylidarterrain`](https://github.com/jborevitz/pylidarterrain)):

| part | area | median grade | % north-facing | demand, radiation only |
|---|---|---|---|---|
| Top farm | 11.6 ha | 14% | 22% | 899 mm/yr |
| Upper forest | 10.2 ha | 48% | 1.5% | **653 mm/yr** |
| Northerly woodland | 11.8 ha | 30% | 85% | **1018 mm/yr** |
| Lower forest flats | 12.7 ha | 29% | 5.5% | 809 mm/yr |

**365 mm/yr separates the thirstiest part from the least thirsty — 68% of annual
rainfall — and the two extremes share a boundary.** For a land manager that is the
whole argument for mapping rather than averaging: these are not gradations of one
site, they are four different places, and a plan written for their mean fits none of
them. The figures above are radiation only; adding wind widens the spread further
rather than evening it out, which is shown in
[`pyterrainwind`](https://github.com/jborevitz/pyterrainwind).

A flat unobstructed cell at 1040 m comes out at 953 mm/yr. SILO's own products at
this location bracket it — `evap_comb` 972, `et_short_crop` 1043 — and ours sits
lower in the way a higher, cooler cell should.

![A year of evaporative demand](docs/avonside_et_year.gif)

*The same days as demand. This one uses a square-root colour ramp, which has since
been retired for flattening exactly the summer days that move the water; the current
linear 0-10 mm/day movies live in
[`pyterrainwind`](https://github.com/jborevitz/pyterrainwind). Daily version:
[avonside_et_year.mp4](https://raw.githubusercontent.com/jborevitz/pyterrainsol/main/docs/avonside_et_year.mp4).*

## Data

[`data/`](data/) holds twenty years of daily SILO climate for this site, with the
script that regenerates it and the traps worth knowing — including the one where
SILO's published vapour pressure deficit is not the FAO-56 quantity of the same name.
See [`data/README.md`](data/README.md).

## Methods

The conventions this package enforces, the traps it exists to prevent, and what it
deliberately does not do are in **[`docs/methods.md`](docs/methods.md)**.

Methods and results are kept apart on purpose. The methods travel to any site and
will be re-checked against new ground and simulated data; the findings above are one
46 ha property in the Snowy Monaro, some of them specific to it and some likely to
generalise. Either can then be argued with without disturbing the other.
