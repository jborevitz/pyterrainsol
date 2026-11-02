# Wind forcing — BARRA-C2 at the Avonside grid cell

A small real-world test dataset for the wind side of the demand model.

## `barra_avonside_20y_20min.npz` (3.3 MB)

Twenty years of near-surface wind at one BARRA-C2 grid cell, 20-minute steps.

| key | meaning |
|---|---|
| `time` | POSIX seconds, **UTC**, 525,960 steps at exactly 1200 s |
| `uas`, `vas` | eastward and northward wind at 10 m, m s⁻¹ |
| `ua50m`, `va50m` | the same at 50 m, m s⁻¹ |
| `lat`, `lon` | −36.45, 148.74 — the grid cell, ~1.7 km from the farm centroid |

2006-01-01 to 2025-12-31 inclusive. 100% finite, no gaps, no duplicates.

```python
import numpy as np
d = np.load("barra_avonside_20y_20min.npz")
speed = np.hypot(d["uas"], d["vas"])                              # m/s at 10 m
drct  = (270 - np.degrees(np.arctan2(d["vas"], d["uas"]))) % 360  # deg the wind comes FROM
```

### Three things that will bite you

**Times are UTC.** Add 10 hours for AEST before any diurnal analysis, or the
afternoon wind maximum lands at dawn. AEST rather than AEDT avoids a
daylight-saving discontinuity mid-year, at the cost of the first 10 hours of
1 January belonging to the previous year's file.

**Direction is circular.** Never average degrees — average the `uas`/`vas`
components, or use circular statistics. Plot with a *cyclic* colormap. Note too
that the vector-mean speed is not the scalar-mean speed; their ratio is the
directional constancy, which here runs from 0.06 in summer to 0.43 in winter.

**It is a 4.4 km grid cell, not your paddock.** It resolves valley-scale
channelling but knows nothing about a ridge or a gully. Terrain redistribution
at farm scale is a separate problem.

### Having two heights is the useful part

10 m and 50 m together invert the neutral log profile for an effective
roughness length, `z0 = exp[(ln z2 − r·ln z1)/(1 − r)]` with `r = u50/u10`.
Here the median ratio is 1.515, giving **z0 ≈ 0.44 m** — more than ten times
rougher than short grass. Restricting to near-neutral windy midday steps gives
1.498 and 0.40 m.

Fit that only on near-neutral steps. The ratio is 1.473 at midday and 1.573 at
night across all speeds, which would invert to 0.33 m and 0.60 m — the ground
does not get rougher after dark, that is the stable nocturnal boundary layer
decoupling the surface from the air above.

Two caveats worth stating in anything published. These are BARRA's own
diagnostic levels, so the inversion recovers *the model's* effective roughness,
a property of its land-surface configuration rather than an independent
measurement. And a 4.4 km cell aggregates vegetation with sub-grid orography,
so the number cannot be attributed to trees.

## `silo_avonside_20y_daily.npz` (0.19 MB)

The same twenty years as the wind, daily, from SILO's 0.05 degree (~5 km) grid:
the atmospheric forcing the wind data needs to be useful.

| key | meaning |
|---|---|
| `date` | days since 1970-01-01, 7,305 of them, 2006-01-01 to 2025-12-31 |
| `daily_rain`, `max_temp`, `min_temp`, `radiation` | mm, C, C, MJ m-2 d-1 |
| `vp`, `vp_deficit` | vapour pressure and **vapour pressure deficit**, hPa |
| `rh_tmax`, `rh_tmin` | relative humidity at the temperature extremes, % |
| `et_short_crop`, `et_tall_crop` | FAO-56 and ASCE reference ET, mm |
| `et_morton_actual`, `et_morton_potential`, `et_morton_wet` | Morton suite, mm |
| `evap_pan`, `evap_syn`, `evap_comb`, `evap_morton_lake` | evaporation estimates, mm |
| `mslp` | mean sea-level pressure, hPa |
| `lat`, `lon`, `elevation_m` | -36.4355, 148.7257, 861.1 m |

All eighteen variables are 100% finite over the full record.

**`vp_deficit` is the variable that unlocks Penman-Monteith.** The aerodynamic
term goes as `u x VPD`, so without humidity a wind speed multiplies an unknown
and contributes nothing - which is why the wind archive alongside this one
could not move an ET number on its own.

### What it says about this place

Rainfall averages **537 mm/yr** over the twenty years, ranging 320 to 784.
Reference ET for a short crop averages **1043 mm/yr**, so the mean annual
deficit is around **506 mm**. That is a drier site than the montane setting
suggests, and it reframes the terrain result in the main README: a 472 mm/yr
spread in demand between north- and south-facing slopes is not merely "most of
a year's rainfall", it is **88% of it**.

VPD runs 3.4 hPa in June to 15.4 in January, annual mean 0.84 kPa.

### Two cautions

**The cell is 861 m; the farm is 888-1197 m**, mean about 1040. SILO's grid
point sits roughly 180 m below the farm's mean elevation, so its temperatures
run warm and its ET high for this block. Lapse-correct before using it as a
farm-scale forcing.

**It is one 5 km cell**, so like the wind it is a time axis, not a map. Every
bit of spatial structure still comes from the 3 m terrain.

As a cross-check on `pyterrainsol.demand`: our independently computed
Priestley-Taylor estimate for a flat unobstructed cell at 1040 m is 953 mm/yr,
which sits between SILO's `evap_comb` at 972 and `et_short_crop` at 1043 -
close enough to suggest the radiation chain is sound, and lower in the way a
higher, cooler cell should be.

`fetch_silo_point.py` regenerates this file. SILO is free but requires a
registered email as the `username` parameter; register your own at
<https://www.longpaddock.qld.gov.au/silo/>. Data is supplied by the Queensland
Department of Environment, Science and Innovation; see that page for licence
and citation. Nothing here is endorsed by them.

## Supporting tables

- `barra_wind_2025_20min.csv` — 2025 only, 26,250 rows, local AEST timestamps
  with speed, direction and the raw components. Convenient for plotting.
- `wind_weibull_2025.csv` — Weibull shape `k` and scale `c` by season and
  6-hour block, with mean speed and the modelled fraction of hours above
  30 km/h. Fitted on 2025 alone; the tail behaviour needs the full 20 years
  before anyone leans on it.

## Provenance and attribution

Derived from the **Bureau of Meteorology Atmospheric high-resolution Regional
Reanalysis for Australia, version 2 (BARRA2)**, model BARRA-C2 on the AUST-04
grid, obtained from the National Computing Infrastructure (NCI) THREDDS server.

- Dataset DOI and licence terms: <https://doi.org/10.25914/1x6g-2v48>
- Data governance contact: `data_governance_office@bom.gov.au`
- NCI project: `ob53`

This file is a **single-point extract** redistributed for convenience and
testing. It is not a Bureau of Meteorology product and has not been endorsed by
the Bureau. Anyone using it in publication should go to the DOI above for the
licence and citation, and should consider pulling their own extract with
`fetch_barra_point.py` — it takes about half an hour.
