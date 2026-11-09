# Climate forcing — SILO daily at the Avonside grid cell

A small real-world test dataset for the demand model: the time axis that the 3 m
terrain supplies the map for.

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

## Wind

The twenty-year BARRA-C2 wind archive for the same cell, the regional 30 m DEM for
upwind fetch, and the fitted gust and harmonic tables now live with the code that
uses them, in
[`pyterrainwind`](https://github.com/jborevitz/pyterrainwind/tree/main/data).
Nothing in this package reads them.
