# Phase 1 notes: dataset construction

These notes record how the Phase 1 dataset was built, which rules were applied, and what each rule removed. Every figure is copied from a file written by the code; the source file is named next to it. Summary files are in `data/interim/summaries/`. The notes describe processing only and do not interpret differences between periods.

## Run order and outputs

The modules run in this order, each from the repository root with `PYTHONPATH=src`:

1. `grid` builds the 500 m grid, the LGA layer and the WorldCover 2021 masks.
2. `population` adds WorldPop totals.
3. `buildings` adds Open Buildings and GHSL built-up shares.
4. `temperature` builds the Landsat land surface temperature composites.
5. `vegetation` downloads the per-date Sentinel-2 data.
6. `s2_haze` runs the haze screen.
7. `vegetation` runs again to build the composites: the 80th-percentile versions, the screened median without the excluded dates, and the unscreened median. All dates are stored by then, so nothing is downloaded. Its last step (`finalize`, also runnable alone with `--finalize`) reads the written composites and computes the vegetated and bare shares, the sand-filling flags and the reference-forest values.
8. `vegetation_check` runs the NDVI diagnostic.
9. `covariates` adds water distances, roads and industrial land.
10. `assemble` builds the final tables and maps.

Row counts of the tables in `data/interim/`:

| File | Rows | Unit |
|---|---|---|
| grid_periods.csv | 14,890 | cell (98 columns, `assemble.json`) |
| lga_periods.csv | 20 | LGA |
| buildings_cells.csv, vegetation_cells.csv, temperature_cells.csv, population_cells.csv, covariates_cells.csv | 14,890 each | cell |
| temperature_cells_landsat8_only.csv | 14,890 | cell (copy of the Landsat 8-only run, kept for comparison) |
| population_lga.csv | 20 | LGA |
| ob_tiles.csv | 116 | Open Buildings tile and year |
| s2_items.csv, s2_offsets.csv | 191 each | Sentinel-2 scene |
| s2_ndvi_check.csv, s2_haze_screen.csv | 50 each | Sentinel-2 date |
| landsat_scenes.csv | 84 | Landsat date key, plus scenes dropped for cloud cover |

## Grid and water

The grid has 500 m cells in EPSG:32631. 14,890 cells intersect Lagos State, whose area is 3,550.8 km², split into 20 LGAs (`grid.json`). A cell is flagged as water when more than half of it is WorldCover 2021 water (class 80). This rule flags 2,452 cells, leaving 12,438 land cells (`grid.json`). Water-flagged cells stay in every table but have no vegetation values and are left out of the LGA means.

For the distance covariates, water is split into ocean, lagoon and other water on a 50 m WorldCover grid:

- Ocean is water connected to the southern edge after an erosion of 1,500 m.
- Lagoon is any remaining water body of at least 5 km².

Within the extended covariate extent this gives 3,965.6 km² of ocean, 704.5 km² of lagoon in 2 bodies, and 29.5 km² of other water (`covariates.json`). The harbour channel joins the lagoon body, so distance to lagoon also measures distance to the harbour.

## Buildings

Building share and height come from Open Buildings 2.5D Temporal for 2016 and 2023. These are 58 tiles per year, read at 4 m (`buildings.json`). Neighbouring tiles overlap and carry identical values in the overlap, so for each cell the single tile with the most valid pixels is used. In 54 cells no single tile covers the whole cell; their tile values are summed and the cells are flagged `ob_overlap_approx_*`. 37 of these 54 cells are land cells (`buildings.json`, `assemble.json`). Open Buildings coverage is complete in all cells for both years, with a minimum of 1.0 (`buildings.json`).

GHSL built-up surface for 2015 and 2020 is kept as a cross-check. Pearson correlations over the 12,438 land cells are 0.957 (Open Buildings 2016 against GHSL 2015) and 0.958 (Open Buildings 2023 against GHSL 2020) (`buildings.json`). GHSL 2025 is a model projection and is not used.

## Sentinel-2 vegetation

### Composites and offset

Sentinel-2 L2A scenes are read at 20 m.

- **Cloud and pixel masks:** scenes with 60% cloud cover or more are dropped. Pixels are kept only where the scene classification layer (SCL) is class 2, 4, 5 or 7.
- **Per-date data:** NDVI uses bands B08 and B04. The bare soil index (BSI) uses B02, B04, B08 and B11. Each date's NDVI and BSI are stored, and per-pixel medians are taken over the dates kept in each period.
- **Relative NDVI:** each date's NDVI minus that date's median over valid Lagos land pixels. The per-pixel median of this value is then taken over the period, which is the same construction as relative LST.

The reflectance offset is decided per scene, not by acquisition date. Each scene's `BOA_ADD_OFFSET` values are read from its product metadata (`MTD_MSIL2A.xml`), and the processing baseline gives a second, independent value. All 118 late-period scenes carry an offset of −1000 in their metadata. None of the 73 early-period scenes carries an offset (`vegetation.json`, `s2_offsets.csv`). The two sources agree for every scene (0 disagreements, `vegetation.json`). The early scenes kept after the cloud filter have baselines 02.12 (49) and 03.00 (3). The late ones have 05.09 (7), 05.10 (45) and 05.11 (39) (`vegetation.json`). No early scene on Planetary Computer has a reprocessed baseline that would carry the offset.

For the 10 m bands (B02, B04, B08), odc-stac reads each file's internal 2× overview rather than averaging the 10 m pixels; B11 is native 20 m. A check against direct reads of the same overview matched pixel for pixel. Against a plain 2×2 mean of the 10 m pixels, the overview differs per pixel but with no systematic shift, and NDVI medians are unchanged. That check was run once by hand and its numbers were not written to a file, so they are not reported here.

The rasterio warning "Dataset has no geotransform, gcps, or rpcs" appears once near the start of some runs. It is raised inside `rasterio.warp.reproject`, which odc-stac calls whenever a band is resampled, and the target grid is always passed explicitly in that call. Several checks found no effect on stored data:

- Re-loading a stored date gave identical NDVI.
- SCL and band loads matched direct reads.
- No stored scene shows the signature of an empty read.

The highest shares of pixels masked for non-positive reflectance (`s2_items.csv`) come from negative blue reflectance (B02) in heavily aerosol-laden scenes, not from failed reads.

### NDVI diagnostic

The median land NDVI of a single date ranges from 0.167 to 0.590 in the early period and from 0.196 to 0.621 in the late period (`s2_ndvi_check.csv`). Surfaces that should be stable show the same date-to-date range:

- **Dense built core:** median NDVI ranges from −0.021 to 0.107 early and from −0.057 to 0.098 late.
- **Lekki Conservation Centre trees:** median NDVI ranges from 0.073 to 0.790 early and from 0.313 to 0.764 late.

The clearest December dates in both periods give near-identical values (`s2_ndvi_check.csv`):

| Date | Land | Built core (GHSL ≥ 0.6) | Lekki CC trees | Unnamed reserve trees (OSM way 797720562) |
|---|---|---|---|---|
| 2016-12-12 | 0.590 | 0.105 | 0.680 | 0.711 |
| 2023-12-16 | 0.584 | 0.090 | 0.694 | 0.699 |
| 2023-12-21 | 0.594 | 0.098 | 0.689 | 0.712 |
| 2024-12-10 | 0.570 | 0.085 | 0.683 | 0.708 |

Medians over all processed dates, early versus late (`vegetation_check.json`):

| Reference | Early | Late |
|---|---|---|
| Land | 0.3524 | 0.3987 |
| Built core, GHSL ≥ 0.6 | 0.0261 | 0.0384 |
| Omu creek trees | 0.396 | 0.5509 |
| Lekki Conservation Centre trees | 0.3546 | 0.5689 |
| Unnamed reserve trees | 0.4934 | 0.5874 |

The difference between the period-level medians therefore depends on which dates fall in each period and on the atmosphere on those dates; the offset is applied correctly. The analysis therefore uses relative NDVI alongside absolute NDVI. Relative NDVI removes a shift shared by all of Lagos on a date. It does not remove effects that vary across the city on the same date, such as patchy haze.

Deviation from the requested reference: the requested built-core rule, GHSL built-up share above 0.9 in both 2015 and 2020, selects 0 cells, because the highest cell-level GHSL share is 0.759 (`vegetation_check.json`, `buildings.json`). Thresholds of 0.6 (228 cells) and 0.7 (32 cells) are used instead. The tree references use 20 m pixels where all four WorldCover 2021 10 m sub-pixels are tree cover, inside OSM protected areas. These are 1,407 pixels in Lekki Conservation Centre, 618 in Omu creek and 40,335 in the unnamed reserve (`vegetation_check.json`). Lekki Conservation Centre is not tagged as a nature reserve in OSM and was added by a name query.

### Primary vegetation measure

The following was decided before any cell-level results were inspected. The primary vegetation measure is the per-pixel 80th percentile of NDVI over all dates of a period (columns `*_p80`). Haze biases NDVI downward, so the 80th percentile favours each pixel's clearer observations while keeping all dates for coverage. A relative version (`ndvi_rel_*_p80`) takes each date's NDVI minus that date's median over Lagos land pixels, then the 80th percentile per pixel.

Two sensitivity versions are kept:

- The screened median: the per-pixel median over the dates kept by the reference-forest screen below (3 early, 9 late), in the unsuffixed columns.
- The unscreened median: the per-pixel median over all dates (`*_noscreen`).

A second decision followed the reference-forest values of the composites (see Observation counts), still without inspecting cell-level results. The absolute 80th-percentile NDVI of the reference forest is 0.5933 early and 0.6785 late, while the screened median gives 0.7099 and 0.6900 (`vegetation.json`). The absolute 80th percentile therefore still carries a haze difference between periods. As a result:

1. The primary measure of vegetation change is the change in relative 80th-percentile NDVI (`ndvi_rel_change_p80`, late minus early). Absolute 80th-percentile NDVI (`ndvi_*_p80`) is used only for levels within a period, and its change is not written.
2. The sensitivity measures of change are the screened median change (`ndvi_change`), restricted to the 746 cells with at least 3 screened observations in both periods (`ndvi_screened_ok`), and the change in the unscreened relative median (`ndvi_rel_change_noscreen`).
3. The vegetated share (`veg_share_*`, share of 20 m pixels with NDVI above 0.4) is computed from the 80th-percentile NDVI within each period. The bare share (`bare_share_*`, share of pixels with BSI above 0) and the sand-filling flag come from the unscreened median BSI. The versions from the screened medians are kept as `*_screened` columns.

### Haze screen

The screen judges haze from NDVI on one stable reference forest, never from the analysed cells. The rule was set before its effect was inspected:

- **Reference forest:** OSM way 797720562, a protected forest of 1,649 ha. It contains 40,335 20 m pixels where all four WorldCover 2021 10 m sub-pixels are tree cover (`s2_haze.json`).
- **Per-date value:** the median NDVI over the reference-forest pixels that are valid that date.
- **Reference level:** the median of the 3 highest date medians, pooled over both periods. These are 2016-12-02 (0.7195), 2023-12-21 (0.7125) and 2016-12-12 (0.7109), giving a level of 0.7125 and a threshold of 0.6125 (`s2_haze.json`).
- **Exclusion:** a date is excluded when its median is more than 0.10 below the reference level. A date with fewer than 200 valid reference pixels (8 ha) cannot be assessed and is also excluded. That second part was decided after the first run of the screen had kept such dates.

Result (`s2_haze.json`, per-date values in `s2_haze_screen.csv`):

| Period | Dates | Not assessable | Below threshold | Kept |
|---|---|---|---|---|
| Early | 19 | 2 | 14 | 3 (2016-12-02, 2016-12-12, 2018-01-01) |
| Late | 31 | 3 | 19 | 9 (2023-12-06, 2023-12-11, 2023-12-16, 2023-12-21, 2023-12-26, 2024-01-25, 2024-12-10, 2024-12-15, 2024-12-20) |

### Failed attempt: AOT screen

A first screen used the Sen2Cor aerosol optical thickness at 550 nm (`AOT` asset, DN / 1000), with the median taken per date over clear Lagos land pixels on a 60 m grid. It would have excluded a date whose median was more than 2 sample SDs above its period's mean, one-sided, in a single pass. It flagged 2018-01-01 early (threshold 0.517) and 2024-01-25, 2024-02-04 and 2024-02-09 late (threshold 1.048) (`s2_haze.json`, `aot_attempt`). It is no longer applied, because it does not measure haze in the early period:

- **Fallback values:** 11 of the 19 early dates have the same date median, 0.204, and on 7 early dates the 90th percentile equals the median. This pattern points to a constant fallback in the older Sen2Cor processing, not a retrieval. The late period also has repeated values: 2 dates share a median of 1.188.
- **Correlation sign:** across dates, AOT correlates negatively with land NDVI in the late period, as haze would produce (Spearman −0.398). In the early period the correlation is positive (0.539).
- **Effect:** the single early date it flagged, 2018-01-01, has a reference-forest NDVI of 0.709 (`s2_haze_screen.csv`). That date passes the reference-forest screen.

The AOT values and the dates that screen would have excluded remain in `s2_haze_screen.csv` (`excluded_aot_attempt`).

### Observation counts

A cell's observation count is the per-pixel number of valid dates, averaged over the cell's valid pixels. The primary version and the unscreened median use the same dates, so they share one count (`s2_nobs_*_noscreen`). A land cell is flagged low-observation when this count is below 3 (`grid_periods.csv`, `assemble.json`):

| Version | Dates early / late | Land cells below 3 obs, early | Land cells below 3 obs, late |
|---|---|---|---|
| 80th percentile (primary) and unscreened median | 19 / 31 | 104 | 53 |
| Screened median | 3 / 9 | 11,692 | 105 |

Over land cells, the all-date count ranges from 1.02 to 17.27 early (median 11.29) and from 1.12 to 28.0 late (median 17.61). The screened count ranges from 1.0 to 3.0 early (median 2.23) and from 1.0 to 9.0 late (median 6.42) (computed from `grid_periods.csv`). With only 3 early dates, most screened early cells have at least one pixel without all 3 observations.

Reference-forest NDVI of each composite (median over the 40,335 reference pixels, `vegetation.json`):

| Version | Early | Late | Late minus early |
|---|---|---|---|
| Relative 80th percentile (`ndvi_rel_*_p80`, primary for change) | 0.1838 | 0.2089 | 0.0251 |
| Absolute 80th percentile (`ndvi_*_p80`, levels only) | 0.5933 | 0.6785 | 0.0852 |
| Screened median | 0.7099 | 0.6900 | −0.0199 |
| Screened relative median | 0.1990 | 0.1418 | −0.0572 |
| Unscreened median (`_noscreen`) | 0.4660 | 0.5752 | 0.1092 |
| Unscreened relative median (`ndvi_rel_*_noscreen`) | 0.1367 | 0.1368 | 0.0001 |

The differences in the last column are computed from the two values in `vegetation.json`.

Land cells, out of 12,438, with an NDVI value in both periods, and with a value and at least 3 observations in both periods (`assemble.json`, `ndvi_versions_land_cells`):

| Version | Value in both periods | Value and at least 3 obs in both periods |
|---|---|---|
| 80th percentile (`_p80`) | 12,438 (100.0%) | 12,327 (99.1%) |
| Screened median | 12,437 (99.99%) | 746 (6.0%) |
| Unscreened median (`_noscreen`) | 12,438 (100.0%) | 12,327 (99.1%) |

The relative 80th-percentile version (`ndvi_rel_*_p80`, `ndvi_rel_change_p80`) was computed within the available memory. `veg_share_*` comes from the 80th-percentile NDVI, and `bare_share_*` and `sandfill_candidate` from the unscreened median BSI (`bsi_*_noscreen`). `bsi_*`, `veg_share_*_screened`, `bare_share_*_screened` and `sandfill_candidate_screened` come from the screened composites. Recomputing the screened shares in the `finalize` step reproduced the earlier values exactly (maximum difference 0).

## Landsat land surface temperature

LST comes from Landsat Collection 2 Level-2 surface temperature (ST_B10 band, scaled as DN × 0.00341802 + 149.0 K), read at 30 m on a grid aligned to the Landsat pixels.

- **Sensors:** Landsat 8 is used in both periods. Landsat 9 is added to the late period only, because it launched in 2021.
- **Masks and coverage:** scenes with 80% cloud cover or more are dropped. QA_PIXEL bits 0, 1, 2, 3, 4, 5 and 7 are masked. A date is kept only if at least 10% of Lagos land pixels are valid.
- **Date keys:** scenes of the same WRS path and date are merged. Landsat 9 date keys carry an `_l9` suffix.

Two rules were set before their effect was inspected:

1. **Valid range:** pixels below 15 °C or above 65 °C are masked as residual cloud or error.
2. **Date outliers:** after rule 1, a date is excluded when its median land LST is more than 2 sample SDs from the mean of its period's date medians. This is a single pass, and the period statistics are recomputed after Landsat 9 is added.

Results (`temperature.json`):

| | Early | Late |
|---|---|---|
| Scenes found / below cloud limit | 46 / 32 | 84 / 72 |
| Date keys / with enough coverage | 19 / 10 | 39 / 24 |
| Mean / SD of date medians (°C) | 32.80 / 2.78 | 32.67 / 3.63 |
| Dates excluded by rule 2 | none | 2024-01-15_p191 (24.26 °C), 2024-02-01_p190_l9 (21.67 °C) |
| Dates used | 10 (Landsat 8) | 22 (9 Landsat 8, 13 Landsat 9) |
| Land pixel observations below 15 °C, all dates / used dates | 29,287 / 29,287 | 207,064 / 83,240 |
| Land pixel observations above 65 °C, all dates / used dates | 71 / 71 | 144 / 144 |
| Land cells with fewer than 3 valid observations | 2,698 | 23 |

`lst_ok` marks cells with at least 3 valid observations (per-pixel mean over the cell) in both periods. No cell is dropped. With Landsat 9, 10,106 cells are `lst_ok`, including 9,739 land cells (`temperature.json`). The Landsat 8-only run gave 4,253 and 4,218 (`temperature_landsat8_only.json`). In that run, 7,385 land cells had fewer than 3 late observations, and 2024-01-15_p191 was the only excluded date.

Adding Landsat 9 makes observation counts unequal: 10 dates early against 22 late. Landsat 8 and 9 TIRS calibrations differ slightly. Relative LST removes differences shared across Lagos on a date, but not differences in how many dates enter each median.

## Population

Population comes from WorldPop Global2 R2025A, constrained, 100 m, for 2016 and 2023 (`population.json`). The LGA totals sum to 10,703,709 (2016) and 12,457,091 (2023) (`assemble.json`). Growth from 2016 to 2023 ranges from 15.7% to 18.7% across the 20 LGAs, with a median of 16.5% (`assemble.json`). Such uniform growth suggests the two years share one spatial model scaled to different totals, so differences in population change between places carry little information.

## Covariates

Major roads and industrial land come from OpenStreetMap via Overpass:

- **Roads:** classes motorway, trunk and primary, 1,419 ways (`covariates.json`). They give each cell its distance to the nearest major road.
- **Industrial land:** 286 polygons covering 75.7 km² (`covariates.json`).

Both were queried on 2026-11-01 and reflect OSM at that time, not 2016 or 2023.

## Combined flags

`analysis_ok` marks cells that are not water-flagged, have Open Buildings coverage of at least 0.5 in both years, and have at least 3 observations in both periods for both LST and the primary vegetation measure. The vegetation part uses the all-date count of the 80th-percentile version, not the screened count. No cell is dropped; flags only mark cells (`assemble.json`):

- Water-flagged cells: 2,452. Land cells: 12,438.
- Land cells with low Open Buildings coverage: 0.
- Land cells with low LST observations in either period: 2,699 (2,698 early, 23 late). Land cells with low vegetation observations (80th-percentile version) in either period but not low LST observations: 9.
- Land cells flagged low-observation on any of these: 2,708 (21.8% of land cells).
- `analysis_ok` cells: 9,730. If the screened-median observation counts were used instead, 634 cells would pass.
- `lst_ok` cells: 10,106, of which 9,739 are land cells.
- `ndvi_screened_ok` cells (not water-flagged, at least 3 screened observations in both periods): 746. This is the sample for the screened-median sensitivity check.
- Land cells flagged `sandfill_candidate` (share of 20 m pixels with BSI above 0 rises by at least 0.25 between periods, unscreened median BSI): 1,465. With the screened median BSI (`sandfill_candidate_screened`): 822 (`vegetation.json`, `assemble.json`). The flag is approximate: BSI, like NDVI, is computed from reflectances that haze alters, and the unscreened median includes the dates the screen marks as hazy. It marks candidates for inspection, not confirmed sand filling.

Quick-look maps are in `figures/phase1/`. `relative_ndvi_change.png` shows the primary change measure (relative 80th percentile). `ndvi_change_screened.png` (screened median, `ndvi_screened_ok` cells only) and `relative_ndvi_change_noscreen.png` (unscreened relative median) show the sensitivity measures. No map of absolute 80th-percentile change is drawn.

## Known mismatches and limits

- **Different dates per source:** buildings refer to 2016 and 2023. Vegetation and LST use the dry seasons December 2016 to February 2018 and December 2023 to February 2025. GHSL is available for 2015 and 2020 only.
- **Haze:** dry-season Harmattan haze varies strongly between dates. The 80th percentile favours each pixel's clearer observations but its absolute value still differs between periods on the reference forest, so change is measured on the relative 80th percentile. The reference-forest screen retains 3 early and 9 late dates for the sensitivity version. The bare share and sand-filling flag use BSI from all dates and are approximate. The AOT screen failed for the early period (see above). Relative NDVI and relative LST reduce the effect of a shift shared across Lagos on a date, not haze that is patchy within a date.
- **Partial Landsat coverage:** the two WRS paths each cover only part of Lagos, so the Lagos-wide median used as the relative LST reference changes with each date's footprint.
- **WorldCover:** WorldCover 2021 serves only as a baseline mask (water, tree references). WorldCover 2020 and 2021 come from different algorithm versions and are not compared.
- **Sensor mix:** the late LST period combines Landsat 8 and 9; the early period is Landsat 8 only.
