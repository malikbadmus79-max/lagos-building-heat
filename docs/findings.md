# Findings: building density, vegetation loss and surface heat in Lagos, 2016–2025

## Summary

Across 9,157 cells of 500 m in Lagos State, cells where building share grew more between 2016 and 2023 warmed more relative to the rest of the city between the dry seasons of 2016–18 and 2023–25. Relative land surface temperature (LST) change was 0.579 °C higher per 10 percentage points of building-share increase (95% interval 0.447 to 0.711), after controls and LGA fixed effects. The association kept its sign in all 13 sensitivity runs. Building growth was also associated with vegetation loss: relative NDVI change was 0.039 lower per 10 percentage points. Adding vegetation change to the heat model reduced the building coefficient by about half (share 0.522). The association was weaker in cells already built up in 2016 (0.190 °C per 10 percentage points) than in mostly vegetated cells (0.731). 329 hotspot cells, in the top fifth of both building change and relative LST change, hold 896,165 residents by WorldPop 2023 estimates, 8.5% of the main-sample total. These are associations in surface temperature, not measured causes or air temperature. Model residuals are strongly spatially autocorrelated (Moran's I 0.53 to 0.65), so the intervals may be too narrow.

## Scope and sources

This document reports the analysis fixed in `docs/phase2_analysis_plan.md`. Deviations from the plan are listed at the end of that file. It is the final version of `docs/findings_draft.md`, with the area names from Phase 3, a summary, the priority areas of Phase 4 and a complete Limits section. Every number below is copied from a file written by the project code, and the file is named next to it. `data/processed/findings_trace.csv` records the source of each number and whether it matches. The results describe associations between changes across 500 m cells, not causes. Intervals are 95% intervals: confidence intervals with standard errors clustered by 5 × 5-cell blocks, or block-bootstrap intervals where stated.

## Sample

The base sample has 9,730 `analysis_ok` cells. 573 cells with more than 10% water are removed, leaving a main sample of 9,157 cells in 621 blocks. 1,097 sand-filling candidates fall in the main sample; without them, 8,060 cells remain (`phase2.json`, `phase2_models.csv`).

| Run | Cells | Blocks |
|---|---|---|
| Main sample | 9,157 | 621 |
| Without sand-filling candidates | 8,060 | 615 |
| Without cells within 1 km of ocean or lagoon | 7,462 | 522 |
| Without cells with industrial land | 8,561 | 618 |
| Without cells within 250 m of a major road | 8,187 | 615 |
| Landsat 8-only LST, cells with at least 3 Landsat 8 observations per period | 4,081 | 451 |
| At least 5 LST observations per period | 2,465 | 320 |
| Terciles of early relative LST | 3,053; 3,052; 3,052 | 365; 463; 332 |
| Screened-median NDVI check (`ndvi_screened_ok`) | 746 | 216 |

Source: `phase2_models.csv`.

In the main sample, building share changed by 0 percentage points at the 25th percentile, 0.6 at the median and 4.8 at the 75th percentile. The interquartile range is therefore 4.8 percentage points, against the 10-point step used to report effects (`phase2.json`, `bshare_change_main`).

## H1: building change and relative LST change

Mean relative LST change rises from 0.005 °C (interval −0.123 to 0.140) in the decile with the largest building loss (D1, mean −0.7 percentage points) to 0.436 °C (0.320 to 0.550) in the decile with the largest gain (D10, mean +11.7 percentage points). The means of deciles D2 to D9 lie between 0.211 and 0.288 °C. The intervals are block-bootstrap intervals (`phase2_deciles.csv`, `figures/phase2/deciles_lst.png`).

In M1, relative LST change is 0.579 °C higher per 10 percentage points of building-share increase (interval 0.447 to 0.711; 9,157 cells; `phase2_models.csv`, run `H1_primary`). By the plan's rule, the association is positive.

The coefficient is positive with an interval excluding zero in all 13 H1 runs. By the plan's rule, H1 is robust (`phase2.json`, `H1`; `figures/phase2/coefficients.png`):

| Run | °C per 10 pp | 95% interval |
|---|---|---|
| Primary (M1) | 0.579 | 0.447 to 0.711 |
| Without sand-filling candidates | 0.504 | 0.368 to 0.640 |
| Without cells within 1 km of ocean or lagoon | 0.549 | 0.412 to 0.685 |
| Without cells with industrial land | 0.620 | 0.474 to 0.766 |
| Without cells within 250 m of a major road | 0.595 | 0.468 to 0.722 |
| Landsat 8-only LST | 0.311 | 0.119 to 0.504 |
| LST observation counts as controls | 0.536 | 0.400 to 0.673 |
| At least 5 LST observations per period | 0.436 | 0.259 to 0.613 |
| Without early relative LST as control | 0.608 | 0.470 to 0.747 |
| Lowest tercile of early relative LST | 1.193 | 0.586 to 1.800 |
| Middle tercile of early relative LST | 0.935 | 0.795 to 1.075 |
| Highest tercile of early relative LST | 0.139 | 0.035 to 0.243 |
| GHSL built-up change 2015–2020 instead of Open Buildings | 0.547 | 0.326 to 0.767 |

Source: `phase2_models.csv`.

## H2: building change and relative NDVI change

Mean relative 80th-percentile NDVI change falls from 0.024 (D1) to −0.033 (D10) across the deciles of building change (`phase2_deciles.csv`, `figures/phase2/deciles_ndvi.png`).

In M2, relative NDVI change is 0.039 lower per 10 percentage points of building-share increase (interval −0.046 to −0.032; `phase2_models.csv`, run `H2_primary`). By the plan's rule, the association is negative.

All four H2 runs are negative with intervals excluding zero, so H2 is robust by the plan's rule (`phase2.json`, `H2`):

| Run | Relative NDVI per 10 pp | 95% interval |
|---|---|---|
| Primary (M2) | −0.039 | −0.046 to −0.032 |
| Without sand-filling candidates | −0.034 | −0.040 to −0.028 |
| Unscreened relative median | −0.040 | −0.048 to −0.033 |
| With bare-share change as control | −0.034 | −0.040 to −0.028 |

The 80th-percentile and unscreened versions agree in sign, as the plan requires. The descriptive check on the 746 `ndvi_screened_ok` cells, which uses the screened median, no controls and no fixed effects, is also negative: −0.082 per 10 percentage points (−0.104 to −0.059), with a Pearson correlation of −0.406 (`phase2_models.csv`, run `H2_screened_descriptive`; `phase2.json`).

## H3: building–LST association with and without vegetation change

Fitting M1 and M3 on the same 9,157 cells gives the following (`phase2_h3.csv`; bootstrap with 2,000 block resamples):

| NDVI version | Building coefficient, M1 (°C per 10 pp) | Building coefficient, M3 | Difference | Share of M1 coefficient |
|---|---|---|---|---|
| Relative 80th percentile | 0.579 (0.447 to 0.711) | 0.277 (0.141 to 0.413) | 0.302 (0.249 to 0.365) | 0.522 (0.403 to 0.698) |
| Relative unscreened median | 0.580 (0.449 to 0.710) | 0.242 (0.097 to 0.387) | 0.338 (0.275 to 0.408) | 0.583 (0.446 to 0.802) |

When NDVI change is added, the building coefficient falls by about half in both versions and remains positive. The NDVI-change coefficient in M3 is −8.29 °C per unit of relative NDVI (−9.32 to −7.25) for the 80th-percentile version and −9.12 (−10.11 to −8.14) for the unscreened version (`phase2_h3.csv`).

The difference has the same sign in both versions (`phase2.json`, `H3`). This is a statistical decomposition. NDVI change is measured over the same period as LST change, and the comparison does not establish that vegetation loss carries the building–LST association.

## H4: built-up versus mostly vegetated cells in 2016

With the main group rules (`phase2_h4_groups.csv`, definition `main`):

| Group | Cells | Median building change | 90th percentile of building change | °C per 10 pp, joint model with group slopes | 95% interval |
|---|---|---|---|---|---|
| Built-up 2016 | 878 | 0.065 | 0.125 | 0.190 | 0.021 to 0.359 |
| Mostly vegetated | 4,018 | 0.000 | 0.017 | 0.731 | 0.297 to 1.166 |
| Neither | 4,261 | 0.034 | 0.094 | 0.653 | 0.508 to 0.799 |

Building change is a share of cell area, so 0.065 means 6.5 percentage points. No main-sample cell meets both group rules.

The difference between the built-up and vegetated slopes is −0.542 °C per 10 percentage points (−1.008 to −0.075). By the plan's rule, it is negative. M1 refitted separately within each group gives −0.036 (−0.184 to 0.113) for the built-up group, which is not distinguishable from zero, and 0.791 (0.357 to 1.226) for the vegetated group (`phase2_models.csv`, runs `H4_main_*_separate`).

Across the six group definitions, the difference has the same sign each time, so the sign is robust by the plan's rule. Its interval excludes zero in four of the six definitions (`phase2_h4_groups.csv`, `phase2.json`):

| Definition | Difference, °C per 10 pp | 95% interval |
|---|---|---|
| Main (built ≥ 0.30, tree ≥ 0.50) | −0.542 | −1.008 to −0.075 |
| Built threshold 0.20 | −0.547 | −1.000 to −0.094 |
| Built threshold 0.40 | −0.390 | −0.909 to 0.129 |
| Tree threshold 0.40 | −0.756 | −1.155 to −0.356 |
| Tree threshold 0.60 | −0.241 | −0.828 to 0.345 |
| Vegetated by early vegetated share ≥ 0.50 | −1.002 | −1.276 to −0.728 |

In the vegetated group, half or more of the cells have no building change, and 90% have less than 1.7 percentage points. Its slope therefore rests on a minority of cells with some building gain. The plan anticipated this.

## Rival explanations

All checks below come from `phase2_models.csv`; the H1 runs are listed in the H1 table above.

- Without the sand-filling candidates, the H1 coefficient is 0.504 and the H2 coefficient −0.034 per 10 percentage points. Both keep their sign and have intervals excluding zero. The candidate flag is approximate because of haze.
- Without cells within 1 km of the coast or the lagoon, H1 is 0.549 (0.412 to 0.685).
- Without industrial cells, H1 is 0.620 (0.474 to 0.766). Without cells within 250 m of a major road, it is 0.595 (0.468 to 0.722).
- With Landsat 8-only LST, which removes the mix of sensors between periods, H1 is 0.311 (0.119 to 0.504) on 4,081 cells. It keeps its sign, but it is smaller than the primary estimate and is estimated on a different, smaller set of cells.
- The unscreened NDVI versions agree in sign with the 80th-percentile versions for H2 and H3, and the screened-median descriptive check is also negative. Haze in the Sentinel-2 composites therefore does not account for the H2 result.
- With LST observation counts as controls, H1 is 0.536. With at least 5 observations per period, it is 0.436 on 2,465 cells.
- Without early relative LST as a control, H1 is 0.608. Within terciles of early relative LST, it is 1.193, 0.935 and 0.139. All are positive, but the coefficient is much smaller in the tercile that was warmest relative to Lagos in the early period, which is consistent with some regression to the mean.
- With GHSL built-up change (2015–2020) as the building measure, H1 is 0.547 (0.326 to 0.767).

## Spatial dependence of residuals

Moran's I of the residuals, with rook contiguity and 999 permutations (`phase2_moran.csv`):

| Model | Moran's I | Permutation p-value | Cells (cells without neighbours) |
|---|---|---|---|
| M1 | 0.639 | 0.001 | 9,105 (52) |
| M2 | 0.532 | 0.001 | 9,105 (52) |
| M3 (80th percentile) | 0.645 | 0.001 | 9,105 (52) |

The residuals remain strongly autocorrelated. As the plan states, this is reported as a limit of the block-clustered standard errors, and no spatial model was added. The intervals above may be too narrow if the dependence extends beyond the 2.5 km blocks.

## Exposure

Hotspot cells are cells in the top quintile of both relative LST change (at least 1.139 °C) and building-share change (at least 0.058) within the main sample. There are 329 hotspot cells (`phase2.json`, `exposure`).

- By WorldPop 2023, 896,165 residents live in hotspot cells, 8.5% of the main-sample total of 10,579,120.
- By GHS-POP 2020, used as a cross-check, 872,254 residents live in hotspot cells, 9.1% of the main-sample total of 9,616,655.

The two sources give similar shares (`phase2.json`, `phase2_exposure_lga.csv`).

LGAs with the largest hotspot populations by WorldPop 2023 (`phase2_exposure_lga.csv`):

| LGA | Hotspot cells | WorldPop 2023 | GHS-POP 2020 |
|---|---|---|---|
| Alimosho | 74 | 353,967 | 309,957 |
| Kosofe | 20 | 93,322 | 118,230 |
| Oshodi/Isolo | 24 | 77,077 | 93,445 |
| Ikorodu | 56 | 71,299 | 43,041 |
| Ifako/Ijaye | 11 | 58,598 | 54,660 |

These are counts of residents in cells that meet the hotspot rule. They do not describe individual exposure or health effects. The hotspot map is `figures/phase2/hotspots.png`.

## Case studies

The hotspot cells form 155 rook-contiguous clusters. Under the plan's rule, the three largest are chosen; the third has 10 cells, at least the required 4 (`phase2_case_studies.csv`, `phase2.json`):

| Case | Area (LGA) | Cells | Mean building change | Mean relative LST change (°C) | Mean relative NDVI change (p80) | WorldPop 2023 | Centre (lon, lat) |
|---|---|---|---|---|---|---|---|
| 1 | Lekki Peninsula II (Chevron area), Eti-Osa | 19 | 0.245 | 2.667 | −0.084 | 18,949 | 3.551, 6.438 |
| 2 | Abule Egba / Iyana Ipaja, Alimosho (cells also in Ifako/Ijaye) | 13 | 0.107 | 1.628 | −0.011 | 65,196 | 3.293, 6.643 |
| 3 | Oke, Ikorodu | 10 | 0.090 | 1.522 | −0.036 | 11,855 | 3.571, 6.627 |

Area names are OpenStreetMap names near each cluster centre (Nominatim reverse geocoding), confirmed by the author (`docs/phase3_validation.md`, `phase3_case_names.csv`). Case 1 lies in Lekki Peninsula II, around the Chevron area of Eti-Osa; case 2 around Abule Egba and Iyana Ipaja in Alimosho; case 3 around Oke in Ikorodu. The before/after maps of buildings, NDVI and LST for these areas are `figures/phase2/case1_before_after.png`, `case2_before_after.png` and `case3_before_after.png` (deviation 7 in the plan).

## Priority areas

Two lists of priority areas were ranked with rules fixed in `docs/phase4_priority_rules.md` before ranking (`data/processed/priority_areas_lga.csv`, `data/processed/priority_cells.csv`, `figures/phase4/`).

- "Cool now" ranks the 17 LGAs and 155 hotspot clusters that contain hotspot cells. The ranks of three measures get equal weight: residents in hotspot cells, mean relative LST change and mean relative NDVI loss. The first-ranked LGA is Eti Osa, and the first-ranked cluster is case 1 in Lekki Peninsula II (`phase4.json`). The 2021 Climate Change Risk Assessment, as cited in the C40 request for proposals, named Lagos Island, Lagos Mainland, Agege, Kosofe and Oshodi-Isolo as the LGAs most affected by heat, which describes where heat levels are highest, whereas this study measures where surface heat grew fastest relative to the rest of the city, which is concentrated on the expanding fringe (for example Eti-Osa, Epe and Alimosho); the two views are complementary, not contradictory.
- "Protect before it is lost" lists 104 cells, in 10 LGAs, that were still mostly vegetated in 2023–25 (vegetated share at least 0.50) but where building share grew by at least the hotspot cut of 0.058. Ikorodu has the most such cells, 33 (`phase4.json`).

The lists rank places by the measured changes and by modelled population. They do not measure need, land tenure or planning status.

## Limits

- The study measures land surface temperature from Landsat, not air temperature. Surface and air temperature differ, especially by day, and heat stress depends on air temperature, humidity and exposure.
- The results describe associations between changes across cells. They do not establish that building growth or vegetation loss causes warming, although the sensitivity runs rule out several other explanations.
- Building-change values in cells that were already densely built in 2016 are less reliable. In the imagery spot-check, the disagreements between the dataset and the imagery are concentrated in such cells, where the dataset records gains of several percentage points but no clear new building is visible (`docs/phase3_validation.md`). This is consistent with the weak association in the built-up group under H4.
- The residuals of the models are strongly spatially autocorrelated (Moran's I 0.53 to 0.65; `phase2_moran.csv`). The block-clustered intervals may therefore be too narrow.
- Population comes from WorldPop, a modelled distribution of census-based totals, not a count. Population change between 2016 and 2023 is nearly uniform across LGAs in this source (`docs/phase1_notes.md`). Resident counts in hotspot cells are therefore estimates.
- The satellite data end in February 2025 (`src/lagos_heat/common.py`, `PERIODS`). Buildings refer to mid-2016 and mid-2023, while vegetation and LST refer to the dry seasons of 2016–18 and 2023–25. Sensors and haze differ between periods (`docs/phase1_notes.md`). The Landsat 8-only estimate is about half the primary estimate.
- The first late dry season overlaps the 2023/24 El Niño. NOAA's Relative Oceanic Niño Index was above 0.5 °C from July–September 2023 to January–March 2024, reaching 1.4 in October–December and November–January (NOAA Climate Prediction Center, https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/, accessed 2026-11-03). Relative LST subtracts each date's Lagos-wide median, which removes city-wide shifts such as this one. Shifts that vary within the city are not removed.
- The estimate varies strongly across terciles of early relative LST, and building change in the vegetated group is concentrated in few cells.
- An AI coding agent (Claude Code) wrote the processing and analysis code under the author's direction, and carried out the validation recompute and the findings trace (`docs/phase3_validation.md`). The imagery spot-check of 40 cells was a visual check by an AI assistant, blinded to the data. A human reviewer confirmed its classification for 10 of the 40 cells, rather than coding them independently. The validation checks consistency and plausibility; it is not an independent field survey.
