# Phase 2 analysis plan

This plan was written after the Phase 1 dataset was built and before any model was run or any relationship between the analysis variables was inspected. Every choice below, including thresholds, sample rules and decision rules, is fixed in advance. Any later change to the plan is recorded in a dated "Deviations" section at the end of this file, with its reason, and the original wording stays in place.

Column names refer to `data/interim/grid_periods.csv`. Counts quoted from Phase 1 come from `data/interim/summaries/assemble.json` and `vegetation.json`; all other counts are produced by the Phase 2 code.

## 1. Sample

The unit is the 500 m grid cell.

- **Base sample:** cells with `analysis_ok` (9,730 cells, `assemble.json`). These are cells that are not water-flagged, have Open Buildings coverage of at least 0.5 in both years, and have at least 3 observations in both periods for LST and for the 80th-percentile NDVI.
- **Main sample:** the base sample minus cells with more than 10% water (`water_frac > 0.10`). Mixed land–water cells are excluded because water lowers LST and NDVI in a way that does not depend on building change. The code reports the number of cells removed by this rule and the size of the main sample.
- **Sand-filling sensitivity:** the main sample minus the cells flagged `sandfill_candidate` (1,465 land cells, `vegetation.json`; the code reports how many of them fall in the main sample).
- Further restricted samples used for rival explanations are listed in section 5.

Every analysis uses complete cases on the variables of its model. Each output table reports its number of cells.

## 2. Variables

All changes are late minus early.

| Role | Primary | Sensitivity or cross-check |
|---|---|---|
| Building change | `bshare_change`: Open Buildings building presence share, 2023 minus 2016 (share of cell area, 0 to 1) | `ghsl_share_change`: GHSL built-up share, 2020 minus 2015 (cross-check; it covers a different period) |
| Vegetation change | `ndvi_rel_change_p80`: relative 80th-percentile NDVI | `ndvi_rel_change_noscreen`: relative unscreened median NDVI. A vegetation finding is reported as a finding only if its sign agrees in both versions. |
| Heat change | `lst_rel_change`: relative LST, °C (Landsat 8 early; Landsat 8 and 9 late) | Landsat 8-only relative LST change from `data/interim/temperature_cells_landsat8_only.csv` |

The vegetation and heat periods are the dry seasons December 2016 to February 2018 and December 2023 to February 2025. The building years are 2016 and 2023. Absolute 80th-percentile NDVI is used only as a level within a period, never as a change (`docs/phase1_notes.md`).

Controls and grouping variables:

- **Baseline levels:**
  - 2016 building share: `bshare_2016`
  - early relative 80th-percentile NDVI: `ndvi_rel_early_p80`; `ndvi_rel_early_noscreen` in the unscreened sensitivity
  - early relative LST: `lst_rel_early`
- **Location:**
  - distance to ocean: `dist_ocean_m`
  - distance to lagoon: `dist_lagoon_m`
  - distance to the nearest major road (OSM motorway, trunk or primary): `dist_major_road_m`
  - industrial share: `industrial_share`
- **Bare-soil change:** `bare_share_change`, from the unscreened median BSI. It is approximate because of haze.
- **LGA:** `lga`, the LGA with the largest overlap with the cell. 20 LGAs, used as fixed effects.

Distances enter in kilometres.

## 3. Questions and expected directions

- **H1:** cells with a larger increase in building share show a larger increase in relative LST. Expected sign of the building coefficient: positive.
- **H2:** a larger increase in building share goes with a larger loss of relative NDVI. Expected sign: negative.
- **H3:** how much of the building–LST association is accounted for by vegetation loss. This is measured by comparing the building coefficient with and without NDVI change in the model. The comparison is described as a statistical decomposition only, not as causal mediation.
- **H4:** whether the building–LST association differs between two groups of cells:
  - **Built-up in 2016:** `bshare_2016 >= 0.30` and `tree_frac_wc2021 < 0.50`.
  - **Mostly vegetated:** `tree_frac_wc2021 >= 0.50` and `bshare_2016 < 0.30`.
  - Cells meeting neither rule, or both, belong to neither group. Their count is reported.
  - The thresholds were set without inspecting the distributions in the sample. Sensitivity thresholds are 0.20 and 0.40 for building share, and 0.40 and 0.60 for tree share, each varied with the other held at its main value.
  - WorldCover tree cover is from 2021, after the 2016 baseline, so cells cleared between 2016 and 2021 do not count as vegetated. A second definition of the vegetated group therefore uses the early-period vegetated share, `veg_share_early >= 0.50` (share of 20 m pixels with 80th-percentile NDVI above 0.4, dry seasons 2016–18), with `bshare_2016 < 0.30`.

Decision rules for reporting:

- A coefficient is described as positive or negative when its 95% confidence interval, with block-clustered standard errors, excludes zero. Otherwise it is described as not distinguishable from zero.
- A result is described as robust when its sign is the same in the primary model and in every sensitivity run in section 5. Sign changes are reported as they occur.
- No sensitivity run replaces the primary result.

## 4. Methods

### 4.1 Main chart

- Cells in the main sample are split into deciles of `bshare_change`, by rank so that ties are split evenly.
- For each decile, the chart shows the mean `lst_rel_change` with a 95% confidence interval. The interval comes from a spatial block bootstrap: 2,000 resamples of the 5 × 5-cell blocks defined in 4.3, percentile intervals.
- The x-axis shows each decile's mean building change in percentage points.
- A companion chart shows the same layout for `ndvi_rel_change_p80` (H2).
- No fitted line is drawn on the main chart.

### 4.2 Regression models

All models are ordinary least squares at cell level.

- **M1 (H1):**

  `lst_rel_change ~ bshare_change + bshare_2016 + ndvi_rel_early_p80 + lst_rel_early + dist_ocean_km + dist_lagoon_km + dist_major_road_km + industrial_share + bare_share_change + C(lga)`

- **M2 (H2):** `ndvi_rel_change_p80` on `bshare_change` with the same controls as M1, except two:
  - `lst_rel_early` is left out.
  - `bare_share_change` is left out, because BSI and NDVI share the red and near-infrared bands. A sensitivity run adds it back.
- **M3 (H3):** M1 plus `ndvi_rel_change_p80`, fitted on exactly the same cells as M1.
  - Reported: the building coefficient in M1 and in M3, their difference, and the difference as a share of the M1 coefficient.
  - A 95% interval for the difference and for the share comes from the block bootstrap in 4.1, refitting both models in every resample.
  - The share is not reported if the M1 coefficient's confidence interval includes zero.
- **M4 (H4):** M1 with `bshare_change` interacted with group (built-up, vegetated, neither), fitted on the main sample.
  - Reported: the building coefficient within each group, and the difference between the built-up and vegetated groups with its confidence interval.
  - M1 is also refitted separately within each of the two groups.
  - The vegetated group is expected to have little room for decline in building share and possibly little building change. The number of cells and the distribution of `bshare_change` in each group are reported alongside.

### 4.3 Uncertainty and spatial dependence

- **Standard errors:** clustered by spatial blocks of 5 × 5 cells (2.5 km × 2.5 km). The block id is `(row // 5, col // 5)` on the Phase 1 grid. The number of blocks in each model is reported.
- **Residual autocorrelation:** for M1, M2 and M3, Moran's I of the residuals uses rook contiguity between sample cells (shared edge on the grid), row-standardised. Significance comes from a permutation test with 999 permutations.
  - libpysal is not installed, so Moran's I is implemented in the project code with scipy sparse matrices. A unit test checks it against a hand-computed small case.
  - If residual Moran's I stays clearly positive, this is reported as a limit of the clustered standard errors. No spatial lag or error model is added after seeing results.

### 4.4 Effect in plain terms

- The building effect is reported as the change in relative LST, in °C, per 10-percentage-point increase in building share. This is the `bshare_change` coefficient × 0.10, with its 95% interval.
- The same scaling is used for H2, in NDVI units per 10 percentage points.
- The text also gives the interquartile range of `bshare_change` in the main sample, so the 10-point step can be compared with the observed spread.

## 5. Rival explanations and how each is tested

Each test reruns M1 (and M2 where stated) on the stated sample or with the stated change, and reports the building coefficient per 10 percentage points next to the primary value.

| Rival explanation | Test |
|---|---|
| Sand filling (bare sand raises LST and lowers NDVI without buildings) | Main sample excluding `sandfill_candidate` cells (M1, M2). `bare_share_change` is already a control in M1. The flag is approximate because of haze. |
| Coast and lagoon (sea and lagoon breezes, shoreline change) | Distances to ocean and lagoon are controls. Sensitivity: exclude cells within 1 km of the ocean or the lagoon. The water rule (more than 10% water) already removes mixed cells. |
| Industry and roads (heat sources and paved surfaces not in the building layer) | Industrial share and road distance are controls. Sensitivity: exclude cells with `industrial_share > 0`. A further sensitivity excludes cells within 250 m of a major road. |
| Sensor mix (Landsat 8 early, Landsat 8 and 9 late) | Rerun M1 with the Landsat 8-only relative LST change on cells with at least 3 Landsat 8-only observations in both periods (`lst_ok` in the Landsat 8-only file) that are in the main sample. The cell count of this sample is reported; it is expected to be much smaller. |
| Haze (NDVI and BSI biased by atmospheric haze that varies between dates) | Primary NDVI change is relative and based on the 80th percentile. Rerun M2 and M3 with `ndvi_rel_change_noscreen`; signs must agree (section 2). As a descriptive check only, `ndvi_change` (screened median) is compared with `bshare_change` on the 746 `ndvi_screened_ok` cells, with no controls. |
| Observation counts (10 LST dates early, 22 late; unequal per-cell counts) | Add `lst_nobs_early` and `lst_nobs_late` as controls. Sensitivity: restrict to cells with at least 5 LST observations in both periods. |
| Regression to the mean (cells with extreme early values move back toward the mean) | Early relative LST is a control in M1. Sensitivity: refit M1 without `lst_rel_early`, and report both. Also report the building coefficient within terciles of `lst_rel_early`. |
| Building measure (Open Buildings model error) | Replace `bshare_change` with `ghsl_share_change` (2015–2020) in M1, as a cross-check of direction only, since the periods differ. |

## 6. Exposure

- **Hotspot cells:** defined within the main sample as cells in the top quintile of `lst_rel_change` and also in the top quintile of `bshare_change`. Quintile cut points are computed on the main sample and written to the output file.
- **Population:** the WorldPop 2023 population (`pop_2023`, constrained, 100 m) of hotspot cells is summed for Lagos and by LGA. It is also reported as a share of the main-sample population.
- **Cross-check with GHSL GHS-POP 2020:** the 100 m product, epoch 2020, release R2023A. It is a small download of the tile(s) covering Lagos, aggregated to cells by area. The same totals are computed from it.
  - GHS-POP refers to 2020, not 2023, so the comparison concerns the share of population in hotspot cells, not the absolute numbers.
  - Each source's totals for the main sample are reported alongside.
- **Wording:** population counts are counts of residents in cells that meet the hotspot rule. They do not describe individual exposure or health effects.

## 7. Case studies

- **Clusters:** hotspot cells (section 6) are grouped into clusters by rook contiguity.
- **Ranking:** clusters are ranked by number of cells, with ties broken by the sum of `bshare_change`.
- **Selection:**
  - The two largest clusters are chosen.
  - The third is added only if it has at least 4 cells (1 km²).
  - If fewer than two clusters have 2 or more cells, this is reported, and the case studies use the largest clusters available.
- **Maps:** for each case study, a map window covers the cluster's bounding box plus 1 km on every side. Each window has early/late map pairs of:
  - building presence (Open Buildings 2.5D Temporal, re-read at 4 m for the window from the same source as Phase 1)
  - 80th-percentile NDVI (20 m composites in `data/processed/`)
  - relative LST (30 m composites in `data/processed/`)
- **Text:** each case study gives the LGA, the number of cells and the cell-level values of the three change variables. The text is descriptive.

## 8. Outputs

- **Tables:** `data/processed/phase2_*.csv`:
  - `phase2_sample.csv`: sample flags per cell
  - `phase2_deciles.csv`
  - `phase2_models.csv`: coefficients, intervals, number of cells and number of blocks for every model and sensitivity run
  - `phase2_moran.csv`
  - `phase2_h4_groups.csv`
  - `phase2_exposure_lga.csv`
  - `phase2_case_studies.csv`
- **Summary:** `data/interim/summaries/phase2.json`.
- **Figures:** `figures/phase2/`: the main decile chart, the H2 decile chart, a coefficient chart across sensitivity runs, an exposure map, and the case-study map panels.
- **Draft:** `docs/findings_draft.md`.
- **Rules:**
  - Every number in the findings comes from a file the code writes, and the file is named next to it.
  - Findings describe associations between changes across cells, not causes.
  - Results that contradict the expected directions are reported in the same way as those that agree.

## Deviations

### 2026-11-02

1. **Decile ties (4.1).** The first implementation ranked tied values in file order. File order follows the grid rows from north to south, so the cells with zero building change were assigned to deciles 2 to 4 by location. This was noticed in the first decile chart, before the findings were written. Ties are now ordered at random with a fixed seed (`SEED` in `phase2.py`), which matches the plan's intent that ties are split evenly. Only the decile table and charts are affected; no model uses the deciles.
2. **Decile chart axis (4.1).** Several deciles have a mean building change at or near 0 percentage points, so points placed at their mean change overlapped. The chart instead places the deciles at equal spacing and labels each tick with the decile's mean change in percentage points.
3. **Constant controls.** In the run that excludes industrial cells, `industrial_share` is 0 for every cell and is collinear with the intercept, so it is left out of that model. Any control dropped this way is listed in `constant_controls_dropped` in `phase2_models.csv`.
4. **H4 sensitivity thresholds (3), interpretation.** When the building-share threshold is varied (0.20, 0.40), the same threshold is used in both group rules: built-up means at least the threshold, and vegetated means below it. For the second vegetated definition (`veg_share_early >= 0.50`), the built-up group keeps its main definition.
5. **Haze check (5), interpretation.** The screened-median comparison on the 746 `ndvi_screened_ok` cells is described in the plan as a descriptive check only. It is therefore reported separately and not counted in the robustness rule for H2.
6. **Hotspot cut points (6), interpretation.** "Top quintile" means at or above the 80th percentile of the main sample.
7. **Case-study maps (7), produced in a separate run.** The Open Buildings window reads for the maps were first stopped by the Claude Code memory reaper while system memory was low. The case-study table (`phase2_case_studies.csv`) was not affected. The maps were then drawn in a separate run of `python -m lagos_heat.phase2 --maps-only` and are `figures/phase2/case1_before_after.png` to `case3_before_after.png` (noted 2026-11-03).
