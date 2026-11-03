# Phase 3: validation

This phase checks the Phase 1 dataset and the Phase 2 findings draft in four ways:

1. An independent recompute of 40 sampled cells.
2. A trace of every number in `docs/findings_draft.md`.
3. A spot-check of the same 40 cells against historical imagery: a blinded AI visual check, with a human check of 10 cells.
4. Place names for the case-study areas.

**Steps 1 and 2 were done by an AI coding agent (Claude Code).** The agent wrote and ran the code in `src/lagos_heat/validate.py` and `src/lagos_heat/trace_findings.py`. The agent also chose the sample rules and tolerances, which were fixed in the code before the recompute ran. These checks therefore test the consistency of the code and outputs, not their correctness against the ground. **Step 3 was also done by an AI assistant (Claude),** as a visual check of imagery blinded to the dataset values. A person then reviewed 10 of its cells. The reviewer confirmed the AI classification rather than coding the cells independently (section 3).

All figures below come from `data/interim/summaries/phase3.json` and the files named next to them.

## 1. Independent recompute of 40 cells

### Sample

40 cells were drawn at random (seed 20261103) from the 9,157-cell main sample of Phase 2, in three strata (`phase3.json`, `sample`):

| Stratum | Rule | Cells in stratum | Drawn |
|---|---|---|---|
| Hotspot | Hotspot cell in `phase2_sample.csv` | 329 | 15 |
| High change, not hotspot | Not a hotspot, and building-share change at or above the hotspot cut (0.058) | 1,503 | 15 |
| Low change | Absolute building-share change at or below the main-sample median (0.006) | 4,384 | 10 |

### Method

For each cell, six values were recomputed from the source rasters and tiles, then compared with `data/interim/grid_periods.csv`:

- building share in 2016 and 2023
- relative LST in the early and late periods
- relative 80th-percentile NDVI in the early and late periods

The recompute does not use the project's aggregation code (`common.aggregate`, `GridSpec.pixel_cells`, `buildings.read_tile`, `vegetation.cell_mean`):

- **Cell bounds** come from `grid_spec.json`.
- **LST and NDVI:** each processed raster is read as a small window. The value is the mean of the valid pixels whose centre lies inside the cell. For LST, validity is taken from the absolute LST raster, as in Phase 1.
- **Building share:** the Open Buildings 2.5D Temporal tiles were read again over HTTP through GDAL's own `/vsicurl/` driver, at the native 0.5 m pixel size. Building share is the fraction of valid pixels with `building_presence` above 0.5. Overlapping tiles are mosaicked so that each pixel is counted once.
  - Phase 1 read the same tiles at 4 m, through a separate reader based on the `requests` library.
  - The files cached in Phase 1 (`data/raw/_work/ob`) hold per-cell counts, not pixels, so they could not serve as an independent source.

The tolerances were fixed before the run:

| Value | Tolerance |
|---|---|
| Relative LST | 0.005 °C |
| Relative NDVI | 0.0005 |
| Building share | 0.02 (2 percentage points) |

The building-share tolerance is wider because the recompute counts 0.5 m pixels while the dataset sampled them at 4 m.

### Results

All 240 values (40 cells × 6 variables) match within tolerance (`phase3_recompute.csv`):

| Variable | Matches | Largest absolute difference | Mean difference |
|---|---|---|---|
| `bshare_2016` | 40 / 40 | 0.0126 | 0.0004 |
| `bshare_2023` | 40 / 40 | 0.0188 | −0.0020 |
| `lst_rel_early` | 40 / 40 | 0.0000 | 0.0000 |
| `lst_rel_late` | 40 / 40 | 0.0000 | 0.0000 |
| `ndvi_rel_early_p80` | 40 / 40 | 0.0000 | 0.0000 |
| `ndvi_rel_late_p80` | 40 / 40 | 0.0000 | 0.0000 |

Notes on the results:

- **Building share:** the largest differences are in densely built cells. Cell 14943, a hotspot, has a share of 0.563 in 2023 in the dataset and 0.544 at 0.5 m. In these cells the 4 m sample gives a slightly higher share than the full-resolution count. The largest difference, 0.0188, is close to the tolerance of 0.02.
- **Tile coverage:** all 80 building reads had full tile coverage of the cell. 1 to 4 tiles were used per cell.
- **Grid geometry:** the cell ids, rows and columns are consistent with the grid. The recomputed cell centres match the stored longitude and latitude to within 0.001 m.
- **Mismatches:** none.

## 2. Findings trace

`data/processed/phase3_findings_trace.csv` has one row per number in `docs/findings_draft.md`. Each row records the line, the text, the source file, the row and column (or JSON path, code constant or plan text), the source value and whether the two match. The number shown in the draft is compared with the source value rounded to the same number of decimals. Word quantities such as "four of the six" are checked as claims (`phase3.json`, `findings_trace`).

| | Count |
|---|---|
| Numbers in the draft | 343 |
| from a `phase2_*.csv` cell | 254 |
| from a code constant or `grid_spec.json` | 41 |
| from `phase2.json` | 20 |
| derived from the outputs (for example, the range of decile means) | 9 |
| stated in the analysis plan | 6 |
| labels, not quantities (for example "Landsat 8", "Phase 2", "WorldPop 2023") | 13 |
| Numbers without a source | 0 |
| Numbers checked | 330 |
| Numbers that match | 329 |
| Word claims checked / holding | 9 / 9 |

One number does not match:

- **Line 64:** the upper end of the interval for H2 with bare-share change as a control was given as −0.029. The source, `phase2_models.csv` (run `H2_bare_control`, column `per_10pp_ci_hi`), holds −0.02847, which rounds to −0.028.

One statement in the draft was also out of date. It is not a number. The draft said the case-study maps had not been produced, citing deviation 7 of the plan, but the maps were drawn later (`figures/phase2/case1_before_after.png` to `case3_before_after.png`).

Both were corrected on 2026-11-03, in the draft and in deviation 7 of the plan. The draft also gained a sentence on the case-study area names and a limit on building change in dense cells. After the corrections, the trace was run again on the revised draft:

| | Count |
|---|---|
| Numbers | 344 |
| Numbers checked | 331 |
| Numbers that match | 331 |
| Word claims checked / holding | 9 / 9 |
| Numbers without a source | 0 |

## 3. Imagery spot-check

### Method

**This check was carried out by an AI assistant (Claude), not by a human coder.** A human reviewer then confirmed the AI classification for 10 of the 40 cells (see below). The per-cell readings are in `data/processed/phase3_spotcheck_ai.csv`, and the original notes are in `docs/phase3_spotcheck_ai_notes.md`.

- **Cells:** the 40 cells of step 1, listed in `data/processed/phase3_spotcheck.csv`.
- **Imagery:** Esri World Imagery Wayback tiles at zoom 17, about 1.2 m per pixel. Each cell was shown as a 3 × 3 tile window centred on it, with the 500 m cell outlined.
  - The early image is the tile current in the 2017-06-27 release; the late image is the tile current in the 2024-01-18 release. For all 40 cells, these tiles were first published on 2017-03-15 and 2023-12-07.
  - Publication dates are not capture dates; the capture dates are earlier.
  - These are not the releases in the sheet's Wayback links (2016-07-06 and 2023-06-29).
- **Blinding:** the image pairs were shown in a shuffled order under a number only. The cell id, stratum and dataset values were joined after all 40 pairs had been classified.
- **Recorded for each pair:**
  - land cover in each image: bush/farmland, sand/cleared, sparse buildings, dense buildings or water/wetland
  - new buildings seen: none, small or large
  - greenery loss seen: none, small, some or large
- **Comparison rule,** set before unblinding: the dataset's building change is classed as none (below 2 percentage points), small (2 to 10) or large (above 10).
  - The same category as seen counts as agreement.
  - An adjacent category counts as partial agreement.
  - None against large counts as disagreement.

### Results (AI check)

- **Not assessable:** 4 cells (17634, 17997, 3634 and 15408), where cloud covered most of the cell in the 2017 image.
- **Presence or absence of new buildings,** where the dataset shows a change of at least 2 percentage points or not: 31 of the 36 assessable cells agree (86%).
- **Three categories:** 20 cells agree, 15 agree partly and 1 disagrees.

  | Stratum | Agree | Partly | Disagree |
  |---|---|---|---|
  | Hotspot | 10 | 5 | 0 |
  | High change, not hotspot | 3 | 9 | 1 |
  | Low change | 7 | 1 | 0 |

- **Low-change stratum:** all 8 assessable cells show little or no new building. One of them (6178) shows a few new structures in a forested area, against a dataset change of 0.1 percentage points.
- **Dense cells:** the disagreements are concentrated in cells that were already densely built in the early image.
  - Of the 10 assessable cells read as dense buildings in 2017, 3 agree, 6 agree partly and 1 disagrees. Of the other 26 cells, 17 agree, 9 agree partly and none disagree.
  - Four presence/absence disagreements are dense cells: 12344, 7908, 5671 and 14943, of which two are hotspot cells. In these the dataset records a gain of 7 to 10 percentage points, but no clear new building is visible at this resolution.
  - This may reflect infill too small to see, or a change between 2016 and 2023 in how Open Buildings detects buildings in dense areas.
  - It is consistent with the weak building–heat association in cells already built up in 2016 (H4). It also means that building-change values in dense cells should be read with caution. This is now listed among the limits in `docs/findings_draft.md`.
- **Precision:** the size of the change was judged by eye, so the three-category comparison is approximate. The presence/absence comparison is the main result.

### Human check

Ten cells were marked for a human check (`human_check = yes`): the four dense-cell disagreements (5671, 7908, 12344 and 14943) and six agreements (8264, 21649, 10862, 2889, 8358 and 8842). The AI notes describe the six agreements as randomly chosen.

- **Agreement:** the reviewer agreed with the AI classification in 10 of 10 cells. This is the count of rows with `human_agrees_with_ai = yes` in `data/processed/phase3_spotcheck_ai.csv`.
- **How the review was done:** the reviewer looked at the AI classification and confirmed it; the cells were not coded independently. The reviewer did not record separate land-cover categories, so `human_seen_2016` and `human_seen_2023` read "same as AI (confirmed by reviewer)".
- **Recording:** the entries were written into the CSV from the reviewer's confirmation in the Claude Code conversation.
- **What it shows:** the result is a confirmation of the AI readings, including the four dense-cell disagreements with the dataset. It is not an independent agreement rate between two coders, and it does not test the 30 cells that were not reviewed.

## 4. Case-study place names

The centre of each case-study cluster was reverse-geocoded with OpenStreetMap Nominatim at zoom levels 18, 16 and 14. The centre of every cell in each cluster was also queried at zoom level 16.

- 51 requests were sent, at most one per second, with an identifying user agent. One request timed out (a cell in case 1).
- Raw responses are in `data/interim/phase3_nominatim.json`, and the tally of names is in `data/processed/phase3_case_names.csv`.
- Names are from OpenStreetMap: data © OpenStreetMap contributors, ODbL.

In the table, "points" counts the query points (centre queries and cell centres) whose address includes the name. "Centre" marks names returned for the cluster centre.

| Case | Candidates from Nominatim (points) | Name used, confirmed by the author |
|---|---|---|
| 1 | Eleganza (19, centre), Lekki Peninsula II (13, centre), Ajah (8), Lekki County Homes (6), Chevron (2), Ikota (1), Ocean Bay Estate (1) | Lekki Peninsula II (Chevron area), Eti-Osa |
| 2 | Abule Egba (15, centre), Iyana Ipaja (2), Iju (1); Alimosho also appears once, as the LGA name | Abule Egba / Iyana Ipaja, Alimosho |
| 3 | Oke (11, centre), Elesin (1), Parafa (1) | Oke, Ikorodu |

The confirmed names are used in the case-study table of `docs/findings_draft.md`.
