# Phase 3, step 3: imagery spot-check (AI visual check)

This check was carried out by an AI assistant (Claude), not by a human coder. A human check of 10 of the 40 cells is pending (column `human_check` in `data/processed/phase3_spotcheck_ai.csv`).

## Method

- Cells: the 40 cells in `data/processed/phase3_spotcheck.csv` (15 hotspot, 15 high building change but not hotspot, 10 low change).
- Imagery: Esri World Imagery Wayback tiles at zoom 17 (about 1.2 m per pixel), a 3 × 3 tile window centred on each cell with the 500 m cell outlined. The early image is the tile current in the 2017-06-27 release; the late image is the tile current in the 2024-01-18 release. For all 40 cells these tiles were first published on 2017-03-15 and 2023-12-07. Publication dates are not capture dates; capture dates are earlier.
- Blinding: the image pairs were shown in a shuffled order under a number only. The cell ID, stratum and dataset values were joined after all 40 pairs had been classified.
- Recorded for each pair: land cover seen in each image (bush/farmland, sand/cleared, sparse buildings, dense buildings, water/wetland), new buildings seen (none, small, large) and greenery loss seen (none, small, some, large).
- Comparison rule, set before unblinding: dataset building change below 2 percentage points = none, 2 to 10 = small, above 10 = large. Same category = agrees; adjacent category = partly; none against large = no.

## Results

- Not assessable: 4 cells, where cloud covered most of the cell in the 2017 image.
- Of the 36 assessable cells:
  - New buildings present or absent (dataset change of at least 2 percentage points against any new buildings seen): 31 of 36 agree (86%).
  - Three categories: 20 agree, 15 agree partly (adjacent category), 1 disagrees.
  - Low-change stratum: all 9 assessable cells agree that there was little or no new building. One shows a few new structures in a forested area.
- Disagreements are concentrated in cells that were already densely built in 2017. In four such cells (12344, 7908, 5671 and 14943; two are hotspot cells), the dataset records a gain of 7 to 10 percentage points, but no clear new building is visible at this resolution. This may reflect infill too small to see, or a change in how Open Buildings detects buildings in dense areas between 2016 and 2023. It is consistent with the Phase 2 finding (H4) that the building–heat association is weak in cells already built up in 2016, and it means building-change values in dense cells should be read with caution.
- The size of the change was judged by eye, so the three-category comparison is approximate; the presence/absence comparison is the main result.

## Human check (pending)

Ten cells are marked for a human check: the four dense-cell disagreements and six randomly chosen agreements. The person checking fills `human_seen_2016`, `human_seen_2023` and `human_agrees_with_ai`, using the Wayback links in the sheet.
