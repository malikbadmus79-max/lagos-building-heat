# Phase 4: rules for the priority-area lists

These rules were written on 2026-11-03, before any ranking was computed. `src/lagos_heat/phase4.py` applies them. Any departure is listed under Deviations at the end, with its date and reason.

All inputs come from earlier phases:

- `data/processed/phase2_sample.csv`: main sample, hotspot flag, hotspot cluster id
- `data/interim/grid_periods.csv`: cell values
- `data/interim/summaries/phase2.json`: hotspot cut points

Population is WorldPop 2023 (`pop_2023`), which is modelled. Both lists use only cells in the Phase 2 main sample, so the same water, data-quality and coverage rules apply as in the findings.

## Common ranking rule

For each list, every unit (LGA, hotspot cluster or cell) gets a rank on each metric. Rank 1 is the most severe value, and tied values share the lowest rank (method `min`). The ranks are averaged with equal weight, and the units are ordered by mean rank, lowest first. Ties in mean rank are broken by the population metric of the list, higher first, and then by unit id.

## List A: "Cool now"

The units are the hotspot cells defined in Phase 2: cells in the top quintile of both relative LST change and building-share change within the main sample. Two rankings are made over the same three metrics.

- A1, LGAs: LGAs with at least one hotspot cell. The unit is the LGA of each hotspot cell, as recorded in `grid_periods.csv`.
- A2, hotspot clusters: the rook-contiguous clusters of hotspot cells from Phase 2 (`cluster` in `phase2_sample.csv`). Single-cell clusters are included.

| Metric | Definition | Direction |
|---|---|---|
| Residents in hotspot cells | sum of `pop_2023` over the unit's hotspot cells | higher = rank 1 |
| Mean relative LST change | mean of `lst_rel_change` over the unit's hotspot cells, °C | higher = rank 1 |
| Mean relative NDVI loss | minus the mean of `ndvi_rel_change_p80` over the unit's hotspot cells | higher loss = rank 1 |

The tie-break is residents in hotspot cells.

## List B: "Protect before it is lost"

The units are main-sample cells that meet two conditions:

- Still mostly vegetated in the late period: `veg_share_late >= 0.50`. This column is the late-period share of 20 m pixels with NDVI above 0.4, computed from the 80th-percentile NDVI composite (`docs/phase1_notes.md`). It is the measure requested as `veg_share_late_p80`; the dataset stores the p80 version under the unsuffixed name.
- Building share grew above the hotspot cut: `bshare_change` at or above the Phase 2 hotspot cut for building change, `exposure/quintile_cut_bshare_change` in `phase2.json`.

| Metric | Definition | Direction |
|---|---|---|
| Building change | `bshare_change` | higher = rank 1 |
| Residents nearby | sum of `pop_2023` over the cell and its eight neighbours (a 1.5 km × 1.5 km block), including neighbours outside the main sample | higher = rank 1 |

The tie-break is building change, then cell id.

The LGA summary (B, by LGA) gives, for each LGA with at least one qualifying cell:

- the number of qualifying cells
- the mean and maximum building change
- the sum of `pop_2023` in the qualifying cells
- the best (lowest) cell rank

LGAs are listed by number of cells, then by the sum of `pop_2023`.

## Names and outputs

Area names are added for the top 10 hotspot clusters (A2) and the top 10 cells (B):

- Each cluster centre or cell centre is reverse-geocoded with OpenStreetMap Nominatim at zoom level 16, at most one request per second.
- The name recorded is the most local of `neighbourhood`, `quarter`, `suburb`, `residential`, `hamlet`, `village`, `town` and `city_district`.
- For the three case-study clusters, the names confirmed in Phase 3 are used instead.
- Names are labels for orientation and are not verified on the ground.

Outputs:

- `data/processed/priority_areas_lga.csv`: lists A1 and B by LGA, with a `list` column.
- `data/processed/priority_cells.csv`: lists A2 (clusters) and B (cells), with a `list` column.
- `data/interim/summaries/phase4.json`: rules, counts and top 10s.
- `figures/phase4/`: maps of both lists.

## Deviations

### 2026-11-03

1. Area names (labels only; ranks unchanged).
   - Where the zoom-16 lookup returned no local name, the lookup is repeated at zoom 14 and then zoom 12 (wider areas). The first local name found is recorded as "near <name>".
   - A name is also recorded as "near <name>" when OpenStreetMap places it in a different LGA than the boundary data (the LGA in the address's `county` field). The LGA from the boundary data (geoBoundaries) is kept, and the conflict is described in the `name_note` column.
   - Reason: three of the top 10 protect cells had no name at zoom 16. Two names (Ikotun, Glory Garden) are placed by OpenStreetMap in a neighbouring LGA.
