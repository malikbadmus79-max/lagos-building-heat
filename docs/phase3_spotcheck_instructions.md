# Phase 3 spot-check: instructions

This check compares what the dataset says about 40 cells with what can be seen on free historical satellite imagery. The cells are listed in `data/processed/phase3_spotcheck.csv`, one row per cell. Each cell is a 500 m × 500 m square. The sample was drawn at random with a fixed seed from the main analysis sample, in three groups:

- `hotspot`: 15 cells in the top fifth of both building change and relative LST change.
- `high_change_not_hotspot`: 15 cells with large building change that are not hotspots.
- `low_change`: 10 cells with little or no building change.

## What the dataset says

The `data_*` columns hold the dataset values for each cell:

| Column | Meaning |
|---|---|
| `data_building_share_2016`, `data_building_share_2023` | Share of the cell covered by buildings in Open Buildings 2.5D Temporal (mid-2016 and mid-2023), from 0 to 1. A value of 0.25 means about a quarter of the cell is roofs. |
| `data_building_change_pp` | Change in building share, in percentage points. |
| `data_ndvi_change_direction` | Direction of the change in relative vegetation index (80th-percentile NDVI) between the dry seasons 2016–18 and 2023–25. A change smaller than 0.02 in either direction is called "little change". |
| `data_ndvi_rel_change_p80`, `data_lst_rel_change_c` | The vegetation change and the relative land surface temperature change, in °C, for reference. |

## Imagery to use

Two free sources work, and either is enough:

1. **Google Earth Pro** (desktop, free).
   - Open `data/processed/phase3_spotcheck_cells.kml` with File > Open. It draws the 40 cell squares in yellow, each named by its `cell_id`.
   - Double-click a cell name in the Places panel to fly to it.
   - Turn on historical imagery with the clock icon in the toolbar (View > Historical Imagery). Move the time slider to the image closest to **2016** and then to the image closest to **2023**. The date of each image is shown in the slider and at the bottom of the view.
   - Prefer images from 2016 or early 2017 for "2016", and from 2023 or early 2024 for "2023". If the nearest image is more than about a year from those dates, note its date in `notes`.
2. **Esri World Imagery Wayback** (in a web browser).
   - The `wayback_2016_url` and `wayback_2023_url` columns open Wayback at the cell centre, using the release closest to mid-2016 and to mid-2023: release 13240 (2016-07-06) and release 47963 (2023-06-29). These are recorded in `data/interim/summaries/phase3.json`, under `spotcheck`.
   - A Wayback release is a snapshot of the imagery archive on its release date, not the capture date. The image shown for Lagos may be older. Click on the map to see the capture date of the image at that point.
   - If a link does not open the expected release, choose it from the release list on the left.
   - Wayback does not show the cell square. The cell is the 500 m square centred on the link's point. At the zoom level of the links, that is roughly the middle of the view.

`google_maps_url` opens the current satellite view at the cell centre. It is useful for orientation but shows recent imagery only.

## How to fill the sheet

For each row, look at the whole 500 m square in both years. Then fill four columns:

- `seen_2016`, `seen_2023`: the land cover that covers most of the cell. Use one of:
  - `bush/farmland`
  - `sand/cleared`
  - `sparse buildings`
  - `dense buildings`
  - `water/wetland`

  If two cover types are about equal, write both separated by `;`, for example `sand/cleared; sparse buildings`.
- `agrees_with_data`:
  - `yes` if the imagery shows the change the dataset reports: more, fewer or about the same buildings, and greener, less green or about the same vegetation.
  - `partly` if only one of the two agrees, or if the direction agrees but the size clearly does not. For example, the data says +25 percentage points but almost nothing new is visible.
  - `no` if neither agrees.
- `notes`: anything that affects the reading, such as:
  - cloud or haze on the image
  - image dates far from 2016 or 2023
  - construction sites
  - sand filling of wetland
  - large single buildings such as warehouses
  - a cell that is mostly a road or a port

A rough guide to building share: 0.05 is a few scattered roofs, 0.15–0.25 is a typical built-up neighbourhood with gardens, and 0.35 or more is densely packed roofs.

## After the check

When the sheet is filled in, the results can be summarised by stratum: the share of `yes`, `partly` and `no`, and the cells marked `no` with their notes. They then go into the "Manual spot-check" section of `docs/phase3_validation.md`, which is marked pending until then.
