# Phase 0 plan: data feasibility

Purpose: confirm that each core dataset can be accessed programmatically for Lagos State before any analysis is designed around it. Checks use metadata queries and small test windows so that downloads stay small. Results are written by `src/lagos_heat/phase0_check.py` to `docs/phase0_results.md`.

## Definitions

- Study area: Lagos State, Nigeria. A coarse bounding box (approx. 2.70 to 4.40 E, 6.35 to 6.75 N) is used for catalogue searches.
- Test window: a small box over central mainland Lagos (approx. 3.33 to 3.38 E, 6.48 to 6.53 N) is used for pixel reads.
- Fringe window: a small box on the urban fringe in Ibeju-Lekki (3.96 to 3.98 E, 6.44 to 6.46 N) is used where change over time needs to be visible. The central window is almost fully built up and shows little change.
- Dry season: December to February. "Dry season 2017" means 1 Dec 2016 to 28 Feb 2017; "dry season 2025" means 1 Dec 2024 to 28 Feb 2025.
- Low cloud: scene-level `eo:cloud_cover` below 20 percent.

## Checks

a. Lagos LGA boundaries
   - Source: geoBoundaries gbOpen NGA ADM2 (fallback candidate: GRID3 Nigeria LGA boundaries).
   - Method: geoBoundaries API, download ADM1 and ADM2, select ADM2 units falling inside the Lagos ADM1 polygon.
   - Pass condition: 20 LGAs returned.

b. Sentinel-2 L2A
   - Source: Microsoft Planetary Computer STAC, collection `sentinel-2-l2a`.
   - Method: anonymous STAC search with signed asset URLs, no API key.
   - Report: scene count (all and low-cloud) for dry seasons 2017 and 2025, cloud cover summary, test read of one band.

c. Landsat Collection 2 Level-2 surface temperature
   - Source: Planetary Computer STAC, collection `landsat-c2-l2` (ST_B10, asset `lwir11`).
   - Report: scene count (all and low-cloud) for the same seasons, by platform, test read of ST_B10.

d. ESA WorldCover 2020 and 2021
   - Source: Planetary Computer STAC, collection `esa-worldcover`.
   - Report: whether both versions load for the test window, and the land cover classes present.

e. Built-up change
   - e1. GHSL GHS-BUILT-S R2023A, epochs 2015, 2020, 2025: confirm open download of the tile covering Lagos and read the test window; repeat the read for 2015 and 2020 on the fringe window.
   - e2. Google Open Buildings 2.5D Temporal: determine whether access requires Google Earth Engine.
   - e3. Open Buildings 2.5D Temporal for Lagos: identify the GeoTIFF tiles in gs://open-buildings-temporal-data/v1 that cover Lagos State, read the fringe window for 2016 and 2023, and report building presence for both years.

f. WorldPop population
   - Source: WorldPop REST API and data.worldpop.org.
   - Report: which Nigeria 100 m population years exist (Global 2000-2020 and the newer Global2 R2025A release), whether a national raster can be downloaded, and a clip of the test window.

For each check record: source URL, access method, resolution, years available, whether it worked, and any errors.
