# Data sources, licences and attribution

Each source below lists:

- the product used and how it was accessed (from `docs/phase0_results.md` and `docs/phase1_notes.md`)
- the licence as stated by the provider
- the attribution or citation text the provider asks for

Licence pages were opened on 2026-11-03, at the URLs given. Raw downloads are not tracked in the repository (`data/raw/README.md`).

## Google Open Buildings 2.5D Temporal

Use: building presence and height for mid-2016 and mid-2023. Version v1, read from the public bucket `gs://open-buildings-temporal-data/v1`.

Licence: "The data is shared under the Creative Commons Attribution (CC-BY 4.0) license and the Open Data Commons Open Database License (ODbL) v1.0 license. As the user, you can pick which of the two licenses you prefer and use the data under the terms of that license." This project uses CC BY 4.0. Source: https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_Research_open-buildings-temporal_v1

Citation: Wojciech Sirko, Emmanuel Asiedu Brempong, Juliana T. C. Marcos, Abigail Annkah, Abel Korme, Mohammed Alewi Hassen, Krishna Sapkota, Tomer Shekel, Abdoulaye Diack, Sella Nevo, Jason Hickey, John Quinn. High-Resolution Building and Road Detection from Sentinel-2, 2023. https://arxiv.org/abs/2310.11622

## Sentinel-2 Level-2A

Use: NDVI and BSI composites, dry seasons 2016–18 and 2023–25. Accessed through the Microsoft Planetary Computer STAC API (collection `sentinel-2-l2a`).

Licence: Copernicus Sentinel data terms. The Planetary Computer collection links to these terms (https://planetarycomputer.microsoft.com/api/stac/v1/collections/sentinel-2-l2a). Access to the data is free, full and open (European Commission, "Legal notice on the use of Copernicus Sentinel Data and Service Information", https://sentinels.copernicus.eu/documents/247904/690755/Sentinel_Data_Legal_Notice).

Attribution, as required by the legal notice for modified data: "Contains modified Copernicus Sentinel data 2016–2025".

## Landsat Collection 2 Level-2 surface temperature

Use: land surface temperature (ST_B10). Landsat 8 in both periods; Landsat 9 added in the late period. Accessed through the Microsoft Planetary Computer (collection `landsat-c2-l2`).

Licence: public domain. The Planetary Computer collection links to the USGS data policy (https://planetarycomputer.microsoft.com/api/stac/v1/collections/landsat-c2-l2). The USGS citation page could not be opened (timeout).

Citation (DOI given by the collection): Landsat 8-9 OLI/TIRS Collection 2 Level-2, U.S. Geological Survey, https://doi.org/10.5066/P9OGBGM6

## ESA WorldCover 2021 (v200)

Use: baseline water mask and tree-cover reference pixels. Accessed through the Microsoft Planetary Computer (collection `esa-worldcover`).

Licence: CC BY 4.0. "The ESA WorldCover product is provided free of charge, without restriction of use." Source: https://esa-worldcover.org/en/data-access

Attribution: "© ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium".

Citation: Zanaga, D., Van De Kerchove, R., Daems, D., De Keersmaecker, W., Brockmann, C., Kirches, G., Wevers, J., Cartus, O., Santoro, M., Fritz, S., Lesiv, M., Herold, M., Tsendbazar, N.E., Xu, P., Ramoino, F., Arino, O., 2022. ESA WorldCover 10 m 2021 v200. https://doi.org/10.5281/zenodo.7254221

## JRC Global Human Settlement Layer: GHS-BUILT-S and GHS-POP, R2023A

Use:

- GHS-BUILT-S R2023A, 2015 and 2020: built-up surface, used as an alternative building measure
- GHS-POP R2023A, epoch 2020: population, used as a cross-check of exposure

Tile R9_C19, downloaded from the JRC open data server.

Licence: "The GHSL data are provided under the Creative Commons Attribution 4.0 International (CC BY 4.0) license." "Reuse of this data is authorised with proper acknowledgment of the source." Source: https://human-settlement.emergency.copernicus.eu/GHSLhowToCite.php

Citation: Pesaresi, M., Schiavina, M., Politis, P., Freire, S., Krasnodębska, K., Uhl, J. H., … Kemper, T. (2024). Advances on the Global Human Settlement Layer by joint assessment of Earth Observation and population survey data. International Journal of Digital Earth, 17(1). https://doi.org/10.1080/17538947.2024.2390454

Dataset DOIs:

- GHS-BUILT-S R2023A: https://doi.org/10.2905/9f06f36f-4b11-47ec-abb0-4f8b7b1d72ea
- GHS-POP R2023A: https://doi.org/10.2905/2ff68a52-5b5b-4a22-8f40-c41da8332cfe

## WorldPop Global2 R2025A, constrained, 100 m

Use: residents per cell, 2016 and 2023 (modelled estimates).

Licence: Creative Commons Attribution 4.0 International. Source: https://hub.worldpop.org/doi/10.5258/SOTON/WP00839

Citation: Bondarenko M., Priyatikanto R., Tejedor-Garavito N., Zhang W., McKeen T., Cunningham A., Woods T., Hilton J., Cihan D., Nosatiuk B., Brinkhoff T., Tatem A., Sorichetta A.. 2025 Constrained estimates of 2015-2030 total number of people per grid square at a resolution of 3 arc (approximately 100m at the equator) R2025A version v1. Global Demographic Data Project - Funded by The Bill and Melinda Gates Foundation (INV-045237). WorldPop - School of Geography and Environmental Science, University of Southampton. DOI:10.5258/SOTON/WP00839

## OpenStreetMap

Use:

- major roads, industrial land and protected areas, through the Overpass API
- place names, through Nominatim reverse geocoding (Phases 3 and 4)

Licence: Open Data Commons Open Database License (ODbL) 1.0. Source: https://www.openstreetmap.org/copyright

Attribution: the copyright page asks for credit to OpenStreetMap through its attribution notice, a clear statement that the data are available under the ODbL, and optionally a link to the copyright page. This repository uses "© OpenStreetMap contributors, data available under the ODbL (https://www.openstreetmap.org/copyright)".

Nominatim was used under its usage policy (https://operations.osmfoundation.org/policies/nominatim/):

- at most one request per second
- an identifying user agent
- results cached in `data/interim/phase3_nominatim.json` and `data/interim/phase4_nominatim.json`
- single-threaded, small one-off jobs

## geoBoundaries

Use: LGA boundaries of Lagos State (gbOpen NGA ADM2, derived from GRID3).

Licence: CC BY 4.0. "Our license requires an acknowledgement in any products you produce which use this data." Source: https://www.geoboundaries.org/

Citation: Runfola, D. et al. (2020) geoBoundaries: A global database of political administrative boundaries. PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866

## Esri World Imagery Wayback

Use: visual spot-check of 40 cells in Phase 3. The images were viewed only; no tiles are redistributed in this repository.

Attribution, as given by the World Imagery service: "Source: Esri, Vantor, Earthstar Geographics, and the GIS User Community" (https://services.arcgisonline.com/arcgis/rest/services/World_Imagery/MapServer?f=pjson).

Licence: the item states "This work is licensed under the Esri Master License Agreement" (https://www.arcgis.com/sharing/rest/content/items/10df2279f9684e4a9f6a7f08febac2a9?f=json). No Wayback-specific terms page was found.

## NOAA Relative Oceanic Niño Index

Use: context for the El Niño overlap described in the limits.

Source: NOAA Climate Prediction Center, https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/ (accessed 2026-11-03). This is a U.S. government source; no licence statement was recorded.
