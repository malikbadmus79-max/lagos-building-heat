# Building out the breeze: building density, vegetation loss and surface heat in Lagos, 2016-2025

## Question

Lagos has densified and expanded rapidly over the past decade, with new construction replacing vegetated and open land across both the core and the peri-urban fringe. This project asks how changes in building density and vegetation cover between 2016 and 2025 relate to changes in dry-season land surface temperature across the 20 Local Government Areas of Lagos State, and whether the warming associated with densification differs between areas that were already built up and areas that were converted from vegetation during the period.

## Data sources (to be confirmed)

- Administrative boundaries: geoBoundaries gbOpen NGA ADM2 (derived from GRID3 Nigeria LGA boundaries, CC BY 4.0).
- Vegetation and surface reflectance: Sentinel-2 Level-2A, accessed through the Microsoft Planetary Computer STAC API.
- Land surface temperature: Landsat Collection 2 Level-2 surface temperature (ST_B10), Microsoft Planetary Computer.
- Land cover: ESA WorldCover 2020 and 2021, Microsoft Planetary Computer.
- Built-up surface: JRC GHSL GHS-BUILT-S R2023A (2015, 2020, 2025 epochs).
- Building presence and height: Google Open Buildings 2.5D Temporal (2016-2023), via Earth Engine or the public Google Cloud Storage bucket.

Access tests for each source are documented in `docs/phase0_plan.md` and `docs/phase0_results.md`.

## Status

Phase 0: feasibility.

## Repository layout

```
data/raw         raw downloads (not tracked)
data/interim     intermediate products
data/processed   analysis-ready data
docs             plans and results notes
figures          figures
notebooks        exploratory notebooks
src/lagos_heat   project code
tests            tests
```

## Environment

The target environment is Python 3.11, specified in `environment.yml` (conda) and `requirements.txt` (pip). The Phase 0 checks can be reproduced with:

```
conda env create -f environment.yml
conda activate lagos-heat
PYTHONPATH=src python -m lagos_heat.phase0_check
```

## Note on code

Code in this repository is AI-assisted, written with Claude Code.
