# Building out the breeze: building density, vegetation loss and surface heat in Lagos, 2016-2025

## Question

This project examines the relationship between built form, vegetation and surface heat in Lagos State. It asks how changes in building density and vegetation cover between 2016 and 2025 relate to changes in dry-season land surface temperature across the 20 Local Government Areas of Lagos State, and whether the warming associated with densification differs between areas that were already built up and areas that were converted from vegetation during the period.

## Data sources (to be confirmed)

- Administrative boundaries: geoBoundaries gbOpen NGA ADM2 (derived from GRID3 Nigeria LGA boundaries, CC BY 4.0).
- Vegetation and surface reflectance: Sentinel-2 Level-2A, accessed through the Microsoft Planetary Computer STAC API.
- Land surface temperature: Landsat Collection 2 Level-2 surface temperature (ST_B10), Microsoft Planetary Computer. Landsat 8 in both periods; Landsat 9 added in the late period.
- Land cover: ESA WorldCover 2020 and 2021, Microsoft Planetary Computer.
- Built-up surface: JRC GHSL GHS-BUILT-S R2023A (2015, 2020, 2025 epochs).
- Building presence and height: Google Open Buildings 2.5D Temporal (2016-2023), via Earth Engine or the public Google Cloud Storage bucket.
- Population: WorldPop Global2 R2025A, constrained, 100 m (2016 and 2023).
- Roads and industrial land use: OpenStreetMap, via the Overpass API.

Access tests for each source are documented in `docs/phase0_plan.md` and `docs/phase0_results.md`. Phase 1 processing choices are documented in `docs/phase1_notes.md`.

## Status

Phase 1: dataset construction.

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

The Phase 1 dataset is built by running the modules in order:

```
for m in grid population buildings temperature vegetation s2_haze vegetation vegetation_check covariates assemble; do
  PYTHONPATH=src python -m lagos_heat.$m
done
python -m pytest tests
```

The second `vegetation` run rebuilds the composites without the dates excluded by the haze screen (`s2_haze`); it reads stored data and downloads nothing.

## Note on code

Code in this repository is AI-assisted, written with Claude Code.
