"""Step 5: WorldPop population per cell and per LGA, 2016 and 2023.

Product: WorldPop Global2 R2025A, constrained, 100 m (alias G2_CN_POP_R25A_100m).
Cell totals: sum of WorldPop pixels whose centres fall in the cell (pixel centres
reprojected to UTM 31N). LGA totals: WorldPop pixels assigned to LGAs by rasterising the
LGA polygons on the WorldPop grid (pixel-centre rule), independent of the 500 m grid.

Outputs
- data/interim/population_cells.csv
- data/interim/population_lga.csv
- data/interim/summaries/population.json
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from lagos_heat.common import (
    INTERIM, RAW, GridSpec, aggregate, fetch, lagos_lgas, load_grid, write_summary,
)

API = "https://hub.worldpop.org/rest/data/pop"
PRODUCT = "G2_CN_POP_R25A_100m"
YEARS = (2016, 2023)


def product_records():
    data = json.loads(fetch(f"{API}/{PRODUCT}", params={"iso3": "NGA"}))
    return {int(r["popyear"]): r for r in data.get("data", [])}


def main():
    import rasterio
    from rasterio.features import rasterize
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds

    spec = GridSpec.load()
    grid = load_grid()
    ids = grid["cell_id"].to_numpy()
    lgas = lagos_lgas()
    recs = product_records()

    cells_out = pd.DataFrame({"cell_id": ids})
    lga_out = pd.DataFrame({"lga": lgas["lga"].to_numpy()})
    summ = {"product": PRODUCT, "constrained": "_CN_" in PRODUCT, "years": {}}
    for year in YEARS:
        rec = recs[year]
        url = rec["files"][0]
        local = RAW / "worldpop" / url.rsplit("/", 1)[-1]
        if not local.exists():
            fetch(url, local)
        with rasterio.open(local) as src:
            b = transform_bounds("EPSG:32631", src.crs, *spec.bounds, densify_pts=50)
            pad = 0.01
            win = from_bounds(b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad,
                              transform=src.transform).round_offsets().round_lengths()
            arr = src.read(1, window=win).astype("float64")
            transform = src.window_transform(win)
            nodata = src.nodata
            crs = src.crs
        ok = np.isfinite(arr) & (arr != nodata) if nodata is not None else np.isfinite(arr)
        pop = np.where(ok, arr, 0.0)
        cells = spec.pixel_cells_reprojected(transform, arr.shape, crs)
        total, _ = aggregate(cells, values=pop, ncells=spec.ncells)
        cells_out[f"pop_{year}"] = total[ids]

        lg = lgas.to_crs(crs)
        labels = rasterize([(g, i + 1) for i, g in enumerate(lg.geometry)], out_shape=arr.shape,
                           transform=transform, fill=0, dtype="int32")
        lga_tot = np.bincount(labels.ravel(), weights=pop.ravel(), minlength=len(lg) + 1)
        lga_out[f"pop_{year}"] = lga_tot[1:]
        summ["years"][year] = {
            "file": url, "doi": rec.get("doi"), "title": rec.get("title"),
            "description": rec.get("desc"),
            "lagos_total_lga_raster": float(lga_tot[1:].sum()),
            "grid_cells_total": float(total[ids].sum()),
            "negative_pixels": int((pop < 0).sum()),
            "cell_max": float(total[ids].max()),
        }
    cells_out["pop_change"] = cells_out["pop_2023"] - cells_out["pop_2016"]
    lga_out["pop_change"] = lga_out["pop_2023"] - lga_out["pop_2016"]
    cells_out.to_csv(INTERIM / "population_cells.csv", index=False)
    lga_out.to_csv(INTERIM / "population_lga.csv", index=False)
    write_summary("population", summ)
    print("population: done")


if __name__ == "__main__":
    main()
