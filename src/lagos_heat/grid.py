"""Step 1: 500 m analysis grid over Lagos State with LGA tags and a WorldCover water mask.

Outputs
- data/interim/lagos_lgas.gpkg         full-resolution Lagos LGAs (UTM 31N)
- data/interim/grid_spec.json          grid origin, size and cell size
- data/interim/grid.gpkg               one polygon per cell intersecting Lagos State
- data/processed/worldcover2021_10m_utm31.tif   WorldCover 2021 on a 10 m grid aligned to the cells
- data/interim/summaries/grid.json
"""

from __future__ import annotations

import io

import numpy as np

from lagos_heat.common import (
    CELL, INTERIM, LAGOS_BBOX, PROCESSED, RAW, UTM31N, GridSpec, aggregate, configure_gdal,
    fetch, pc_client, write_raster, write_summary,
)

GB_API = "https://www.geoboundaries.org/api/current/gbOpen/NGA/ADM2/"
WATER_CLASS, TREE_CLASS, BUILT_CLASS = 80, 10, 50
WATER_FLAG_THRESHOLD = 0.5


def lagos_lgas_full_resolution():
    """Full-resolution geoBoundaries ADM2 polygons for the 20 Lagos LGAs.

    Lagos LGAs are identified by shapeID from the Phase 0 selection (simplified ADM2
    units inside the Lagos ADM1 polygon); names alone are ambiguous (e.g. Surulere
    also exists in Oyo State).
    """
    import geopandas as gpd
    import requests

    meta = requests.get(GB_API, timeout=60).json()
    local = RAW / "boundaries" / "geoBoundaries-NGA-ADM2.geojson"
    if not local.exists():
        fetch(meta["gjDownloadURL"], local)
    adm2 = gpd.read_file(local)
    simplified = gpd.read_file(RAW / "lagos_lga_geoboundaries_simplified.geojson")
    lgas = adm2[adm2["shapeID"].isin(simplified["shapeID"])].copy()
    lgas = lgas[["shapeName", "shapeID", "geometry"]].rename(columns={"shapeName": "lga"})
    return lgas.to_crs(UTM31N), meta


def build_grid(lgas):
    import geopandas as gpd
    from shapely.geometry import box

    union = lgas.union_all()
    minx, miny, maxx, maxy = union.bounds
    xmin, ymin = np.floor(minx / CELL) * CELL, np.floor(miny / CELL) * CELL
    xmax, ymax = np.ceil(maxx / CELL) * CELL, np.ceil(maxy / CELL) * CELL
    spec = GridSpec(xmin=float(xmin), ymax=float(ymax),
                    ncol=int((xmax - xmin) / CELL), nrow=int((ymax - ymin) / CELL))

    rows, cols = np.divmod(np.arange(spec.ncells), spec.ncol)
    x0 = spec.xmin + cols * CELL
    y1 = spec.ymax - rows * CELL
    cells = gpd.GeoDataFrame(
        {"cell_id": np.arange(spec.ncells), "row": rows, "col": cols},
        geometry=[box(a, b - CELL, a + CELL, b) for a, b in zip(x0, y1)],
        crs=UTM31N,
    )
    cells = cells[cells.intersects(union)].copy()

    # LGA by largest overlap; share of the cell inside Lagos State
    inter = gpd.overlay(cells[["cell_id", "geometry"]], lgas[["lga", "geometry"]],
                        how="intersection", keep_geom_type=True)
    inter["a"] = inter.area
    best = inter.sort_values("a", ascending=False).drop_duplicates("cell_id")
    cells = cells.merge(best[["cell_id", "lga"]], on="cell_id", how="left")
    in_lagos = inter.groupby("cell_id")["a"].sum() / CELL**2
    cells["frac_in_lagos"] = cells["cell_id"].map(in_lagos).fillna(0.0).clip(upper=1.0)
    return spec, cells


def worldcover_fractions(spec):
    """WorldCover 2021 on a 10 m grid aligned to the cells; per-cell class fractions."""
    import odc.stac

    configure_gdal()
    client = pc_client()
    lon0, lat0, lon1, lat1 = LAGOS_BBOX
    items = [
        i for i in client.search(collections=["esa-worldcover"], bbox=LAGOS_BBOX).items()
        if (i.properties.get("start_datetime") or "").startswith("2021")
    ]
    gbox = spec.geobox(10)
    ds = odc.stac.load(items, bands=["map"], geobox=gbox, resampling="nearest",
                       groupby="solar_day", pool=4)
    wc = ds["map"].isel(time=0).values.astype("uint8")
    out = PROCESSED / "worldcover2021_10m_utm31.tif"
    write_raster(out, wc, gbox.transform, nodata=0)

    cells = spec.pixel_cells(gbox.transform, wc.shape)
    _, n_all = aggregate(cells, ncells=spec.ncells)
    fr = {}
    for name, cls in (("water", WATER_CLASS), ("tree", TREE_CLASS), ("built_wc", BUILT_CLASS)):
        n_cls, _ = aggregate(cells, mask=wc == cls, ncells=spec.ncells)
        fr[name] = np.divide(n_cls, n_all, out=np.full(spec.ncells, np.nan), where=n_all > 0)
    n_nodata, _ = aggregate(cells, mask=wc == 0, ncells=spec.ncells)
    fr["wc_nodata"] = np.divide(n_nodata, n_all, out=np.full(spec.ncells, np.nan),
                                where=n_all > 0)
    return fr, [i.id for i in items], out


def main():
    lgas, meta = lagos_lgas_full_resolution()
    INTERIM.mkdir(parents=True, exist_ok=True)
    lgas.to_file(INTERIM / "lagos_lgas.gpkg", driver="GPKG")
    spec, cells = build_grid(lgas)
    spec.save(INTERIM / "grid_spec.json")

    fr, wc_items, wc_path = worldcover_fractions(spec)
    ids = cells["cell_id"].to_numpy()
    cells["water_frac"] = fr["water"][ids]
    cells["tree_frac_wc2021"] = fr["tree"][ids]
    cells["built_frac_wc2021"] = fr["built_wc"][ids]
    cells["wc_nodata_frac"] = fr["wc_nodata"][ids]
    cells["water_flag"] = cells["water_frac"] > WATER_FLAG_THRESHOLD
    centroids = cells.geometry.centroid.to_crs("EPSG:4326")
    cells["lon"], cells["lat"] = centroids.x, centroids.y
    cells.to_file(INTERIM / "grid.gpkg", driver="GPKG")

    write_summary("grid", {
        "boundary_source": meta.get("gjDownloadURL"),
        "boundary_year": meta.get("boundaryYearRepresented"),
        "n_lgas": int(len(lgas)),
        "lagos_area_km2": float(lgas.area.sum() / 1e6),
        "grid_spec": spec.__dict__,
        "n_cells_bbox": spec.ncells,
        "n_cells_intersecting_lagos": int(len(cells)),
        "n_cells_fully_inside": int((cells["frac_in_lagos"] >= 0.999).sum()),
        "n_cells_water_flag": int(cells["water_flag"].sum()),
        "n_cells_land": int((~cells["water_flag"]).sum()),
        "water_flag_rule": f"WorldCover 2021 class {WATER_CLASS} share > {WATER_FLAG_THRESHOLD}",
        "worldcover_items": wc_items,
        "worldcover_raster": wc_path.relative_to(INTERIM.parents[1]),
        "cells_with_worldcover_nodata": int((cells["wc_nodata_frac"] > 0).sum()),
        "cells_missing_lga": int(cells["lga"].isna().sum()),
        "cells_per_lga": cells.groupby("lga").size().to_dict(),
        "water_flag_cells_per_lga": cells.groupby("lga")["water_flag"].sum().to_dict(),
    })
    print(f"grid: {len(cells)} cells, {int(cells['water_flag'].sum())} water-flagged")


if __name__ == "__main__":
    main()
