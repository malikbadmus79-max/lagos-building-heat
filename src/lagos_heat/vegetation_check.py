"""Diagnostic for the NDVI difference between periods (run after vegetation.py, before or
after compositing). Reads the per-date NDVI stored by vegetation.py; downloads only an
OpenStreetMap query for protected areas.

For every processed date it reports the median NDVI over
- all valid Lagos land pixels (non-water cells),
- dense built-up core: pixels in cells with GHSL built-up share > CORE_SHARE in both
  2015 and 2020 (stable impervious surface),
- stable trees: pixels in OSM protected areas where all four WorldCover 2021 10 m
  sub-pixels are tree cover (class 10).
If the period difference appears equally on these stable surfaces, it reflects
acquisition or processing conditions rather than land change.

Outputs
- data/interim/s2_ndvi_check.csv       per date
- data/interim/s2_offsets.csv          baseline and offset for every scene considered
- data/interim/summaries/vegetation_check.json
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from lagos_heat.common import (
    INTERIM, PERIODS, PROCESSED, WORK, DayStack, GridSpec, load_grid, write_summary,
)
from lagos_heat.covariates import osm_bbox, osm_polygons, overpass
from lagos_heat.vegetation import RES

# The requested rule (GHSL built-up > 0.9 in both years) selects no cells: the highest
# cell-level GHSL share is about 0.76. The closest feasible thresholds are used instead
# and the 0.9 count is reported.
CORE_SHARE_REQUESTED = 0.9
CORE_SHARES = (0.6, 0.7)
MIN_AREA_TREE_PX = 200  # 8 ha of 20 m tree pixels
NAMED_SITES = "Lekki Conservation|Omu"


def protected_areas(spec):
    s, w, n, e = osm_bbox(spec)
    q = (f'[out:json][timeout:300];('
         f'way["leisure"="nature_reserve"]({s},{w},{n},{e});'
         f'relation["leisure"="nature_reserve"]({s},{w},{n},{e});'
         f'way["boundary"="protected_area"]({s},{w},{n},{e});'
         f'relation["boundary"="protected_area"]({s},{w},{n},{e});'
         f'way["name"~"{NAMED_SITES}",i]({s},{w},{n},{e});'
         f'relation["name"~"{NAMED_SITES}",i]({s},{w},{n},{e}););out geom tags;')
    gdf = osm_polygons(overpass("protected_areas_v2", q))
    cent = gdf.geometry.centroid.to_crs("EPSG:4326")
    return gdf.assign(lon=cent.x, lat=cent.y)


def tree_mask_20m(gbox):
    """20 m pixels whose four WorldCover 10 m sub-pixels are all tree cover."""
    import rasterio

    with rasterio.open(PROCESSED / "worldcover2021_10m_utm31.tif") as src:
        wc = src.read(1)
        t = src.transform
    assert (t.c, t.f) == (gbox.transform.c, gbox.transform.f), "10 m and 20 m grids not aligned"
    h, w = gbox.shape
    wc = wc[: 2 * h, : 2 * w] == 10
    return wc[0::2, 0::2] & wc[1::2, 0::2] & wc[0::2, 1::2] & wc[1::2, 1::2]


def main():
    from rasterio.features import geometry_mask

    spec = GridSpec.load()
    grid = load_grid()
    ids = grid["cell_id"].to_numpy()
    gbox = spec.geobox(RES)
    cells = spec.pixel_cells(gbox.transform, gbox.shape).astype(np.int32)
    safe = np.maximum(cells, 0)

    landcell = np.zeros(spec.ncells, bool)
    landcell[ids[~grid["water_flag"].to_numpy()]] = True
    masks = {"land": (cells >= 0) & landcell[safe]}

    b = pd.read_csv(INTERIM / "buildings_cells.csv")
    stable = lambda t: (b["ghsl_share_2015"] >= t) & (b["ghsl_share_2020"] >= t)  # noqa: E731
    n_core = {f">{CORE_SHARE_REQUESTED} (requested)": int(
        ((b["ghsl_share_2015"] > CORE_SHARE_REQUESTED) &
         (b["ghsl_share_2020"] > CORE_SHARE_REQUESTED)).sum())}
    for t in CORE_SHARES:
        core_ids = b.loc[stable(t), "cell_id"].to_numpy()
        n_core[f">={t}"] = int(len(core_ids))
        corecell = np.zeros(spec.ncells, bool)
        corecell[core_ids] = True
        masks[f"built_core_ghsl>={t}"] = (cells >= 0) & corecell[safe] & masks["land"]

    trees = tree_mask_20m(gbox)
    areas = protected_areas(spec)
    area_rows = []
    for _, a in areas.iterrows():
        inside = ~geometry_mask([a.geometry], out_shape=gbox.shape, transform=gbox.transform)
        n_tree = int((inside & trees).sum())
        label = a["name"] if isinstance(a["name"], str) else f"unnamed OSM {a.osm_type} {a.osm_id}"
        area_rows.append({"osm_id": a.osm_id, "name": label, "area_ha": a.geometry.area / 1e4,
                          "lon": round(a.lon, 4), "lat": round(a.lat, 4),
                          "tree_px_20m": n_tree})
        if n_tree >= MIN_AREA_TREE_PX:
            masks[f"trees: {label}"] = inside & trees
    area_tbl = pd.DataFrame(area_rows).sort_values("tree_px_20m", ascending=False)

    rows = []
    for period in PERIODS:
        meta = json.loads((WORK / f"s2_{period}_ndvi.json").read_text())
        st = DayStack(f"s2_{period}_ndvi", meta["keys"], tuple(meta["shape"]), 10000)
        for d in st.keys:
            if not st.is_done(d):
                continue
            recs = st.done[d]["records"]
            if not any(r.get("used") for r in recs):
                continue
            sl = st.arr[st.keys.index(d)]
            row = {"period": period, "date": d,
                   "tiles_used": ";".join(r["mgrs_tile"] for r in recs if r.get("used")),
                   "baselines": ";".join(sorted({r["processing_baseline"] for r in recs})),
                   "offsets_applied": ";".join(sorted({str(r.get("boa_add_offset"))
                                                       for r in recs if r.get("used")})),
                   "mean_scene_cloud_cover": float(np.mean([r["cloud_cover"] for r in recs]))}
            for name, m in masks.items():
                v = sl[m]
                v = v[v != DayStack.NODATA]
                row[f"{name} median"] = float(np.median(v)) / 1e4 if v.size else np.nan
                row[f"{name} n_px"] = int(v.size)
            rows.append(row)
    per_date = pd.DataFrame(rows)
    per_date.to_csv(INTERIM / "s2_ndvi_check.csv", index=False)

    offs = []
    for f in sorted((WORK / "s2_metadata").glob("*.json")):
        vals = sorted(set(json.loads(f.read_text()).values()))
        offs.append({"id": f.stem, "metadata_boa_add_offset": vals[0] if vals else None,
                     "metadata_has_offset": bool(vals)})
    offs = pd.DataFrame(offs)
    offs.to_csv(INTERIM / "s2_offsets.csv", index=False)

    med_cols = [c for c in per_date.columns if c.endswith(" median")]
    summ = {
        "core_rule": "GHSL built-up share >= threshold in both 2015 and 2020",
        "n_core_cells": n_core,
        "protected_areas": area_tbl.to_dict("records"),
        "masks_px": {k: int(v.sum()) for k, v in masks.items()},
        "dates_per_period": per_date.groupby("period").size().to_dict(),
        "median_over_dates": {p: g[med_cols].median().round(4).to_dict()
                              for p, g in per_date.groupby("period")},
        "range_over_dates": {p: {c: [float(g[c].min()), float(g[c].max())] for c in med_cols}
                             for p, g in per_date.groupby("period")},
        "scenes_with_metadata_checked": int(len(offs)),
        "scenes_metadata_offset_counts": offs["metadata_boa_add_offset"].value_counts(
            dropna=False).to_dict() if len(offs) else {},
    }
    write_summary("vegetation_check", summ)
    print(per_date[["period", "date", "baselines", "offsets_applied"] + med_cols]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
