"""Step 3b: haze screen for the Sentinel-2 dates (run after vegetation.py downloads, before
its composites are used). Reads the stored per-date NDVI; downloads only the AOT layers of
the documented, unused AOT screen (cached after the first run).

Primary screen (stable reference forest, the same rule for both periods):
1. Reference forest: the WorldCover 2021 tree pixels (all four 10 m sub-pixels tree cover)
   inside OSM way REFERENCE_OSM_ID, a large protected forest.
2. Per date, the median NDVI over the reference-forest pixels that are valid that date.
   A date is assessable when at least MIN_REF_PX pixels are valid; dates that are not
   assessable are excluded.
3. The reference level is the median of the CLEAREST_N highest date medians, pooled over
   both periods. A date is excluded when its median is more than MAX_DROP below that level.
The screen looks at NDVI on the reference forest only, never on the analysed cells.

Documented failed attempt (reported, not applied): Sen2Cor aerosol optical thickness at
550 nm (asset AOT, DN / 1000), median over clear land on a 60 m grid; a date would be
excluded when its median is more than AOT_SD_RULE sample SDs above the period mean. In the
early period (older Sen2Cor processing) most dates carry one repeated fallback value, so
this indicator does not measure haze there.

vegetation.py reads data/interim/s2_haze_screen.csv and leaves dates with excluded_haze
out of its screened composites.

Outputs
- data/interim/s2_haze_screen.csv     one row per date
- data/interim/summaries/s2_haze.json
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from lagos_heat.common import (
    INTERIM, PERIODS, WORK, DayStack, GridSpec, configure_gdal, load_grid, pc_client, resign,
    retry, write_summary,
)
from lagos_heat.vegetation import RES

# Primary screen, set before inspecting its effect
REFERENCE_OSM_ID = 797720562
MIN_REF_PX = 200  # 8 ha of 20 m pixels
CLEAREST_N = 3
MAX_DROP = 0.10
# Failed attempt, kept for documentation
AOT_RES = 60
AOT_SCALE = 1000.0
AOT_SD_RULE = 2.0
CACHE = WORK / "s2_haze.json"


def reference_forest_mask(spec, gbox):
    """WorldCover tree pixels (all four 10 m sub-pixels) inside the reference forest."""
    from rasterio.features import geometry_mask

    from lagos_heat.vegetation_check import protected_areas, tree_mask_20m

    areas = protected_areas(spec)
    poly = areas.loc[(areas["osm_type"] == "way") & (areas["osm_id"] == REFERENCE_OSM_ID)]
    assert len(poly) == 1, f"reference forest way {REFERENCE_OSM_ID} not found"
    inside = ~geometry_mask(list(poly.geometry), out_shape=gbox.shape, transform=gbox.transform)
    return inside & tree_mask_20m(gbox), float(poly.geometry.area.iloc[0] / 1e4)


def reference_screen(medians, n_px, min_px=MIN_REF_PX, clearest_n=CLEAREST_N,
                     max_drop=MAX_DROP):
    """Reference level, threshold and excluded keys for {key: median} pooled over periods.

    Keys with fewer than min_px valid pixels cannot be assessed and are excluded.
    """
    ok = {k: v for k, v in medians.items() if v is not None and n_px[k] >= min_px}
    top = sorted(ok.values(), reverse=True)[:clearest_n]
    level = float(np.median(top))
    threshold = level - max_drop
    excluded = {k for k in medians if k not in ok} | {k for k, v in ok.items() if v < threshold}
    return level, threshold, excluded, set(ok)


def hazy_dates(medians, n_sd=AOT_SD_RULE):
    """Keys whose median exceeds the mean of all medians by more than n_sd sample SDs."""
    vals = np.array([v for v in medians.values() if v is not None], dtype="float64")
    if vals.size < 3:
        return set(), float("nan"), float("nan")
    mean, sd = float(vals.mean()), float(vals.std(ddof=1))
    return {k for k, v in medians.items() if v is not None and v > mean + n_sd * sd}, mean, sd


def load_aot(item, gbox):
    import odc.stac

    def _load():
        ds = odc.stac.load([resign(item)], bands=["AOT"], geobox=gbox, resampling="nearest",
                           pool=2, fail_on_error=True)
        return ds["AOT"].isel(time=0).values.astype("float32")

    return retry(_load)


def aot_medians(st, period, dates, land, gbox, cache, client):
    """Median AOT over clear land per date (cached in WORK/s2_haze.json)."""
    from affine import Affine
    from odc.geo.geobox import GeoBox

    f = AOT_RES // RES
    h, w = gbox.shape[0] // f, gbox.shape[1] // f
    t = gbox.transform
    gbox60 = GeoBox((h, w), Affine(AOT_RES, 0, t.c, 0, -AOT_RES, t.f), gbox.crs)
    centre = (slice(f // 2, f // 2 + f * h, f), slice(f // 2, f // 2 + f * w, f))
    land60 = land[centre]
    for k, d in enumerate(dates, 1):
        ckey = f"{period}/{d}"
        if ckey in cache:
            continue
        t0 = time.time()
        used_ids = [r["id"] for r in st.done[d]["records"] if r.get("used")]
        items = retry(lambda: list(client.search(collections=["sentinel-2-l2a"],
                                                 ids=used_ids).items()))
        items = sorted(items, key=lambda i: used_ids.index(i.id))
        aot = np.full((h, w), np.nan, dtype="float32")
        for item in items:
            a = load_aot(item, gbox60)
            ok = (a > 0) & np.isnan(aot)
            aot[ok] = a[ok] / AOT_SCALE
            del a
        sl = st.arr[st.keys.index(d)]
        clear60 = land60 & (sl[centre] != DayStack.NODATA) & np.isfinite(aot)
        cache[ckey] = {
            "aot_median": float(np.median(aot[clear60])) if clear60.any() else None,
            "aot_p90": float(np.percentile(aot[clear60], 90)) if clear60.any() else None,
            "aot_n_px": int(clear60.sum()), "items": used_ids}
        CACHE.write_text(json.dumps(cache))
        print(f"  AOT {period} {d} ({k}/{len(dates)}): median "
              f"{cache[ckey]['aot_median']} {time.time() - t0:.0f} s", flush=True)
    return cache


def main():
    configure_gdal()
    spec = GridSpec.load()
    grid = load_grid()
    ids = grid["cell_id"].to_numpy()
    gbox = spec.geobox(RES)
    cells = spec.pixel_cells(gbox.transform, gbox.shape)
    landcell = np.zeros(spec.ncells, bool)
    landcell[ids[~grid["water_flag"].to_numpy()]] = True
    land = (cells >= 0) & landcell[np.maximum(cells, 0)]
    del cells
    forest, forest_ha = reference_forest_mask(spec, gbox)

    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    client = None
    rows = []
    for period in PERIODS:
        meta = json.loads((WORK / f"s2_{period}_ndvi.json").read_text())
        st = DayStack(f"s2_{period}_ndvi", meta["keys"], tuple(meta["shape"]), 10000)
        dates = [d for d in st.keys if st.is_done(d)
                 and any(r.get("used") for r in st.done[d]["records"])]
        if any(f"{period}/{d}" not in cache for d in dates):
            client = client or retry(pc_client)
        cache = aot_medians(st, period, dates, land, gbox, cache, client)
        for d in dates:
            fv = st.arr[st.keys.index(d)][forest]
            fv = fv[fv != DayStack.NODATA]
            c = cache[f"{period}/{d}"]
            rows.append({"period": period, "date": d,
                         "forest_ndvi_median": float(np.median(fv)) / 1e4 if fv.size else None,
                         "forest_n_px": int(fv.size),
                         "land_ndvi_median": next(r.get("date_land_ndvi_median")
                                                  for r in st.done[d]["records"]),
                         "aot_median": c["aot_median"], "aot_p90": c["aot_p90"],
                         "aot_n_px": c["aot_n_px"]})
    out = pd.DataFrame(rows)
    key = out["period"] + "/" + out["date"]

    level, threshold, excluded, assessable = reference_screen(
        dict(zip(key, out["forest_ndvi_median"])), dict(zip(key, out["forest_n_px"])))
    out["assessable"] = key.isin(assessable)
    out["excluded_haze"] = key.isin(excluded)
    clearest = out[out["assessable"]].nlargest(CLEAREST_N, "forest_ndvi_median")

    summ = {"primary_screen": {
        "rule": (f"exclude a date if its median NDVI over the reference forest is more than "
                 f"{MAX_DROP} below the median of the {CLEAREST_N} highest date medians, "
                 f"pooled over both periods"),
        "reference_osm_way": REFERENCE_OSM_ID, "reference_area_ha": forest_ha,
        "reference_tree_px_20m": int(forest.sum()), "min_valid_px": MIN_REF_PX,
        "clearest_dates": (clearest["period"] + "/" + clearest["date"]).tolist(),
        "clearest_medians": clearest["forest_ndvi_median"].round(4).tolist(),
        "reference_level": level, "threshold": threshold}}
    summ["aot_attempt"] = {"indicator": "Sen2Cor AOT 550 nm (DN / 1000), median over clear land",
                           "aot_res_m": AOT_RES, "aot_sd_rule": AOT_SD_RULE, "applied": False}
    out["excluded_aot_attempt"] = False
    for period, g in out.groupby("period", sort=False):
        summ["primary_screen"][period] = {
            "dates": len(g), "dates_assessable": int(g["assessable"].sum()),
            "dates_not_assessable_excluded": g.loc[~g["assessable"], "date"].tolist(),
            "dates_below_threshold_excluded": g.loc[g["assessable"] & g["excluded_haze"],
                                                    "date"].tolist(),
            "dates_kept_list": g.loc[~g["excluded_haze"], "date"].tolist(),
            "dates_kept": int((~g["excluded_haze"]).sum())}
        medians = dict(zip(g["date"], g["aot_median"]))
        ex, mean, sd = hazy_dates(medians)
        out.loc[g.index, "excluded_aot_attempt"] = g["date"].isin(ex)
        mode = pd.Series(medians).round(3).value_counts()
        summ["aot_attempt"][period] = {
            "aot_date_median_mean": mean, "aot_date_median_sd": sd,
            "aot_threshold": mean + AOT_SD_RULE * sd, "dates_would_exclude": sorted(ex),
            # repeated identical date medians indicate a Sen2Cor fallback, not a retrieval
            "aot_most_common_date_median": float(mode.index[0]),
            "aot_dates_at_most_common_median": int(mode.iloc[0]),
            "aot_dates_p90_equal_median": int((g["aot_p90"] == g["aot_median"]).sum()),
            "spearman_aot_vs_land_ndvi": float(g["aot_median"].corr(g["land_ndvi_median"],
                                                                    method="spearman")),
            "spearman_aot_vs_forest_ndvi": float(g["aot_median"].corr(
                g["forest_ndvi_median"], method="spearman"))}
    out.to_csv(INTERIM / "s2_haze_screen.csv", index=False)
    write_summary("s2_haze", summ)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
