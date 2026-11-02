"""Step 4: Landsat C2 L2 land surface temperature composites per period (30 m), aggregated to cells.

Sensors: Landsat 8 in both periods; Landsat 9 (launched 2021) in the late period only
(PLATFORMS). Both use TIRS band 10 with the same Collection 2 scaling.

Per acquisition date (WRS rows of the same path and date are merged):
1. ST_B10 (asset lwir11) is scaled to Kelvin (DN * 0.00341802 + 149.0) and converted to
   degrees Celsius; DN 0 is no data.
2. QA_PIXEL masks fill, dilated cloud, cirrus, cloud, cloud shadow, snow and water.
   Dates with fewer than MIN_VALID_FRACTION of Lagos land pixels valid are dropped.
   The per-date LST is stored at this stage, so the rules below can be re-applied
   without downloading (--composite-only).
3. Rule 1: pixels outside VALID_RANGE_C are masked as residual cloud or error.
4. Rule 2: the Lagos-wide median land LST is computed per date (after rule 1); dates
   more than DATE_SD_RULE sample SDs from the mean of the period's date medians are
   excluded (single pass).
5. Relative LST = pixel LST minus that date's Lagos-wide median land LST.
6. Per-pixel medians of absolute and relative LST per period; valid observation counts;
   lst_ok = at least MIN_OBS valid observations in both periods (flag only).

The 30 m grid is aligned to the Landsat pixel grid (edges at 15 m + k * 30 m), so no
resampling shift is introduced.

Outputs
- data/processed/lst_{period}_{abs,rel}_median_30m.tif, lst_{period}_nobs_30m.tif
- data/interim/temperature_cells.csv
- data/interim/landsat_scenes.csv      every item considered, date-level medians, use
- data/interim/summaries/temperature.json
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from lagos_heat.common import (
    INTERIM, LAGOS_BBOX, MIN_OBS, PERIODS, PROCESSED, SUMMARIES, WORK, DayStack, GridSpec,
    aggregate, configure_gdal, lagos_lgas, load_grid, pc_client, resign, retry, safe_div,
    write_raster, write_summary,
)

RES, ANCHOR = 30, 15
ST_SCALE, ST_OFFSET = 0.00341802, 149.0
CLOUD_MAX = 80  # scene-level eo:cloud_cover filter, percent
MIN_VALID_FRACTION = 0.10  # of Lagos land pixels, per date
# QA_PIXEL bits: 0 fill, 1 dilated cloud, 2 cirrus, 3 cloud, 4 cloud shadow, 5 snow, 7 water
QA_MASK_BITS = (0, 1, 2, 3, 4, 5, 7)
PLATFORMS = {"early": ("landsat-8",), "late": ("landsat-8", "landsat-9")}
# Rules set before inspecting their effect:
VALID_RANGE_C = (15.0, 65.0)  # pixels outside: residual cloud or error
DATE_SD_RULE = 2.0  # exclude dates whose Lagos land median is > 2 SD from the period mean


# ----------------------------------------------------------- pure functions
def lst_celsius(dn):
    """Landsat C2 L2 ST_B10 digital numbers to degrees Celsius; DN 0 -> NaN."""
    dn = np.asarray(dn)
    kelvin = dn.astype("float64") * ST_SCALE + ST_OFFSET
    return np.where(dn == 0, np.nan, kelvin - 273.15).astype("float32")


def qa_clear(qa, bits=QA_MASK_BITS):
    """True where none of the given QA_PIXEL bits is set."""
    mask = 0
    for b in bits:
        mask |= 1 << b
    return (np.asarray(qa).astype("uint16") & mask) == 0


def relative_lst(lst, reference_mask):
    """LST minus the median of LST over reference_mask (valid land pixels); (rel, median)."""
    vals = lst[reference_mask & np.isfinite(lst)]
    if vals.size == 0:
        return np.full(lst.shape, np.nan, dtype="float32"), np.nan
    med = float(np.median(vals))
    return (lst - med).astype("float32"), med


# ---------------------------------------------------------------- loading
def search_items(period):
    client = retry(pc_client)
    items = []
    for season in PERIODS[period]:
        items += retry(lambda: list(client.search(
            collections=["landsat-c2-l2"], bbox=LAGOS_BBOX, datetime=season,
            query={"platform": {"in": list(PLATFORMS[period])}}).items()))
    return sorted([i for i in items if "lwir11" in i.assets], key=lambda i: (i.datetime, i.id))


def load_item(item, gbox):
    """Load ST_B10 and QA_PIXEL; network errors are retried, then raised (never skipped)."""
    import odc.stac

    def _load():
        ds = odc.stac.load([resign(item)], bands=["lwir11", "qa_pixel"], geobox=gbox,
                           resampling="nearest", pool=2, fail_on_error=True)
        return ds["lwir11"].isel(time=0).values, ds["qa_pixel"].isel(time=0).values

    return retry(_load)


def date_key(item):
    """Date and WRS path; Landsat 9 keys carry an "_l9" suffix (Landsat 8 keys are unchanged)."""
    key = f"{item.datetime.date().isoformat()}_p{item.properties.get('landsat:wrs_path')}"
    return key + "_l9" if item.properties.get("platform") == "landsat-9" else key


def download_period(period, gbox, land, records):
    """Download and store per-date LST (QA-masked, before the range and date rules)."""
    n_land = int(land.sum())
    items = search_items(period)
    kept = [i for i in items if i.properties.get("eo:cloud_cover", 100) < CLOUD_MAX]
    kept_ids = {i.id for i in kept}
    for i in items:
        if i.id not in kept_ids:
            records.append({"period": period, "date_key": date_key(i), "ids": i.id,
                            "platform": i.properties.get("platform"),
                            "cloud_cover": i.properties.get("eo:cloud_cover"),
                            "used": False, "note": f"scene cloud cover >= {CLOUD_MAX}%"})
    keys = sorted({date_key(i) for i in kept})
    abs_stack = DayStack(f"lst_{period}_abs", keys, gbox.shape, 100)
    for k, key in enumerate(keys, 1):
        if abs_stack.is_done(key):
            if not str(abs_stack.done[key].get("note", "")).startswith("error"):
                continue
            abs_stack.forget(key)  # written by an earlier run that hit read errors
        t0 = time.time()
        day_items = [i for i in kept if date_key(i) == key]
        lst = np.full(gbox.shape, np.nan, dtype="float32")
        rec = {"period": period, "date_key": key, "date": key[:10],
               "platform": day_items[0].properties.get("platform"),
               "ids": ";".join(i.id for i in day_items),
               "wrs": ";".join(f"{i.properties.get('landsat:wrs_path')}/"
                               f"{i.properties.get('landsat:wrs_row')}" for i in day_items),
               "cloud_cover": ";".join(f"{i.properties.get('eo:cloud_cover'):.1f}"
                                       for i in day_items)}
        for item in day_items:
            dn, qa = load_item(item, gbox)
            t = lst_celsius(dn)
            ok = qa_clear(qa) & np.isfinite(t) & np.isnan(lst)
            lst[ok] = t[ok]
            del dn, qa, t, ok
        frac = float((land & np.isfinite(lst)).sum() / n_land)
        rec.update({"valid_land_fraction": frac, "used": frac >= MIN_VALID_FRACTION})
        if not rec["used"]:
            rec["note"] = f"valid Lagos land fraction < {MIN_VALID_FRACTION}"
            lst[:] = np.nan
        abs_stack.write(key, lst, rec)
        print(f"  LST {period} {key} ({k}/{len(keys)}): used={rec['used']} "
              f"valid={frac:.2f} {time.time() - t0:.0f} s", flush=True)
        del lst
    return {"items_found": len(items), "items_below_cloud_max": len(kept)}


def date_outliers(medians, n_sd=None):
    """Keys whose median lies more than n_sd sample SDs from the mean of all medians."""
    n_sd = DATE_SD_RULE if n_sd is None else n_sd
    vals = np.array([v for v in medians.values() if v is not None], dtype="float64")
    if vals.size < 3:
        return set(), float("nan"), float("nan")
    mean, sd = float(vals.mean()), float(vals.std(ddof=1))
    return ({k for k, v in medians.items() if v is not None and abs(v - mean) > n_sd * sd},
            mean, sd)


def main(composite_only=False):
    from rasterio.features import geometry_mask

    configure_gdal()
    spec = GridSpec.load()
    grid = load_grid()
    ids = grid["cell_id"].to_numpy()
    gbox = spec.geobox(RES, anchor=ANCHOR)
    cells = spec.pixel_cells(gbox.transform, gbox.shape)
    inside = ~geometry_mask([lagos_lgas().union_all()], out_shape=gbox.shape,
                            transform=gbox.transform)
    landcell = np.zeros(spec.ncells, bool)
    landcell[ids[~grid["water_flag"].to_numpy()]] = True
    land = inside & (cells >= 0) & landcell[np.maximum(cells, 0)]
    del inside

    out = pd.DataFrame({"cell_id": ids})
    records, prev = [], {}
    if composite_only:
        # cloud-filtered scene rows and item counts from the run that downloaded the data
        old = pd.read_csv(INTERIM / "landsat_scenes.csv")
        records = old[old["note"].astype(str).str.startswith("scene cloud cover")].to_dict("records")
        prev = json.loads((SUMMARIES / "temperature.json").read_text())
    summ = {"resolution_m": RES, "grid_anchor_m": ANCHOR,
            "scale_kelvin": [ST_SCALE, ST_OFFSET], "cloud_max_scene": CLOUD_MAX,
            "qa_mask_bits": QA_MASK_BITS, "min_valid_fraction": MIN_VALID_FRACTION,
            "valid_range_c": VALID_RANGE_C, "date_sd_rule": DATE_SD_RULE,
            "lagos_land_pixels": int(land.sum()), "min_obs": MIN_OBS}
    vmin, vmax = VALID_RANGE_C
    for period in PERIODS:
        if composite_only:
            counts = {k: prev[period][k] for k in ("items_found", "items_below_cloud_max")}
        else:
            counts = download_period(period, gbox, land, records)
        meta = json.loads((WORK / f"lst_{period}_abs.json").read_text())
        abs_stack = DayStack(f"lst_{period}_abs", meta["keys"], tuple(meta["shape"]), 100)
        keys = abs_stack.keys
        coverage_ok = [k for k in keys if abs_stack.done[k].get("used")]

        # Rule 1: valid range (residual cloud or error), counted per rule on land pixels
        outside = abs_stack.count_outside(coverage_ok, land, vmin, vmax)
        # Rule 2: date-level outliers in the Lagos-wide land median (after rule 1)
        medians = abs_stack.slice_medians(coverage_ok, land, vmin, vmax)
        excluded, mean_med, sd_med = date_outliers(medians)
        used = [k for k in coverage_ok if k not in excluded]

        for k in keys:
            rec = dict(abs_stack.done[k])
            rec.pop("lagos_median_c", None)
            below, above = outside.get(k, (None, None))
            rec.update({"lagos_median_c": medians.get(k), "px_below_15c": below,
                        "px_above_65c": above, "excluded_date_outlier": k in excluded,
                        "used_in_composite": k in used})
            records.append(rec)

        a_med, nobs = abs_stack.median(used, vmin=vmin, vmax=vmax)
        r_med, _ = abs_stack.median(used, offsets=medians, vmin=vmin, vmax=vmax)
        write_raster(PROCESSED / f"lst_{period}_abs_median_30m.tif", a_med, gbox.transform)
        write_raster(PROCESSED / f"lst_{period}_rel_median_30m.tif", r_med, gbox.transform)
        write_raster(PROCESSED / f"lst_{period}_nobs_30m.tif", nobs, gbox.transform, nodata=None)

        has = np.isfinite(a_med)
        n_px, _ = aggregate(cells, ncells=spec.ncells)
        a_sum, n_has = aggregate(cells, values=np.where(has, a_med, 0), mask=has,
                                 ncells=spec.ncells)
        r_sum, _ = aggregate(cells, values=np.where(has, r_med, 0), mask=has,
                             ncells=spec.ncells)
        o_sum, _ = aggregate(cells, values=nobs, mask=has, ncells=spec.ncells)
        out[f"lst_{period}"] = safe_div(a_sum, n_has)[ids]
        out[f"lst_rel_{period}"] = safe_div(r_sum, n_has)[ids]
        # mean valid observations per pixel, over pixels with at least one (0 if none)
        out[f"lst_nobs_{period}"] = np.nan_to_num(safe_div(o_sum, n_has)[ids], nan=0.0)
        out[f"lst_px_share_{period}"] = safe_div(n_has, n_px)[ids]
        out[f"lst_low_obs_{period}"] = ~(out[f"lst_nobs_{period}"] >= MIN_OBS)

        land_rows = ~grid["water_flag"].to_numpy()
        lv = a_med[land & has]
        summ[period] = {
            **counts,
            "dates_considered": len(keys), "dates_coverage_ok": len(coverage_ok),
            "dates_excluded_outlier": sorted(excluded),
            "date_median_mean_c": mean_med, "date_median_sd_c": sd_med,
            "date_medians_c": medians,
            "platforms": list(PLATFORMS[period]),
            "dates_used": len(used), "dates_used_list": used,
            "dates_used_by_platform": {
                "landsat-9": sum(k.endswith("_l9") for k in used),
                "landsat-8": sum(not k.endswith("_l9") for k in used)},
            "scene_ids_used": [abs_stack.done[k]["ids"] for k in used],
            "land_px_obs_removed_below_15c": int(sum(v[0] for v in outside.values())),
            "land_px_obs_removed_above_65c": int(sum(v[1] for v in outside.values())),
            "land_px_obs_removed_below_15c_used_dates": int(sum(outside[k][0] for k in used)),
            "land_px_obs_removed_above_65c_used_dates": int(sum(outside[k][1] for k in used)),
            "nobs_pixel_max": int(nobs.max()),
            "land_cells_low_obs": int(out.loc[land_rows, f"lst_low_obs_{period}"].sum()),
            "land_cells": int(land_rows.sum()),
            "composite_land_px_min_c": float(lv.min()) if lv.size else None,
            "composite_land_px_max_c": float(lv.max()) if lv.size else None,
        }
        del a_med, r_med, nobs

    out["lst_change"] = out["lst_late"] - out["lst_early"]
    out["lst_rel_change"] = out["lst_rel_late"] - out["lst_rel_early"]
    out["lst_ok"] = (out["lst_nobs_early"] >= MIN_OBS) & (out["lst_nobs_late"] >= MIN_OBS)
    out.to_csv(INTERIM / "temperature_cells.csv", index=False)
    pd.DataFrame(records).to_csv(INTERIM / "landsat_scenes.csv", index=False)
    land_rows = ~grid["water_flag"].to_numpy()
    summ["cells_lst_ok"] = int(out["lst_ok"].sum())
    summ["land_cells_lst_ok"] = int(out.loc[land_rows, "lst_ok"].sum())
    write_summary("temperature", summ)
    print("temperature: done")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--composite-only", action="store_true",
                   help="rebuild composites from stored per-date LST (no download)")
    main(p.parse_args().composite_only)
