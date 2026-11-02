"""Step 3: Sentinel-2 L2A median NDVI and BSI composites per period (20 m), aggregated to cells.

Processing per acquisition date (tiles of the same date are merged, first valid pixel wins):
1. SCL is read first; pixels are kept only if SCL is in VALID_SCL. If an item has fewer
   than MIN_VALID_FRACTION of Lagos land pixels clear, its reflectance bands are not read.
2. B02, B04, B08, B11 are read at 20 m and converted to reflectance with the
   processing-baseline offset (see s2_reflectance). For the 10 m bands (B02, B04, B08)
   odc-stac reads the COG's internal 2x overview rather than averaging the 10 m pixels;
   B11 is native 20 m.
3. Pixels with any band reflectance <= 0 after the offset are treated as invalid.
4. NDVI and BSI are written to a disk-backed day stack; per-pixel medians and valid
   observation counts are computed per period over the dates kept by the haze screen
   (s2_haze.py; all dates if the screen has not run).
5. Comparison versions over all dates, without the haze screen: unscreened medians of
   NDVI, relative NDVI and BSI (columns *_noscreen), and the per-pixel 80th percentile of
   NDVI and of relative NDVI (columns *_p80). The p80 version is the primary vegetation
   measure for levels within a period; relative p80 NDVI is the primary measure of change
   (docs/phase1_notes.md). The screened and unscreened medians are sensitivity versions.
6. Finalize (also runnable alone with --finalize, from the written rasters): vegetated
   share (NDVI > NDVI_VEGETATED) from the p80 NDVI, bare share (BSI > BSI_BARE) and the
   approximate sand-filling flag from the unscreened median BSI; the same shares from the
   screened medians in *_screened columns; reference-forest values of every NDVI version.

Outputs
- data/processed/s2_{period}_{ndvi,bsi}_median_20m.tif, s2_{period}_nobs_20m.tif
- data/interim/vegetation_cells.csv
- data/interim/s2_items.csv            every item considered, with offset and use
- data/interim/summaries/vegetation.json
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from lagos_heat.common import (
    INTERIM, LAGOS_BBOX, MIN_OBS, PERIODS, PROCESSED, WORK, DayStack, GridSpec, aggregate,
    configure_gdal, fetch, load_grid, pc_client, resign, retry, safe_div, write_raster,
    write_summary,
)

RES = 20
CLOUD_MAX = 60  # scene-level eo:cloud_cover filter, percent (reduces download volume)
# SCL classes kept: 2 dark area pixels, 4 vegetation, 5 not vegetated, 7 unclassified.
# Masked: 0 no data, 1 saturated/defective, 3 cloud shadow, 6 water, 8 cloud medium
# probability, 9 cloud high probability, 10 thin cirrus, 11 snow/ice.
VALID_SCL = (2, 4, 5, 7)
MIN_VALID_FRACTION = 0.01  # of Lagos land pixels, per item
BANDS = ("B02", "B04", "B08", "B11")
NDVI_VEGETATED = 0.4
BSI_BARE = 0.0  # pixels with BSI > 0 counted as bare-dominated
SANDFILL_BARE_SHARE_RISE = 0.25  # candidate sand-filling flag: bare share rises by >= 25 pp
HAZE_SCREEN = INTERIM / "s2_haze_screen.csv"  # written by s2_haze.py


# ----------------------------------------------------------- pure functions
def baseline_offset(processing_baseline):
    """BOA_ADD_OFFSET for a Sentinel-2 L2A processing baseline string such as '04.00'.

    ESA introduced an additive offset of -1000 (DN) with processing baseline 04.00
    (products from 25 January 2022). Earlier baselines have no offset.
    """
    if processing_baseline is None:
        raise ValueError("processing baseline missing")
    return -1000 if float(processing_baseline) >= 4.0 else 0


def s2_reflectance(dn, processing_baseline=None, offset=None):
    """Surface reflectance from L2A digital numbers: (DN + offset) / 10000; DN 0 -> NaN.

    `offset` (BOA_ADD_OFFSET from the product metadata) takes precedence; otherwise it is
    derived from the processing baseline.
    """
    dn = np.asarray(dn)
    if offset is None:
        offset = baseline_offset(processing_baseline)
    refl = (dn.astype("float32") + offset) / 10000.0
    return np.where(dn == 0, np.nan, refl).astype("float32")


def parse_boa_offsets(xml_text):
    """BOA_ADD_OFFSET values ({band_id: offset}) from an L2A MTD_MSIL2A.xml; {} if absent."""
    import re

    return {int(b): int(float(v)) for b, v in re.findall(
        r'<BOA_ADD_OFFSET[^>]*band_id="(\d+)"[^>]*>\s*(-?[\d.]+)\s*</BOA_ADD_OFFSET>', xml_text)}


def resolve_offset(baseline, metadata_offsets):
    """Offset to apply for one scene, with its source; metadata wins over the baseline rule."""
    from_baseline = baseline_offset(baseline)
    if metadata_offsets:
        values = set(metadata_offsets.values())
        if len(values) != 1:
            raise ValueError(f"band-dependent BOA_ADD_OFFSET not handled: {metadata_offsets}")
        applied, source = values.pop(), "metadata BOA_ADD_OFFSET"
    else:
        applied, source = from_baseline, "processing baseline (no BOA_ADD_OFFSET in metadata)"
    return {"offset_from_baseline": from_baseline,
            "offset_from_metadata": (sorted(set(metadata_offsets.values()))[0]
                                     if metadata_offsets else None),
            "boa_add_offset": applied, "offset_source": source,
            "offset_sources_agree": applied == from_baseline}


def scene_offset(item):
    """Resolve the offset for an item from its own product metadata (cached on disk)."""
    import json

    cache = WORK / "s2_metadata" / f"{item.id}.json"
    if cache.exists():
        meta = {int(k): v for k, v in json.loads(cache.read_text()).items()}
    else:
        href = item.assets["product-metadata"].href
        meta = parse_boa_offsets(fetch(href).decode("utf-8", errors="replace"))
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(meta))
    return resolve_offset(item.properties.get("s2:processing_baseline"), meta)


def ndvi(nir, red):
    nir, red = np.asarray(nir, "float32"), np.asarray(red, "float32")
    den = nir + red
    out = np.full(den.shape, np.nan, dtype="float32")
    np.divide(nir - red, den, out=out, where=den > 0)
    return out


def bsi(blue, red, nir, swir):
    """Bare soil index: ((SWIR + red) - (NIR + blue)) / ((SWIR + red) + (NIR + blue))."""
    blue, red, nir, swir = (np.asarray(a, "float32") for a in (blue, red, nir, swir))
    a, b = swir + red, nir + blue
    den = a + b
    out = np.full(den.shape, np.nan, dtype="float32")
    np.divide(a - b, den, out=out, where=den > 0)
    return out


def haze_excluded(period):
    """Dates excluded by the haze screen (s2_haze.py); empty if the screen has not run."""
    if not HAZE_SCREEN.exists():
        return set()
    t = pd.read_csv(HAZE_SCREEN)
    return set(t.loc[(t["period"] == period) & t["excluded_haze"], "date"])


# ---------------------------------------------------------------- loading
def search_items(period):
    client = retry(pc_client)
    items = []
    for season in PERIODS[period]:
        items += retry(lambda: list(client.search(collections=["sentinel-2-l2a"],
                                                  bbox=LAGOS_BBOX, datetime=season).items()))
    return sorted(items, key=lambda i: (i.datetime, i.id))


def load_bands(item, gbox, bands, resampling):
    """Load bands for one item; network errors are retried, then raised (never skipped)."""
    import odc.stac

    def _load():
        ds = odc.stac.load([resign(item)], bands=list(bands), geobox=gbox, resampling=resampling,
                           pool=4, fail_on_error=True)
        return {b: ds[b].isel(time=0).values for b in bands}

    return retry(_load)


def item_window(item, gbox):
    """Row/column slices of gbox covered by the item's footprint (None if no overlap)."""
    from pyproj import Transformer
    from shapely.geometry import shape

    tr = Transformer.from_crs("EPSG:4326", str(gbox.crs), always_xy=True)
    x0, y0, x1, y1 = shape(item.geometry).bounds
    xs, ys = tr.transform([x0, x1, x0, x1], [y0, y0, y1, y1])
    t = gbox.transform
    c0 = max(int(np.floor((min(xs) - t.c) / t.a)), 0)
    c1 = min(int(np.ceil((max(xs) - t.c) / t.a)), gbox.shape[1])
    r0 = max(int(np.floor((t.f - max(ys)) / -t.e)), 0)
    r1 = min(int(np.ceil((t.f - min(ys)) / -t.e)), gbox.shape[0])
    if c0 >= c1 or r0 >= r1:
        return None
    return slice(r0, r1), slice(c0, c1)


def process_date(items, gbox, land):
    """NDVI and BSI for one date (all tiles), plus per-item records.

    Each item is read only over its own footprint window to keep memory low.
    """
    nd = np.full(gbox.shape, np.nan, dtype="float32")
    bs = np.full(gbox.shape, np.nan, dtype="float32")
    n_land = int(land.sum())
    records = []
    for item in items:
        rec = {"id": item.id, "date": item.datetime.date().isoformat(),
               "mgrs_tile": item.properties.get("s2:mgrs_tile"),
               "cloud_cover": item.properties.get("eo:cloud_cover"),
               "processing_baseline": item.properties.get("s2:processing_baseline")}
        rec.update(scene_offset(item))
        fw = item_window(item, gbox)
        if fw is None:
            rec.update({"used": False, "note": "footprint outside grid"})
            records.append(rec)
            continue
        scl = load_bands(item, gbox[fw], ["SCL"], "nearest")["SCL"]
        valid_fw = np.isin(scl, VALID_SCL) & land[fw] & np.isnan(nd[fw])
        rec["clear_land_fraction"] = float(valid_fw.sum() / n_land)
        del scl
        if rec["clear_land_fraction"] < MIN_VALID_FRACTION:
            rec["used"] = False
            rec["note"] = "too few clear land pixels; bands not read"
            records.append(rec)
            continue
        rows, cols = np.where(valid_fw.any(axis=1))[0], np.where(valid_fw.any(axis=0))[0]
        r0, c0 = fw[0].start, fw[1].start
        win = (slice(r0 + rows[0], r0 + rows[-1] + 1), slice(c0 + cols[0], c0 + cols[-1] + 1))
        v = valid_fw[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
        del valid_fw
        raw = load_bands(item, gbox[win], BANDS, {"*": "average"})
        refl = {}
        for b in BANDS:  # convert band by band to limit temporaries
            refl[b] = s2_reflectance(raw.pop(b), offset=rec["boa_add_offset"])
        nonpos = np.zeros(v.shape, bool)
        for b in BANDS:
            nonpos |= ~(refl[b] > 0)
        rec["masked_nonpositive_px"] = int((v & nonpos).sum())
        v &= ~nonpos
        del nonpos
        nd_item = ndvi(refl["B08"], refl["B04"])
        np.copyto(nd[win], nd_item, where=v)
        del nd_item
        bs_item = bsi(refl["B02"], refl["B04"], refl["B08"], refl["B11"])
        np.copyto(bs[win], bs_item, where=v)
        del bs_item, refl
        rec["used"] = True
        rec["valid_land_px"] = int(v.sum())
        records.append(rec)
    return nd, bs, records


def cell_mean(cells, arr, ids, ncells):
    """Mean of the finite pixels of arr per cell, for the cells in ids."""
    has = np.isfinite(arr)
    total, n = aggregate(cells, values=np.where(has, arr, 0), mask=has, ncells=ncells)
    return safe_div(total, n)[ids]


# ------------------------------------------------------------------- main
def main(limit_dates=None):
    configure_gdal()
    spec = GridSpec.load()
    grid = load_grid()
    ids = grid["cell_id"].to_numpy()
    gbox = spec.geobox(RES)
    cells = spec.pixel_cells(gbox.transform, gbox.shape)
    landcell = np.zeros(spec.ncells, bool)
    landcell[ids[~grid["water_flag"].to_numpy()]] = True
    land = (cells >= 0) & landcell[np.maximum(cells, 0)]
    del cells  # rebuilt for the composites; not held during downloads
    out = pd.DataFrame({"cell_id": ids})
    all_records, summ = [], {"resolution_m": RES, "cloud_max_scene": CLOUD_MAX,
                             "valid_scl": VALID_SCL, "min_valid_fraction": MIN_VALID_FRACTION,
                             "periods": PERIODS, "min_obs": MIN_OBS}
    for period in PERIODS:
        items = search_items(period)
        kept = [i for i in items if i.properties.get("eo:cloud_cover", 100) < CLOUD_MAX]
        kept_ids = {i.id for i in kept}
        offsets = {i.id: retry(scene_offset, i) for i in items}  # every scene, logged
        for i in items:
            if i.id not in kept_ids:
                all_records.append({"period": period, "id": i.id,
                                    "date": i.datetime.date().isoformat(),
                                    "mgrs_tile": i.properties.get("s2:mgrs_tile"),
                                    "cloud_cover": i.properties.get("eo:cloud_cover"),
                                    "processing_baseline": i.properties.get("s2:processing_baseline"),
                                    **offsets[i.id],
                                    "used": False, "note": f"scene cloud cover >= {CLOUD_MAX}%"})
        dates = sorted({i.datetime.date().isoformat() for i in kept})
        nd_stack = DayStack(f"s2_{period}_ndvi", dates, gbox.shape, 10000)
        bs_stack = DayStack(f"s2_{period}_bsi", dates, gbox.shape, 10000)
        pending = []
        for d in dates:
            if nd_stack.is_done(d) and bs_stack.is_done(d):
                done_recs = nd_stack.done[d]["records"]
                errors = any(str(r.get("note", "")).startswith("error") for r in done_recs)
                stale = any(r.get("used") and r.get("boa_add_offset") !=
                            offsets[r["id"]]["boa_add_offset"] for r in done_recs)
                if not errors and not stale:
                    for r in done_recs:  # refresh offset fields from the metadata check
                        r.update(offsets[r["id"]])
                    all_records += done_recs
                    continue
                # written by an earlier run that hit read errors or used another offset
                nd_stack.forget(d)
                bs_stack.forget(d)
            pending.append(d)
        if limit_dates is not None:
            # diagnostic mode: the clearest pending dates only, no composites
            mean_cc = {d: np.mean([i.properties["eo:cloud_cover"] for i in kept
                                   if i.datetime.date().isoformat() == d]) for d in pending}
            pending = sorted(sorted(pending, key=mean_cc.get)[:limit_dates])
        for k, d in enumerate(pending, 1):
            t0 = time.time()
            nd, bs, recs = process_date([i for i in kept if i.datetime.date().isoformat() == d],
                                        gbox, land)
            # Date-level median NDVI over valid land pixels: a check on haze or
            # calibration differences between dates (not used in the composites).
            vals = nd[land & np.isfinite(nd)]
            date_med = float(np.median(vals)) if vals.size else None
            recs = [{"period": period, "date_land_ndvi_median": date_med,
                     "date_valid_land_fraction": float(vals.size / land.sum()), **r}
                    for r in recs]
            del vals
            nd_stack.write(d, nd, {"records": recs})
            bs_stack.write(d, bs, {})
            all_records += recs
            print(f"  S2 {period} {d} ({k}/{len(pending)} pending): "
                  f"{sum(r['used'] for r in recs)}/{len(recs)} items used, "
                  f"{time.time() - t0:.0f} s", flush=True)
            del nd, bs
        if limit_dates is not None:
            continue

        all_dates = [d for d in dates
                     if any(r["used"] for r in nd_stack.done[d]["records"])]
        hazy = haze_excluded(period)
        used_dates = [d for d in all_dates if d not in hazy]
        nd_med, nobs = nd_stack.median(used_dates)
        bs_med, _ = bs_stack.median(used_dates)
        # Relative NDVI: each date's NDVI minus that date's median over valid Lagos land
        # pixels, then the per-pixel median over dates (same construction as relative LST).
        land_medians = nd_stack.slice_medians(used_dates, land)
        rel_med, _ = nd_stack.median(used_dates, offsets=land_medians)
        write_raster(PROCESSED / f"s2_{period}_ndvi_median_20m.tif", nd_med, gbox.transform)
        write_raster(PROCESSED / f"s2_{period}_ndvi_rel_median_20m.tif", rel_med, gbox.transform)
        write_raster(PROCESSED / f"s2_{period}_bsi_median_20m.tif", bs_med, gbox.transform)
        write_raster(PROCESSED / f"s2_{period}_nobs_20m.tif", nobs.astype("uint16"),
                     gbox.transform, nodata=None)

        cells = spec.pixel_cells(gbox.transform, gbox.shape)
        has = np.isfinite(nd_med)
        n_px, _ = aggregate(cells, ncells=spec.ncells)
        nd_sum, n_has = aggregate(cells, values=np.where(has, nd_med, 0), mask=has,
                                  ncells=spec.ncells)
        bs_has = np.isfinite(bs_med)
        bs_sum, n_bs = aggregate(cells, values=np.where(bs_has, bs_med, 0), mask=bs_has,
                                 ncells=spec.ncells)
        obs_sum, _ = aggregate(cells, values=nobs, mask=has, ncells=spec.ncells)
        rel_sum, _ = aggregate(cells, values=np.where(has, rel_med, 0), mask=has,
                               ncells=spec.ncells)
        out[f"ndvi_{period}"] = safe_div(nd_sum, n_has)[ids]
        out[f"ndvi_rel_{period}"] = safe_div(rel_sum, n_has)[ids]
        out[f"bsi_{period}"] = safe_div(bs_sum, n_bs)[ids]
        out[f"s2_nobs_{period}"] = safe_div(obs_sum, n_has)[ids]
        out[f"s2_px_share_{period}"] = safe_div(n_has, n_px)[ids]
        out[f"s2_low_obs_{period}"] = ~(out[f"s2_nobs_{period}"] >= MIN_OBS)

        recs = [r for r in all_records if r["period"] == period]
        land_rows = ~grid["water_flag"].to_numpy()
        summ[period] = {
            "items_found": len(items), "items_below_cloud_max": len(kept),
            "items_used": sum(r.get("used", False) for r in recs),
            "dates_considered": len(dates), "dates_used": len(used_dates),
            "dates_used_noscreen_p80": len(all_dates),
            "dates_excluded_haze": sorted(hazy), "haze_screen_applied": HAZE_SCREEN.exists(),
            "items_by_baseline": pd.Series([i.properties.get("s2:processing_baseline")
                                            for i in kept]).value_counts().to_dict(),
            "items_offset_sources_disagree": int(sum(
                r.get("offset_sources_agree") is False for r in recs)),
            "items_offset_from_metadata": int(sum(
                str(r.get("offset_source", "")).startswith("metadata") for r in recs)),
            "offsets_all_items": pd.Series([r.get("boa_add_offset") for r in recs
                                            ]).value_counts().to_dict(),
            "offsets_applied": pd.Series([r.get("boa_add_offset") for r in recs
                                          if r.get("used")]).value_counts().to_dict(),
            "item_errors": [r["note"] for r in recs if str(r.get("note", "")).startswith("error")],
            "nobs_pixel_max": int(nobs.max()),
            "land_cells_low_obs": int(out.loc[land_rows, f"s2_low_obs_{period}"].sum()),
            "land_cells": int(land_rows.sum()),
            "ndvi_cell_min": float(np.nanmin(out[f"ndvi_{period}"])),
            "ndvi_cell_max": float(np.nanmax(out[f"ndvi_{period}"])),
            "ndvi_pixel_min": float(np.nanmin(nd_med)), "ndvi_pixel_max": float(np.nanmax(nd_med)),
            "bsi_pixel_min": float(np.nanmin(bs_med)), "bsi_pixel_max": float(np.nanmax(bs_med)),
        }
        summ[period]["date_land_ndvi_medians"] = land_medians
        del nd_med, bs_med, nobs, rel_med

        # Comparison versions over all dates: unscreened medians (_noscreen) and the
        # 80th-percentile NDVI (_p80), which favours the clearer dates per pixel. Each is
        # computed, written and freed in turn to keep memory low.
        ns_land = nd_stack.slice_medians(all_dates, land)
        versions = (
            ("ndvi_median_noscreen", "ndvi_{}_noscreen",
             lambda: nd_stack.median(all_dates)),
            ("ndvi_rel_median_noscreen", "ndvi_rel_{}_noscreen",
             lambda: nd_stack.median(all_dates, offsets=ns_land)),
            ("bsi_median_noscreen", "bsi_{}_noscreen", lambda: bs_stack.median(all_dates)),
            ("ndvi_p80", "ndvi_{}_p80", lambda: nd_stack.median(all_dates, strip=16, q=80)),
            ("ndvi_rel_p80", "ndvi_rel_{}_p80",
             lambda: nd_stack.median(all_dates, strip=16, offsets=ns_land, q=80)),
        )
        for raster, col, compute in versions:
            arr, n = compute()
            write_raster(PROCESSED / f"s2_{period}_{raster}_20m.tif", arr, gbox.transform)
            out[col.format(period)] = cell_mean(cells, arr, ids, spec.ncells)
            if raster == "ndvi_median_noscreen":
                out[f"s2_nobs_{period}_noscreen"] = cell_mean(
                    cells, np.where(np.isfinite(arr), n, np.nan), ids, spec.ncells)
            del arr, n

    if limit_dates is not None:
        print(f"vegetation: diagnostic run ({limit_dates} dates per period), no composites")
        return
    out["ndvi_change"] = out["ndvi_late"] - out["ndvi_early"]
    out["ndvi_rel_change"] = out["ndvi_rel_late"] - out["ndvi_rel_early"]
    out["bsi_change"] = out["bsi_late"] - out["bsi_early"]
    # Absolute p80 NDVI is a within-period level only; its change is not written.
    out["ndvi_change_noscreen"] = out["ndvi_late_noscreen"] - out["ndvi_early_noscreen"]
    for v in ("noscreen", "p80"):
        out[f"ndvi_rel_change_{v}"] = out[f"ndvi_rel_late_{v}"] - out[f"ndvi_rel_early_{v}"]
    out["bsi_change_noscreen"] = out["bsi_late_noscreen"] - out["bsi_early_noscreen"]
    out["tree_frac_wc2021"] = grid["tree_frac_wc2021"].to_numpy()
    out.to_csv(INTERIM / "vegetation_cells.csv", index=False)
    pd.DataFrame(all_records).to_csv(INTERIM / "s2_items.csv", index=False)
    write_summary("vegetation", summ)
    finalize()
    print("vegetation: done")


# Rasters read by finalize: (summary key, raster name). The forest value of each NDVI
# version is the median over the reference-forest tree pixels (s2_haze.py).
FOREST_VERSIONS = (
    ("p80", "ndvi_p80"), ("rel_p80", "ndvi_rel_p80"),
    ("screened", "ndvi_median"), ("rel_screened", "ndvi_rel_median"),
    ("noscreen", "ndvi_median_noscreen"), ("rel_noscreen", "ndvi_rel_median_noscreen"),
)
# Cell shares: (column stem, column suffix, raster name, pixel threshold)
SHARES = (
    ("veg_share", "", "ndvi_p80", NDVI_VEGETATED),
    ("veg_share", "_screened", "ndvi_median", NDVI_VEGETATED),
    ("bare_share", "", "bsi_median_noscreen", BSI_BARE),
    ("bare_share", "_screened", "bsi_median", BSI_BARE),
)


def read_composite(period, raster):
    import rasterio

    with rasterio.open(PROCESSED / f"s2_{period}_{raster}_20m.tif") as src:
        return src.read(1).astype("float32")


def share_above(cells, arr, threshold, ids, ncells):
    """Per cell, the share of valid pixels whose value exceeds threshold (cells < 0: outside)."""
    has = np.isfinite(arr) & (cells >= 0)
    n_has = np.bincount(cells[has], minlength=ncells).astype("float64")
    n_hit = np.bincount(cells[has & (arr > threshold)], minlength=ncells).astype("float64")
    return safe_div(n_hit, n_has)[ids]


def finalize():
    """Vegetated and bare shares, sand-filling flags and reference-forest values, from the
    composites written by main(). Updates vegetation_cells.csv and vegetation.json."""
    import json

    from lagos_heat.s2_haze import reference_forest_mask

    spec = GridSpec.load()
    grid = load_grid()
    ids = grid["cell_id"].to_numpy()
    gbox = spec.geobox(RES)
    cells = spec.pixel_cells(gbox.transform, gbox.shape)
    forest, _ = reference_forest_mask(spec, gbox)
    # ndvi_change_p80 was written by earlier runs; absolute p80 is not a change measure.
    out = pd.read_csv(INTERIM / "vegetation_cells.csv").drop(columns="ndvi_change_p80",
                                                               errors="ignore")
    assert (out["cell_id"].to_numpy() == ids).all()
    forest_ndvi = {}
    for period in PERIODS:
        forest_ndvi[period] = {}
        for key, raster in FOREST_VERSIONS:
            arr = read_composite(period, raster)
            forest_ndvi[period][key] = float(np.nanmedian(arr[forest]))
            del arr
        for stem, suffix, raster, threshold in SHARES:
            arr = read_composite(period, raster)
            out[f"{stem}_{period}{suffix}"] = share_above(cells, arr, threshold, ids,
                                                          spec.ncells)
            del arr
    for stem in ("veg_share", "bare_share"):
        for suffix in ("", "_screened"):
            out[f"{stem}_change{suffix}"] = (out[f"{stem}_late{suffix}"]
                                             - out[f"{stem}_early{suffix}"])
    for suffix in ("", "_screened"):
        out[f"sandfill_candidate{suffix}"] = (out[f"bare_share_change{suffix}"]
                                              >= SANDFILL_BARE_SHARE_RISE)
    out.to_csv(INTERIM / "vegetation_cells.csv", index=False)

    land = ~grid["water_flag"].to_numpy()
    summ = json.loads((INTERIM / "summaries" / "vegetation.json").read_text())
    summ["reference_forest_ndvi_of_composites"] = {
        "note": "median over the reference-forest tree pixels (s2_haze.py) of each composite",
        **forest_ndvi}
    summ["veg_share_rule"] = (f"share of valid 20 m pixels with NDVI > {NDVI_VEGETATED}; "
                              "veg_share_* from the p80 NDVI, *_screened from the screened "
                              "median")
    summ["sandfill_candidate_rule"] = (
        f"share of pixels with BSI > {BSI_BARE} rises by >= {SANDFILL_BARE_SHARE_RISE} "
        "between periods; sandfill_candidate from the unscreened median BSI, "
        "sandfill_candidate_screened from the screened median BSI; approximate (haze)")
    summ["sandfill_candidate_land_cells"] = int(out.loc[land, "sandfill_candidate"].sum())
    summ["sandfill_candidate_screened_land_cells"] = int(
        out.loc[land, "sandfill_candidate_screened"].sum())
    write_summary("vegetation", summ)
    print("vegetation: finalize done")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--limit-dates", type=int, default=None,
                   help="diagnostic: process only the N clearest pending dates per period")
    p.add_argument("--finalize", action="store_true",
                   help="only recompute shares, flags and forest values from the rasters")
    a = p.parse_args()
    if a.finalize:
        finalize()
    else:
        main(a.limit_dates)
