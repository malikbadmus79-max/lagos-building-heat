"""Phase 3 validation: independent recompute of sampled cells, a ground-truth spot-check sheet,
and candidate place names for the case-study areas.

Step 1 (recompute) draws a fixed-seed stratified sample of 40 main-sample cells and recomputes
six cell values straight from the source rasters and tiles. It does not use the project's
aggregation code (common.aggregate, GridSpec.pixel_cells, buildings.read_tile,
vegetation.cell_mean): cell bounds come from grid_spec.json, rasters are read as small windows
with rasterio, pixel centres are tested against the cell bounds here, and the Open Buildings
tiles are re-read over HTTP through GDAL's own /vsicurl/ driver at the native 0.5 m pixel
size (the project read them at 4 m through a separate requests-based reader). The cached
files in data/raw/_work/ob hold per-cell counts, not pixels, so they cannot serve as an
independent source.

Outputs
- data/processed/phase3_recompute.csv        one row per cell and variable
- data/processed/phase3_spotcheck.csv        sheet for manual checking against imagery
- data/processed/phase3_spotcheck_cells.kml  the 40 cell squares, for Google Earth Pro
- data/processed/phase3_case_names.csv       Nominatim names near the case-study areas
- data/interim/phase3_nominatim.json         raw Nominatim responses
- data/interim/summaries/phase3.json         summary (sample, tolerances, match rates)

Run steps one at a time:
    python -m lagos_heat.validate --sample      (no network)
    python -m lagos_heat.validate --recompute   (HTTP window reads of Open Buildings tiles)
    python -m lagos_heat.validate --spotcheck   (one request for the Esri Wayback release list)
    python -m lagos_heat.validate --names       (Nominatim, one request per second)
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
SUMMARY = INTERIM / "summaries" / "phase3.json"
WORK = ROOT / "data" / "raw" / "_work"

SEED = 20261103
STRATA = {"hotspot": 15, "high_change_not_hotspot": 15, "low_change": 10}
# Tolerances, fixed before the recompute was run. LST and NDVI are re-read from the same
# rasters with the same rule (mean over valid pixels whose centre lies in the cell), so they
# should agree to float32 rounding. Building share is recomputed from the 0.5 m pixels, while
# the dataset sampled the tiles at 4 m, so a sampling difference of up to 2 percentage
# points is allowed.
TOL = {"bshare": 0.02, "lst_rel": 0.005, "ndvi_rel_p80": 0.0005}
VARIABLES = {  # dataset column -> (kind, period)
    "bshare_2016": ("bshare", 2016), "bshare_2023": ("bshare", 2023),
    "lst_rel_early": ("lst_rel", "early"), "lst_rel_late": ("lst_rel", "late"),
    "ndvi_rel_early_p80": ("ndvi_rel_p80", "early"), "ndvi_rel_late_p80": ("ndvi_rel_p80", "late"),
}
OB_PRESENCE_MIN = 0.5  # building_presence above this counts as built (as in Phase 1)
OB_NODATA = -99.0
OB_PIXEL = 0.5
NDVI_FLAT = 0.02  # |relative NDVI change| below this is described as "little change"
USER_AGENT = "lagos-building-heat-validation/0.1 (research script; low-volume reverse geocoding)"


# ------------------------------------------------------------------ pure helpers
def load_spec():
    return json.loads((INTERIM / "grid_spec.json").read_text())


def cell_bounds(cell_id, spec):
    """(x0, y0, x1, y1) in UTM 31N of a cell; cell id = row * ncol + col, row 0 at the top."""
    row, col = divmod(int(cell_id), int(spec["ncol"]))
    size = float(spec["cell"])
    x0 = spec["xmin"] + col * size
    y1 = spec["ymax"] - row * size
    return x0, y1 - size, x0 + size, y1


def centres_in_cell(transform, shape, bounds):
    """Boolean mask of pixels whose centre lies in the cell, with the same edge rule as the
    grid (x0 <= x < x1, y0 < y <= y1)."""
    x0, y0, x1, y1 = bounds
    nr, nc = shape
    xs = transform.c + (np.arange(nc) + 0.5) * transform.a
    ys = transform.f + (np.arange(nr) + 0.5) * transform.e
    return ((ys > y0) & (ys <= y1))[:, None] & ((xs >= x0) & (xs < x1))[None, :]


def stratified_sample(df, hot_cut_b, median_b, rng):
    """Fixed-seed draw of STRATA cells. df needs cell_id, hotspot, bshare_change (main sample)."""
    strata = {
        "hotspot": df["hotspot"],
        "high_change_not_hotspot": ~df["hotspot"] & (df["bshare_change"] >= hot_cut_b),
        "low_change": df["bshare_change"].abs() <= median_b,
    }
    out = []
    for name, k in STRATA.items():
        ids = np.sort(df.loc[strata[name], "cell_id"].to_numpy())
        pick = rng.choice(ids, size=k, replace=False)
        out.append(pd.DataFrame({"cell_id": np.sort(pick), "stratum": name,
                                 "stratum_size": len(ids)}))
    return pd.concat(out, ignore_index=True)


def compare(dataset, recomputed, tol):
    diff = recomputed - dataset
    both_nan = np.isnan(dataset) & np.isnan(recomputed)
    return diff, both_nan | (np.abs(diff) <= tol)


def ndvi_direction(change):
    if not np.isfinite(change):
        return "no data"
    if abs(change) < NDVI_FLAT:
        return "little change"
    return "decrease" if change < 0 else "increase"


def wayback_pick(releases, target):
    """Release (num, date) whose date is closest to target (YYYY-MM-DD)."""
    t = pd.Timestamp(target)
    return min(releases, key=lambda r: (abs(pd.Timestamp(r[1]) - t), r[1]))


MANUAL_COLUMNS = ["seen_2016", "seen_2023", "agrees_with_data", "notes"]


def keep_manual_columns(new, old):
    """Carry hand-entered columns (MANUAL_COLUMNS and any human_*) from an existing sheet
    into a regenerated one, matched by cell_id, so that rerunning never erases them."""
    new = new.copy()
    manual = list(MANUAL_COLUMNS)
    if old is not None:
        manual += [c for c in old.columns if c.startswith("human_") and c not in manual]
        old = old.assign(cell_id=old["cell_id"].astype(int)).set_index("cell_id")
    for c in manual:
        if old is not None and c in old.columns:
            new[c] = new["cell_id"].map(old[c]).fillna("")
        else:
            new[c] = ""
    return new


def google_maps_url(lat, lon, zoom=17):
    return (f"https://www.google.com/maps/@?api=1&map_action=map&center={lat:.6f},{lon:.6f}"
            f"&zoom={zoom}&basemap=satellite")


def wayback_url(num, lat, lon, zoom=17):
    return (f"https://livingatlas.arcgis.com/wayback/#active={num}"
            f"&mapCenter={lon:.6f}%2C{lat:.6f}%2C{zoom}&mode=explore")


def name_candidates(responses):
    """Tally neighbourhood-level names from Nominatim reverse responses (dicts)."""
    keys = ("neighbourhood", "quarter", "suburb", "hamlet", "village", "town",
            "city_district", "residential")
    rows = []
    for r in responses:
        addr = (r.get("response") or {}).get("address", {})
        for k in keys:
            if addr.get(k):
                rows.append({"case": r["case"], "level": k, "name": addr[k],
                             "query": r["query"]})
    if not rows:
        return pd.DataFrame(columns=["case", "level", "name", "n_points", "from_centre"])
    t = pd.DataFrame(rows)
    g = (t.groupby(["case", "level", "name"])
         .agg(n_points=("query", "size"),
              from_centre=("query", lambda q: bool(q.str.startswith("centre").any())))
         .reset_index())
    return g.sort_values(["case", "from_centre", "n_points"], ascending=[True, False, False])


# ------------------------------------------------------------------ step 1: sample
def write_summary(key, value):
    d = json.loads(SUMMARY.read_text()) if SUMMARY.exists() else {}
    d[key] = value
    SUMMARY.write_text(json.dumps(d, indent=2, default=float))


def sample():
    p2 = json.loads((INTERIM / "summaries" / "phase2.json").read_text())
    s = pd.read_csv(PROCESSED / "phase2_sample.csv", usecols=["cell_id", "main", "hotspot"])
    g = pd.read_csv(INTERIM / "grid_periods.csv", usecols=["cell_id", "bshare_change"])
    df = s[s["main"]].merge(g, on="cell_id")
    hot_cut = p2["exposure"]["quintile_cut_bshare_change"]
    med = p2["bshare_change_main"]["median"]
    smp = stratified_sample(df, hot_cut, med, np.random.default_rng(SEED))
    smp.to_csv(WORK / "phase3_sample.csv", index=False)
    write_summary("sample", {
        "seed": SEED, "population": "phase2 main sample", "n_main": int(len(df)),
        "strata": {k: {"n_drawn": v, "rule": r, "stratum_size": int(
            smp.loc[smp["stratum"] == k, "stratum_size"].iloc[0])}
            for (k, v), r in zip(STRATA.items(), [
                "hotspot cell (phase2_sample.csv)",
                f"not hotspot and bshare_change >= {hot_cut:.6f} (hotspot quintile cut)",
                f"|bshare_change| <= {med:.6f} (main-sample median)"])},
    })
    print(smp.groupby("stratum").size())
    return smp


# ------------------------------------------------------------------ step 1: recompute
def raster_cell_mean(path, bounds, valid_path=None):
    """Mean of finite pixels with centre in the cell; validity optionally from another raster."""
    import rasterio
    from rasterio.windows import from_bounds

    x0, y0, x1, y1 = bounds
    with rasterio.open(path) as src:
        pad = 2 * abs(src.res[0])
        win = from_bounds(x0 - pad, y0 - pad, x1 + pad, y1 + pad, transform=src.transform)
        win = win.round_offsets().round_lengths()
        arr = src.read(1, window=win, boundless=True, fill_value=np.nan).astype("float64")
        tr = src.window_transform(win)
    ok = np.isfinite(arr)
    if valid_path is not None:
        with rasterio.open(valid_path) as src:
            ok &= np.isfinite(src.read(1, window=win, boundless=True, fill_value=np.nan))
    inside = centres_in_cell(tr, arr.shape, bounds)
    sel = ok & inside
    return (float(arr[sel].mean()) if sel.any() else np.nan), int(sel.sum()), int(inside.sum())


def gdal_env():
    import rasterio

    return rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                        GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="3",
                        GDAL_HTTP_TIMEOUT="60", GDAL_HTTP_CONNECTTIMEOUT="30")


def vsicurl(uri):
    return "/vsicurl/" + uri.replace("gs://", "https://storage.googleapis.com/")


def tile_bounds(year):
    """{uri: (x0, y0, x1, y1)} for the Open Buildings tiles read in Phase 1, from tile headers."""
    import rasterio

    cache = WORK / f"phase3_ob_tile_bounds_{year}.json"
    if cache.exists():
        return {k: tuple(v) for k, v in json.loads(cache.read_text()).items()}
    uris = pd.read_csv(INTERIM / "ob_tiles.csv").query("year == @year")["uri"].tolist()
    out = {}
    with gdal_env():
        for k, uri in enumerate(uris, 1):
            with rasterio.open(vsicurl(uri)) as src:
                assert src.crs.to_epsg() == 32631, (uri, src.crs)
                out[uri] = tuple(src.bounds)
            print(f"  header {year} {k}/{len(uris)}", flush=True)
    cache.write_text(json.dumps(out))
    return out


def ob_share(bounds, tiles):
    """Built share of valid 0.5 m pixels in the cell, mosaicking all tiles (first valid wins)."""
    import rasterio
    from rasterio.windows import Window, from_bounds

    x0, y0, x1, y1 = bounds
    n = int(round((x1 - x0) / OB_PIXEL))
    mosaic = np.full((n, n), np.nan, dtype="float32")
    used = 0
    for uri, (tx0, ty0, tx1, ty1) in sorted(tiles.items()):
        if tx1 <= x0 or tx0 >= x1 or ty1 <= y0 or ty0 >= y1:
            continue
        with rasterio.open(vsicurl(uri)) as src:
            band = list(src.descriptions).index("building_presence") + 1
            win = from_bounds(x0, y0, x1, y1, transform=src.transform)
            win = win.round_offsets().round_lengths().intersection(
                Window(0, 0, src.width, src.height))
            arr = src.read(band, window=win)
            tr = src.window_transform(win)
        used += 1
        xs = tr.c + (np.arange(arr.shape[1]) + 0.5) * tr.a
        ys = tr.f + (np.arange(arr.shape[0]) + 0.5) * tr.e
        ci = np.floor((xs - x0) / OB_PIXEL).astype(int)
        ri = np.floor((y1 - ys) / OB_PIXEL).astype(int)
        okc, okr = (ci >= 0) & (ci < n), (ri >= 0) & (ri < n)
        sub = arr[np.ix_(okr, okc)]
        tgt = mosaic[np.ix_(ri[okr], ci[okc])]
        fill = np.isnan(tgt) & (sub != OB_NODATA)
        tgt[fill] = sub[fill]
        mosaic[np.ix_(ri[okr], ci[okc])] = tgt
        del arr, sub
    valid = np.isfinite(mosaic)
    share = float((mosaic[valid] > OB_PRESENCE_MIN).mean()) if valid.any() else np.nan
    return share, float(valid.mean()), used


def recompute():
    from pyproj import Transformer

    spec = load_spec()
    smp = pd.read_csv(WORK / "phase3_sample.csv")
    g = pd.read_csv(INTERIM / "grid_periods.csv",
                    usecols=["cell_id", "row", "col", "lga", "lon", "lat", "ob_overlap_approx_2016",
                             "ob_overlap_approx_2023", *VARIABLES])
    smp = smp.merge(g, on="cell_id", how="left")
    to_ll = Transformer.from_crs("EPSG:32631", "EPSG:4326", always_xy=True)
    tiles = {}
    with gdal_env():
        for y in (2016, 2023):
            tiles[y] = tile_bounds(y)
    rows, cells = [], []
    for k, r in enumerate(smp.itertuples(index=False), 1):
        b = cell_bounds(r.cell_id, spec)
        lon, lat = to_ll.transform((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
        cells.append({"cell_id": r.cell_id, "row_col_ok": (r.cell_id // spec["ncol"], r.cell_id
                      % spec["ncol"]) == (r.row, r.col),
                      "lon_recomputed": lon, "lat_recomputed": lat,
                      "centre_diff_m": float(np.hypot((lon - r.lon) * 110574 * np.cos(np.radians(lat)),
                                                      (lat - r.lat) * 110574))})
        for col, (kind, per) in VARIABLES.items():
            extra = {}
            cache = WORK / "phase3_ob" / f"{r.cell_id}_{per}.json"
            if kind == "bshare" and cache.exists():  # finished in an earlier run
                val, cov, used = json.loads(cache.read_text())
                extra = {"ob_valid_share": cov, "ob_tiles_used": used,
                         "ob_overlap_approx": bool(getattr(r, f"ob_overlap_approx_{per}"))}
            elif kind == "bshare":
                for attempt in range(4):
                    try:
                        with gdal_env():
                            val, cov, used = ob_share(b, tiles[per])
                        break
                    except Exception as e:  # network hiccup: wait and retry
                        print(f"    retry {attempt + 1}: {e}", flush=True)
                        time.sleep(5 * (attempt + 1))
                else:
                    raise RuntimeError(f"Open Buildings read failed for cell {r.cell_id}")
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_text(json.dumps([val, cov, used]))
                extra = {"ob_valid_share": cov, "ob_tiles_used": used,
                         "ob_overlap_approx": bool(getattr(r, f"ob_overlap_approx_{per}"))}
            elif kind == "lst_rel":
                val, npx, nin = raster_cell_mean(PROCESSED / f"lst_{per}_rel_median_30m.tif", b,
                                                 PROCESSED / f"lst_{per}_abs_median_30m.tif")
                extra = {"pixels_used": npx, "pixels_in_cell": nin}
            else:
                val, npx, nin = raster_cell_mean(PROCESSED / f"s2_{per}_ndvi_rel_p80_20m.tif", b)
                extra = {"pixels_used": npx, "pixels_in_cell": nin}
            ds = float(getattr(r, col))
            diff, ok = compare(np.array([ds]), np.array([val]), TOL[kind])
            rows.append({"cell_id": r.cell_id, "stratum": r.stratum, "lga": r.lga,
                         "variable": col, "dataset": ds, "recomputed": val,
                         "difference": float(diff[0]), "tolerance": TOL[kind],
                         "match": bool(ok[0]), **extra})
        print(f"  cell {k}/{len(smp)} {r.cell_id} done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(PROCESSED / "phase3_recompute.csv", index=False)
    cdf = pd.DataFrame(cells)
    by_var = (out.groupby("variable")
              .agg(n=("match", "size"), matches=("match", "sum"),
                   max_abs_difference=("difference", lambda d: float(np.nanmax(np.abs(d)))),
                   mean_difference=("difference", "mean"))
              .reset_index())
    mism = out[~out["match"]]
    write_summary("recompute", {
        "tolerances": TOL, "n_cells": int(smp["cell_id"].nunique()), "n_values": int(len(out)),
        "n_match": int(out["match"].sum()), "match_rate": float(out["match"].mean()),
        "by_variable": by_var.to_dict(orient="records"),
        "cell_ids_row_col_ok": bool(cdf["row_col_ok"].all()),
        "max_centre_diff_m": float(cdf["centre_diff_m"].max()),
        "mismatches": mism[["cell_id", "variable", "dataset", "recomputed",
                            "difference"]].to_dict(orient="records"),
    })
    print(by_var.to_string(index=False))
    print(f"match {out['match'].sum()}/{len(out)}; centre max diff {cdf['centre_diff_m'].max():.3f} m")


# ------------------------------------------------------------------ step 3: spot-check sheet
def wayback_releases():
    from lagos_heat.common import fetch

    cfg = json.loads(fetch("https://s3-us-west-2.amazonaws.com/config.maptiles.arcgis.com/"
                           "waybackconfig.json"))
    rel = []
    for num, item in cfg.items():
        title = item.get("itemTitle", "")
        if "Wayback" in title:
            rel.append((int(num), title.split("Wayback")[-1].strip(" )")))
    return rel


def kml(cells, spec):
    from pyproj import Transformer

    to_ll = Transformer.from_crs("EPSG:32631", "EPSG:4326", always_xy=True)
    marks = []
    for r in cells.itertuples(index=False):
        x0, y0, x1, y1 = cell_bounds(r.cell_id, spec)
        ring = [to_ll.transform(x, y) for x, y in
                ((x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0))]
        coords = " ".join(f"{lo:.6f},{la:.6f},0" for lo, la in ring)
        marks.append(
            f"<Placemark><name>{r.cell_id}</name><description>{r.stratum}; {r.lga}</description>"
            "<styleUrl>#cell</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>"
            f"{coords}</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>")
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">'
            "<Document><name>Phase 3 spot-check cells (500 m)</name>"
            '<Style id="cell"><LineStyle><color>ff00ffff</color><width>2</width></LineStyle>'
            "<PolyStyle><fill>0</fill></PolyStyle></Style>" + "".join(marks) +
            "</Document></kml>\n")


def spotcheck():
    from pyproj import Transformer

    spec = load_spec()
    smp = pd.read_csv(WORK / "phase3_sample.csv")
    g = pd.read_csv(INTERIM / "grid_periods.csv",
                    usecols=["cell_id", "lga", "bshare_2016", "bshare_2023", "bshare_change",
                             "ndvi_rel_change_p80", "lst_rel_change"])
    df = smp.merge(g, on="cell_id", how="left")
    to_ll = Transformer.from_crs("EPSG:32631", "EPSG:4326", always_xy=True)
    ll = [to_ll.transform((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
          for b in (cell_bounds(c, spec) for c in df["cell_id"])]
    df["centre_lon"] = [round(p[0], 6) for p in ll]
    df["centre_lat"] = [round(p[1], 6) for p in ll]
    rel = wayback_releases()
    w16, w23 = wayback_pick(rel, "2016-06-30"), wayback_pick(rel, "2023-06-30")
    df["google_maps_url"] = [google_maps_url(la, lo) for la, lo in zip(df.centre_lat, df.centre_lon)]
    df["wayback_2016_url"] = [wayback_url(w16[0], la, lo) for la, lo in
                              zip(df.centre_lat, df.centre_lon)]
    df["wayback_2023_url"] = [wayback_url(w23[0], la, lo) for la, lo in
                              zip(df.centre_lat, df.centre_lon)]
    df["data_building_share_2016"] = df["bshare_2016"].round(3)
    df["data_building_share_2023"] = df["bshare_2023"].round(3)
    df["data_building_change_pp"] = (100 * df["bshare_change"]).round(1)
    df["data_ndvi_change_direction"] = df["ndvi_rel_change_p80"].map(ndvi_direction)
    df["data_ndvi_rel_change_p80"] = df["ndvi_rel_change_p80"].round(3)
    df["data_lst_rel_change_c"] = df["lst_rel_change"].round(2)
    sheet = PROCESSED / "phase3_spotcheck.csv"
    df = keep_manual_columns(df, pd.read_csv(sheet, dtype=str) if sheet.exists() else None)
    cols = ["cell_id", "stratum", "lga", "centre_lat", "centre_lon", "google_maps_url",
            "wayback_2016_url", "wayback_2023_url", "data_building_share_2016",
            "data_building_share_2023", "data_building_change_pp", "data_ndvi_change_direction",
            "data_ndvi_rel_change_p80", "data_lst_rel_change_c", *MANUAL_COLUMNS]
    cols += [c for c in df.columns if c.startswith("human_") and c not in cols]
    df[cols].to_csv(sheet, index=False)
    (PROCESSED / "phase3_spotcheck_cells.kml").write_text(kml(df, spec), encoding="utf-8")
    write_summary("spotcheck", {
        "n_cells": int(len(df)), "wayback_release_2016": {"num": w16[0], "date": w16[1]},
        "wayback_release_2023": {"num": w23[0], "date": w23[1]},
        "ndvi_flat_threshold": NDVI_FLAT,
        "ndvi_direction_counts": df["data_ndvi_change_direction"].value_counts().to_dict(),
    })
    print(f"wayback {w16} {w23}")


# ------------------------------------------------------------------ step 4: names
def nominatim(lat, lon, zoom):
    import requests

    r = requests.get("https://nominatim.openstreetmap.org/reverse",
                     params={"lat": f"{lat:.6f}", "lon": f"{lon:.6f}", "zoom": zoom,
                             "format": "jsonv2", "addressdetails": 1},
                     headers={"User-Agent": USER_AGENT, "Accept-Language": "en"}, timeout=60)
    r.raise_for_status()
    return r.json()


def names():
    from pyproj import Transformer

    spec = load_spec()
    cases = pd.read_csv(PROCESSED / "phase2_case_studies.csv").query("chosen")
    smp = pd.read_csv(PROCESSED / "phase2_sample.csv", usecols=["cell_id", "cluster"])
    to_ll = Transformer.from_crs("EPSG:32631", "EPSG:4326", always_xy=True)
    queries = []
    for c in cases.itertuples(index=False):
        for z in (18, 16, 14):
            queries.append((int(c.rank), f"centre_z{z}", c.lat, c.lon, z))
        for cid in smp.loc[smp["cluster"] == c.cluster, "cell_id"].sort_values():
            b = cell_bounds(cid, spec)
            lon, lat = to_ll.transform((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)
            queries.append((int(c.rank), f"cell_{cid}_z16", lat, lon, 16))
    out, last = [], 0.0
    for k, (case, q, lat, lon, z) in enumerate(queries, 1):
        wait = 1.1 - (time.time() - last)
        if wait > 0:
            time.sleep(wait)
        last = time.time()
        try:
            resp = nominatim(lat, lon, z)
        except Exception as e:  # record and continue; one request per second is kept
            resp = {"error": str(e)}
        out.append({"case": case, "query": q, "lat": lat, "lon": lon, "zoom": z,
                    "response": resp})
        print(f"  {k}/{len(queries)} case {case} {q}", flush=True)
    (INTERIM / "phase3_nominatim.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    tab = name_candidates(out)
    centre = {(o["case"]): o["response"].get("display_name", "") for o in out
              if o["query"] == "centre_z18" and "response" in o}
    tab["centre_display_name_z18"] = tab["case"].map(centre)
    tab["confirmed_by_author"] = ""
    tab.to_csv(PROCESSED / "phase3_case_names.csv", index=False)
    write_summary("names", {"n_requests": len(queries),
                            "n_errors": sum("error" in o["response"] for o in out),
                            "attribution": "Data (c) OpenStreetMap contributors, ODbL",
                            "user_agent": USER_AGENT})
    print(tab.to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = ap.add_mutually_exclusive_group(required=True)
    for s in ("sample", "recompute", "spotcheck", "names"):
        g.add_argument(f"--{s}", action="store_true")
    a = ap.parse_args()
    {"sample": sample, "recompute": recompute, "spotcheck": spotcheck,
     "names": names}[next(k for k, v in vars(a).items() if v)]()


if __name__ == "__main__":
    main()
