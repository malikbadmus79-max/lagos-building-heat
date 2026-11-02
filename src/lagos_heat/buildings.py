"""Step 2: building presence and height per cell (Open Buildings 2.5D Temporal, 2016 and 2023),
with GHSL GHS-BUILT-S (2015, 2020) as a cross-check.

Open Buildings GeoTIFFs (0.5 m pixels, 4 m effective resolution) are read one tile at a
time at 4 m through HTTP Range requests (GDAL uses the internal overviews), and reduced
to per-cell accumulators immediately, so memory stays low.

Outputs
- data/interim/buildings_cells.csv
- data/interim/ob_tiles.csv            tiles read, per year
- data/interim/summaries/buildings.json
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from lagos_heat.common import (
    INTERIM, RAW, UTM31N, WORK, GridSpec, HTTPRangeFile, aggregate, fetch, http_opener,
    lagos_lgas,
    load_grid, retry, safe_div, write_summary,
)

OB_BUCKET = "open-buildings-temporal-data"
OB_HTTP = f"https://storage.googleapis.com/{OB_BUCKET}"
OB_YEARS = (2016, 2023)
OB_RES = 4.0  # metres, read resolution
PRESENCE_THRESHOLD = 0.5
OB_NODATA = -99.0

GHSL_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A"
GHSL_TILE = "R9_C19"
GHSL_EPOCHS = (2015, 2020)


def list_manifests(year):
    api = f"https://storage.googleapis.com/storage/v1/b/{OB_BUCKET}/o"
    names, token = [], None
    while True:
        params = {"prefix": "v1/manifests/", "fields": "items(name),nextPageToken",
                  "maxResults": 1000}
        if token:
            params["pageToken"] = token
        page = json.loads(fetch(api, params=params))
        names += [i["name"] for i in page.get("items", [])]
        token = page.get("nextPageToken")
        if not token:
            break
    epsg = UTM31N.split(":")[1]
    return sorted(n for n in names if f"_EPSG_{epsg}_{year}_06_30" in n)


def tile_footprints(manifest):
    from rasterio.transform import Affine
    from shapely.geometry import Polygon

    out = []
    for tileset in manifest["tilesets"]:
        for source in tileset["sources"]:
            a = source["affineTransform"]
            t = Affine.translation(a["translateX"], a["translateY"]) * Affine.scale(
                a["scaleX"], a["scaleY"])
            w, h = source["dimensions"]["width"], source["dimensions"]["height"]
            out.append((manifest["uriPrefix"] + source["uris"][0],
                        Polygon([t * c for c in [(0, 0), (w, 0), (w, h), (0, h)]])))
    return out


def lagos_tiles(year, lagos_union):
    tiles = {}
    for name in list_manifests(year):
        manifest = json.loads(fetch(f"{OB_HTTP}/{name}"))
        for uri, poly in tile_footprints(manifest):
            if poly.intersects(lagos_union):
                tiles[uri] = poly
    return tiles


def read_tile(url, spec):
    """Per-cell accumulators for one tile at OB_RES: valid px, built px, height sum over built px."""
    import rasterio
    from rasterio.enums import Resampling

    with rasterio.open(url, opener=http_opener) as src:
        names = list(src.descriptions)
        b_pres = names.index("building_presence") + 1
        b_height = names.index("building_height") + 1
        factor = OB_RES / abs(src.res[0])
        shape = (int(round(src.height / factor)), int(round(src.width / factor)))
        arr = src.read([b_pres, b_height], out_shape=(2, *shape),
                       resampling=Resampling.nearest)
        transform = src.transform * src.transform.scale(src.width / shape[1],
                                                        src.height / shape[0])
    pres, height = arr[0], arr[1]
    cells = spec.pixel_cells(transform, pres.shape)
    valid = pres != OB_NODATA
    built = valid & (pres > PRESENCE_THRESHOLD)
    n_valid, _ = aggregate(cells, mask=valid, ncells=spec.ncells)
    n_built, _ = aggregate(cells, mask=built, ncells=spec.ncells)
    h_ok = built & (height != OB_NODATA) & np.isfinite(height)
    h_sum, h_n = aggregate(cells, values=np.where(h_ok, height, 0), mask=h_ok,
                           ncells=spec.ncells)
    stats = {"valid_px": int(valid.sum()), "built_px": int(built.sum()),
             "height_min": float(height[h_ok].min()) if h_ok.any() else None,
             "height_max": float(height[h_ok].max()) if h_ok.any() else None}
    return n_valid, n_built, h_sum, h_n, stats


def open_buildings(spec, lagos_union, limit=None):
    acc, tile_rows, info = {}, [], {}
    for year in OB_YEARS:
        tiles = lagos_tiles(year, lagos_union)
        uris = sorted(tiles)[:limit] if limit else sorted(tiles)
        n_valid = np.zeros(spec.ncells)
        n_built = np.zeros(spec.ncells)
        h_sum = np.zeros(spec.ncells)
        h_n = np.zeros(spec.ncells)
        # Tiles of neighbouring S2-cell folders overlap and carry identical values in the
        # overlap, so summing tiles double-counts pixels. Keep, per cell, the single tile
        # with the most valid pixels; where it covers the whole cell it is used alone.
        best = np.zeros((4, spec.ncells))  # nv, nb, hs, hn of the best single tile
        before, t0 = HTTPRangeFile.bytes_fetched, time.time()
        for k, uri in enumerate(uris, 1):
            url = uri.replace("gs://", "https://storage.googleapis.com/")
            b0 = HTTPRangeFile.bytes_fetched
            cache = WORK / "ob" / f"{'_'.join(uri.split('/')[-2:])}.npz"
            if cache.exists():  # tile finished in an earlier run
                z = np.load(cache, allow_pickle=True)
                nv, nb, hs, hn = (np.zeros(spec.ncells) for _ in range(4))
                for arr, key in ((nv, "nv"), (nb, "nb"), (hs, "hs"), (hn, "hn")):
                    arr[z["idx"]] = z[key]
                st = z["stats"].item()
                st["mb_fetched_original"] = float(z["mb"])
            else:
                nv, nb, hs, hn, st = retry(read_tile, url, spec)
                HTTPRangeFile.clear_cache()
                idx = np.nonzero(nv)[0]
                cache.parent.mkdir(parents=True, exist_ok=True)
                np.savez(cache, idx=idx, nv=nv[idx], nb=nb[idx], hs=hs[idx], hn=hn[idx],
                         stats=np.array(st, dtype=object),
                         mb=(HTTPRangeFile.bytes_fetched - b0) / 1e6)
            n_valid += nv
            n_built += nb
            h_sum += hs
            h_n += hn
            better = nv > best[0]
            best[:, better] = np.vstack([nv, nb, hs, hn])[:, better]
            mb = st.pop("mb_fetched_original", (HTTPRangeFile.bytes_fetched - b0) / 1e6)
            tile_rows.append({"year": year, "uri": uri, "mb_fetched": mb, **st})
            print(f"  OB {year} tile {k}/{len(uris)} {mb:.1f} MB", flush=True)
        full = (spec.cell / OB_RES) ** 2
        single = best[0] >= full
        summed = np.vstack([n_valid, n_built, h_sum, h_n])
        dedup = np.where(single, best, summed)
        # cells where tiles overlap but none covers the cell alone: duplicates remain
        approx = ~single & (n_valid > full * 1.01)
        acc[year] = (*dedup, approx)
        info[year] = {"n_tiles": len(uris), "n_tiles_intersecting_lagos": len(tiles),
                      "mb_fetched": float(sum(r["mb_fetched"] for r in tile_rows
                                              if r["year"] == year)),
                      "seconds_this_run": round(time.time() - t0)}
    return acc, pd.DataFrame(tile_rows), info


def ghsl(spec):
    """Mean built-up fraction of GHSL 100 m pixels whose centres fall in each cell."""
    import rasterio
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds

    out = {}
    for epoch in GHSL_EPOCHS:
        stem = f"GHS_BUILT_S_E{epoch}_GLOBE_R2023A_54009_100"
        name = f"{stem}_V1_0_{GHSL_TILE}"
        local = RAW / "ghsl" / f"{name}.zip"
        if not local.exists():
            fetch(f"{GHSL_BASE}/{stem}/V1-0/tiles/{name}.zip", local)
        with rasterio.open(f"/vsizip/{local.as_posix()}/{name}.tif") as src:
            b = transform_bounds(UTM31N, src.crs, *spec.bounds, densify_pts=50)
            win = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
            arr = src.read(1, window=win, boundless=True, fill_value=src.nodata or 0)
            transform = src.window_transform(win)
            cells = spec.pixel_cells_reprojected(transform, arr.shape, src.crs)
            ok = arr != src.nodata if src.nodata is not None else np.ones(arr.shape, bool)
        total, count = aggregate(cells, values=np.where(ok, arr, 0) / 10000.0, mask=ok,
                                 ncells=spec.ncells)
        out[epoch] = (safe_div(total, count), count)
    return out


def main(limit=None):
    spec = GridSpec.load()
    grid = load_grid()
    lagos_union = lagos_lgas().union_all()
    acc, tiles, info = open_buildings(spec, lagos_union, limit=limit)
    ids = grid["cell_id"].to_numpy()
    df = pd.DataFrame({"cell_id": ids})
    px_area = OB_RES**2
    for year, (nv, nb, hs, hn, approx) in acc.items():
        df[f"bshare_{year}"] = safe_div(nb, nv)[ids]
        df[f"bheight_{year}"] = safe_div(hs, hn)[ids]
        df[f"ob_coverage_{year}"] = np.minimum(nv * px_area / spec.cell**2, 1.0)[ids]
        df[f"ob_overlap_approx_{year}"] = approx[ids]
    for epoch, (frac, count) in ghsl(spec).items():
        df[f"ghsl_share_{epoch}"] = frac[ids]
        df[f"ghsl_npix_{epoch}"] = count[ids]
    df.to_csv(INTERIM / "buildings_cells.csv", index=False)
    tiles.to_csv(INTERIM / "ob_tiles.csv", index=False)

    land = ~grid["water_flag"].to_numpy()
    summ = {"open_buildings": info,
            "read_resolution_m": OB_RES, "presence_threshold": PRESENCE_THRESHOLD,
            "total_mb_fetched": float(tiles["mb_fetched"].sum()) if len(tiles) else 0.0,
            "limit": limit}
    for y in OB_YEARS:
        cov = df[f"ob_coverage_{y}"]
        summ[f"ob_{y}"] = {
            "cells_no_data": int((cov == 0).sum()),
            "coverage_min": float(cov.min()),
            "cells_coverage_below_1": int((cov < 0.999).sum()),
            "cells_overlap_approx": int(df[f"ob_overlap_approx_{y}"].sum()),
            "bshare_mean_land_cells": float(df.loc[land, f"bshare_{y}"].mean()),
            "bheight_min": float(df[f"bheight_{y}"].min()),
            "bheight_max": float(df[f"bheight_{y}"].max()),
            "tile_height_max": float(tiles.loc[tiles.year == y, "height_max"].max()),
        }
    for e in GHSL_EPOCHS:
        s = df[f"ghsl_share_{e}"]
        summ[f"ghsl_{e}"] = {"min": float(s.min()), "max": float(s.max()),
                             "mean_land_cells": float(s[land].mean()),
                             "npix_min": float(df[f"ghsl_npix_{e}"].min()),
                             "npix_max": float(df[f"ghsl_npix_{e}"].max())}
    valid = df[["bshare_2016", "ghsl_share_2015", "bshare_2023", "ghsl_share_2020"]][land].dropna()
    summ["cross_check_pearson_land_cells"] = {
        "ob2016_vs_ghsl2015": float(valid["bshare_2016"].corr(valid["ghsl_share_2015"])),
        "ob2023_vs_ghsl2020": float(valid["bshare_2023"].corr(valid["ghsl_share_2020"])),
        "n_cells": int(len(valid)),
    }
    write_summary("buildings", summ)
    print(f"buildings: {summ['total_mb_fetched']:.0f} MB fetched")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=None, help="tiles per year (testing)")
    main(p.parse_args().limit)
