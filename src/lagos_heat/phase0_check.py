"""Phase 0 feasibility checks for the Lagos building density / surface heat project.

Run from the repository root:

    python -m lagos_heat.phase0_check        (with src/ on PYTHONPATH)

Each check returns a dict with: name, source_url, access, resolution,
years_available, worked, details (list of str), errors (list of str).
Results are written to docs/phase0_results.md. Catalogue searches use the
Lagos bounding box; pixel reads use a small test window so downloads stay small.
"""

from __future__ import annotations

import datetime as dt
import io
import re
import statistics
import traceback
from pathlib import Path

import numpy as np

from lagos_heat.common import (
    LAGOS_BBOX, PC_STAC, RAW, ROOT, UTM31N, HTTPRangeFile, fetch, http_opener, pc_client,
)

RESULTS = ROOT / "docs" / "phase0_results.md"

TEST_BBOX = (3.33, 6.48, 3.38, 6.53)  # central mainland Lagos, about 5.5 x 5.5 km
# Urban fringe, Ibeju-Lekki, about 2.2 x 2.2 km. Placed at 3.96-3.98 E (not centred on
# 3.95 E) so it falls inside a single Open Buildings Temporal GeoTIFF tile.
FRINGE_BBOX = (3.96, 6.44, 3.98, 6.46)
SEASONS = {
    "2017": "2016-12-01/2017-02-28",
    "2025": "2024-12-01/2025-02-28",
}
CLOUD_MAX = 20  # percent, scene level


def new_result(name, source_url, access):
    return {
        "name": name,
        "source_url": source_url,
        "access": access,
        "resolution": "",
        "years_available": "",
        "worked": False,
        "details": [],
        "errors": [],
    }


def run_safely(fn):
    """Run one check; any uncaught exception is recorded rather than raised."""
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001
        res = new_result(fn.__name__, "", "")
        res["errors"].append(f"{type(exc).__name__}: {exc}")
        res["errors"].append(traceback.format_exc(limit=3))
        return res


def collection_extent(client, collection_id):
    col = client.get_collection(collection_id)
    start, end = col.extent.temporal.intervals[0]
    fmt = lambda d: d.date().isoformat() if d else "ongoing"  # noqa: E731
    return f"{fmt(start)} to {fmt(end)}"


def summarise_cloud(values):
    if not values:
        return "n/a"
    return (
        f"min {min(values):.1f}%, median {statistics.median(values):.1f}%, "
        f"max {max(values):.1f}%"
    )


def intersects_test(item):
    from shapely.geometry import box, shape

    return shape(item.geometry).intersects(box(*TEST_BBOX))


# ---------------------------------------------------------------- a. LGAs
def check_lga_boundaries():
    import geopandas as gpd
    import requests

    api = "https://www.geoboundaries.org/api/current/gbOpen/NGA/{}/"
    res = new_result(
        "a. Lagos LGA boundaries (geoBoundaries NGA ADM2)",
        api.format("ADM2"),
        "geoBoundaries REST API -> GeoJSON (simplified geometry), no key",
    )
    meta = {lvl: requests.get(api.format(lvl), timeout=60).json() for lvl in ("ADM1", "ADM2")}
    m2 = meta["ADM2"]
    res["years_available"] = f"boundaryYearRepresented = {m2.get('boundaryYearRepresented')}"
    res["resolution"] = "vector polygons (simplified geometry used for this check)"
    res["details"].append(
        f"ADM2 upstream source: {m2.get('boundarySource')}; licence: {m2.get('boundaryLicense')}; "
        f"national unit count: {m2.get('admUnitCount')}"
    )
    res["details"].append(f"Full-resolution GeoJSON: {m2.get('gjDownloadURL')}")

    # Fetched with requests rather than GDAL's curl: on some Windows setups GDAL's
    # schannel backend fails certificate revocation checks.
    adm1 = gpd.read_file(io.BytesIO(fetch(meta["ADM1"]["simplifiedGeometryGeoJSON"])))
    adm2 = gpd.read_file(io.BytesIO(fetch(m2["simplifiedGeometryGeoJSON"])))
    lagos = adm1[adm1["shapeName"].str.lower() == "lagos"]
    if len(lagos) != 1:
        res["errors"].append(f"Expected one ADM1 named Lagos, found {len(lagos)}")
        return res
    lagos_geom = lagos.geometry.iloc[0]
    pts = adm2.geometry.representative_point()
    lgas = adm2[pts.within(lagos_geom)].copy()
    n = len(lgas)
    res["details"].append(f"ADM2 units inside Lagos ADM1 (representative point test): {n}")
    res["details"].append("Names: " + ", ".join(sorted(lgas["shapeName"])))
    RAW.mkdir(parents=True, exist_ok=True)
    out = RAW / "lagos_lga_geoboundaries_simplified.geojson"
    lgas.to_file(out, driver="GeoJSON")
    res["details"].append(f"Saved to {out.relative_to(ROOT).as_posix()}")
    bounds = lgas.total_bounds
    res["details"].append(
        "Lagos LGA bounds (lon/lat): " + ", ".join(f"{b:.3f}" for b in bounds)
    )
    res["worked"] = n == 20
    if n != 20:
        res["errors"].append(f"Expected 20 LGAs, found {n}")
    return res


# ---------------------------------------------------------- b. Sentinel-2
def check_sentinel2():
    import odc.stac

    res = new_result(
        "b. Sentinel-2 L2A (Planetary Computer)",
        f"{PC_STAC}/collections/sentinel-2-l2a",
        "STAC search (pystac-client) + signed COG reads (odc-stac), no key",
    )
    client = pc_client()
    res["years_available"] = collection_extent(client, "sentinel-2-l2a")
    best = None
    for label, period in SEASONS.items():
        items = list(
            client.search(collections=["sentinel-2-l2a"], bbox=LAGOS_BBOX, datetime=period).items()
        )
        cc = [it.properties.get("eo:cloud_cover") for it in items]
        cc = [c for c in cc if c is not None]
        low = [it for it in items if it.properties.get("eo:cloud_cover", 100) < CLOUD_MAX]
        tiles = sorted({it.properties.get("s2:mgrs_tile") for it in items})
        low_dates = sorted({it.datetime.date().isoformat() for it in low})
        res["details"].append(
            f"Dry season {label} ({period}): {len(items)} scenes over Lagos bbox; "
            f"{len(low)} with cloud < {CLOUD_MAX}%; cloud cover {summarise_cloud(cc)}; "
            f"MGRS tiles {', '.join(t for t in tiles if t)}; "
            f"distinct low-cloud acquisition dates: {len(low_dates)}"
        )
        cand = sorted(
            (it for it in low if intersects_test(it)),
            key=lambda it: it.properties["eo:cloud_cover"],
        )
        if label == "2025" and cand:
            best = cand[0]
        if not items:
            res["errors"].append(f"No scenes returned for {label}")

    if best is not None:
        ds = odc.stac.load([best], bands=["B04"], bbox=TEST_BBOX, crs=UTM31N, resolution=10)
        arr = ds["B04"].isel(time=0).values
        valid = float(np.mean(arr > 0))
        res["resolution"] = "10 m (B04 test read); 20 m and 60 m bands also exist"
        res["details"].append(
            f"Test read: {best.id} (cloud {best.properties['eo:cloud_cover']:.1f}%), "
            f"B04 over test window -> array {arr.shape}, valid fraction {valid:.3f}"
        )
        res["worked"] = not res["errors"]
    else:
        res["errors"].append("No low-cloud 2025 scene intersecting the test window for a test read")
    return res


# ------------------------------------------------------------- c. Landsat
def check_landsat_st():
    import odc.stac

    res = new_result(
        "c. Landsat C2 L2 surface temperature ST_B10 (Planetary Computer)",
        f"{PC_STAC}/collections/landsat-c2-l2",
        "STAC search + signed COG reads (asset 'lwir11' = ST_B10), no key",
    )
    client = pc_client()
    res["years_available"] = collection_extent(client, "landsat-c2-l2")
    best = None
    for label, period in SEASONS.items():
        items = list(
            client.search(collections=["landsat-c2-l2"], bbox=LAGOS_BBOX, datetime=period).items()
        )
        with_st = [it for it in items if "lwir11" in it.assets]
        by_platform = {}
        for it in with_st:
            p = it.properties.get("platform", "unknown")
            by_platform[p] = by_platform.get(p, 0) + 1
        cc = [it.properties.get("eo:cloud_cover") for it in with_st]
        cc = [c for c in cc if c is not None]
        low = [it for it in with_st if it.properties.get("eo:cloud_cover", 100) < CLOUD_MAX]
        paths = sorted(
            {
                f"{it.properties.get('landsat:wrs_path')}/{it.properties.get('landsat:wrs_row')}"
                for it in with_st
            }
        )
        res["details"].append(
            f"Dry season {label} ({period}): {len(items)} scenes, {len(with_st)} with ST_B10 "
            f"(by platform: {by_platform}); {len(low)} with cloud < {CLOUD_MAX}%; "
            f"cloud cover {summarise_cloud(cc)}; WRS path/row {', '.join(paths)}"
        )
        cand = sorted(
            (it for it in low if intersects_test(it)),
            key=lambda it: it.properties["eo:cloud_cover"],
        )
        if label == "2025" and cand:
            best = cand[0]
        if not with_st:
            res["errors"].append(f"No ST_B10 scenes returned for {label}")

    if best is not None:
        ds = odc.stac.load([best], bands=["lwir11"], bbox=TEST_BBOX, crs=UTM31N, resolution=30)
        raw = ds["lwir11"].isel(time=0).values.astype("float64")
        raw[raw == 0] = np.nan
        # Collection 2 Level-2 ST scale/offset (USGS): K = DN * 0.00341802 + 149.0
        celsius = raw * 0.00341802 + 149.0 - 273.15
        res["resolution"] = "30 m grid (thermal band resampled by USGS; native TIRS 100 m)"
        res["details"].append(
            f"Test read: {best.id} (cloud {best.properties['eo:cloud_cover']:.1f}%), "
            f"lwir11 over test window -> array {raw.shape}; LST "
            f"min {np.nanmin(celsius):.1f}, median {np.nanmedian(celsius):.1f}, "
            f"max {np.nanmax(celsius):.1f} deg C (no cloud masking applied)"
        )
        res["worked"] = not res["errors"]
    else:
        res["errors"].append("No low-cloud 2025 scene intersecting the test window for a test read")
    return res


# -------------------------------------------------------- d. WorldCover
def check_worldcover():
    import rasterio
    from rasterio.windows import from_bounds

    res = new_result(
        "d. ESA WorldCover 2020 / 2021 (Planetary Computer)",
        f"{PC_STAC}/collections/esa-worldcover",
        "STAC search + signed COG window read (asset 'map'), no key",
    )
    client = pc_client()
    res["years_available"] = collection_extent(client, "esa-worldcover")
    items = list(client.search(collections=["esa-worldcover"], bbox=TEST_BBOX).items())
    found = set()
    for it in sorted(items, key=lambda i: i.properties.get("start_datetime", "")):
        version = it.properties.get("esa_worldcover:product_version")
        year = (it.properties.get("start_datetime") or str(it.datetime))[:4]
        with rasterio.open(it.assets["map"].href) as src:
            win = from_bounds(*TEST_BBOX, transform=src.transform)
            arr = src.read(1, window=win)
            res["resolution"] = f"{abs(src.res[0]):.7f} deg pixel (about 10 m), CRS {src.crs}"
        vals, counts = np.unique(arr, return_counts=True)
        share = ", ".join(f"{v}: {c / arr.size:.1%}" for v, c in zip(vals, counts))
        res["details"].append(
            f"{year} (product {version}), item {it.id}: test window array {arr.shape}; "
            f"class shares {share}"
        )
        found.add(year)
    res["details"].append(
        "Class codes: 10 tree, 20 shrub, 30 grass, 40 crop, 50 built-up, 60 bare, "
        "70 snow/ice, 80 water, 90 herbaceous wetland, 95 mangrove, 100 moss/lichen"
    )
    missing = {"2020", "2021"} - found
    if missing:
        res["errors"].append(f"Missing years for test window: {sorted(missing)}")
    res["worked"] = not missing
    return res


# ------------------------------------------------------------ e. GHSL
GHSL_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A"
GHSL_TILE = "R9_C19"  # Mollweide 1000 km tile grid; covers Lagos


def ghsl_url(epoch):
    stem = f"GHS_BUILT_S_E{epoch}_GLOBE_R2023A_54009_100"
    name = f"{stem}_V1_0_{GHSL_TILE}"
    return f"{GHSL_BASE}/{stem}/V1-0/tiles/{name}.zip", f"{name}.tif"


def check_ghsl():
    import rasterio
    import requests
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds

    res = new_result(
        "e1. GHSL GHS-BUILT-S R2023A (JRC open download)",
        f"{GHSL_BASE}/",
        f"HTTPS tile zip ({GHSL_TILE}, Mollweide ESRI:54009, 100 m), downloaded to "
        "data/raw/ghsl/ and read locally for the test window",
    )
    res["details"].append(
        "Remote read via GDAL /vsizip//vsicurl/ failed on this machine (schannel revocation "
        "check in GDAL's curl), so whole tile zips are downloaded with requests instead."
    )
    def read_window(src, bbox):
        bounds = transform_bounds("EPSG:4326", src.crs, *bbox)
        arr = src.read(1, window=from_bounds(*bounds, transform=src.transform))
        valid = arr[arr != src.nodata] if src.nodata is not None else arr
        return arr.shape, valid

    windows = {"central test window": (TEST_BBOX, (2015, 2020, 2025)),
               "fringe window (Ibeju-Lekki)": (FRINGE_BBOX, (2015, 2020))}
    ok = []
    for epoch in (2015, 2020, 2025):
        url, tif = ghsl_url(epoch)
        local = RAW / "ghsl" / url.rsplit("/", 1)[-1]
        try:
            if local.exists():
                # Downloaded on an earlier run; the JRC server can be slow to answer HEAD.
                status = f"local copy from earlier download, tile zip {local.stat().st_size / 1e6:.1f} MB"
            else:
                head = requests.head(url, timeout=60, allow_redirects=True)
                if head.status_code != 200:
                    res["errors"].append(f"E{epoch}: HTTP {head.status_code} for {url}")
                    continue
                size = head.headers.get("Content-Length")
                size_txt = f"{int(size) / 1e6:.1f} MB" if size else "unknown size"
                status = f"HTTP 200, tile zip {size_txt}"
                fetch(url, local)
            parts = []
            with rasterio.open(f"/vsizip/{local.as_posix()}/{tif}") as src:
                res["resolution"] = f"{abs(src.res[0]):.0f} m, CRS {src.crs.to_string()}"
                for label, (bbox, epochs) in windows.items():
                    if epoch not in epochs:
                        continue
                    shape, valid = read_window(src, bbox)
                    parts.append(
                        f"{label} {bbox}: array {shape}, built-up surface sum "
                        f"{valid.sum() / 1e4:,.1f} ha, cells with any built-up "
                        f"{np.mean(valid > 0):.1%}"
                    )
            res["details"].append(f"E{epoch}: {status}; " + "; ".join(parts))
            ok.append(epoch)
        except Exception as exc:  # noqa: BLE001
            res["errors"].append(f"E{epoch}: {type(exc).__name__}: {exc}")
    res["details"].append("Cell values are m2 of built-up surface per 100 m cell (max 10,000).")
    res["years_available"] = (
        "R2023A epochs 1975-2030 in 5-year steps (2025/2030 are model projections per JRC); "
        f"tested: {', '.join(str(e) for e in ok) or 'none'}"
    )
    res["worked"] = len(ok) == 3
    return res


def check_open_buildings_temporal():
    import requests

    catalog = (
        "https://developers.google.com/earth-engine/datasets/catalog/"
        "GOOGLE_Research_open-buildings-temporal_v1"
    )
    bucket = "open-buildings-temporal-data"
    res = new_result(
        "e2. Google Open Buildings 2.5D Temporal",
        catalog,
        'Earth Engine ImageCollection "GOOGLE/Research/open-buildings-temporal/v1"; '
        f"also public GCS bucket gs://{bucket}/v1",
    )
    page = requests.get(catalog, timeout=60).text
    m = re.search(r"(\d{4})-\d\d-\d\dT[\d:]+Z\s*(?:&ndash;|–|-)\s*(\d{4})-\d\d-\d\dT", page)
    res["years_available"] = (
        f"{m.group(1)} to {m.group(2)} (annual, from EE catalog page)" if m else "not parsed"
    )
    rm = re.search(r"(\d+)\s*m effective spatial resolution", page)
    res["resolution"] = (
        f"{rm.group(1)} m effective (from EE catalog page)" if rm else "not parsed"
    )
    api = f"https://storage.googleapis.com/storage/v1/b/{bucket}/o"
    r = requests.get(api, params={"prefix": "v1/", "delimiter": "/"}, timeout=60)
    if r.status_code == 200:
        prefixes = r.json().get("prefixes", [])
        res["details"].append(
            f"Anonymous GCS listing of gs://{bucket}/v1/: HTTP 200, prefixes {prefixes}"
        )
        r2 = requests.get(
            api, params={"prefix": "v1/geotiffs/", "delimiter": "/", "maxResults": 50}, timeout=60
        )
        sub = r2.json().get("prefixes", []) if r2.status_code == 200 else []
        res["details"].append(
            f"First-level entries under v1/geotiffs/ (up to 50): {sub[:10]}"
            + (" ..." if len(sub) > 10 else "")
        )
        res["details"].append(
            "Earth Engine is not strictly required: GeoTIFFs are in a public bucket, and Google "
            "provides a Colab for regional download "
            "(google-research/building_detection/open_buildings_temporal_download_region_geotiffs"
            ".ipynb). Lagos tiles are located and read in check e3."
        )
        res["worked"] = True
    else:
        res["errors"].append(f"GCS listing failed: HTTP {r.status_code}")
        res["details"].append("Access appears to require Earth Engine.")
    return res


OB_BUCKET = "open-buildings-temporal-data"
OB_HTTP = f"https://storage.googleapis.com/{OB_BUCKET}"


def ob_tile_footprints(manifest):
    """(gs:// uri, footprint polygon) for each GeoTIFF in an Open Buildings manifest.

    Footprints follow Google's download Colab: affine transform + pixel dimensions.
    """
    from rasterio.transform import Affine
    from shapely.geometry import Polygon

    out = []
    for tileset in manifest["tilesets"]:
        for source in tileset["sources"]:
            a = source["affineTransform"]
            t = Affine.translation(a["translateX"], a["translateY"]) * Affine.scale(
                a["scaleX"], a["scaleY"]
            )
            w, h = source["dimensions"]["width"], source["dimensions"]["height"]
            poly = Polygon([t * c for c in [(0, 0), (w, 0), (w, h), (0, h)]])
            out.append((manifest["uriPrefix"] + source["uris"][0], poly))
    return out


def check_open_buildings_lagos():
    import geopandas as gpd
    import rasterio
    import requests
    from rasterio.enums import Resampling
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds
    from shapely.geometry import box

    res = new_result(
        "e3. Open Buildings 2.5D Temporal: Lagos tiles and fringe read",
        f"gs://{OB_BUCKET}/v1 ({OB_HTTP}/v1/)",
        "Anonymous GCS JSON API listing + manifest JSON; windowed GeoTIFF reads via HTTP "
        "Range requests (no Earth Engine, no key)",
    )
    api = f"https://storage.googleapis.com/storage/v1/b/{OB_BUCKET}/o"
    names, token = [], None
    while True:
        params = {"prefix": "v1/manifests/", "fields": "items(name),nextPageToken",
                  "maxResults": 1000}
        if token:
            params["pageToken"] = token
        page = requests.get(api, params=params, timeout=60).json()
        names += [i["name"] for i in page.get("items", [])]
        token = page.get("nextPageToken")
        if not token:
            break
    # Manifest names: {S2 level-2 token}_EPSG_{utm}_{yyyy}_06_30.json. Lagos is in UTM 31N.
    utm = sorted(n for n in names if f"_EPSG_{UTM31N.split(':')[1]}_" in n)
    years_found = sorted({int(n.rsplit("_", 3)[-3]) for n in utm})
    res["years_available"] = f"{years_found[0]} to {years_found[-1]} annual (manifests for UTM 31N)"
    res["details"].append(f"Manifests in bucket: {len(names)}; for {UTM31N}: {len(utm)}")

    lga_path = RAW / "lagos_lga_geoboundaries_simplified.geojson"
    if lga_path.exists():
        lagos = gpd.read_file(lga_path).to_crs(UTM31N).union_all()
        area_note = "Lagos LGA polygons (union)"
    else:
        lagos = gpd.GeoSeries([box(*LAGOS_BBOX)], crs="EPSG:4326").to_crs(UTM31N).iloc[0]
        area_note = "Lagos bbox (LGA file missing)"
    fringe = box(*transform_bounds("EPSG:4326", UTM31N, *FRINGE_BBOX))

    presence = {}
    for year in (2016, 2023):
        manifests = [n for n in utm if f"_{year}_06_30" in n]
        lagos_tiles, fringe_tiles = set(), []
        for name in manifests:
            manifest = requests.get(f"{OB_HTTP}/{name}", timeout=120).json()
            for uri, poly in ob_tile_footprints(manifest):
                if poly.intersects(lagos):
                    lagos_tiles.add(uri)
                if poly.contains(fringe):
                    fringe_tiles.append(uri)
        folders = sorted({u.split("/")[-2] for u in lagos_tiles})
        res["details"].append(
            f"{year}: {len(lagos_tiles)} GeoTIFF tiles intersect {area_note}, in "
            f"{len(folders)} S2-cell folders ({', '.join(folders)}); manifests used: "
            f"{', '.join(m.split('/')[-1] for m in manifests)}"
        )
        # Tiles from neighbouring S2 cells overlap; each carries data only inside its own
        # cell (-99 elsewhere), so the window is merged across all tiles that contain it.
        merged, before = None, HTTPRangeFile.bytes_fetched
        for uri in fringe_tiles:
            with rasterio.open(uri.replace("gs://", "https://storage.googleapis.com/"),
                               opener=http_opener) as src:
                band = src.descriptions.index("building_presence") + 1
                win = from_bounds(*fringe.bounds, transform=src.transform)
                # 0.5 m pixels, 4 m effective resolution: read at 4 m (factor 8) via overviews
                shape = (round(win.height / 8), round(win.width / 8))
                arr = src.read(band, window=win, out_shape=shape, resampling=Resampling.nearest)
                nodata = src.nodata
            merged = arr if merged is None else np.where(merged == nodata, arr, merged)
        if merged is None:
            res["errors"].append(f"{year}: no tile fully contains the fringe window")
            continue
        valid = merged[merged != nodata]
        presence[year] = valid
        res["resolution"] = "GeoTIFF pixels 0.5 m (4 m effective); fringe read at 4 m"
        res["details"].append(
            f"{year} fringe window {FRINGE_BBOX}: tiles {[u.split('/', 4)[-1] for u in fringe_tiles]}; "
            f"array {merged.shape} at 4 m, valid pixels {valid.size}; "
            f"building_presence mean {valid.mean():.4f}, pixels with presence > 0.5: "
            f"{np.mean(valid > 0.5):.2%}; bytes fetched "
            f"{(HTTPRangeFile.bytes_fetched - before) / 1e6:.1f} MB"
        )

    # Context for the fringe window: WorldCover 2021 class shares
    try:
        client = pc_client()
        items = [i for i in client.search(collections=["esa-worldcover"], bbox=FRINGE_BBOX).items()
                 if (i.properties.get("start_datetime") or "").startswith("2021")]
        with rasterio.open(items[0].assets["map"].href) as src:
            wc = src.read(1, window=from_bounds(*FRINGE_BBOX, transform=src.transform))
        vals, counts = np.unique(wc, return_counts=True)
        res["details"].append(
            "Context, WorldCover 2021 in fringe window, class shares: "
            + ", ".join(f"{v}: {c / wc.size:.1%}" for v, c in zip(vals, counts))
        )
    except Exception as exc:  # noqa: BLE001
        res["errors"].append(f"WorldCover context read: {type(exc).__name__}: {exc}")

    res["worked"] = set(presence) == {2016, 2023} and not any(
        e.startswith(("2016", "2023")) for e in res["errors"]
    )
    return res


# ------------------------------------------------------------ f. WorldPop
WORLDPOP_API = "https://hub.worldpop.org/rest/data/pop"
WORLDPOP_PRODUCTS = {
    "wpgp": "Global 2000-2020, unconstrained 100 m",
    "G2_CN_POP_R25A_100m": "Global2 R2025A 2015-2030, constrained 100 m",
}
WORLDPOP_TEST = ("G2_CN_POP_R25A_100m", 2016)


def check_worldpop():
    import rasterio
    import requests
    from rasterio.windows import from_bounds

    res = new_result(
        "f. WorldPop Nigeria 100 m population",
        f"{WORLDPOP_API}/<product>?iso3=NGA (files on https://data.worldpop.org)",
        "WorldPop REST API for file listing; HTTPS download of national GeoTIFF, no key",
    )
    files, years_txt = {}, []
    for alias, label in WORLDPOP_PRODUCTS.items():
        data = requests.get(f"{WORLDPOP_API}/{alias}", params={"iso3": "NGA"}, timeout=120).json()
        recs = data.get("data", [])
        yrs = sorted(int(r["popyear"]) for r in recs)
        doi = recs[0].get("doi") if recs else None
        files[alias] = {int(r["popyear"]): r["files"][0] for r in recs if r.get("files")}
        years_txt.append(f"{alias}: {yrs[0]}-{yrs[-1]}" if yrs else f"{alias}: none")
        res["details"].append(
            f"{alias} ({label}): {len(recs)} Nigeria records, years "
            f"{', '.join(map(str, yrs))}; DOI {doi}"
        )
    res["years_available"] = "; ".join(years_txt)

    alias, year = WORLDPOP_TEST
    url = files[alias][year]
    head = requests.head(url, timeout=60, allow_redirects=True)
    size_mb = int(head.headers.get("Content-Length", 0)) / 1e6
    res["details"].append(f"Test file {url}: HTTP {head.status_code}, {size_mb:.1f} MB")
    with requests.get(url, headers={"Range": "bytes=0-1023"}, stream=True, timeout=60) as r:
        ranged = r.status_code == 206
    res["details"].append(
        "Byte-range request honoured (windowed remote read possible)" if ranged else
        "Server ignores byte-range requests (returns HTTP 200 and the full file), so a "
        "clip requires downloading the national raster"
    )
    local = RAW / "worldpop" / url.rsplit("/", 1)[-1]
    if not local.exists():
        fetch(url, local)
    with rasterio.open(local) as src:
        res["resolution"] = (
            f"{abs(src.res[0]):.7f} deg pixel (about 100 m), CRS {src.crs.to_string()}"
        )
        for label, bbox in (("central test window", TEST_BBOX),
                            ("fringe window", FRINGE_BBOX)):
            win = from_bounds(*bbox, transform=src.transform)
            arr = src.read(1, window=win)
            valid = arr[arr != src.nodata] if src.nodata is not None else arr
            res["details"].append(
                f"{year} clip, {label} {bbox}: array {arr.shape}, populated cells "
                f"{valid.size}, population sum {valid.sum():,.0f}"
            )
            if label == "central test window":
                profile = src.profile | {
                    "height": arr.shape[0], "width": arr.shape[1],
                    "transform": src.window_transform(win), "compress": "deflate",
                }
                clip = RAW / "worldpop" / f"{local.stem}_test_clip.tif"
                with rasterio.open(clip, "w", **profile) as dst:
                    dst.write(arr, 1)
                res["details"].append(f"Clip saved to {clip.relative_to(ROOT).as_posix()}")
    res["details"].append(
        f"National file kept at {local.relative_to(ROOT).as_posix()} (gitignored)"
    )
    res["worked"] = head.status_code == 200
    return res


# ------------------------------------------------------------- report
DESIGN_DECISIONS = """\
## Design decisions

Recorded after Phase 0. These set the analysis design for Phase 1 and later.

1. Building change is measured with Google Open Buildings 2.5D Temporal, 2016 to 2023. \
GHSL GHS-BUILT-S 2015 and 2020 are used as an independent cross-check. GHSL 2025 is a \
model projection and is not used as observed data.
2. Vegetation change is measured with Sentinel-2 L2A NDVI composites.
3. ESA WorldCover is used only as a baseline tree mask. It is not used to measure change, \
because the 2020 (v100) and 2021 (v200) maps were produced with different algorithm versions.
4. Dry-season composites are pooled over two periods: early, December 2016 to February 2018 \
(dry seasons 2016/17 and 2017/18); late, December 2023 to February 2025 (dry seasons \
2023/24 and 2024/25). Pooling two seasons per period compensates for the small number of \
low-cloud scenes in a single dry season.
5. Surface temperature comes from Landsat 8 only in both periods, so that sensor \
differences between Landsat 8 and Landsat 9 do not enter the change estimate.
"""


def write_report(results):
    lines = [
        "# Phase 0 results: data feasibility",
        "",
        f"Generated by `src/lagos_heat/phase0_check.py` on "
        f"{dt.datetime.now().strftime('%Y-%m-%d %H:%M')} (local time).",
        "",
        f"Lagos search bbox (lon/lat): {LAGOS_BBOX}. Test window: {TEST_BBOX}. "
        f"Fringe window (Ibeju-Lekki): {FRINGE_BBOX}. "
        f"Dry seasons: {SEASONS}. Low cloud: scene cloud cover < {CLOUD_MAX}%.",
        "",
        "## Summary",
        "",
        "| Check | Worked | Access | Resolution | Years available |",
        "|---|---|---|---|---|",
    ]
    for r in results:
        cell = lambda s: str(s).replace("|", "/")  # noqa: E731
        lines.append(
            f"| {cell(r['name'])} | {'yes' if r['worked'] else 'no'} | {cell(r['access'])} | "
            f"{cell(r['resolution'])} | {cell(r['years_available'])} |"
        )
    for r in results:
        lines += ["", f"## {r['name']}", "", f"- Source: {r['source_url']}",
                  f"- Access: {r['access']}", f"- Resolution: {r['resolution']}",
                  f"- Years available: {r['years_available']}",
                  f"- Worked: {'yes' if r['worked'] else 'no'}"]
        lines += [f"- {d}" for d in r["details"]]
        if r["errors"]:
            lines += ["- Errors:", "", "```"] + r["errors"] + ["```"]
        else:
            lines.append("- Errors: none")
    lines += ["", DESIGN_DECISIONS]
    RESULTS.write_text("\n".join(lines), encoding="utf-8")


def main():
    checks = [
        check_lga_boundaries,
        check_sentinel2,
        check_landsat_st,
        check_worldcover,
        check_ghsl,
        check_open_buildings_temporal,
        check_open_buildings_lagos,
        check_worldpop,
    ]
    results = []
    for fn in checks:
        print(f"Running {fn.__name__} ...", flush=True)
        r = run_safely(fn)
        print(f"  worked={r['worked']} errors={len(r['errors'])}", flush=True)
        results.append(r)
    write_report(results)
    print(f"Wrote {RESULTS}")


if __name__ == "__main__":
    main()
