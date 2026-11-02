"""Step 6: covariates per cell: distance to ocean and lagoon, distance to major roads,
share of industrial land use.

Water bodies (ESA WorldCover 2021, class 80, 50 m, grid extent + 15 km margin):
- ocean: water connected to the southern edge of the extent after an erosion of
  OCEAN_EROSION_M (which cuts narrow channels such as the Lagos harbour entrance),
  dilated back by the same distance and intersected with water;
- lagoon: remaining connected water bodies of at least LAGOON_MIN_KM2.
Distances are Euclidean distances on the 50 m raster from each cell centre.

OpenStreetMap (Overpass API, Lagos bbox + margin): ways tagged highway=motorway, trunk or
primary (link roads excluded); ways and multipolygon relations tagged landuse=industrial.

Outputs
- data/processed/water_bodies_50m.tif  (0 land, 1 ocean, 2 lagoon, 3 other water)
- data/raw/osm/overpass_*.json         raw Overpass responses
- data/interim/covariates_cells.csv
- data/interim/summaries/covariates.json
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from lagos_heat.common import (
    CELL, INTERIM, PROCESSED, RAW, UTM31N, GridSpec, configure_gdal, load_grid, pc_client,
    retry, write_raster, write_summary,
)

WC_RES = 50
MARGIN_M = 15000
OCEAN_EROSION_M = 1500
LAGOON_MIN_KM2 = 5.0
OVERPASS = "https://overpass-api.de/api/interpreter"
OSM_MARGIN_DEG = 0.05
ROAD_CLASSES = ("motorway", "trunk", "primary")


# ----------------------------------------------------------------- water
def water_bodies(spec):
    import odc.stac
    from affine import Affine
    from odc.geo.geobox import GeoBox
    from scipy import ndimage

    configure_gdal()
    x0, y0, x1, y1 = spec.bounds
    x0, y0, x1, y1 = x0 - MARGIN_M, y0 - MARGIN_M, x1 + MARGIN_M, y1 + MARGIN_M
    shape = (int((y1 - y0) / WC_RES), int((x1 - x0) / WC_RES))
    gbox = GeoBox(shape, Affine(WC_RES, 0, x0, 0, -WC_RES, y1), UTM31N)
    lonlat = gbox.extent.to_crs("EPSG:4326").boundingbox
    found = retry(lambda: list(pc_client().search(collections=["esa-worldcover"],
                                                  bbox=tuple(lonlat)).items()))
    items = [i for i in found if (i.properties.get("start_datetime") or "").startswith("2021")]
    ds = odc.stac.load(items, bands=["map"], geobox=gbox, resampling="mode",
                       groupby="solar_day", pool=4)
    wc = ds["map"].isel(time=0).values
    nodata = wc == 0
    water = (wc == 80) | nodata  # no-data treated as water (open sea beyond tile coverage)

    r = OCEAN_EROSION_M / WC_RES
    core = ndimage.distance_transform_edt(water) > r
    lab, _ = ndimage.label(core)
    south = np.unique(lab[-1][lab[-1] > 0])
    ocean_core = np.isin(lab, south)
    ocean = water & (ndimage.distance_transform_edt(~ocean_core) <= r + 1)

    rest = water & ~ocean
    lab2, n2 = ndimage.label(rest)
    sizes = np.bincount(lab2.ravel(), minlength=n2 + 1) * WC_RES**2 / 1e6
    big = np.where(sizes >= LAGOON_MIN_KM2)[0]
    big = big[big > 0]
    lagoon = np.isin(lab2, big)

    classes = np.zeros(wc.shape, "uint8")
    classes[ocean] = 1
    classes[lagoon] = 2
    classes[rest & ~lagoon] = 3
    write_raster(PROCESSED / "water_bodies_50m.tif", classes, gbox.transform, nodata=255)

    d_ocean = ndimage.distance_transform_edt(~ocean) * WC_RES
    d_lagoon = ndimage.distance_transform_edt(~lagoon) * WC_RES
    info = {"items": [i.id for i in items], "extent_utm": [x0, y0, x1, y1],
            "nodata_pixels": int(nodata.sum()),
            "ocean_km2": float(ocean.sum() * WC_RES**2 / 1e6),
            "lagoon_km2": float(lagoon.sum() * WC_RES**2 / 1e6),
            "n_lagoon_bodies": int(len(big)),
            "other_water_km2": float((rest & ~lagoon).sum() * WC_RES**2 / 1e6)}
    return gbox, d_ocean, d_lagoon, info


def sample_at(gbox, raster, x, y):
    t = gbox.transform
    col = np.floor((x - t.c) / t.a).astype(int)
    row = np.floor((y - t.f) / t.e).astype(int)
    return raster[row, col]


# ------------------------------------------------------------------- OSM
def overpass(name, query):
    import requests

    path = RAW / "osm" / f"overpass_{name}.json"
    if not path.exists():
        def _post():
            r = requests.post(OVERPASS, data={"data": query}, timeout=600,
                              headers={"User-Agent": "lagos-heat-research/0.1"})
            r.raise_for_status()
            return r

        r = retry(_post)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(r.content)
    return json.loads(path.read_text(encoding="utf-8"))


def osm_bbox(spec):
    from pyproj import Transformer

    tr = Transformer.from_crs(UTM31N, "EPSG:4326", always_xy=True)
    x0, y0, x1, y1 = spec.bounds
    lon0, lat0 = tr.transform(x0, y0)
    lon1, lat1 = tr.transform(x1, y1)
    m = OSM_MARGIN_DEG
    return lat0 - m, lon0 - m, lat1 + m, lon1 + m  # south, west, north, east


def osm_roads(spec):
    import geopandas as gpd
    from shapely.geometry import LineString

    s, w, n, e = osm_bbox(spec)
    q = (f'[out:json][timeout:300];way["highway"~"^({"|".join(ROAD_CLASSES)})$"]'
         f"({s},{w},{n},{e});out geom;")
    data = overpass("roads", q)
    rows = [{"osm_id": el["id"], "highway": el["tags"].get("highway"),
             "geometry": LineString([(p["lon"], p["lat"]) for p in el["geometry"]])}
            for el in data["elements"] if el["type"] == "way" and len(el.get("geometry", [])) > 1]
    gdf = gpd.GeoDataFrame(rows, crs="EPSG:4326").to_crs(UTM31N)
    return gdf, data.get("osm3s", {}).get("timestamp_osm_base")


def osm_polygons(data):
    """Polygons (UTM 31N) from an Overpass 'out geom' response: closed ways and
    multipolygon relations (outer rings minus inner rings). Keeps OSM id and name."""
    import geopandas as gpd
    from shapely.geometry import LineString, Polygon
    from shapely.ops import linemerge, polygonize, unary_union

    rows = []
    for el in data["elements"]:
        geom = None
        if el["type"] == "way":
            pts = [(p["lon"], p["lat"]) for p in el.get("geometry", [])]
            if len(pts) >= 4 and pts[0] == pts[-1]:
                geom = Polygon(pts)
        elif el["type"] == "relation":
            outer, inner = [], []
            for m in el.get("members", []):
                if m.get("type") == "way" and len(m.get("geometry", [])) > 1:
                    line = LineString([(p["lon"], p["lat"]) for p in m["geometry"]])
                    (inner if m.get("role") == "inner" else outer).append(line)
            shell = unary_union(list(polygonize(linemerge(outer)))) if outer else None
            if shell is not None and not shell.is_empty:
                if inner:
                    shell = shell.difference(unary_union(list(polygonize(linemerge(inner)))))
                geom = shell
        if geom is not None:
            rows.append({"osm_type": el["type"], "osm_id": el["id"],
                         "name": el.get("tags", {}).get("name"), "geometry": geom.buffer(0)})
    if not rows:
        return gpd.GeoDataFrame({"osm_type": [], "osm_id": [], "name": []},
                                geometry=[], crs=UTM31N)
    return gpd.GeoDataFrame(rows, crs="EPSG:4326").to_crs(UTM31N)


def osm_industrial(spec):
    s, w, n, e = osm_bbox(spec)
    q = (f'[out:json][timeout:300];(way["landuse"="industrial"]({s},{w},{n},{e});'
         f'relation["landuse"="industrial"]({s},{w},{n},{e}););out geom;')
    data = overpass("industrial", q)
    return osm_polygons(data), data.get("osm3s", {}).get("timestamp_osm_base")


def main():
    import geopandas as gpd

    spec = GridSpec.load()
    grid = load_grid()
    cent = grid.geometry.centroid
    out = pd.DataFrame({"cell_id": grid["cell_id"].to_numpy()})

    gbox, d_ocean, d_lagoon, winfo = water_bodies(spec)
    out["dist_ocean_m"] = sample_at(gbox, d_ocean, cent.x.to_numpy(), cent.y.to_numpy())
    out["dist_lagoon_m"] = sample_at(gbox, d_lagoon, cent.x.to_numpy(), cent.y.to_numpy())
    out["dist_large_water_m"] = np.minimum(out["dist_ocean_m"], out["dist_lagoon_m"])

    roads, ts_roads = osm_roads(spec)
    pts = gpd.GeoDataFrame({"cell_id": out["cell_id"]}, geometry=cent.to_numpy(), crs=UTM31N)
    near = gpd.sjoin_nearest(pts, roads[["highway", "geometry"]], distance_col="d")
    near = near.drop_duplicates("cell_id").set_index("cell_id")
    out["dist_major_road_m"] = out["cell_id"].map(near["d"])

    ind, ts_ind = osm_industrial(spec)
    if len(ind):
        ind_union = ind.union_all()
        inter = grid.geometry.intersection(ind_union).area
        out["industrial_share"] = (inter / CELL**2).clip(upper=1.0).to_numpy()
    else:
        out["industrial_share"] = 0.0
    out.to_csv(INTERIM / "covariates_cells.csv", index=False)

    write_summary("covariates", {
        "water": winfo, "water_rules": {"erosion_m": OCEAN_EROSION_M,
                                        "lagoon_min_km2": LAGOON_MIN_KM2, "res_m": WC_RES},
        "osm_source": OVERPASS, "osm_timestamp_roads": ts_roads,
        "osm_timestamp_industrial": ts_ind,
        "road_classes": ROAD_CLASSES, "n_road_ways": int(len(roads)),
        "road_km_by_class": (roads.assign(km=roads.length / 1000)
                             .groupby("highway")["km"].sum().round(1).to_dict()),
        "n_industrial_polygons": int(len(ind)),
        "industrial_km2": float(ind.union_all().area / 1e6) if len(ind) else 0.0,
        "cells_with_industrial": int((out["industrial_share"] > 0).sum()),
        "dist_ocean_m_max": float(out["dist_ocean_m"].max()),
        "dist_lagoon_m_max": float(out["dist_lagoon_m"].max()),
        "dist_major_road_m_max": float(out["dist_major_road_m"].max()),
    })
    print("covariates: done")


if __name__ == "__main__":
    main()
