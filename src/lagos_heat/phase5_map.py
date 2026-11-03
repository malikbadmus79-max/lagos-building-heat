"""Phase 5: data build for the public web map (docs/map/index.html, GitHub Pages).

Builds one self-contained HTML page (Leaflet from a pinned CDN, data inlined) from the
project files. Every number on the page is read here from a project file:

- data/processed/phase2_sample.csv            main sample, hotspot flag, clusters
- data/interim/grid_periods.csv                relative LST change, LGA, WorldPop 2023
- data/processed/priority_areas_lga.csv        "Cool now" ranks, "Protect" counts by LGA
- data/processed/priority_cells.csv            "Protect" cells
- data/processed/phase2_case_studies.csv       case-study clusters
- data/processed/phase2_models.csv             effect per 10 percentage points (H1)
- data/interim/summaries/phase2.json           hotspot count and residents

Outputs
- docs/map/index.html
- data/processed/map_lga_stats.csv             the figures shown in the LGA pop-ups
- data/interim/summaries/phase5_map.json       counts and page size

    python -m lagos_heat.phase5_map
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from lagos_heat.common import INTERIM, PROCESSED, ROOT, SUMMARIES, GridSpec, lagos_lgas

OUT_HTML = ROOT / "docs" / "map" / "index.html"
LEAFLET = "1.9.4"
LEAFLET_JS_SRI = "sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo="
LEAFLET_CSS_SRI = "sha256-p4NxAoJBhIIN+hmNHrzRCf9tD/miZyoHS5obTRR9BMY="
REPO = "https://github.com/malikbadmus79-max/lagos-building-heat"
BRIEF = REPO + "/blob/main/docs/policy_brief.md"
# Heat-change classes (°C, relative to Lagos): upper bounds of the lower six classes.
BINS = [-2.0, -1.0, -0.25, 0.25, 1.0, 2.0]
# Diverging blue - grey - red (the Phase 2 palette; adjacent pairs pass the CVD check).
COLOURS = ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f3a5a1", "#e34948", "#9c2526"]
BIN_LABELS = ["2 °C or more cooler", "1 to 2 °C cooler", "0.25 to 1 °C cooler",
              "About the same (within 0.25 °C)", "0.25 to 1 °C hotter", "1 to 2 °C hotter",
              "2 °C or more hotter"]
SIMPLIFY_M = 60.0  # LGA boundary simplification tolerance, metres
DP = 5             # decimal places kept for coordinates (about 1 m)
CASE_LABELS = {1: "Lekki Peninsula II / Chevron area", 2: "Abule Egba / Iyana Ipaja",
               3: "Oke, Ikorodu"}


# ------------------------------------------------------------------ pure helpers
def heat_class(v):
    """Index (0-6) of the heat-change class for each value."""
    return np.searchsorted(np.asarray(BINS), np.asarray(v, float), side="right")


def encode_nodes(lat, lon):
    """Grid-node coordinates as integers of 1e-5 degrees from the minimum (smaller JSON)."""
    lat0, lon0 = float(np.min(lat)), float(np.min(lon))
    return {"lat0": round(lat0, DP), "lon0": round(lon0, DP),
            "lat": np.rint((np.ravel(lat) - round(lat0, DP)) * 10 ** DP).astype(int).tolist(),
            "lon": np.rint((np.ravel(lon) - round(lon0, DP)) * 10 ** DP).astype(int).tolist()}


def decode_node(enc, i):
    return (enc["lat0"] + enc["lat"][i] / 10 ** DP, enc["lon0"] + enc["lon"][i] / 10 ** DP)


def lga_stats(cells, cool, protect):
    """Per-LGA figures for the pop-ups. cells: main-sample cells with lga, hotspot,
    pop_2023, lst_rel_change. cool: A1 rows (lga, priority_rank). protect: B-by-LGA rows."""
    g = cells.groupby("lga")
    out = pd.DataFrame({
        "analysed_squares": g.size(),
        "hotspot_squares": g["hotspot"].sum().astype(int),
        "residents_in_hotspots": g.apply(lambda d: d.loc[d["hotspot"], "pop_2023"].sum(),
                                         include_groups=False),
        "mean_heat_change_all_squares": g["lst_rel_change"].mean(),
    }).reset_index()
    out = out.merge(cool[["lga", "priority_rank", "mean_lst_rel_change"]].rename(
        columns={"priority_rank": "cool_now_rank",
                 "mean_lst_rel_change": "mean_heat_change_hotspot_squares"}),
        on="lga", how="left")
    out = out.merge(protect[["lga", "n_cells"]].rename(columns={"n_cells": "protect_squares"}),
                    on="lga", how="left")
    out["protect_squares"] = out["protect_squares"].fillna(0).astype(int)
    out["cool_now_of"] = int(len(cool))
    return out


def round_coords(geom, dp=DP):
    from shapely.geometry import mapping

    def rnd(c):
        if isinstance(c[0], (int, float)):
            return [round(c[0], dp), round(c[1], dp)]
        return [rnd(x) for x in c]

    m = mapping(geom)
    return {"type": m["type"], "coordinates": rnd(m["coordinates"])}


# ------------------------------------------------------------------ build
def build():
    import shapely
    from pyproj import Transformer
    from shapely.geometry import box
    from shapely.ops import unary_union

    spec = GridSpec.load()
    to_ll = Transformer.from_crs("EPSG:32631", "EPSG:4326", always_xy=True)
    s = pd.read_csv(PROCESSED / "phase2_sample.csv", usecols=["cell_id", "main", "hotspot",
                                                               "cluster"])
    g = pd.read_csv(INTERIM / "grid_periods.csv",
                    usecols=["cell_id", "row", "col", "lga", "pop_2023", "lst_rel_change"])
    cells = s[s["main"]].merge(g, on="cell_id").sort_values("cell_id").reset_index(drop=True)
    pl = pd.read_csv(PROCESSED / "priority_areas_lga.csv")
    cool = pl[pl["list"] == "A1_cool_now_lga"]
    prot_lga = pl[pl["list"] == "B_protect_by_lga"]
    pc = pd.read_csv(PROCESSED / "priority_cells.csv")
    prot_ids = set(pc.loc[pc["list"] == "B_protect_cell", "unit_id"].astype(int))
    cells["protect"] = cells["cell_id"].isin(prot_ids)

    stats = lga_stats(cells, cool, prot_lga)
    stats.to_csv(PROCESSED / "map_lga_stats.csv", index=False)

    # grid nodes (cell corners) in lon/lat
    rr, cc = np.meshgrid(np.arange(spec.nrow + 1), np.arange(spec.ncol + 1), indexing="ij")
    lon, lat = to_ll.transform(spec.xmin + cc * spec.cell, spec.ymax - rr * spec.cell)
    nodes = encode_nodes(lat, lon)

    lgas = sorted(cells["lga"].unique())
    lidx = {n: i for i, n in enumerate(lgas)}
    flags = cells["hotspot"].astype(int) + 2 * cells["protect"].astype(int)
    cdata = {"r": cells["row"].tolist(), "c": cells["col"].tolist(),
             "v": np.rint(cells["lst_rel_change"] * 100).astype(int).tolist(),
             "f": flags.tolist(), "l": cells["lga"].map(lidx).tolist()}

    # LGA polygons, simplified
    geo = lagos_lgas()
    geo["geometry"] = geo.geometry.simplify(SIMPLIFY_M, preserve_topology=True)
    geo = geo.to_crs("EPSG:4326")
    st = stats.set_index("lga")
    feats = []
    for _, r in geo.iterrows():
        p = {"name": r["lga"]}
        if r["lga"] in st.index:
            x = st.loc[r["lga"]]
            p.update({"res": int(round(x["residents_in_hotspots"])),
                      "hot": int(x["hotspot_squares"]),
                      "heat": round(float(x["mean_heat_change_all_squares"]), 2),
                      "hheat": None if pd.isna(x["mean_heat_change_hotspot_squares"])
                      else round(float(x["mean_heat_change_hotspot_squares"]), 2),
                      "rank": None if pd.isna(x["cool_now_rank"]) else int(x["cool_now_rank"]),
                      "of": int(x["cool_now_of"]), "prot": int(x["protect_squares"]),
                      "n": int(x["analysed_squares"])})
        feats.append({"type": "Feature", "properties": p, "geometry": round_coords(r.geometry)})

    # case-study areas: outline of the cluster's squares
    cases = pd.read_csv(PROCESSED / "phase2_case_studies.csv").query("chosen").sort_values("rank")
    case_out = []
    for c in cases.itertuples(index=False):
        sq = g.merge(s.loc[s["cluster"] == c.cluster, ["cell_id"]], on="cell_id")
        poly = unary_union([box(spec.xmin + col * spec.cell, spec.ymax - (row + 1) * spec.cell,
                                spec.xmin + (col + 1) * spec.cell, spec.ymax - row * spec.cell)
                            for row, col in zip(sq["row"], sq["col"])])
        ll = shapely.transform(poly, lambda xy: np.column_stack(to_ll.transform(xy[:, 0],
                                                                                 xy[:, 1])))
        case_out.append({"rank": int(c.rank), "name": CASE_LABELS[int(c.rank)],
                         "lga": c.lga, "squares": int(c.n_cells),
                         "geom": round_coords(ll)})

    # headline numbers for the side panel
    p2 = json.loads((SUMMARIES / "phase2.json").read_text())
    m = pd.read_csv(PROCESSED / "phase2_models.csv").set_index("run")
    text = {
        "per10": round(float(m.loc["H1_primary", "per_10pp"]), 1),
        "hotspots": int(p2["exposure"]["hotspot_cells"]),
        "residents_k": int(round(p2["exposure"]["worldpop2023_hotspot"], -3)),
        "protect": int(len(prot_ids)),
        "squares": int(len(cells)),
    }
    data = {"nodes": nodes, "ncol1": spec.ncol + 1, "cells": cdata, "lgas": lgas,
            "lgaGeo": {"type": "FeatureCollection", "features": feats}, "cases": case_out,
            "bins": BINS, "colours": COLOURS, "binLabels": BIN_LABELS, "text": text}
    return data, stats


def page(data):
    t = data["text"]
    html = TEMPLATE
    subs = {
        "__LEAFLET__": LEAFLET, "__JS_SRI__": LEAFLET_JS_SRI, "__CSS_SRI__": LEAFLET_CSS_SRI,
        "__REPO__": REPO, "__BRIEF__": BRIEF,
        "__PER10__": f"{t['per10']:.1f}", "__HOTSPOTS__": f"{t['hotspots']:,}",
        "__RESIDENTS__": f"{t['residents_k']:,}", "__PROTECT__": f"{t['protect']:,}",
        "__SQUARES__": f"{t['squares']:,}",
        "__DATA__": json.dumps(data, separators=(",", ":")),
    }
    for k, v in subs.items():
        html = html.replace(k, v)
    assert "__" not in html.replace("__proto__", ""), "unfilled placeholder"
    return html


def main():
    data, stats = build()
    html = page(data)
    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(html, encoding="utf-8")
    size = OUT_HTML.stat().st_size
    summ = {"page": "docs/map/index.html", "bytes": size, "squares": data["text"]["squares"],
            "hotspot_squares": int(sum(f & 1 for f in data["cells"]["f"])),
            "protect_squares": int(sum((f & 2) // 2 for f in data["cells"]["f"])),
            "lgas": len(data["lgaGeo"]["features"]), "case_areas": len(data["cases"]),
            "text": data["text"], "leaflet": LEAFLET}
    (SUMMARIES / "phase5_map.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ))


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Is Lagos building itself hotter?</title>
<meta name="description" content="Map of where the ground in Lagos got hotter than the rest of the city as buildings spread, 2016-2025.">
<link rel="stylesheet" href="https://unpkg.com/leaflet@__LEAFLET__/dist/leaflet.css" integrity="__CSS_SRI__" crossorigin="">
<script src="https://unpkg.com/leaflet@__LEAFLET__/dist/leaflet.js" integrity="__JS_SRI__" crossorigin=""></script>
<style>
:root{--ink:#0b0b0b;--ink2:#52514e;--line:#d9d7d0;--panel:#fcfcfb;--hot:#c0392b;--protect:#1e7b34}
*{box-sizing:border-box}
html,body{margin:0;height:100%;font:15px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;color:var(--ink);background:var(--panel)}
#app{display:flex;height:100%}
#panel{width:340px;max-width:100%;overflow-y:auto;padding:16px 18px;border-right:1px solid var(--line);background:var(--panel)}
#map{flex:1;min-height:0}
h1{font-size:21px;line-height:1.2;margin:0 0 10px}
p{margin:0 0 10px}
.small{font-size:13px;color:var(--ink2)}
.legend{margin:12px 0}
.legend h2{font-size:14px;margin:0 0 6px}
.row{display:flex;align-items:center;gap:8px;font-size:13px;margin:3px 0}
.sw{width:18px;height:14px;border-radius:2px;flex:none;border:1px solid rgba(0,0,0,.12)}
.sw.out{background:transparent;border:2px solid var(--hot)}
.sw.pro{background:var(--protect);border-color:var(--protect)}
.sw.case{background:transparent;border:2px dashed var(--ink)}
a{color:#1d5fb3}
#toggle{display:none}
.leaflet-popup-content{font-size:14px;line-height:1.4;margin:10px 12px}
.leaflet-popup-content b{font-size:15px}
.case-label{background:rgba(255,255,255,.9);border:1px solid var(--ink);border-radius:3px;padding:1px 5px;font-size:12px;font-weight:600;white-space:nowrap;box-shadow:none}
.case-label:before{display:none}
@media (max-width:760px){
  #app{flex-direction:column-reverse}
  #panel{width:100%;max-height:42%;border-right:0;border-top:1px solid var(--line);padding:12px 16px}
  #panel.collapsed{max-height:52px;overflow:hidden}
  #toggle{display:block;width:100%;margin:0 0 8px;padding:6px;border:1px solid var(--line);border-radius:6px;background:#fff;font:inherit;font-size:14px}
  h1{font-size:18px}
}
</style>
</head>
<body>
<div id="app">
<aside id="panel">
<button id="toggle" aria-expanded="true">Hide details</button>
<h1>Is Lagos building itself hotter?</h1>
<p>Each square on the map is 500 m across. Its colour shows how much hotter or cooler the ground got than the rest of Lagos, comparing the dry seasons of 2016&ndash;18 with 2023&ndash;25.</p>
<p>Across __SQUARES__ squares, places where more buildings went up generally warmed more: about __PER10__ &deg;C extra for every extra 10% of the land covered by buildings.</p>
<p>__HOTSPOTS__ hotspot squares, among the top fifth for both new building and warming, are home to about __RESIDENTS__ people; green squares (__PROTECT__) are still mostly green but are being built on fast.</p>
<p class="small">Data years: buildings 2016&ndash;2023; heat and greenery: dry seasons 2016&ndash;18 vs 2023&ndash;25.</p>
<p class="small">Limit: satellite ground temperature, not air temperature; shows patterns, not proof of cause.</p>
<div class="legend">
<h2>How much hotter the ground got than the rest of Lagos</h2>
<div id="heatlegend"></div>
<div class="row"><span class="sw out"></span>Hotspot square</div>
<div class="row"><span class="sw pro"></span>Protect: still green, being built on</div>
<div class="row"><span class="sw case"></span>Case-study area</div>
</div>
<p class="small">Tap a square or an area for details. Use the layer button (top right) to switch layers on and off.</p>
<p class="small"><a href="__REPO__">Project on GitHub</a> &middot; <a href="__BRIEF__">Policy brief</a></p>
</aside>
<div id="map" role="region" aria-label="Map of Lagos"></div>
</div>
<script>
const D = __DATA__;
const map = L.map('map', {preferCanvas: true, zoomControl: true, zoomSnap: 0.25}).setView([6.55, 3.45], 10);
L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
  maxZoom: 18,
  attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
}).addTo(map);
map.attributionControl.addAttribution('Data: Google Open Buildings, Copernicus Sentinel-2, USGS Landsat, WorldPop, geoBoundaries');

const N = D.nodes, W = D.ncol1;
function node(i){ return [N.lat0 + N.lat[i] / 1e5, N.lon0 + N.lon[i] / 1e5]; }
function square(r, c){ return [node(r*W+c), node(r*W+c+1), node((r+1)*W+c+1), node((r+1)*W+c)]; }
function cls(v){ let k = 0; while (k < D.bins.length && v >= D.bins[k]) k++; return k; }
function fmt(v){ const s = Math.abs(v).toFixed(1); return v >= 0.05 ? '+' + s + ' &deg;C hotter' : (v <= -0.05 ? s + ' &deg;C cooler' : 'about the same'); }

const heat = L.layerGroup(), hot = L.layerGroup(), pro = L.layerGroup();
const C = D.cells, renderer = L.canvas({padding: 0.3});
for (let i = 0; i < C.r.length; i++) {
  const v = C.v[i] / 100, sq = square(C.r[i], C.c[i]);
  const tip = '<b>' + fmt(v) + '</b><br>than the rest of Lagos (ground temperature, 2016&ndash;18 to 2023&ndash;25)<br>' + D.lgas[C.l[i]]
    + ((C.f[i] & 1) ? '<br><span style="color:#c0392b">Hotspot square</span>' : '')
    + ((C.f[i] & 2) ? '<br><span style="color:#1e7b34">Protect: still green, being built on</span>' : '');
  L.polygon(sq, {renderer, stroke: false, fillColor: D.colours[cls(v)], fillOpacity: 0.75}).bindPopup(tip).addTo(heat);
  if (C.f[i] & 1) L.polygon(sq, {renderer, color: '#c0392b', weight: 2, fill: false, interactive: false}).addTo(hot);
  if (C.f[i] & 2) L.polygon(sq, {renderer, color: '#1e7b34', weight: 1, fillColor: '#1e7b34', fillOpacity: 0.8}).bindPopup(tip).addTo(pro);
}

function lgaPopup(p){
  if (p.n === undefined) return '<b>' + p.name + '</b><br>No analysed squares.';
  const rank = p.rank ? p.rank + ' of ' + p.of : 'not ranked (no hotspot squares)';
  return '<b>' + p.name + '</b>'
    + '<br>People living in hotspot squares: <b>' + p.res.toLocaleString('en-GB') + '</b> (2023 estimate)'
    + '<br>Hotspot squares: <b>' + p.hot + '</b>'
    + '<br>Average extra warming in its hotspot squares: <b>' + (p.hheat === null ? 'no hotspot squares' : (p.hheat >= 0 ? '+' : '') + p.hheat.toFixed(2) + ' &deg;C') + '</b>'
    + '<br>Average extra warming across all its squares: <b>' + (p.heat >= 0 ? '+' : '') + p.heat.toFixed(2) + ' &deg;C</b> (' + p.n.toLocaleString('en-GB') + ' squares)'
    + '<br><span style="color:#52514e;font-size:13px">Warming is compared with Lagos as a whole. Cooler squares elsewhere in an area can pull its all-squares average down.</span>'
    + '<br>&ldquo;Cool now&rdquo; priority rank: <b>' + rank + '</b>'
    + '<br>&ldquo;Protect&rdquo; squares: <b>' + p.prot + '</b>';
}
const lga = L.geoJSON(D.lgaGeo, {
  style: {color: '#3a3a38', weight: 1.2, fill: true, fillOpacity: 0, opacity: 0.8},
  onEachFeature: (f, layer) => layer.bindPopup(lgaPopup(f.properties))
});

const cases = L.layerGroup(), caseLayers = [];
D.cases.forEach(c => {
  const g = L.geoJSON(c.geom, {style: {color: '#0b0b0b', weight: 2.5, dashArray: '5 4', fill: false}});
  g.bindTooltip('Case ' + c.rank, {permanent: true, direction: 'top', className: 'case-label'});
  caseLayers.push([g, c]);
  g.bindPopup('<b>Case ' + c.rank + ': ' + c.name + '</b><br>' + c.lga + '<br>' + c.squares + ' hotspot squares');
  g.addTo(cases);
});

heat.addTo(map); lga.addTo(map); hot.addTo(map); cases.addTo(map);
L.control.layers(null, {
  'Ground heat change': heat, 'Local government areas': lga, 'Hotspots': hot,
  'Protect: green areas being built on': pro, 'Case-study areas': cases
}, {collapsed: true}).addTo(map);
// narrow screens open on the built-up core; wider screens on the whole state
if (window.innerWidth < 760) map.fitBounds([[6.42, 3.18], [6.72, 3.66]]);
else map.fitBounds(lga.getBounds(), {padding: [10, 10]});
// short case labels when zoomed out, full names when zoomed in
function caseLabels(){
  const full = map.getZoom() >= 11;
  caseLayers.forEach(([g, c]) => g.setTooltipContent('Case ' + c.rank + (full ? ': ' + c.name : '')));
}
map.on('zoomend', caseLabels); caseLabels();

const leg = document.getElementById('heatlegend');
for (let k = D.colours.length - 1; k >= 0; k--) {
  leg.insertAdjacentHTML('beforeend', '<div class="row"><span class="sw" style="background:' + D.colours[k] + '"></span>' + D.binLabels[k] + '</div>');
}
const panel = document.getElementById('panel'), tg = document.getElementById('toggle');
tg.addEventListener('click', () => {
  const c = panel.classList.toggle('collapsed');
  tg.textContent = c ? 'Show details' : 'Hide details';
  tg.setAttribute('aria-expanded', String(!c));
  setTimeout(() => map.invalidateSize(), 50);
});
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
