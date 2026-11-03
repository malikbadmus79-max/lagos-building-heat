"""Phase 4: priority-area lists ("Cool now" and "Protect before it is lost").

Rules are fixed in docs/phase4_priority_rules.md (written before ranking).

Outputs
- data/processed/priority_areas_lga.csv   lists A1 (LGAs) and B summarised by LGA
- data/processed/priority_cells.csv       lists A2 (hotspot clusters) and B (cells)
- data/interim/phase4_nominatim.json      cached Nominatim responses for the top 10s
- data/interim/summaries/phase4.json      rules, counts and top 10s
- figures/phase4/cool_now.png, figures/phase4/protect.png

    python -m lagos_heat.phase4            (Nominatim requests only for names not cached)
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from lagos_heat.common import FIGURES, INTERIM, PROCESSED, SUMMARIES, GridSpec, lagos_lgas

VEG_LATE_MIN = 0.50
TOP = 10
NAME_KEYS = ("neighbourhood", "quarter", "suburb", "residential", "hamlet", "village", "town",
             "city_district")
CONFIRMED_CASE_NAMES = {  # Phase 3, confirmed by the author (docs/phase3_validation.md)
    1: "Lekki Peninsula II (Chevron area)",
    2: "Abule Egba / Iyana Ipaja",
    3: "Oke",
}
NOMINATIM_CACHE = INTERIM / "phase4_nominatim.json"
ZOOMS = (16, 14, 12)  # Nominatim zoom: neighbourhood first, then wider areas
OUT_FIG = FIGURES / "phase4"


# ------------------------------------------------------------------ pure helpers
def rank_table(df, metrics, tiebreak, id_col):
    """Equal-weight mean of per-metric ranks (1 = highest value, ties share the lowest rank),
    ordered by mean rank, then tiebreak (higher first), then id."""
    out = df.copy()
    for m in metrics:
        out[f"rank_{m}"] = out[m].rank(ascending=False, method="min")
    out["mean_rank"] = out[[f"rank_{m}" for m in metrics]].mean(axis=1)
    out = out.sort_values(["mean_rank", tiebreak, id_col], ascending=[True, False, True],
                          kind="mergesort").reset_index(drop=True)
    out["priority_rank"] = np.arange(1, len(out) + 1)
    return out


def neighbourhood_sum(row, col, values, nrow, ncol):
    """Sum of values over each cell and its 8 neighbours (3 x 3), for cells given by row/col.
    `values` is a full (nrow, ncol) array with NaN treated as 0."""
    grid = np.nan_to_num(np.asarray(values, float))
    pad = np.pad(grid, 1)
    tot = np.zeros_like(grid)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            tot += pad[1 + dr:1 + dr + nrow, 1 + dc:1 + dc + ncol]
    return tot[np.asarray(row), np.asarray(col)]


def local_name(address):
    for k in NAME_KEYS:
        if address.get(k):
            return address[k]
    return ""


def cool_now(hot):
    """hot: hotspot cells with lga, cluster, pop_2023, lst_rel_change, ndvi_rel_change_p80."""
    agg = dict(n_hotspot_cells=("cell_id", "size"), residents_hotspot=("pop_2023", "sum"),
               mean_lst_rel_change=("lst_rel_change", "mean"),
               mean_ndvi_rel_change=("ndvi_rel_change_p80", "mean"))
    metrics = ["residents_hotspot", "mean_lst_rel_change", "mean_ndvi_loss"]
    lga = hot.groupby("lga").agg(**agg).reset_index()
    lga["mean_ndvi_loss"] = -lga["mean_ndvi_rel_change"]
    lga = rank_table(lga, metrics, "residents_hotspot", "lga")
    cl = hot.groupby("cluster").agg(**agg, lgas=("lga", lambda s: "; ".join(
        s.value_counts().index))).reset_index()
    cl["lga"] = cl["lgas"].str.split("; ").str[0]
    cl["mean_ndvi_loss"] = -cl["mean_ndvi_rel_change"]
    cl = rank_table(cl, metrics, "residents_hotspot", "cluster")
    return lga, cl


def protect(cells, cut):
    """cells: main-sample cells with veg_share_late, bshare_change, residents_nearby, ..."""
    sel = cells[(cells["veg_share_late"] >= VEG_LATE_MIN) & (cells["bshare_change"] >= cut)]
    ranked = rank_table(sel, ["bshare_change", "residents_nearby"], "bshare_change", "cell_id")
    by_lga = (ranked.groupby("lga")
              .agg(n_cells=("cell_id", "size"), mean_bshare_change=("bshare_change", "mean"),
                   max_bshare_change=("bshare_change", "max"),
                   residents_in_cells=("pop_2023", "sum"),
                   best_cell_rank=("priority_rank", "min"))
              .reset_index()
              .sort_values(["n_cells", "residents_in_cells"], ascending=False, kind="mergesort")
              .reset_index(drop=True))
    by_lga["priority_rank"] = np.arange(1, len(by_lga) + 1)
    return ranked, by_lga


# ------------------------------------------------------------------ names
def _norm(s):
    return "".join(ch for ch in str(s).lower() if ch.isalpha())


def label_name(name, osm_lga, lga, widened):
    """(label, note) for a place name. The name is marked "near" when it came from a wider
    lookup, or when OpenStreetMap places it in a different LGA than the boundary data
    (geoBoundaries), whose LGA is kept."""
    if not name:
        return "", "no named place found"
    notes = []
    if widened:
        notes.append("nearest named place from a wider Nominatim lookup")
    if osm_lga and _norm(osm_lga) != _norm(lga):
        notes.append(f"OpenStreetMap places {name} in {osm_lga}; the cell centre lies in "
                     f"{lga} in the boundary data")
    return (f"near {name}" if notes else name), "; ".join(notes)


def names_for(points):
    """points: [(key, lat, lon)]. Nominatim at zoom 16; where that gives no local name, zoom 14
    and then 12 (wider areas). Cached; at most 1 request per second.
    Returns {key: (name, osm_lga, widened)}."""
    from lagos_heat.validate import nominatim

    cache = json.loads(NOMINATIM_CACHE.read_text(encoding="utf-8")) if NOMINATIM_CACHE.exists() \
        else {}
    last, out = 0.0, {}
    for key, lat, lon in points:
        for zoom in ZOOMS:
            ck = key if zoom == ZOOMS[0] else f"{key}_z{zoom}"
            if ck not in cache or "error" in cache[ck]:
                wait = 1.1 - (time.time() - last)
                if wait > 0:
                    time.sleep(wait)
                last = time.time()
                try:
                    cache[ck] = nominatim(lat, lon, zoom)
                except Exception as e:  # recorded; the name is left blank
                    cache[ck] = {"error": str(e)}
                print(f"  nominatim {ck}", flush=True)
            addr = cache[ck].get("address", {})
            name = local_name(addr)
            if name:
                break
        out[key] = (name, addr.get("county", ""), bool(name) and zoom != ZOOMS[0])
    NOMINATIM_CACHE.write_text(json.dumps(cache, indent=1), encoding="utf-8")
    return out


# ------------------------------------------------------------------ data
def load():
    s = pd.read_csv(PROCESSED / "phase2_sample.csv",
                    usecols=["cell_id", "main", "hotspot", "cluster"])
    g = pd.read_csv(INTERIM / "grid_periods.csv",
                    usecols=["cell_id", "row", "col", "lga", "lon", "lat", "pop_2023",
                             "lst_rel_change", "ndvi_rel_change_p80", "veg_share_late",
                             "bshare_2016", "bshare_2023", "bshare_change"])
    return g.merge(s, on="cell_id", how="left")


# ------------------------------------------------------------------ maps
def _base(ax, spec):
    from lagos_heat.phase2 import SURFACE

    x0, y0, x1, y1 = lagos_lgas().total_bounds
    ax.set_facecolor(SURFACE)
    ax.set_xlim(x0 - 1000, x1 + 1000)
    ax.set_ylim(y0 - 1000, y1 + 1000)
    ax.set_axis_off()


def _rank_labels(ax, t):
    import matplotlib.patheffects as pe

    from lagos_heat.phase2 import TEXT_PRIMARY

    for _, c in t.iterrows():
        ax.annotate(str(int(c["priority_rank"])), (c["x"], c["y"]), xytext=(3, 3),
                    textcoords="offset points", fontsize=7.5, color=TEXT_PRIMARY,
                    fontweight="bold", path_effects=[pe.withStroke(linewidth=2.2,
                                                                   foreground="white")])


def _side(fig, lines, x=0.80):
    from lagos_heat.phase2 import TEXT_PRIMARY

    fig.text(x, 0.93, "\n".join(lines), fontsize=7.2, va="top", color=TEXT_PRIMARY,
             family="DejaVu Sans", linespacing=1.45)


def map_cool_now(df, lga_rank, clusters, spec, path):
    from lagos_heat.phase2 import HOT, OUTLINE, SURFACE, TEXT_PRIMARY, TEXT_SECONDARY, _plt

    plt = _plt()
    import matplotlib.colors as mcolors

    lgas = lagos_lgas().merge(lga_rank[["lga", "priority_rank"]], left_on="lga", right_on="lga",
                              how="left")
    fig = plt.figure(figsize=(13, 4.6), dpi=150, facecolor=SURFACE)
    ax = fig.add_axes([0.0, 0.08, 0.78, 0.84])
    _base(ax, spec)
    n = int(lga_rank["priority_rank"].max())
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "r", ["#9c2526", "#e34948", "#f3a5a1", "#f6e3e1"])
    lgas[lgas["priority_rank"].isna()].plot(ax=ax, color="#eeede8", edgecolor=OUTLINE,
                                            linewidth=0.4)
    lgas[lgas["priority_rank"].notna()].plot(ax=ax, column="priority_rank", cmap=cmap, vmin=1,
                                             vmax=n, edgecolor=OUTLINE, linewidth=0.4, alpha=0.55)
    code = np.full((spec.nrow, spec.ncol), np.nan)
    h = df[df["hotspot"] == True]  # noqa: E712
    code[h["row"], h["col"]] = 1
    x0, y0, x1, y1 = spec.bounds
    ax.imshow(code, cmap=mcolors.ListedColormap([HOT]), extent=(x0, x1, y0, y1),
              interpolation="nearest")
    _rank_labels(ax, clusters.head(TOP))
    lines = ["Top 10 LGAs", *[f"{int(r.priority_rank):>2}. {r.lga}" for r in
                              lga_rank.head(TOP).itertuples()], "",
             "Top 10 hotspot clusters (map numbers)",
             *[f"{int(r.priority_rank):>2}. {r.area_name or '(no name)'}, {r.lga}" for r in
               clusters.head(TOP).itertuples()]]
    _side(fig, lines)
    fig.text(0.01, 0.95, '"Cool now": LGA rank (shading), hotspot cells (red) and the top 10 '
             "hotspot clusters (numbers)", fontsize=10, color=TEXT_PRIMARY)
    fig.text(0.01, 0.03, "LGA shading: darker = higher priority (equal-weight ranks of residents "
             "in hotspot cells, mean relative LST change, mean relative NDVI loss). "
             "Unshaded LGAs have no hotspot cells. Red cells: hotspot cells. 500 m cells, "
             "EPSG:32631.", fontsize=7, color=TEXT_SECONDARY)
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def map_protect(df, ranked, by_lga, spec, path):
    from lagos_heat.phase2 import OUTLINE, SEQ_GREEN, SURFACE, TEXT_PRIMARY, TEXT_SECONDARY, _plt

    plt = _plt()
    import matplotlib.colors as mcolors

    fig = plt.figure(figsize=(15, 4.6), dpi=150, facecolor=SURFACE)
    ax = fig.add_axes([0.0, 0.08, 0.66, 0.84])
    inset = fig.add_axes([0.665, 0.12, 0.16, 0.70])
    veg = np.full((spec.nrow, spec.ncol), np.nan)
    m = df["main"] == True  # noqa: E712
    veg[df.loc[m, "row"], df.loc[m, "col"]] = df.loc[m, "veg_share_late"]
    code = np.full((spec.nrow, spec.ncol), np.nan)
    code[ranked["row"], ranked["col"]] = 1
    x0, y0, x1, y1 = spec.bounds
    top = ranked.head(TOP)
    # inset: the top-10 cells within 8 km of the first-ranked cell, padded by 1.5 km
    near = top[np.hypot(top["x"] - top["x"].iloc[0], top["y"] - top["y"].iloc[0]) <= 8000]
    ib = (near["x"].min() - 1500, near["y"].min() - 1500, near["x"].max() + 1500,
          near["y"].max() + 1500)
    for a in (ax, inset):
        _base(a, spec)
        a.imshow(veg, cmap=mcolors.LinearSegmentedColormap.from_list("g", SEQ_GREEN[:5]),
                 vmin=0, vmax=1, extent=(x0, x1, y0, y1), interpolation="nearest", alpha=0.8)
        lagos_lgas().boundary.plot(ax=a, color=OUTLINE, linewidth=0.4)
        a.imshow(code, cmap=mcolors.ListedColormap(["#6b3fa0"]), extent=(x0, x1, y0, y1),
                 interpolation="nearest")
    inset.set_xlim(ib[0], ib[2])
    inset.set_ylim(ib[1], ib[3])
    inset.set_axis_on()
    inset.set_xticks([])
    inset.set_yticks([])
    inset.set_xlabel("")
    inset.set_ylabel("")
    for sp in inset.spines.values():
        sp.set_edgecolor(TEXT_SECONDARY)
    inset.set_title("Inset", fontsize=7.5, color=TEXT_SECONDARY, loc="left")
    ax.add_patch(plt.Rectangle((ib[0], ib[1]), ib[2] - ib[0], ib[3] - ib[1], fill=False,
                               edgecolor=TEXT_SECONDARY, linewidth=0.8))
    _rank_labels(ax, top[~top.index.isin(near.index)])
    _rank_labels(inset, near)
    lines = ["Top 10 cells (map numbers)",
             *[f"{int(r.priority_rank):>2}. {r.area_name if isinstance(r.area_name, str) and r.area_name else '(no name)'}, {r.lga}"
               for r in ranked.head(TOP).itertuples()], "", "Qualifying cells by LGA",
             *[f"{r.lga}: {int(r.n_cells)}" for r in by_lga.itertuples()]]
    _side(fig, lines, x=0.84)
    fig.text(0.01, 0.95, '"Protect before it is lost": cells still mostly vegetated in 2023-25 '
             "with building growth above the hotspot cut", fontsize=10, color=TEXT_PRIMARY)
    fig.text(0.01, 0.03, "Purple: qualifying cells; numbers: top 10 by equal-weight ranks of "
             "building change and residents within 1.5 km x 1.5 km. Green: late-period vegetated "
             "share (main sample). 500 m cells, EPSG:32631.", fontsize=7, color=TEXT_SECONDARY)
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


# ------------------------------------------------------------------ main
def main():
    spec = GridSpec.load()
    p2 = json.loads((SUMMARIES / "phase2.json").read_text())
    cut = p2["exposure"]["quintile_cut_bshare_change"]
    df = load()
    df["x"] = spec.xmin + (df["col"] + 0.5) * spec.cell
    df["y"] = spec.ymax - (df["row"] + 0.5) * spec.cell
    main_df = df[df["main"] == True].copy()  # noqa: E712

    hot = main_df[main_df["hotspot"] == True]  # noqa: E712
    a1, a2 = cool_now(hot)
    cases = pd.read_csv(PROCESSED / "phase2_case_studies.csv")
    a2 = a2.merge(cases[["cluster", "lon", "lat", "xmin", "ymin", "xmax", "ymax", "rank",
                         "chosen"]].rename(columns={"rank": "size_rank"}), on="cluster",
                  how="left")
    a2["x"] = (a2["xmin"] + a2["xmax"]) / 2
    a2["y"] = (a2["ymin"] + a2["ymax"]) / 2

    pop = np.full((spec.nrow, spec.ncol), np.nan)
    pop[df["row"], df["col"]] = df["pop_2023"]
    main_df["residents_nearby"] = neighbourhood_sum(main_df["row"], main_df["col"], pop,
                                                    spec.nrow, spec.ncol)
    b, b_lga = protect(main_df, cut)

    # names for the top 10s
    pts = [(f"cluster_{int(c.cluster)}", c.lat, c.lon) for c in a2.head(TOP).itertuples()]
    pts += [(f"cell_{int(c.cell_id)}", c.lat, c.lon) for c in b.head(TOP).itertuples()]
    nm = names_for(pts)
    a2["area_name"], a2["name_note"], a2["osm_lga"] = "", "", ""
    for i, c in a2.head(TOP).iterrows():
        case_rank = int(c["size_rank"]) if bool(c["chosen"]) else None
        name, osm_lga, widened = nm[f"cluster_{int(c.cluster)}"]
        a2.loc[i, "osm_lga"] = osm_lga
        if case_rank in CONFIRMED_CASE_NAMES:
            a2.loc[i, "area_name"] = CONFIRMED_CASE_NAMES[case_rank]
            a2.loc[i, "name_source"] = "Phase 3 confirmed"
            continue
        # a cluster may span LGAs: no conflict if OSM's LGA is any of them
        lgas = c["lgas"].split("; ")
        ref = next((x for x in lgas if _norm(x) == _norm(osm_lga)), c["lga"])
        a2.loc[i, "area_name"], a2.loc[i, "name_note"] = label_name(name, osm_lga, ref, widened)
        a2.loc[i, "name_source"] = "Nominatim (OpenStreetMap)"
    b["area_name"], b["name_note"], b["osm_lga"] = "", "", ""
    for i, c in b.head(TOP).iterrows():
        name, osm_lga, widened = nm[f"cell_{int(c.cell_id)}"]
        b.loc[i, "osm_lga"] = osm_lga
        b.loc[i, "area_name"], b.loc[i, "name_note"] = label_name(name, osm_lga, c["lga"], widened)
        b.loc[i, "name_source"] = "Nominatim (OpenStreetMap)"

    # write tables
    a1_out = a1.assign(list="A1_cool_now_lga")
    b_lga_out = b_lga.assign(list="B_protect_by_lga")
    pd.concat([a1_out, b_lga_out], ignore_index=True).to_csv(
        PROCESSED / "priority_areas_lga.csv", index=False)
    a2_out = a2.assign(list="A2_cool_now_cluster", unit="cluster").rename(
        columns={"cluster": "unit_id"})
    b_out = b.assign(list="B_protect_cell", unit="cell").rename(columns={"cell_id": "unit_id"})
    keep_a2 = ["list", "unit", "unit_id", "priority_rank", "mean_rank", "lga", "lgas",
               "area_name", "name_source", "name_note", "osm_lga", "n_hotspot_cells", "residents_hotspot",
               "mean_lst_rel_change", "mean_ndvi_rel_change", "mean_ndvi_loss",
               "rank_residents_hotspot", "rank_mean_lst_rel_change", "rank_mean_ndvi_loss",
               "lon", "lat"]
    keep_b = ["list", "unit", "unit_id", "priority_rank", "mean_rank", "lga", "area_name", "name_note", "osm_lga",
              "name_source", "veg_share_late", "bshare_2016", "bshare_2023", "bshare_change",
              "residents_nearby", "pop_2023", "lst_rel_change", "ndvi_rel_change_p80",
              "rank_bshare_change", "rank_residents_nearby", "lon", "lat"]
    pd.concat([a2_out[keep_a2], b_out[keep_b]], ignore_index=True).to_csv(
        PROCESSED / "priority_cells.csv", index=False)

    OUT_FIG.mkdir(parents=True, exist_ok=True)
    map_cool_now(df, a1, a2, spec, OUT_FIG / "cool_now.png")
    map_protect(df, b, b_lga, spec, OUT_FIG / "protect.png")

    def top(t, cols):
        return t.head(TOP)[cols].to_dict(orient="records")

    summ = {
        "rules": "docs/phase4_priority_rules.md",
        "hotspot_cells": int(len(hot)), "hotspot_lgas": int(len(a1)),
        "hotspot_clusters": int(len(a2)),
        "protect_rule": {"veg_share_late_min": VEG_LATE_MIN, "bshare_change_min": cut},
        "protect_cells": int(len(b)), "protect_lgas": int(len(b_lga)),
        "protect_residents_in_cells": float(b["pop_2023"].sum()),
        "A1_top10": top(a1, ["priority_rank", "lga", "n_hotspot_cells", "residents_hotspot",
                             "mean_lst_rel_change", "mean_ndvi_rel_change", "mean_rank"]),
        "A2_top10": top(a2, ["priority_rank", "cluster", "lga", "lgas", "area_name", "name_note",
                             "n_hotspot_cells", "residents_hotspot", "mean_lst_rel_change",
                             "mean_ndvi_rel_change", "mean_rank"]),
        "B_top10": top(b, ["priority_rank", "cell_id", "lga", "area_name", "name_note", "veg_share_late",
                           "bshare_change", "residents_nearby", "mean_rank"]),
        "B_by_lga": b_lga.to_dict(orient="records"),
    }
    (SUMMARIES / "phase4.json").write_text(json.dumps(summ, indent=2, default=float))
    print(json.dumps({k: summ[k] for k in ("hotspot_cells", "hotspot_lgas", "hotspot_clusters",
                                           "protect_cells", "protect_lgas")}))


if __name__ == "__main__":
    main()
