"""Step 7: assemble the cell- and LGA-level analysis tables, quick-look maps and checks.

Inputs: the per-cell CSVs written by grid.py, buildings.py, vegetation.py,
temperature.py, population.py and covariates.py.

Outputs
- data/interim/grid_periods.csv   one row per cell
- data/interim/lga_periods.csv    LGA area-weighted means and population totals
- figures/phase1/*.png            quick-look maps
- data/interim/summaries/assemble.json   row counts, flags, range checks
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from lagos_heat.common import (
    FIGURES, INTERIM, MIN_OBS, GridSpec, lagos_lgas, load_grid, write_summary,
)

MIN_OB_COVERAGE = 0.5  # share of a cell observed by Open Buildings, each year

# Plausible ranges for the range check (values outside are counted, not removed).
RANGES = {
    "share": (0.0, 1.0), "index": (-1.0, 1.0), "lst_c": (15.0, 65.0),
    "height_m": (0.0, 150.0), "pop": (0.0, np.inf), "dist_m": (0.0, np.inf),
}
RANGE_COLUMNS = {
    "share": ["bshare_2016", "bshare_2023", "veg_share_early", "veg_share_late",
              "veg_share_early_screened", "veg_share_late_screened", "bare_share_early",
              "bare_share_late", "bare_share_early_screened", "bare_share_late_screened",
              "ghsl_share_2015", "ghsl_share_2020",
              "tree_frac_wc2021", "water_frac", "industrial_share", "frac_in_lagos"],
    "index": ["ndvi_early_p80", "ndvi_late_p80", "ndvi_early", "ndvi_late", "ndvi_early_noscreen",
              "ndvi_late_noscreen", "bsi_early", "bsi_late", "bsi_early_noscreen",
              "bsi_late_noscreen"],
    "lst_c": ["lst_early", "lst_late"],
    "height_m": ["bheight_2016", "bheight_2023"],
    "pop": ["pop_2016", "pop_2023"],
    "dist_m": ["dist_ocean_m", "dist_lagoon_m", "dist_major_road_m"],
}


def build_table():
    grid = load_grid()
    base = pd.DataFrame(grid.drop(columns="geometry"))
    keep = ["cell_id", "row", "col", "lga", "lon", "lat", "frac_in_lagos", "water_frac",
            "water_flag", "tree_frac_wc2021", "built_frac_wc2021"]
    df = base[keep]
    for name in ("buildings", "vegetation", "temperature", "population", "covariates"):
        part = pd.read_csv(INTERIM / f"{name}_cells.csv")
        part = part.drop(columns=[c for c in part.columns if c in df.columns and c != "cell_id"])
        df = df.merge(part, on="cell_id", how="left", validate="one_to_one")

    df["bshare_change"] = df["bshare_2023"] - df["bshare_2016"]
    df["bheight_change"] = df["bheight_2023"] - df["bheight_2016"]
    df["ghsl_share_change"] = df["ghsl_share_2020"] - df["ghsl_share_2015"]
    df["ob_low_coverage"] = ~((df["ob_coverage_2016"] >= MIN_OB_COVERAGE) &
                              (df["ob_coverage_2023"] >= MIN_OB_COVERAGE))
    # The primary vegetation measure (_p80) uses all dates, so its observation count is the
    # all-date count (s2_nobs_*_noscreen). s2_low_obs_{period} refers to the screened composite.
    for p in ("early", "late"):
        df[f"s2_low_obs_{p}_p80"] = ~(df[f"s2_nobs_{p}_noscreen"] >= MIN_OBS)
    df["low_obs_any"] = (df["s2_low_obs_early_p80"] | df["s2_low_obs_late_p80"] |
                         df["lst_low_obs_early"] | df["lst_low_obs_late"])
    df["analysis_ok"] = ~df["water_flag"] & ~df["low_obs_any"] & ~df["ob_low_coverage"]
    # Sample for the screened-median sensitivity check of vegetation change.
    df["ndvi_screened_ok"] = (~df["water_flag"] & ~df["s2_low_obs_early"] &
                              ~df["s2_low_obs_late"])
    lead = keep + ["analysis_ok", "lst_ok", "ndvi_screened_ok", "low_obs_any", "ob_low_coverage"]
    return df[lead + [c for c in df.columns if c not in lead]], grid


def lga_table(df):
    """Area-weighted LGA means over non-water cells; weight = cell area in Lagos x land share."""
    land = df[~df["water_flag"]].copy()
    land["w"] = land["frac_in_lagos"] * (1 - land["water_frac"])
    value_cols = ["bshare_2016", "bshare_2023", "bshare_change", "bheight_2016",
                  "bheight_2023", "ghsl_share_2015", "ghsl_share_2020", "ndvi_early_p80",
                  "ndvi_late_p80", "ndvi_rel_early_p80", "ndvi_rel_late_p80",
                  "ndvi_rel_change_p80", "ndvi_rel_change_noscreen", "ndvi_early", "ndvi_late",
                  "ndvi_change", "ndvi_rel_early", "ndvi_rel_late", "ndvi_rel_change",
                  "veg_share_early", "veg_share_late",
                  "veg_share_change", "bare_share_early", "bare_share_late",
                  "bare_share_change", "bsi_early_noscreen", "bsi_late_noscreen",
                  "bsi_change_noscreen", "bsi_early", "bsi_late", "bsi_change", "lst_early",
                  "lst_late", "lst_change", "lst_rel_early", "lst_rel_late", "lst_rel_change",
                  "tree_frac_wc2021", "industrial_share", "dist_ocean_m", "dist_lagoon_m",
                  "dist_major_road_m"]
    rows = []
    for lga, g in land.groupby("lga"):
        r = {"lga": lga, "n_cells": int((df["lga"] == lga).sum()), "n_land_cells": len(g),
             "n_analysis_ok": int(g["analysis_ok"].sum()),
             "land_area_km2": float((g["w"] * 0.25).sum())}
        for c in value_cols:
            ok = g[c].notna() & (g["w"] > 0)
            r[c] = float(np.average(g.loc[ok, c], weights=g.loc[ok, "w"])) if ok.any() else np.nan
        rows.append(r)
    out = pd.DataFrame(rows)
    pop = pd.read_csv(INTERIM / "population_lga.csv")
    return out.merge(pop, on="lga", how="left")


def range_checks(df):
    res = {}
    for kind, cols in RANGE_COLUMNS.items():
        lo, hi = RANGES[kind]
        for c in cols:
            v = df[c]
            res[c] = {"min": float(v.min()), "max": float(v.max()),
                      "n_missing": int(v.isna().sum()),
                      "n_outside": int(((v < lo) | (v > hi)).sum()), "range": [lo, hi]}
    return res


NDVI_VERSIONS = {  # version: (column suffix, observation-count column suffix)
    "p80": ("_p80", "_noscreen"), "screened": ("", ""), "noscreen": ("_noscreen", "_noscreen")}


def ndvi_coverage(df, land):
    """Per NDVI version, land cells with a value in both periods, and with >= MIN_OBS too."""
    n = int(land.sum())
    res = {}
    for v, (suf, osuf) in NDVI_VERSIONS.items():
        both = land & df[f"ndvi_early{suf}"].notna() & df[f"ndvi_late{suf}"].notna()
        obs = both & (df[f"s2_nobs_early{osuf}"] >= MIN_OBS) & (df[f"s2_nobs_late{osuf}"] >= MIN_OBS)
        res[v] = {"land_cells": n, "valid_both": int(both.sum()),
                  "valid_both_share": float(both.sum() / n),
                  "valid_both_min_obs": int(obs.sum()),
                  "valid_both_min_obs_share": float(obs.sum() / n)}
    return res


# ---------------------------------------------------------------- figures
SURFACE = "#fcfcfb"
NODATA_FILL = "#e6e5e0"
TEXT_PRIMARY, TEXT_SECONDARY, OUTLINE = "#0b0b0b", "#52514e", "#8a8984"
DIVERGING = ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f3a5a1", "#e34948", "#9c2526"]
SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def cell_raster(df, spec, col):
    arr = np.full((spec.nrow, spec.ncol), np.nan)
    v = df[col].astype(float).where(~df["water_flag"])
    arr[df["row"].to_numpy(), df["col"].to_numpy()] = v.to_numpy()
    return arr


def quicklook(df, spec, col, title, label, path, diverging=True, log=False):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    arr = cell_raster(df, spec, col)
    finite = arr[np.isfinite(arr)]
    if diverging:
        lim = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 1.0
        cmap = mcolors.LinearSegmentedColormap.from_list("div", DIVERGING)
        norm = mcolors.TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
        note = f"Colour limits: +/-{lim:.3g} (98th percentile of |value|); extremes clipped."
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list("seq", SEQUENTIAL)
        pos = finite[finite > 0]
        norm = (mcolors.LogNorm(vmin=max(float(np.percentile(pos, 2)), 1.0),
                                vmax=float(np.percentile(pos, 99.5)))
                if log else None)
        note = "Log colour scale (2nd to 99.5th percentile of cells with population > 0)."
    cmap.set_bad(NODATA_FILL)

    x0, y0, x1, y1 = spec.bounds
    fig, ax = plt.subplots(figsize=(12, 3.6), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    im = ax.imshow(arr, cmap=cmap, norm=norm, extent=(x0, x1, y0, y1), interpolation="nearest")
    lagos_lgas().boundary.plot(ax=ax, color=OUTLINE, linewidth=0.4)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_axis_off()
    ax.set_title(title, loc="left", color=TEXT_PRIMARY, fontsize=11)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label(label, color=TEXT_SECONDARY, fontsize=8)
    cb.ax.tick_params(labelsize=7, colors=TEXT_SECONDARY)
    cb.outline.set_visible(False)
    fig.text(0.01, 0.02, f"500 m cells, EPSG:32631. Grey: water-flagged cells or no data. "
             f"{note}", fontsize=7, color=TEXT_SECONDARY)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main():
    df, grid = build_table()
    spec = GridSpec.load()
    df.to_csv(INTERIM / "grid_periods.csv", index=False)
    lga = lga_table(df)
    lga.to_csv(INTERIM / "lga_periods.csv", index=False)

    fig_dir = FIGURES / "phase1"
    quicklook(df, spec, "bshare_change", "Change in building share, 2016 to 2023 "
              "(Open Buildings 2.5D Temporal)", "share of cell area, 2023 minus 2016",
              fig_dir / "building_share_change.png")
    for stale in ("ndvi_change.png", "ndvi_change_noscreen.png"):  # earlier versions
        (fig_dir / stale).unlink(missing_ok=True)
    quicklook(df, spec, "ndvi_rel_change_p80", "Change in relative NDVI, early to late "
              "(80th percentile of NDVI minus Lagos land median per date, all dates)",
              "relative NDVI, late minus early", fig_dir / "relative_ndvi_change.png")
    quicklook(df.assign(v=df["ndvi_change"].where(df["ndvi_screened_ok"])), spec, "v",
              "Change in median dry-season NDVI, early to late (reference-forest screened "
              "dates; cells with at least 3 observations in both periods)",
              "NDVI, late minus early", fig_dir / "ndvi_change_screened.png")
    quicklook(df, spec, "ndvi_rel_change_noscreen", "Change in relative NDVI, early to late "
              "(median of NDVI minus Lagos land median per date, all dates)",
              "relative NDVI, late minus early", fig_dir / "relative_ndvi_change_noscreen.png")
    quicklook(df, spec, "lst_rel_change", "Change in relative land surface temperature, "
              "early to late (Landsat 8; late period also Landsat 9)",
              "deg C relative to Lagos median, late minus early",
              fig_dir / "relative_lst_change.png")
    quicklook(df, spec, "pop_2023", "Population 2023 (WorldPop R2025A, constrained)",
              "people per 500 m cell", fig_dir / "population_2023.png", diverging=False,
              log=True)

    land = ~df["water_flag"]
    write_summary("assemble", {
        "grid_periods_rows": int(len(df)), "grid_periods_columns": int(df.shape[1]),
        "lga_periods_rows": int(len(lga)),
        "cells_water_flag": int(df["water_flag"].sum()), "cells_land": int(land.sum()),
        "land_cells_s2_low_obs_early_p80": int(df.loc[land, "s2_low_obs_early_p80"].sum()),
        "land_cells_s2_low_obs_late_p80": int(df.loc[land, "s2_low_obs_late_p80"].sum()),
        "land_cells_s2_low_obs_early_screened": int(df.loc[land, "s2_low_obs_early"].sum()),
        "land_cells_s2_low_obs_late_screened": int(df.loc[land, "s2_low_obs_late"].sum()),
        "ndvi_versions_land_cells": ndvi_coverage(df, land),
        "land_cells_lst_low_obs_early": int(df.loc[land, "lst_low_obs_early"].sum()),
        "land_cells_lst_low_obs_late": int(df.loc[land, "lst_low_obs_late"].sum()),
        "land_cells_low_obs_any": int(df.loc[land, "low_obs_any"].sum()),
        "land_cells_low_obs_any_share": float(df.loc[land, "low_obs_any"].mean()),
        "land_cells_ob_low_coverage": int(df.loc[land, "ob_low_coverage"].sum()),
        "cells_analysis_ok": int(df["analysis_ok"].sum()),
        "cells_ndvi_screened_ok": int(df["ndvi_screened_ok"].sum()),
        "cells_analysis_ok_if_screened_obs": int((
            df["analysis_ok"] & ~df["s2_low_obs_early"] & ~df["s2_low_obs_late"]).sum()),
        "cells_lst_ok": int(df["lst_ok"].sum()),
        "land_cells_lst_ok": int(df.loc[land, "lst_ok"].sum()),
        "land_cells_ob_overlap_approx": int(df.loc[land, "ob_overlap_approx_2016"].sum()),
        "land_cells_sandfill_candidate": int(df.loc[land, "sandfill_candidate"].sum()),
        "land_cells_sandfill_candidate_screened": int(
            df.loc[land, "sandfill_candidate_screened"].sum()),
        "min_obs_threshold": MIN_OBS, "min_ob_coverage": MIN_OB_COVERAGE,
        "range_checks": range_checks(df),
        "lagos_population_lga_sum": {str(y): float(lga[f"pop_{y}"].sum()) for y in (2016, 2023)},
        "lga_pop_growth_2016_2023": {
            k: float(getattr(lga["pop_2023"] / lga["pop_2016"] - 1, k)())
            for k in ("min", "median", "max")},
        "figures": sorted(p.name for p in fig_dir.glob("*.png")),
    })
    print(f"assemble: {len(df)} cells, {len(lga)} LGAs")


if __name__ == "__main__":
    main()
