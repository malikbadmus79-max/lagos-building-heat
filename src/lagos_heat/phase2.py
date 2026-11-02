"""Phase 2: analysis as fixed in docs/phase2_analysis_plan.md (section numbers below refer
to that plan). Works on the cell table; the only raster reads are small windows for the
case-study maps, plus one GHSL GHS-POP tile for the exposure cross-check.

Outputs
- data/processed/phase2_sample.csv        per-cell sample flags, block, group, hotspot, cluster
- data/processed/phase2_deciles.csv       decile means with block-bootstrap 95% intervals
- data/processed/phase2_models.csv        every model and sensitivity run
- data/processed/phase2_h3.csv            building coefficient with and without NDVI change
- data/processed/phase2_moran.csv         Moran's I of residuals
- data/processed/phase2_h4_groups.csv     H4 groups: sizes, building-change spread, slopes
- data/processed/phase2_exposure_lga.csv  population in hotspot cells by LGA
- data/processed/phase2_case_studies.csv  hotspot clusters, ranked, with the chosen cases
- data/interim/summaries/phase2.json
- figures/phase2/*.png
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from lagos_heat.common import (
    FIGURES, INTERIM, PROCESSED, RAW, UTM31N, GridSpec, aggregate, fetch, lagos_lgas,
    write_summary,
)

# ---------------------------------------------------------------- fixed in the plan
WATER_MAX = 0.10            # 1. main sample: water_frac <= 0.10
BLOCK = 5                   # 4.3 cluster blocks of 5 x 5 cells
N_BOOT = 2000               # 4.1, 4.2 block bootstrap resamples
N_PERM = 999                # 4.3 Moran's I permutations
SEED = 20261102
BUILT_T, TREE_T, VEG_T = 0.30, 0.50, 0.50   # 3. H4 groups
H4_BUILT_ALT, H4_TREE_ALT = (0.20, 0.40), (0.40, 0.60)
COAST_M, ROAD_M = 1000.0, 250.0             # 5. exclusion distances
LST_OBS_STRICT = 5                          # 5. observation-count sensitivity
HOT_Q = 0.8                                 # 6. top quintile
CASE_THIRD_MIN, CASE_PAD_M = 4, 1000.0      # 7. case-study rule and map padding
STEP = 0.10                                 # 4.4 effect per 10 percentage points

DIST = ["dist_ocean_km", "dist_lagoon_km", "dist_major_road_km"]
M1_CONTROLS = ["bshare_2016", "ndvi_rel_early_p80", "lst_rel_early", *DIST,
               "industrial_share", "bare_share_change"]
M2_CONTROLS = ["bshare_2016", "ndvi_rel_early_p80", *DIST, "industrial_share"]

GHSPOP_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A"
GHSPOP_STEM = "GHS_POP_E2020_GLOBE_R2023A_54009_100"
GHSL_TILE = "R9_C19"
OUT_FIG = FIGURES / "phase2"


# ---------------------------------------------------------------- pure functions
def block_ids(row, col, size=BLOCK):
    """Integer id of the size x size block of each cell (row // size, col // size)."""
    row, col = np.asarray(row), np.asarray(col)
    return (row // size) * 100000 + (col // size)


def rank_deciles(x, n=10, rng=None):
    """Deciles 1..n by rank, so groups have equal size +/- 1. Tied values are ordered at
    random (fixed seed), not by file order, which follows the grid rows (north to south)."""
    rng = rng or np.random.default_rng(SEED)
    x = np.asarray(x, float)
    order = np.lexsort((rng.random(len(x)), x))
    r = np.empty(len(x))
    r[order] = np.arange(1, len(x) + 1)
    return np.asarray(pd.qcut(r, n, labels=list(range(1, n + 1)))).astype(int)


def bootstrap_weights(blocks, n_boot, rng):
    """(n_boot x n_cells) weights: how often each cell's block is drawn in a resample of
    all blocks with replacement."""
    codes, uniq = pd.factorize(np.asarray(blocks))
    g = len(uniq)
    draws = rng.integers(0, g, size=(n_boot, g))
    counts = np.stack([np.bincount(d, minlength=g) for d in draws])
    return counts[:, codes]


def group_means_ci(y, groups, blocks, n_boot=N_BOOT, rng=None, level=0.95):
    """Mean of y per group with block-bootstrap percentile intervals."""
    rng = rng or np.random.default_rng(SEED)
    y, groups = np.asarray(y, float), np.asarray(groups)
    labels = np.unique(groups)
    codes, uniq = pd.factorize(np.asarray(blocks))
    g = len(uniq)
    s = np.zeros((g, len(labels)))
    c = np.zeros((g, len(labels)))
    for j, lab in enumerate(labels):
        m = groups == lab
        s[:, j] = np.bincount(codes[m], weights=y[m], minlength=g)
        c[:, j] = np.bincount(codes[m], minlength=g)
    draws = rng.integers(0, g, size=(n_boot, g))
    k = np.stack([np.bincount(d, minlength=g) for d in draws])
    with np.errstate(invalid="ignore", divide="ignore"):
        boot = (k @ s) / (k @ c)
    a = (1 - level) / 2
    return pd.DataFrame({"group": labels, "n": c.sum(0).astype(int),
                         "mean": s.sum(0) / c.sum(0),
                         "ci_lo": np.nanquantile(boot, a, axis=0),
                         "ci_hi": np.nanquantile(boot, 1 - a, axis=0)})


def rook_weights(row, col):
    """Row-standardised rook-contiguity weights (scipy CSR) between the given grid cells,
    and a mask of cells that have at least one neighbour."""
    from scipy import sparse

    row, col = np.asarray(row), np.asarray(col)
    index = {(r, c): i for i, (r, c) in enumerate(zip(row.tolist(), col.tolist()))}
    i_idx, j_idx = [], []
    for i, (r, c) in enumerate(zip(row.tolist(), col.tolist())):
        for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            j = index.get((r + dr, c + dc))
            if j is not None:
                i_idx.append(i)
                j_idx.append(j)
    n = len(row)
    w = sparse.csr_matrix((np.ones(len(i_idx)), (i_idx, j_idx)), shape=(n, n))
    deg = np.asarray(w.sum(1)).ravel()
    inv = np.divide(1.0, deg, out=np.zeros(n), where=deg > 0)
    return sparse.diags(inv) @ w, deg > 0


def moran_i(values, row, col, n_perm=N_PERM, rng=None):
    """Moran's I with rook weights (row-standardised); cells without neighbours are left
    out. One-sided permutation p-value for positive autocorrelation."""
    rng = rng or np.random.default_rng(SEED)
    values, row, col = np.asarray(values, float), np.asarray(row), np.asarray(col)
    _, has = rook_weights(row, col)
    w, _ = rook_weights(row[has], col[has])
    z = values[has] - values[has].mean()
    n, s0 = len(z), w.sum()
    zz = z @ z

    def stat(v):
        return n / s0 * (v @ (w @ v)) / zz

    i_obs = stat(z)
    perm = np.array([stat(rng.permutation(z)) for _ in range(n_perm)])
    return {"moran_i": float(i_obs), "expected": -1.0 / (n - 1),
            "p_perm": float((1 + (perm >= i_obs).sum()) / (n_perm + 1)),
            "n": int(n), "n_islands": int((~has).sum()), "n_perm": n_perm}


def weighted_coefs(x, y, w):
    """Weighted least-squares coefficients for each weight vector in w (n_boot x n);
    rank-deficient designs use the pseudo-inverse."""
    out = np.empty((w.shape[0], x.shape[1]))
    for b, wb in enumerate(w):
        xtw = x.T * wb
        out[b] = np.linalg.pinv(xtw @ x) @ (xtw @ y)
    return out


def top_quintile_both(a, b, q=HOT_Q):
    """Cells at or above the q-quantile of a and of b, with the two cut points."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ca, cb = float(np.quantile(a, q)), float(np.quantile(b, q))
    return (a >= ca) & (b >= cb), ca, cb


def label_clusters(row, col, flag, nrow, ncol):
    """Rook-contiguous clusters of flagged cells; 0 = not flagged."""
    from scipy import ndimage

    grid = np.zeros((nrow, ncol), bool)
    grid[np.asarray(row)[flag], np.asarray(col)[flag]] = True
    lab, _ = ndimage.label(grid)  # default structure = rook (4-neighbour)
    return lab[np.asarray(row), np.asarray(col)] * np.asarray(flag)


def choose_cases(table, third_min=CASE_THIRD_MIN):
    """Rank clusters by cells, then summed building change; the two largest, plus the third
    if it has at least third_min cells. table: cluster, n_cells, sum_bshare_change."""
    t = table.sort_values(["n_cells", "sum_bshare_change"], ascending=False)
    t = t.reset_index(drop=True)
    t["rank"] = np.arange(1, len(t) + 1)
    chosen = t["rank"] <= 2
    if len(t) >= 3 and t.loc[2, "n_cells"] >= third_min:
        chosen |= t["rank"] == 3
    t["chosen"] = chosen
    return t


def sign_label(lo, hi):
    if lo > 0:
        return "positive"
    if hi < 0:
        return "negative"
    return "not distinguishable from zero"


# ---------------------------------------------------------------- data
def load_table():
    df = pd.read_csv(INTERIM / "grid_periods.csv")
    l8 = pd.read_csv(INTERIM / "temperature_cells_landsat8_only.csv")
    l8 = l8[["cell_id", "lst_rel_early", "lst_rel_change", "lst_ok"]].rename(
        columns=lambda c: c if c == "cell_id" else f"{c}_l8")
    df = df.merge(l8, on="cell_id", how="left", validate="one_to_one")
    for c in ("ocean", "lagoon", "major_road"):
        df[f"dist_{c}_km"] = df[f"dist_{c}_m"] / 1000.0
    df["block"] = block_ids(df["row"], df["col"])
    df["base"] = df["analysis_ok"].astype(bool)
    df["main"] = df["base"] & (df["water_frac"] <= WATER_MAX)
    df["main_nosand"] = df["main"] & ~df["sandfill_candidate"].astype(bool)
    return df


# ---------------------------------------------------------------- models
def fit(df, mask, outcome, xvar, controls, fe=True):
    """OLS with block-clustered SEs (t with G-1 df) on complete cases."""
    import statsmodels.formula.api as smf

    cols = [outcome, xvar, *controls, "block", "lga"]
    d = df.loc[mask, cols + ["row", "col", "cell_id"]].dropna(subset=cols)
    # A control that is constant in the sample (e.g. industrial_share after excluding
    # industrial cells) is collinear with the intercept and is left out.
    used = [c for c in controls if d[c].nunique() > 1]
    rhs = " + ".join([xvar, *used] + (["C(lga)"] if fe else []))
    res = smf.ols(f"{outcome} ~ {rhs}", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(d["block"])[0]}, use_t=True)
    res.dropped_controls = [c for c in controls if c not in used]
    return res, d


def run_row(run, hyp, res, d, xvar, outcome, note=""):
    lo, hi = res.conf_int().loc[xvar]
    b = res.params[xvar]
    return {"run": run, "hypothesis": hyp, "outcome": outcome, "term": xvar,
            "coef": b, "se": res.bse[xvar], "ci_lo": lo, "ci_hi": hi, "p": res.pvalues[xvar],
            "per_10pp": b * STEP, "per_10pp_ci_lo": lo * STEP, "per_10pp_ci_hi": hi * STEP,
            "sign": sign_label(lo, hi), "n": int(res.nobs),
            "n_blocks": int(d["block"].nunique()), "r2": res.rsquared, "note": note,
            "constant_controls_dropped": "; ".join(getattr(res, "dropped_controls", []))}


def h1_runs(df):
    """M1 and every section-5 rerun of it: (run, mask, outcome, xvar, controls, note)."""
    main = df["main"]
    m1 = M1_CONTROLS
    terc = pd.qcut(df.loc[main, "lst_rel_early"], 3, labels=[1, 2, 3])
    df["lst_early_tercile"] = terc.reindex(df.index)
    l8_controls = [c if c != "lst_rel_early" else "lst_rel_early_l8" for c in m1]
    runs = [
        ("H1_primary", main, "lst_rel_change", "bshare_change", m1, "M1"),
        ("H1_excl_sandfill", df["main_nosand"], "lst_rel_change", "bshare_change", m1,
         "sand filling"),
        ("H1_excl_coast_lagoon_1km", main & (df["dist_ocean_m"] >= COAST_M)
         & (df["dist_lagoon_m"] >= COAST_M), "lst_rel_change", "bshare_change", m1,
         "coast and lagoon"),
        ("H1_excl_industrial", main & (df["industrial_share"] == 0), "lst_rel_change",
         "bshare_change", m1, "industry"),
        ("H1_excl_road_250m", main & (df["dist_major_road_m"] >= ROAD_M), "lst_rel_change",
         "bshare_change", m1, "roads"),
        ("H1_landsat8_only", main & (df["lst_ok_l8"] == True), "lst_rel_change_l8",  # noqa: E712
         "bshare_change", l8_controls, "sensor mix"),
        ("H1_obs_controls", main, "lst_rel_change", "bshare_change",
         m1 + ["lst_nobs_early", "lst_nobs_late"], "observation counts"),
        ("H1_obs_min5", main & (df["lst_nobs_early"] >= LST_OBS_STRICT)
         & (df["lst_nobs_late"] >= LST_OBS_STRICT), "lst_rel_change", "bshare_change", m1,
         "observation counts"),
        ("H1_no_baseline_lst", main, "lst_rel_change", "bshare_change",
         [c for c in m1 if c != "lst_rel_early"], "regression to the mean"),
    ]
    for t in (1, 2, 3):
        runs.append((f"H1_lst_early_tercile{t}", main & (df["lst_early_tercile"] == t),
                     "lst_rel_change", "bshare_change", m1, "regression to the mean"))
    runs.append(("H1_ghsl_2015_2020", main, "lst_rel_change", "ghsl_share_change", m1,
                 "building measure"))
    return runs


def h2_runs(df):
    main = df["main"]
    ns = [c if c != "ndvi_rel_early_p80" else "ndvi_rel_early_noscreen" for c in M2_CONTROLS]
    return [
        ("H2_primary", main, "ndvi_rel_change_p80", "bshare_change", M2_CONTROLS, "M2"),
        ("H2_excl_sandfill", df["main_nosand"], "ndvi_rel_change_p80", "bshare_change",
         M2_CONTROLS, "sand filling"),
        ("H2_noscreen", main, "ndvi_rel_change_noscreen", "bshare_change", ns, "haze"),
        ("H2_bare_control", main, "ndvi_rel_change_p80", "bshare_change",
         M2_CONTROLS + ["bare_share_change"], "bare-share control"),
    ]


def h3(df, version, rng):
    """M1 and M3 (M1 + NDVI change) on the same cells; bootstrap of the difference."""
    from patsy import dmatrices

    if version == "p80":
        ctrl, nd = M1_CONTROLS, "ndvi_rel_change_p80"
    else:
        ctrl = [c if c != "ndvi_rel_early_p80" else "ndvi_rel_early_noscreen"
                for c in M1_CONTROLS]
        nd = "ndvi_rel_change_noscreen"
    cols = ["lst_rel_change", "bshare_change", nd, *ctrl, "block", "lga"]
    common = df["main"] & df[cols].notna().all(axis=1)
    r1, d1 = fit(df, common, "lst_rel_change", "bshare_change", ctrl)
    r3, d3 = fit(df, common, "lst_rel_change", "bshare_change", ctrl + [nd])
    assert len(d1) == len(d3)
    b1, b3 = r1.params["bshare_change"], r3.params["bshare_change"]
    rhs = " + ".join(["bshare_change", *ctrl, "C(lga)"])
    assert (d1.index == d3.index).all()
    y, x1 = dmatrices(f"lst_rel_change ~ {rhs}", d3, return_type="dataframe")
    _, x3 = dmatrices(f"lst_rel_change ~ {rhs} + {nd}", d3, return_type="dataframe")
    w = bootstrap_weights(d1["block"].to_numpy(), N_BOOT, rng)
    j1, j3 = list(x1.columns).index("bshare_change"), list(x3.columns).index("bshare_change")
    bb1 = weighted_coefs(x1.to_numpy(), y.to_numpy().ravel(), w)[:, j1]
    bb3 = weighted_coefs(x3.to_numpy(), y.to_numpy().ravel(), w)[:, j3]
    diff, share = bb1 - bb3, (bb1 - bb3) / bb1
    lo1, hi1 = r1.conf_int().loc["bshare_change"]
    share_ok = lo1 > 0 or hi1 < 0
    row = {"version": version, "ndvi_change_var": nd, "n": len(d1),
           "n_blocks": int(d1["block"].nunique()),
           "coef_m1": b1, "coef_m1_ci_lo": lo1, "coef_m1_ci_hi": hi1,
           "coef_m3": b3, "coef_m3_ci_lo": r3.conf_int().loc["bshare_change"][0],
           "coef_m3_ci_hi": r3.conf_int().loc["bshare_change"][1],
           "ndvi_change_coef_m3": r3.params[nd],
           "ndvi_change_ci_lo_m3": r3.conf_int().loc[nd][0],
           "ndvi_change_ci_hi_m3": r3.conf_int().loc[nd][1],
           "difference": b1 - b3, "difference_ci_lo": np.quantile(diff, 0.025),
           "difference_ci_hi": np.quantile(diff, 0.975),
           "share_reported": share_ok,
           "share": (b1 - b3) / b1 if share_ok else np.nan,
           "share_ci_lo": np.quantile(share, 0.025) if share_ok else np.nan,
           "share_ci_hi": np.quantile(share, 0.975) if share_ok else np.nan,
           "n_boot": N_BOOT}
    for k in ("coef_m1", "coef_m1_ci_lo", "coef_m1_ci_hi", "coef_m3", "coef_m3_ci_lo",
              "coef_m3_ci_hi", "difference", "difference_ci_lo", "difference_ci_hi"):
        row[f"{k}_per_10pp"] = row[k] * STEP
    return row, (r1, d1), (r3, d3)


def h4_groups(df, built_t=BUILT_T, tree_t=TREE_T, veg_alt=False):
    b = df["bshare_2016"] >= built_t
    if veg_alt:
        built = b & (df["tree_frac_wc2021"] < tree_t)
        veg = (df["veg_share_early"] >= VEG_T) & ~b
    else:
        t = df["tree_frac_wc2021"] >= tree_t
        built, veg = b & ~t, t & ~b
    if veg_alt:
        both = b & (df["veg_share_early"] >= VEG_T)
    else:
        both = b & (df["tree_frac_wc2021"] >= tree_t)
    g = np.where(built, "built", np.where(veg, "vegetated", "neither"))
    return pd.Series(g, index=df.index), int((df["main"] & both).sum())


def h4(df, name, groups, n_both):
    import statsmodels.formula.api as smf

    d = df.assign(h4=groups)
    for g in ("built", "vegetated", "neither"):
        d[f"bx_{g}"] = d["bshare_change"] * (d["h4"] == g)
    cols = ["lst_rel_change", "bshare_change", *M1_CONTROLS, "block", "lga"]
    d = d.loc[d["main"]].dropna(subset=cols)
    rhs = " + ".join(["bx_built", "bx_vegetated", "bx_neither", "C(h4)", *M1_CONTROLS,
                      "C(lga)"])
    res = smf.ols(f"lst_rel_change ~ {rhs}", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(d["block"])[0]}, use_t=True)
    ci = res.conf_int()
    rows = []
    for g in ("built", "vegetated", "neither"):
        x = d.loc[d["h4"] == g, "bshare_change"]
        lo, hi = ci.loc[f"bx_{g}"]
        rows.append({"definition": name, "group": g, "n": len(x),
                     "bshare_change_p10": x.quantile(0.1), "bshare_change_p25": x.quantile(0.25),
                     "bshare_change_p50": x.median(), "bshare_change_p75": x.quantile(0.75),
                     "bshare_change_p90": x.quantile(0.9),
                     "slope": res.params[f"bx_{g}"], "slope_ci_lo": lo, "slope_ci_hi": hi,
                     "per_10pp": res.params[f"bx_{g}"] * STEP,
                     "per_10pp_ci_lo": lo * STEP, "per_10pp_ci_hi": hi * STEP,
                     "n_cells_meeting_both_rules": n_both})
    tt = res.t_test("bx_built - bx_vegetated = 0")
    lo, hi = np.asarray(tt.conf_int()).ravel()
    diff = {"definition": name, "difference_built_minus_vegetated": float(tt.effect[0]),
            "ci_lo": float(lo), "ci_hi": float(hi), "p": float(tt.pvalue),
            "sign": sign_label(lo, hi), "n": int(res.nobs),
            "n_blocks": int(d["block"].nunique()),
            "difference_per_10pp": float(tt.effect[0]) * STEP,
            "difference_per_10pp_ci_lo": float(lo) * STEP,
            "difference_per_10pp_ci_hi": float(hi) * STEP}
    return rows, diff, d


# ---------------------------------------------------------------- exposure
def ghs_pop_cells(spec):
    """GHS-POP 2020 (R2023A, 100 m) summed per cell over pixels whose centres fall in it."""
    import rasterio
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds

    name = f"{GHSPOP_STEM}_V1_0_{GHSL_TILE}"
    local = RAW / "ghsl" / f"{name}.zip"
    if not local.exists():
        fetch(f"{GHSPOP_BASE}/{GHSPOP_STEM}/V1-0/tiles/{name}.zip", local)
    with rasterio.open(f"/vsizip/{local.as_posix()}/{name}.tif") as src:
        b = transform_bounds(UTM31N, src.crs, *spec.bounds, densify_pts=50)
        win = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
        arr = src.read(1, window=win, boundless=True, fill_value=src.nodata or 0)
        cells = spec.pixel_cells_reprojected(src.window_transform(win), arr.shape, src.crs)
        ok = (arr != src.nodata) if src.nodata is not None else np.ones(arr.shape, bool)
        ok &= arr > 0
    total, _ = aggregate(cells, values=np.where(ok, arr, 0), mask=ok, ncells=spec.ncells)
    return total, name


def exposure(df, hot):
    main = df["main"]
    rows = []
    for lga, g in df[main].groupby("lga"):
        h = hot[g.index]
        rows.append({"lga": lga, "main_cells": len(g), "hotspot_cells": int(h.sum()),
                     "worldpop2023_main": g["pop_2023"].sum(),
                     "worldpop2023_hotspot": g.loc[h, "pop_2023"].sum(),
                     "ghspop2020_main": g["ghs_pop_2020"].sum(),
                     "ghspop2020_hotspot": g.loc[h, "ghs_pop_2020"].sum()})
    t = pd.DataFrame(rows)
    tot = t.drop(columns="lga").sum()
    t = pd.concat([t.sort_values("worldpop2023_hotspot", ascending=False),
                   pd.DataFrame([{"lga": "All LGAs", **tot.to_dict()}])], ignore_index=True)
    t[["main_cells", "hotspot_cells"]] = t[["main_cells", "hotspot_cells"]].astype(int)
    t["worldpop2023_hotspot_share"] = t["worldpop2023_hotspot"] / t["worldpop2023_main"]
    t["ghspop2020_hotspot_share"] = t["ghspop2020_hotspot"] / t["ghspop2020_main"]
    return t


# ---------------------------------------------------------------- figures
SURFACE = "#fcfcfb"
NODATA_FILL = "#e6e5e0"
TEXT_PRIMARY, TEXT_SECONDARY, OUTLINE, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
ACCENT = "#256abf"
HOT = "#c0392b"
DIVERGING = ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f3a5a1", "#e34948", "#9c2526"]
SEQ_GREEN = ["#f3f1e7", "#d5e4b8", "#a9cc86", "#74ad5a", "#468b3a", "#256a26", "#0f4a17"]


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": OUTLINE,
                         "axes.labelcolor": TEXT_SECONDARY, "xtick.color": TEXT_SECONDARY,
                         "ytick.color": TEXT_SECONDARY, "axes.titlecolor": TEXT_PRIMARY})
    return plt


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def decile_chart(dec, var, ylabel, title, path, n_main):
    plt = _plt()
    d = dec[dec["outcome"] == var]
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=150, facecolor=SURFACE)
    _style(ax)
    ax.axhline(0, color=OUTLINE, linewidth=0.8)
    # one position per decile; tick labels give the decile's mean building change
    xs = d["decile"].to_numpy()
    ax.errorbar(xs, d["mean"], yerr=[d["mean"] - d["ci_lo"], d["ci_hi"] - d["mean"]],
                fmt="o", color=ACCENT, ecolor=ACCENT, elinewidth=1.5, capsize=0, markersize=6,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.set_xticks(xs, [f"D{k}\n{v:+.1f}" for k, v in zip(xs, d["x_mean_pp"])], fontsize=7)
    ax.tick_params(axis="x", length=0)
    ax.set_xlabel("Decile of building-share change, with its mean change "
                  "(percentage points, 2023 minus 2016)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=10)
    fig.text(0.01, -0.09, f"{n_main:,} main-sample cells in deciles of equal size (tied values "
             f"ordered at random). Points: decile means.\nBars: 95% block-bootstrap intervals "
             f"({N_BOOT:,} resamples of {BLOCK} x {BLOCK}-cell blocks).", fontsize=7,
             color=TEXT_SECONDARY)
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def coef_chart(models, path):
    plt = _plt()
    panels = [("H1", "°C of relative LST per 10 pp building share"),
              ("H2", "relative NDVI per 10 pp building share")]
    sizes = [int((models["hypothesis"] == h).sum()) for h, _ in panels]
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 0.32 * sum(sizes) + 2.0), dpi=150,
                             facecolor=SURFACE, gridspec_kw={"height_ratios": sizes})
    for ax, (h, lab) in zip(axes, panels):
        d = models[models["hypothesis"] == h].iloc[::-1]
        y = np.arange(len(d))
        ax.set_facecolor(SURFACE)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="x", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.axvline(0, color=OUTLINE, linewidth=0.8)
        prim = d["run"].str.endswith("_primary").to_numpy()
        ax.hlines(y, d["per_10pp_ci_lo"], d["per_10pp_ci_hi"], color=ACCENT, linewidth=2)
        ax.scatter(d["per_10pp"], y, s=np.where(prim, 46, 26), color=ACCENT, zorder=3,
                   edgecolor=SURFACE, linewidth=1.2)
        ax.set_yticks(y, [f"{r}  (n={n:,})" for r, n in zip(d["run"], d["n"])], fontsize=7)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel(lab)
        ax.set_title(f"{h}: building-change coefficient across runs (95% intervals, "
                     f"block-clustered)", loc="left", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def hotspot_map(df, spec, cases, path):
    plt = _plt()
    import matplotlib.colors as mcolors

    code = np.full((spec.nrow, spec.ncol), np.nan)
    land = ~df["water_flag"]
    for val, m in ((0, land & ~df["main"]), (1, df["main"]), (2, df["hotspot"])):
        code[df.loc[m, "row"], df.loc[m, "col"]] = val
    cmap = mcolors.ListedColormap(["#d9d7d0", "#eef3fa", HOT])
    cmap.set_bad(SURFACE)
    x0, y0, x1, y1 = spec.bounds
    fig, ax = plt.subplots(figsize=(12, 3.8), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.imshow(code, cmap=cmap, vmin=-0.5, vmax=2.5, extent=(x0, x1, y0, y1),
              interpolation="nearest")
    lagos_lgas().boundary.plot(ax=ax, color=OUTLINE, linewidth=0.4)
    for _, c in cases[cases["chosen"]].iterrows():
        ax.add_patch(plt.Rectangle((c["xmin"], c["ymin"]), c["xmax"] - c["xmin"],
                                   c["ymax"] - c["ymin"], fill=False, edgecolor=TEXT_PRIMARY,
                                   linewidth=0.9))
        ax.text(c["xmin"], c["ymax"] + 300, f"Case {int(c['rank'])}", fontsize=7,
                color=TEXT_PRIMARY)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_axis_off()
    ax.set_title("Hotspot cells: top quintile of relative LST change and of building-share "
                 "change (main sample)", loc="left", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (HOT, "#eef3fa", "#d9d7d0")]
    ax.legend(handles, ["Hotspot cell", "Main sample, not hotspot", "Land, outside main sample"],
              loc="upper left", bbox_to_anchor=(0.0, 0.02), ncol=3, fontsize=7,
              frameon=False)
    fig.text(0.01, 0.02, "500 m cells, EPSG:32631. White: water-flagged cells. Boxes: "
             "case-study clusters.", fontsize=7, color=TEXT_SECONDARY)
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def read_window(path, bounds):
    import rasterio
    from rasterio.windows import from_bounds

    with rasterio.open(path) as src:
        win = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
        arr = src.read(1, window=win, boundless=True, fill_value=np.nan).astype("float32")
        if src.nodata is not None and not np.isnan(src.nodata):
            arr[arr == src.nodata] = np.nan
        t = src.window_transform(win)
    ext = (t.c, t.c + t.a * arr.shape[1], t.f + t.e * arr.shape[0], t.f)
    return arr, ext


def ob_presence_window(year, bounds, tiles, res=4.0):
    """Open Buildings presence (> 0.5) on a res-metre grid covering bounds (window reads of
    the tiles {uri: footprint} that intersect bounds). Cached in WORK."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import Window, from_bounds
    from shapely.geometry import box

    from lagos_heat.buildings import OB_NODATA, PRESENCE_THRESHOLD
    from lagos_heat.common import WORK, http_opener, retry

    x0, y0, x1, y1 = bounds
    cache = WORK / f"phase2_ob_{year}_{int(x0)}_{int(y0)}_{int(x1)}_{int(y1)}.npy"
    if cache.exists():
        return np.load(cache), (x0, x1, y0, y1)
    ny, nx = int(round((y1 - y0) / res)), int(round((x1 - x0) / res))
    out = np.full((ny, nx), np.nan, dtype="float32")
    win_box = box(*bounds)
    for uri in sorted(u for u, poly in tiles.items() if poly.intersects(win_box)):
        url = uri.replace("gs://", "https://storage.googleapis.com/")

        def _read():
            with rasterio.open(url, opener=http_opener) as src:
                band = list(src.descriptions).index("building_presence") + 1
                win = from_bounds(*bounds, transform=src.transform)
                win = win.intersection(Window(0, 0, src.width, src.height))
                win = win.round_offsets().round_lengths()
                left, bottom, right, top = src.window_bounds(win)
                c0, r0 = int(round((left - x0) / res)), int(round((y1 - top) / res))
                nc, nr = int(round((right - left) / res)), int(round((top - bottom) / res))
                if nc <= 0 or nr <= 0:
                    return None
                arr = src.read(band, window=win, out_shape=(nr, nc),
                               resampling=Resampling.nearest)
                return r0, c0, arr

        got = retry(_read, attempts=4, base=5, cap=30)
        print(f"    OB {year} {uri.split('/')[-1]} read", flush=True)
        if got is None:
            continue
        r0, c0, arr = got
        rr = slice(max(r0, 0), min(r0 + arr.shape[0], ny))
        cc = slice(max(c0, 0), min(c0 + arr.shape[1], nx))
        sub = arr[rr.start - r0:rr.stop - r0, cc.start - c0:cc.stop - c0]
        val = np.where(sub == OB_NODATA, np.nan, (sub > PRESENCE_THRESHOLD).astype("float32"))
        tgt = out[rr, cc]
        fill = np.isnan(tgt) & np.isfinite(val)
        tgt[fill] = val[fill]
    np.save(cache, out)
    return out, (x0, x1, y0, y1)


def case_map(case, cells, tiles, path):
    plt = _plt()
    import matplotlib.colors as mcolors

    bounds = (case["xmin"] - CASE_PAD_M, case["ymin"] - CASE_PAD_M,
              case["xmax"] + CASE_PAD_M, case["ymax"] + CASE_PAD_M)
    rows = [("Building presence (Open Buildings, 4 m)", "building", None),
            ("80th-percentile NDVI (Sentinel-2, 20 m)", "ndvi", None),
            ("Relative LST, °C (Landsat, 30 m)", "lst", None)]
    data = {}
    for year, per in ((2016, "early"), (2023, "late")):
        data[("building", per)] = ob_presence_window(year, bounds, tiles[year])
        data[("ndvi", per)] = read_window(PROCESSED / f"s2_{per}_ndvi_p80_20m.tif", bounds)
        data[("lst", per)] = read_window(PROCESSED / f"lst_{per}_rel_median_30m.tif", bounds)
    lst_all = np.concatenate([data[("lst", p)][0][np.isfinite(data[("lst", p)][0])]
                              for p in ("early", "late")])
    lim = float(np.percentile(np.abs(lst_all), 98)) if lst_all.size else 1.0
    styles = {
        "building": dict(cmap=mcolors.ListedColormap(["#f0efec", "#3a3a38"]), vmin=0, vmax=1),
        "ndvi": dict(cmap=mcolors.LinearSegmentedColormap.from_list("g", SEQ_GREEN),
                     vmin=0.0, vmax=0.8),
        "lst": dict(cmap=mcolors.LinearSegmentedColormap.from_list("d", DIVERGING),
                    norm=mcolors.TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)),
    }
    fig, axes = plt.subplots(3, 2, figsize=(8.5, 11), dpi=150, facecolor=SURFACE)
    for i, (label, key, _) in enumerate(rows):
        for j, (per, ttl) in enumerate((("early", "Early"), ("late", "Late"))):
            ax = axes[i, j]
            arr, ext = data[(key, per)]
            st = dict(styles[key])
            cmap = st.pop("cmap").copy()
            cmap.set_bad(NODATA_FILL)
            im = ax.imshow(np.ma.masked_invalid(arr), extent=ext, cmap=cmap,
                           interpolation="nearest", **st)
            for _, c in cells.iterrows():
                ax.add_patch(plt.Rectangle((c["x0"], c["y0"]), 500, 500, fill=False,
                                           edgecolor=HOT, linewidth=0.6))
            ax.set_xlim(bounds[0], bounds[2])
            ax.set_ylim(bounds[1], bounds[3])
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_color(OUTLINE)
                s.set_linewidth(0.5)
            yr = {"building": ("2016", "2023"), "ndvi": ("Dec 2016–Feb 2018", "Dec 2023–Feb 2025"),
                  "lst": ("Dec 2016–Feb 2018", "Dec 2023–Feb 2025")}[key][j]
            ax.set_title(f"{label}\n{ttl}: {yr}", loc="left", fontsize=8)
        cb = fig.colorbar(im, ax=axes[i, :].tolist(), fraction=0.03, pad=0.02)
        cb.ax.tick_params(labelsize=7, colors=TEXT_SECONDARY)
        cb.outline.set_visible(False)
        if key == "building":
            cb.set_ticks([0.25, 0.75], labels=["no building", "building"])
    fig.suptitle(f"Case {int(case['rank'])}: {case['lga']} ({int(case['n_cells'])} hotspot "
                 f"cells, outlined)", x=0.02, ha="left", fontsize=11, color=TEXT_PRIMARY)
    fig.text(0.02, 0.06, "Window: cluster extent plus 1 km. Grey: no data. Absolute "
             "80th-percentile NDVI is a level within each period; its early-late difference "
             "includes a haze difference (docs/phase1_notes.md). LST colour limits: ±"
             f"{lim:.2g} °C (98th percentile of |value| in the window).", fontsize=7,
             color=TEXT_SECONDARY, wrap=True)
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main(skip_maps=False):
    rng = np.random.default_rng(SEED)
    spec = GridSpec.load()
    df = load_table()
    PROCESSED.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    summ = {"plan": "docs/phase2_analysis_plan.md", "seed": SEED}

    # 1. samples
    main_ = df["main"]
    summ["samples"] = {
        "base_analysis_ok": int(df["base"].sum()),
        "removed_water_frac_gt_0_10": int((df["base"] & ~main_).sum()),
        "main": int(main_.sum()),
        "sandfill_candidates_in_main": int((main_ & df["sandfill_candidate"]).sum()),
        "main_excl_sandfill": int(df["main_nosand"].sum()),
        "ndvi_screened_ok": int(df["ndvi_screened_ok"].sum()),
    }
    iqr = df.loc[main_, "bshare_change"].quantile([0.25, 0.5, 0.75])
    summ["bshare_change_main"] = {"p25": iqr[0.25], "median": iqr[0.5], "p75": iqr[0.75],
                                  "iqr": iqr[0.75] - iqr[0.25],
                                  "mean": df.loc[main_, "bshare_change"].mean()}

    # 4.1 deciles
    dec_rows = []
    for var in ("lst_rel_change", "ndvi_rel_change_p80"):
        d = df.loc[main_, ["bshare_change", var, "block"]].dropna()
        d["decile"] = rank_deciles(d["bshare_change"])
        g = group_means_ci(d[var], d["decile"], d["block"], rng=rng)
        g["x_mean_pp"] = d.groupby("decile")["bshare_change"].mean().to_numpy() * 100
        g["x_min_pp"] = d.groupby("decile")["bshare_change"].min().to_numpy() * 100
        g["x_max_pp"] = d.groupby("decile")["bshare_change"].max().to_numpy() * 100
        dec_rows.append(g.rename(columns={"group": "decile"}).assign(outcome=var))
    dec = pd.concat(dec_rows, ignore_index=True)
    dec.to_csv(PROCESSED / "phase2_deciles.csv", index=False)

    # 4.2 M1, M2 and section-5 runs
    rows, resid = [], {}
    for hyp, runs in (("H1", h1_runs(df)), ("H2", h2_runs(df))):
        for run, mask, outcome, xvar, ctrl, note in runs:
            res, d = fit(df, mask, outcome, xvar, ctrl)
            rows.append(run_row(run, hyp, res, d, xvar, outcome, note))
            if run in ("H1_primary", "H2_primary"):
                resid[run] = (res.resid, d)
    # haze, descriptive only: screened median change on ndvi_screened_ok cells, no controls
    scr = df["ndvi_screened_ok"].astype(bool)
    res, d = fit(df, scr, "ndvi_change", "bshare_change", [], fe=False)
    desc = run_row("H2_screened_descriptive", "H2-descriptive", res, d, "bshare_change",
                   "ndvi_change", "haze, descriptive only (no controls)")
    desc["pearson_r"] = float(d["ndvi_change"].corr(d["bshare_change"]))
    rows.append(desc)

    # H3
    h3_rows = []
    for version in ("p80", "noscreen"):
        row, (r1, d1), (r3, d3) = h3(df, version, rng)
        h3_rows.append(row)
        rows.append(run_row(f"H3_M1_{version}", "H3", r1, d1, "bshare_change",
                            "lst_rel_change", "M1 on the M3 cells"))
        rows.append(run_row(f"H3_M3_{version}", "H3", r3, d3, "bshare_change",
                            "lst_rel_change", f"M1 + {row['ndvi_change_var']}"))
        if version == "p80":
            resid["H3_M3_p80"] = (r3.resid, d3)
    h3t = pd.DataFrame(h3_rows)
    h3t.to_csv(PROCESSED / "phase2_h3.csv", index=False)

    # H4
    h4_rows, h4_diffs = [], []
    defs = [("main", dict()), ("built_0.20", dict(built_t=H4_BUILT_ALT[0])),
            ("built_0.40", dict(built_t=H4_BUILT_ALT[1])),
            ("tree_0.40", dict(tree_t=H4_TREE_ALT[0])),
            ("tree_0.60", dict(tree_t=H4_TREE_ALT[1])),
            ("vegetated_by_veg_share_early", dict(veg_alt=True))]
    for name, kw in defs:
        groups, n_both = h4_groups(df, **kw)
        if name == "main":
            df["h4_group"] = groups.where(main_)
        r, diff, d = h4(df, name, groups, n_both)
        h4_rows += r
        h4_diffs.append(diff)
        for g in ("built", "vegetated"):
            m = main_ & (groups == g)
            res, dd = fit(df, m, "lst_rel_change", "bshare_change", M1_CONTROLS)
            rows.append(run_row(f"H4_{name}_{g}_separate", "H4", res, dd, "bshare_change",
                                "lst_rel_change", f"M1 within {g} group ({name})"))
    h4t = pd.DataFrame(h4_rows).merge(pd.DataFrame(h4_diffs), on="definition",
                                      suffixes=("", "_diff"))
    h4t.to_csv(PROCESSED / "phase2_h4_groups.csv", index=False)
    models = pd.DataFrame(rows)
    models.to_csv(PROCESSED / "phase2_models.csv", index=False)

    # 4.3 Moran's I
    mor = []
    for run, (r, d) in resid.items():
        mor.append({"run": run, **moran_i(r.to_numpy(), d["row"], d["col"], rng=rng)})
    mor = pd.DataFrame(mor)
    mor.to_csv(PROCESSED / "phase2_moran.csv", index=False)

    # robustness per the plan's decision rules
    def robust(sub):
        s = np.sign(sub["coef"])
        return bool((s == s.iloc[0]).all()), sub["run"][s != s.iloc[0]].tolist()

    h1m = models[models["hypothesis"] == "H1"]
    h2m = models[models["hypothesis"] == "H2"]
    summ["H1"] = {"primary": models.loc[models["run"] == "H1_primary"].iloc[0].to_dict(),
                  "robust_sign_all_runs": robust(h1m)[0], "runs_with_other_sign": robust(h1m)[1],
                  "runs_ci_excludes_zero": int((h1m["sign"] != "not distinguishable from zero").sum()),
                  "n_runs": len(h1m)}
    summ["H2"] = {"primary": models.loc[models["run"] == "H2_primary"].iloc[0].to_dict(),
                  "robust_sign_all_runs": robust(h2m)[0], "runs_with_other_sign": robust(h2m)[1],
                  "p80_and_noscreen_agree": bool(np.sign(h2m.set_index("run").loc["H2_primary", "coef"])
                                                 == np.sign(h2m.set_index("run").loc["H2_noscreen", "coef"])),
                  "screened_descriptive": desc, "n_runs": len(h2m)}
    summ["H3"] = {r["version"]: r for r in h3_rows}
    summ["H3"]["difference_sign_agrees"] = bool(np.sign(h3_rows[0]["difference"])
                                                == np.sign(h3_rows[1]["difference"]))
    hd = pd.DataFrame(h4_diffs)
    summ["H4"] = {"differences": h4_diffs,
                  "robust_sign_all_definitions": bool(
                      (np.sign(hd["difference_built_minus_vegetated"]) ==
                       np.sign(hd["difference_built_minus_vegetated"].iloc[0])).all())}
    summ["moran"] = mor.to_dict("records")

    # 6. exposure
    m = df.loc[main_]
    hot, c_lst, c_b = top_quintile_both(m["lst_rel_change"], m["bshare_change"])
    df["hotspot"] = False
    df.loc[m.index, "hotspot"] = hot
    ghs, ghs_name = ghs_pop_cells(spec)
    df["ghs_pop_2020"] = ghs[df["cell_id"].to_numpy()]
    expo = exposure(df, df["hotspot"])
    expo.to_csv(PROCESSED / "phase2_exposure_lga.csv", index=False)
    tot = expo[expo["lga"] == "All LGAs"].iloc[0]
    summ["exposure"] = {
        "quintile_cut_lst_rel_change": c_lst, "quintile_cut_bshare_change": c_b,
        "hotspot_cells": int(df["hotspot"].sum()),
        "worldpop2023_hotspot": tot["worldpop2023_hotspot"],
        "worldpop2023_main": tot["worldpop2023_main"],
        "worldpop2023_hotspot_share": tot["worldpop2023_hotspot_share"],
        "ghspop2020_hotspot": tot["ghspop2020_hotspot"], "ghspop2020_main": tot["ghspop2020_main"],
        "ghspop2020_hotspot_share": tot["ghspop2020_hotspot_share"], "ghspop_file": ghs_name,
        "top_lgas_worldpop": expo[expo["lga"] != "All LGAs"].head(5)[
            ["lga", "hotspot_cells", "worldpop2023_hotspot", "ghspop2020_hotspot"]
        ].to_dict("records")}

    # 7. case studies
    df["cluster"] = label_clusters(df["row"], df["col"], df["hotspot"].to_numpy(),
                                   spec.nrow, spec.ncol)
    hc = df[df["cluster"] > 0].copy()
    hc["x0"] = spec.xmin + hc["col"] * spec.cell
    hc["y0"] = spec.ymax - (hc["row"] + 1) * spec.cell
    agg = hc.groupby("cluster").agg(
        n_cells=("cell_id", "size"), sum_bshare_change=("bshare_change", "sum"),
        lga=("lga", lambda s: s.value_counts().index[0]),
        lgas=("lga", lambda s: "; ".join(sorted(s.unique()))),
        mean_bshare_change=("bshare_change", "mean"),
        mean_lst_rel_change=("lst_rel_change", "mean"),
        mean_ndvi_rel_change_p80=("ndvi_rel_change_p80", "mean"),
        worldpop2023=("pop_2023", "sum"),
        xmin=("x0", "min"), ymin=("y0", "min"),
        xmax=("x0", lambda s: s.max() + spec.cell), ymax=("y0", lambda s: s.max() + spec.cell),
        lon=("lon", "mean"), lat=("lat", "mean")).reset_index()
    cases = choose_cases(agg)
    cases.to_csv(PROCESSED / "phase2_case_studies.csv", index=False)
    summ["case_studies"] = {
        "clusters": len(cases), "clusters_with_2plus_cells": int((cases["n_cells"] >= 2).sum()),
        "chosen": cases[cases["chosen"]].drop(columns=["xmin", "ymin", "xmax", "ymax"]
                                              ).to_dict("records")}

    df[["cell_id", "row", "col", "lga", "block", "base", "main", "main_nosand", "h4_group",
        "hotspot", "cluster", "ghs_pop_2020"]].to_csv(PROCESSED / "phase2_sample.csv",
                                                      index=False)

    # figures
    n_main = int(main_.sum())
    decile_chart(dec, "lst_rel_change", "Mean change in relative LST (°C, late minus early)",
                 "Relative LST change by decile of building-share change",
                 OUT_FIG / "deciles_lst.png", n_main)
    decile_chart(dec, "ndvi_rel_change_p80", "Mean change in relative NDVI (80th percentile)",
                 "Relative NDVI change by decile of building-share change",
                 OUT_FIG / "deciles_ndvi.png", n_main)
    coef_chart(models[models["hypothesis"].isin(["H1", "H2"])], OUT_FIG / "coefficients.png")
    hotspot_map(df, spec, cases, OUT_FIG / "hotspots.png")
    summ["figures"] = sorted(p.name for p in OUT_FIG.glob("*.png"))
    write_summary("phase2", summ)
    print("phase2: analysis done")
    if not skip_maps:
        case_maps()


def case_maps():
    """Case-study before/after maps from the written case and sample tables (window reads
    only; Open Buildings tiles are listed once per year)."""
    import json

    from shapely.geometry import box
    from shapely.ops import unary_union

    from lagos_heat.buildings import lagos_tiles
    from lagos_heat.common import SUMMARIES, retry

    spec = GridSpec.load()
    cases = pd.read_csv(PROCESSED / "phase2_case_studies.csv")
    cases = cases[cases["chosen"]]
    sample = pd.read_csv(PROCESSED / "phase2_sample.csv")
    boxes = unary_union([box(c.xmin - CASE_PAD_M, c.ymin - CASE_PAD_M, c.xmax + CASE_PAD_M,
                             c.ymax + CASE_PAD_M) for c in cases.itertuples()])
    tiles = {y: retry(lagos_tiles, y, boxes) for y in (2016, 2023)}
    for _, c in cases.iterrows():
        cells = sample[sample["cluster"] == c["cluster"]].copy()
        cells["x0"] = spec.xmin + cells["col"] * spec.cell
        cells["y0"] = spec.ymax - (cells["row"] + 1) * spec.cell
        print(f"  case {int(c['rank'])}: {len(cells)} cells", flush=True)
        case_map(c, cells, tiles, OUT_FIG / f"case{int(c['rank'])}_before_after.png")
    path = SUMMARIES / "phase2.json"
    summ = json.loads(path.read_text(encoding="utf-8"))
    summ["figures"] = sorted(p.name for p in OUT_FIG.glob("*.png"))
    write_summary("phase2", summ)
    print("phase2: case maps done")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--skip-maps", action="store_true", help="skip the case-study map reads")
    p.add_argument("--maps-only", action="store_true",
                   help="only draw the case-study maps from the written tables")
    a = p.parse_args()
    if a.maps_only:
        case_maps()
    else:
        main(a.skip_maps)
