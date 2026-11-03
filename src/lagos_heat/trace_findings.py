"""Phase 3, step 2: trace every number in docs/findings_draft.md to the file that holds it.

Each line of the draft that contains numbers has an entry in SPEC: a unique piece of the
line's text and, in order, one source per number on that line. The numbers are read from
the line itself (inline code spans removed), so a number that has no source, or a source
with no number, stops the run. Sources:

- ("csv", file, {column: value, ...}, column, transform): one cell of a phase2_*.csv file
- ("json", "a/b/c", transform): a value in data/interim/summaries/phase2.json
- ("const", "module.NAME", transform): a constant in the project code (or grid_spec.json)
- ("plan", "text"): a value stated in docs/phase2_analysis_plan.md (the text must occur there)
- ("derived", description, function): computed here from the phase2 output files
- ("label", reason): not a quantity (a phase number, a sensor name such as Landsat 8, a
  data-set year in a name)

Transforms: "x" as is, "pp" and "pct" times 100, "neg" negated, "km" metres to km.
A number matches when the source value, transformed and rounded to the number of decimals
shown, equals it. Word quantities (for example "four of the six") are checked in CLAIMS.

Output: data/processed/phase3_findings_trace.csv; summary key "findings_trace" in
data/interim/summaries/phase3.json.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DRAFT = ROOT / "docs" / "findings_draft.md"
PLAN = ROOT / "docs" / "phase2_analysis_plan.md"
PROCESSED = ROOT / "data" / "processed"
P2JSON = ROOT / "data" / "interim" / "summaries" / "phase2.json"
OUT = PROCESSED / "phase3_findings_trace.csv"

NUM = re.compile(r"(?<![\w.])[−-]?\d{1,3}(?:,\d{3})+(?:\.\d+)?(?![\w])"
                 r"|(?<![\w.])[−-]?\d+(?:\.\d+)?(?![\w])")
CODE = re.compile(r"`[^`]*`")
URL = re.compile(r"https?://\S+")
ISO_DATE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")  # access dates: not quantities
ALL_OUT = PROCESSED / "findings_trace.csv"
TRANSFORMS = {"x": lambda v: v, "pp": lambda v: 100 * v, "pct": lambda v: 100 * v,
              "neg": lambda v: -v, "km": lambda v: v / 1000}


# ------------------------------------------------------------------ parsing and matching
def numbers(line):
    """Numeric tokens of a markdown line, inline code removed: [(text, value, decimals)]."""
    out = []
    line = ISO_DATE.sub("", URL.sub("", CODE.sub("", line)))
    for t in NUM.findall(line):
        s = t.replace("−", "-").replace(",", "")
        dec = len(s.split(".")[1]) if "." in s else 0
        out.append((t, float(s), dec))
    return out


def matches(value, decimals, source):
    if source is None or not np.isfinite(source):
        return False
    return abs(value - source) <= 0.5 * 10 ** -decimals + 1e-9 * max(1.0, abs(source))


# ------------------------------------------------------------------ sources
@lru_cache(maxsize=None)
def table(name):
    return pd.read_csv(PROCESSED / name)


@lru_cache(maxsize=None)
def p2json():
    return json.loads(P2JSON.read_text())


def csv_value(name, where, col):
    t = table(name)
    m = np.ones(len(t), bool)
    for k, v in where.items():
        m &= (t[k].astype(str) == str(v)).to_numpy()
    if m.sum() != 1:
        raise KeyError(f"{name} {where}: {m.sum()} rows")
    return float(t.loc[m, col].iloc[0])


def json_value(path):
    o = p2json()
    for k in path.split("/"):
        o = o[int(k)] if isinstance(o, list) else o[k]
    return float(o)


def const_value(name):
    import importlib

    if name == "grid_spec.cell":
        return float(json.loads((ROOT / "data/interim/grid_spec.json").read_text())["cell"])
    mod, attr = name.rsplit(".", 1)
    idx = None
    if attr.endswith("]"):
        attr, idx = attr[:-1].split("[")
    v = getattr(importlib.import_module(f"lagos_heat.{mod}"), attr)
    return float(v[int(idx)] if idx is not None else v)


def resolve(src):
    """(value or None, file, locator) for one source entry."""
    kind = src[0]
    if kind == "csv":
        _, name, where, col, tr = src
        loc = ", ".join(f"{k}={v}" for k, v in where.items()) + f"; column {col}"
        return TRANSFORMS[tr](csv_value(name, where, col)), f"data/processed/{name}", loc + _tr(tr)
    if kind == "json":
        _, path, tr = src
        return TRANSFORMS[tr](json_value(path)), "data/interim/summaries/phase2.json", path + _tr(tr)
    if kind == "const":
        _, name, tr = src
        f = "data/interim/grid_spec.json" if name.startswith("grid_spec") else \
            f"src/lagos_heat/{name.split('.')[0]}.py"
        return TRANSFORMS[tr](const_value(name)), f, name + _tr(tr)
    if kind == "plan":
        text = src[1]
        found = text in PLAN.read_text(encoding="utf-8")
        val = float(NUM.findall(text)[0].replace(",", "")) if found else None
        return val, "docs/phase2_analysis_plan.md", f'text "{text}"'
    if kind == "derived":
        _, desc, fn, files = src
        return float(fn()), files, desc
    if kind == "jsonf":
        _, name, path, tr = src
        o = json.loads((ROOT / "data/interim/summaries" / name).read_text())
        for k in path.split("/"):
            o = o[int(k)] if isinstance(o, list) else o[k]
        return TRANSFORMS[tr](float(o)), f"data/interim/summaries/{name}", path + _tr(tr)
    if kind == "cited":
        _, url, phrase = src
        return None, url, f'external source; text: "{phrase}"'
    if kind == "label":
        return None, "", src[1]
    raise ValueError(kind)


def _tr(tr):
    return {"x": "", "pp": " (x100, percentage points)", "pct": " (x100, percent)",
            "neg": " (sign reversed in the text)", "km": " (metres to km)"}[tr]


# ------------------------------------------------------------------ helpers for the spec
M, H3, H4, DEC, MOR, LGA, CASE = ("phase2_models.csv", "phase2_h3.csv", "phase2_h4_groups.csv",
                                  "phase2_deciles.csv", "phase2_moran.csv",
                                  "phase2_exposure_lga.csv", "phase2_case_studies.csv")


def run(r, col="per_10pp", tr="x"):
    return ("csv", M, {"run": r}, col, tr)


def ci(r, tr="x"):
    return [run(r, "per_10pp", tr), run(r, "per_10pp_ci_lo", tr), run(r, "per_10pp_ci_hi", tr)]


def ci_only(r):
    return [run(r, "per_10pp_ci_lo"), run(r, "per_10pp_ci_hi")]


def nb(r):
    return [run(r, "n"), run(r, "n_blocks")]


def h3(v, col):
    return ("csv", H3, {"version": v}, col, "x")


def h3_triplet(v, stem):
    return [h3(v, f"{stem}_per_10pp"), h3(v, f"{stem}_ci_lo_per_10pp"),
            h3(v, f"{stem}_ci_hi_per_10pp")]


def h4(defn, group, col, tr="x"):
    return ("csv", H4, {"definition": defn, "group": group}, col, tr)


def h4_diff(defn):
    return [h4(defn, "built", c) for c in
            ("difference_per_10pp", "difference_per_10pp_ci_lo", "difference_per_10pp_ci_hi")]


def h4_row(group):
    return [h4("main", group, "n"), h4("main", group, "bshare_change_p50"),
            h4("main", group, "bshare_change_p90"), h4("main", group, "per_10pp"),
            h4("main", group, "per_10pp_ci_lo"), h4("main", group, "per_10pp_ci_hi")]


def dec(outcome, d, col):
    return ("csv", DEC, {"outcome": outcome, "decile": d}, col, "x")


def moran(r):
    return [("csv", MOR, {"run": r}, c, "x") for c in ("moran_i", "p_perm", "n", "n_islands")]


def lga(name):
    return [("csv", LGA, {"lga": name}, c, "x") for c in
            ("hotspot_cells", "worldpop2023_hotspot", "ghspop2020_hotspot")]


def case(rank):
    return [("csv", CASE, {"rank": rank}, c, "x") for c in
            ("rank", "n_cells", "mean_bshare_change", "mean_lst_rel_change",
             "mean_ndvi_rel_change_p80", "worldpop2023", "lon", "lat")]


STEP = ("const", "phase2.STEP", "pp")
P95 = ("plan", "95% confidence interval")
L8 = ("label", "sensor name (Landsat 8)")


def _decile_range(fn):
    def f():
        t = table(DEC)
        return fn(t.loc[(t["outcome"] == "lst_rel_change") & t["decile"].between(2, 9), "mean"])
    return f


def _moran(fn):
    return lambda: round(fn(table(MOR)["moran_i"]), 2)


def _period_year(period, i, j):
    def f():
        from lagos_heat.common import PERIODS
        return int(PERIODS[period][i].split("/")[j][:4]) % (10000 if j == 0 else 100)
    return f


SPEC = [
    ("# Findings draft (Phase 2)", [("label", "phase number")]),
    ("not causes. Intervals are 95% intervals", [
        ("const", "grid_spec.cell", "x"), P95, ("const", "phase2.BLOCK", "x"),
        ("const", "phase2.BLOCK", "x")]),
    ("The base sample has 9,730", [
        ("json", "samples/base_analysis_ok", "x"), ("json", "samples/removed_water_frac_gt_0_10", "x"),
        ("const", "phase2.WATER_MAX", "pct"), ("json", "samples/main", "x"),
        run("H1_primary", "n_blocks"), ("json", "samples/sandfill_candidates_in_main", "x"),
        ("json", "samples/main_excl_sandfill", "x")]),
    ("| Main sample | 9,157", nb("H1_primary")),
    ("| Without sand-filling candidates | 8,060", nb("H1_excl_sandfill")),
    ("| Without cells within 1 km of ocean or lagoon | 7,462",
     [("const", "phase2.COAST_M", "km"), *nb("H1_excl_coast_lagoon_1km")]),
    ("| Without cells with industrial land | 8,561", nb("H1_excl_industrial")),
    ("| Without cells within 250 m of a major road | 8,187",
     [("const", "phase2.ROAD_M", "x"), *nb("H1_excl_road_250m")]),
    ("| Landsat 8-only LST, cells with at least 3", [
        L8, ("const", "common.MIN_OBS", "x"), L8, *nb("H1_landsat8_only")]),
    ("| At least 5 LST observations per period | 2,465",
     [("const", "phase2.LST_OBS_STRICT", "x"), *nb("H1_obs_min5")]),
    ("| Terciles of early relative LST |", [
        run("H1_lst_early_tercile1", "n"), run("H1_lst_early_tercile2", "n"),
        run("H1_lst_early_tercile3", "n"), run("H1_lst_early_tercile1", "n_blocks"),
        run("H1_lst_early_tercile2", "n_blocks"), run("H1_lst_early_tercile3", "n_blocks")]),
    ("| Screened-median NDVI check", nb("H2_screened_descriptive")),
    ("In the main sample, building share changed by 0", [
        ("json", "bshare_change_main/p25", "pp"), ("json", "bshare_change_main/median", "pp"),
        ("json", "bshare_change_main/p75", "pp"), ("json", "bshare_change_main/iqr", "pp"), STEP]),
    ("Mean relative LST change rises from", [
        dec("lst_rel_change", 1, "mean"), dec("lst_rel_change", 1, "ci_lo"),
        dec("lst_rel_change", 1, "ci_hi"), dec("lst_rel_change", 1, "x_mean_pp"),
        dec("lst_rel_change", 10, "mean"), dec("lst_rel_change", 10, "ci_lo"),
        dec("lst_rel_change", 10, "ci_hi"), dec("lst_rel_change", 10, "x_mean_pp"),
        ("derived", "lowest mean of deciles 2-9 (lst_rel_change)", _decile_range(np.min),
         "data/processed/phase2_deciles.csv"),
        ("derived", "highest mean of deciles 2-9 (lst_rel_change)", _decile_range(np.max),
         "data/processed/phase2_deciles.csv")]),
    ("In M1, relative LST change is 0.579", [
        run("H1_primary"), STEP, *ci_only("H1_primary"), run("H1_primary", "n")]),
    ("The coefficient is positive with an interval excluding zero in all 13",
     [("json", "H1/n_runs", "x")]),
    ("| Run | °C per 10 pp | 95% interval |", [STEP, P95]),
    ("| Primary (M1) | 0.579", ci("H1_primary")),
    ("| Without sand-filling candidates | 0.504", ci("H1_excl_sandfill")),
    ("| Without cells within 1 km of ocean or lagoon | 0.549",
     [("const", "phase2.COAST_M", "km"), *ci("H1_excl_coast_lagoon_1km")]),
    ("| Without cells with industrial land | 0.620", ci("H1_excl_industrial")),
    ("| Without cells within 250 m of a major road | 0.595",
     [("const", "phase2.ROAD_M", "x"), *ci("H1_excl_road_250m")]),
    ("| Landsat 8-only LST | 0.311", [L8, *ci("H1_landsat8_only")]),
    ("| LST observation counts as controls |", ci("H1_obs_controls")),
    ("| At least 5 LST observations per period | 0.436",
     [("const", "phase2.LST_OBS_STRICT", "x"), *ci("H1_obs_min5")]),
    ("| Without early relative LST as control |", ci("H1_no_baseline_lst")),
    ("| Lowest tercile of early relative LST |", ci("H1_lst_early_tercile1")),
    ("| Middle tercile of early relative LST |", ci("H1_lst_early_tercile2")),
    ("| Highest tercile of early relative LST |", ci("H1_lst_early_tercile3")),
    ("| GHSL built-up change 2015–2020 instead", [
        ("const", "buildings.GHSL_EPOCHS[0]", "x"), ("const", "buildings.GHSL_EPOCHS[1]", "x"),
        *ci("H1_ghsl_2015_2020")]),
    ("Mean relative 80th-percentile NDVI change falls from", [
        dec("ndvi_rel_change_p80", 1, "mean"), dec("ndvi_rel_change_p80", 10, "mean")]),
    ("In M2, relative NDVI change is 0.039 lower", [
        run("H2_primary", tr="neg"), STEP, *ci_only("H2_primary")]),
    ("| Run | Relative NDVI per 10 pp | 95% interval |", [STEP, P95]),
    ("| Primary (M2) |", ci("H2_primary")),
    ("| Without sand-filling candidates | −0.034", ci("H2_excl_sandfill")),
    ("| Unscreened relative median |", ci("H2_noscreen")),
    ("| With bare-share change as control |", ci("H2_bare_control")),
    ("The descriptive check on the 746", [
        run("H2_screened_descriptive", "n"), run("H2_screened_descriptive"), STEP,
        *ci_only("H2_screened_descriptive"), run("H2_screened_descriptive", "pearson_r")]),
    ("Fitting M1 and M3 on the same 9,157", [h3("p80", "n"), h3("p80", "n_boot")]),
    ("| NDVI version | Building coefficient, M1", [STEP]),
    ("| Relative 80th percentile | 0.579", [
        *h3_triplet("p80", "coef_m1"), *h3_triplet("p80", "coef_m3"),
        *h3_triplet("p80", "difference"), h3("p80", "share"), h3("p80", "share_ci_lo"),
        h3("p80", "share_ci_hi")]),
    ("| Relative unscreened median | 0.580", [
        *h3_triplet("noscreen", "coef_m1"), *h3_triplet("noscreen", "coef_m3"),
        *h3_triplet("noscreen", "difference"), h3("noscreen", "share"),
        h3("noscreen", "share_ci_lo"), h3("noscreen", "share_ci_hi")]),
    ("The NDVI-change coefficient in M3 is", [
        h3("p80", "ndvi_change_coef_m3"), h3("p80", "ndvi_change_ci_lo_m3"),
        h3("p80", "ndvi_change_ci_hi_m3"), h3("noscreen", "ndvi_change_coef_m3"),
        h3("noscreen", "ndvi_change_ci_lo_m3"), h3("noscreen", "ndvi_change_ci_hi_m3")]),
    ("## H4: built-up versus mostly vegetated cells in 2016",
     [("const", "buildings.OB_YEARS[0]", "x")]),
    ("| Group | Cells | Median building change |", [STEP, P95]),
    ("| Built-up 2016 | 878", [("const", "buildings.OB_YEARS[0]", "x"), *h4_row("built")]),
    ("| Mostly vegetated | 4,018", h4_row("vegetated")),
    ("| Neither | 4,261", h4_row("neither")),
    ("Building change is a share of cell area, so 0.065", [
        h4("main", "built", "bshare_change_p50"), h4("main", "built", "bshare_change_p50", "pp")]),
    ("The difference between the built-up and vegetated slopes is", [
        h4("main", "built", "difference_per_10pp"), STEP,
        h4("main", "built", "difference_per_10pp_ci_lo"),
        h4("main", "built", "difference_per_10pp_ci_hi"), *ci("H4_main_built_separate"),
        *ci("H4_main_vegetated_separate")]),
    ("| Definition | Difference, °C per 10 pp |", [STEP, P95]),
    ("| Main (built ≥ 0.30, tree ≥ 0.50) |", [
        ("const", "phase2.BUILT_T", "x"), ("const", "phase2.TREE_T", "x"), *h4_diff("main")]),
    ("| Built threshold 0.20 |", [("const", "phase2.H4_BUILT_ALT[0]", "x"), *h4_diff("built_0.20")]),
    ("| Built threshold 0.40 |", [("const", "phase2.H4_BUILT_ALT[1]", "x"), *h4_diff("built_0.40")]),
    ("| Tree threshold 0.40 |", [("const", "phase2.H4_TREE_ALT[0]", "x"), *h4_diff("tree_0.40")]),
    ("| Tree threshold 0.60 |", [("const", "phase2.H4_TREE_ALT[1]", "x"), *h4_diff("tree_0.60")]),
    ("| Vegetated by early vegetated share ≥ 0.50 |", [
        ("const", "phase2.VEG_T", "x"), *h4_diff("vegetated_by_veg_share_early")]),
    ("In the vegetated group, half or more of the cells", [
        ("label", "percentile named by column bshare_change_p90"),
        h4("main", "vegetated", "bshare_change_p90", "pp")]),
    ("- **Sand filling:** without the sand-filling candidates", [
        run("H1_excl_sandfill"), run("H2_excl_sandfill"), STEP]),
    ("- **Coast and lagoon:**", [("const", "phase2.COAST_M", "km"), *ci("H1_excl_coast_lagoon_1km")]),
    ("- **Industry and roads:**", [
        *ci("H1_excl_industrial"), ("const", "phase2.ROAD_M", "x"), *ci("H1_excl_road_250m")]),
    ("- **Sensor mix:** with Landsat 8-only LST", [
        L8, run("H1_landsat8_only", "n"), *ci("H1_landsat8_only")]),
    ("- **Observation counts:**", [
        run("H1_obs_controls"), ("const", "phase2.LST_OBS_STRICT", "x"),
        run("H1_obs_min5", "n"), run("H1_obs_min5")]),
    ("- **Regression to the mean:** without", [
        run("H1_no_baseline_lst"), run("H1_lst_early_tercile1"), run("H1_lst_early_tercile2"),
        run("H1_lst_early_tercile3")]),
    ("- **Building measure:**", [
        ("const", "buildings.GHSL_EPOCHS[0]", "x"), ("const", "buildings.GHSL_EPOCHS[1]", "x"),
        *ci("H1_ghsl_2015_2020")]),
    ("Moran's I of the residuals, with rook contiguity and 999",
     [("csv", MOR, {"run": "H1_primary"}, "n_perm", "x")]),
    ("| M1 | 0.639", moran("H1_primary")),
    ("| M2 | 0.532", moran("H2_primary")),
    ("| M3 (80th percentile) | 0.645", moran("H3_M3_p80")),
    ("beyond the 2.5 km blocks", [
        ("derived", "phase2.BLOCK x grid cell size, in km",
         lambda: const_value("phase2.BLOCK") * const_value("grid_spec.cell") / 1000,
         "src/lagos_heat/phase2.py; data/interim/grid_spec.json")]),
    ("Hotspot cells are cells in the top quintile", [
        ("json", "exposure/quintile_cut_lst_rel_change", "x"),
        ("json", "exposure/quintile_cut_bshare_change", "x"), ("json", "exposure/hotspot_cells", "x")]),
    ("- **WorldPop 2023:** 896,165", [
        ("label", "data-set name (WorldPop 2023)"), ("json", "exposure/worldpop2023_hotspot", "x"),
        ("json", "exposure/worldpop2023_hotspot_share", "pct"),
        ("json", "exposure/worldpop2023_main", "x")]),
    ("- **GHS-POP 2020 cross-check:**", [
        ("label", "data-set name (GHS-POP 2020)"), ("json", "exposure/ghspop2020_hotspot", "x"),
        ("json", "exposure/ghspop2020_hotspot_share", "pct"),
        ("json", "exposure/ghspop2020_main", "x")]),
    ("LGAs with the largest hotspot populations by WorldPop 2023",
     [("label", "data-set name (WorldPop 2023)")]),
    ("| LGA | Hotspot cells | WorldPop 2023 | GHS-POP 2020 |", [
        ("label", "data-set name (WorldPop 2023)"), ("label", "data-set name (GHS-POP 2020)")]),
    ("| Alimosho | 74", lga("Alimosho")),
    ("| Kosofe |", lga("Kosofe")),
    ("| Oshodi/Isolo |", lga("Oshodi/Isolo")),
    ("| Ikorodu | 56", lga("Ikorodu")),
    ("| Ifako/Ijaye |", lga("Ifako/Ijaye")),
    ("The hotspot cells form 155", [
        ("json", "case_studies/clusters", "x"), ("csv", CASE, {"rank": 3}, "n_cells", "x"),
        ("const", "phase2.CASE_THIRD_MIN", "x")]),
    ("| Case | Area (LGA) |", [("label", "data-set name (WorldPop 2023)")]),
    ("| 1 | Lekki Peninsula II", case(1)),
    ("| 2 | Abule Egba / Iyana Ipaja", case(2)),
    ("| 3 | Oke, Ikorodu |", case(3)),
    ("The before/after maps of buildings, NDVI and LST", [
        ("plan", "7. **Case-study maps (7)")]),
    ("- **Spatial dependence:** the residuals", [
        ("derived", "lowest Moran's I, 2 decimals", _moran(np.min), "data/processed/phase2_moran.csv"),
        ("derived", "highest Moran's I, 2 decimals", _moran(np.max), "data/processed/phase2_moran.csv")]),
    ("- **Sensor mix:** the Landsat 8-only estimate is about half", [L8]),
    ("- **Building change in dense cells:**", [("const", "buildings.OB_YEARS[0]", "x")]),
    ("- **Dates and sensors:**", [
        ("const", "buildings.OB_YEARS[0]", "x"), ("const", "buildings.OB_YEARS[1]", "x"),
        ("derived", "first year of the early period (common.PERIODS)", _period_year("early", 0, 0),
         "src/lagos_heat/common.py"),
        ("derived", "last year of the early period, two digits", _period_year("early", -1, 1),
         "src/lagos_heat/common.py"),
        ("derived", "first year of the late period", _period_year("late", 0, 0),
         "src/lagos_heat/common.py"),
        ("derived", "last year of the late period, two digits", _period_year("late", -1, 1),
         "src/lagos_heat/common.py")]),
]


def _h4_defs():
    return table(H4).drop_duplicates("definition")


CLAIMS = [  # (unique text on the line, phrase, description, function returning True if it holds)
    ("All four H2 runs are negative", "All four H2 runs",
     "H2/n_runs == 4 and all runs negative (phase2.json)",
     lambda: p2json()["H2"]["n_runs"] == 4 and p2json()["H2"]["robust_sign_all_runs"]),
    ("Across the six group definitions", "six group definitions",
     "6 definitions in phase2_h4_groups.csv", lambda: len(_h4_defs()) == 6),
    ("Across the six group definitions", "excludes zero in four of the six",
     "4 definitions with sign negative or positive (phase2_h4_groups.csv)",
     lambda: int((_h4_defs()["sign"] != "not distinguishable from zero").sum()) == 4),
    ("No main-sample cell meets both group rules", "No main-sample cell meets both",
     "n_cells_meeting_both_rules == 0 for definition main (phase2_h4_groups.csv)",
     lambda: csv_value(H4, {"definition": "main", "group": "built"},
                       "n_cells_meeting_both_rules") == 0),
    ("In the vegetated group, half or more", "half or more of the cells have no building change",
     "median building change in the vegetated group == 0 (phase2_h4_groups.csv)",
     lambda: csv_value(H4, {"definition": "main", "group": "vegetated"}, "bshare_change_p50") == 0),
    ("the building coefficient falls by about half", "falls by about half in both versions",
     "share of M1 coefficient between 0.4 and 0.6 for both versions (phase2_h3.csv)",
     lambda: table(H3)["share"].between(0.4, 0.6).all()),
    ("The coefficient is positive with an interval excluding zero in all 13",
     "interval excluding zero in all 13 H1 runs",
     "H1/runs_ci_excludes_zero == H1/n_runs (phase2.json)",
     lambda: p2json()["H1"]["runs_ci_excludes_zero"] == p2json()["H1"]["n_runs"]),
    ("the three largest are chosen", "the three largest are chosen",
     "3 rows with chosen == True, ranks 1-3 (phase2_case_studies.csv)",
     lambda: table(CASE)["chosen"].sum() == 3),
    ("- **Sensor mix:** the Landsat 8-only estimate is about half", "about half the primary",
     "ratio H1_landsat8_only / H1_primary between 0.4 and 0.6 (phase2_models.csv)",
     lambda: 0.4 <= csv_value(M, {"run": "H1_landsat8_only"}, "per_10pp")
     / csv_value(M, {"run": "H1_primary"}, "per_10pp") <= 0.6),
]


def find_line(lines, key):
    hits = [i for i, l in enumerate(lines) if key in l]
    if len(hits) != 1:
        raise KeyError(f"spec key {key!r} found on {len(hits)} lines")
    return hits[0]


def trace(lines=None, spec=None, claims=None):
    lines = lines if lines is not None else DRAFT.read_text(encoding="utf-8").splitlines()
    spec = SPEC if spec is None else spec
    claims = CLAIMS if claims is None else claims
    fence, in_code = False, set()  # lines inside fenced code blocks (commands) are skipped
    for i, l in enumerate(lines):
        if l.strip().startswith("```"):
            fence = not fence
            in_code.add(i)
        elif fence:
            in_code.add(i)
    rows, used = [], set(in_code)
    for key, sources in spec:
        i = find_line(lines, key)
        used.add(i)
        toks = numbers(lines[i])
        if len(toks) != len(sources):
            raise ValueError(f"line {i + 1}: {len(toks)} numbers, {len(sources)} sources")
        for k, ((text, val, d), src) in enumerate(zip(toks, sources), 1):
            sval, f, loc = resolve(src)
            kind = src[0]
            ok = {"label": "n/a", "cited": "cited"}.get(kind) or                 ("yes" if matches(val, d, sval) else "no")
            rows.append({"line": i + 1, "position": k, "text": text, "value": val,
                         "kind": kind, "source_file": f, "locator": loc,
                         "source_value": sval, "match": ok})
    untraced = [(i + 1, t[0]) for i, l in enumerate(lines) if i not in used for t in numbers(l)]
    for key, phrase, desc, fn in claims:
        i = find_line(lines, key)
        if phrase not in lines[i]:
            raise KeyError(f"claim {phrase!r} not on line {i + 1}")
        rows.append({"line": i + 1, "position": 0, "text": phrase, "value": np.nan,
                     "kind": "word claim", "source_file": "", "locator": desc,
                     "source_value": np.nan, "match": "yes" if fn() else "no"})
    out = pd.DataFrame(rows).sort_values(["line", "position"], kind="stable")
    return out, untraced


def summarise(out, untraced, doc):
    num = out[out["kind"] != "word claim"]
    return {
        "document": doc, "numbers": int(len(num)),
        "by_kind": num["kind"].value_counts().to_dict(),
        "numbers_checked": int(num["match"].isin(["yes", "no"]).sum()),
        "numbers_cited_external": int((num["match"] == "cited").sum()),
        "numbers_matched": int((num["match"] == "yes").sum()),
        "labels_not_checked": int((num["match"] == "n/a").sum()),
        "word_claims": int((out["kind"] == "word claim").sum()),
        "word_claims_hold": int(((out["kind"] == "word claim") & (out["match"] == "yes")).sum()),
        "untraced_numbers": untraced,
        "mismatches": out.loc[out["match"] == "no", ["line", "text", "locator",
                                                       "source_value"]].to_dict(orient="records"),
    }


def main():
    from lagos_heat.trace_docs import DOCS
    from lagos_heat.validate import write_summary

    out, untraced = trace()
    out.to_csv(OUT, index=False)  # Phase 3 trace of the draft, kept as before
    summ = summarise(out, untraced, "docs/findings_draft.md")
    summ["draft"] = summ.pop("document")
    write_summary("findings_trace", summ)
    frames, per_doc = [out.assign(document="docs/findings_draft.md")], [summ]
    for doc, (spec, claims) in DOCS.items():
        lines = (ROOT / doc).read_text(encoding="utf-8").splitlines()
        o, u = trace(lines, spec, claims)
        frames.append(o.assign(document=doc))
        per_doc.append(summarise(o, u, doc))
    allout = pd.concat(frames, ignore_index=True)
    allout.insert(0, "document", allout.pop("document"))
    allout.to_csv(ALL_OUT, index=False)
    p4 = ROOT / "data/interim/summaries/phase4.json"
    d = json.loads(p4.read_text())
    d["findings_trace"] = per_doc[1:]
    p4.write_text(json.dumps(d, indent=2, default=float))
    for s_ in per_doc:
        print(json.dumps({k: s_[k] for k in ("numbers", "numbers_checked", "numbers_matched",
                                             "numbers_cited_external", "labels_not_checked",
                                             "word_claims", "word_claims_hold", "untraced_numbers",
                                             "mismatches") if k in s_} | {"doc": s_.get("document", s_.get("draft"))},
                         default=str))


if __name__ == "__main__":
    main()
