"""Phase 4: trace specifications for the final documents (findings, brief, README, policy
context). Same format as trace_findings.SPEC; see that module for the source kinds.

Numbers inside fenced code blocks (the run commands in README.md), URLs and ISO access dates
are not quantities and are skipped by the parser. Facts from external sources are recorded
with kind "cited" (URL and the source's wording); they are not machine-checked.
"""

from __future__ import annotations

import numpy as np

from lagos_heat import trace_findings as tf
from lagos_heat.trace_findings import (CASE, H3, H4, M, MOR, P95, STEP, L8, _period_year, ci,
                                       ci_only, csv_value, h3, h4, moran, nb, run, table)

LGA4, CELLS4 = "priority_areas_lga.csv", "priority_cells.csv"

# common sources
CELL = ("const", "grid_spec.cell", "x")
OB16, OB23 = ("const", "buildings.OB_YEARS[0]", "x"), ("const", "buildings.OB_YEARS[1]", "x")
POP16, POP23 = ("const", "population.YEARS[0]", "x"), ("const", "population.YEARS[1]", "x")
E0 = ("derived", "first year of the early period", _period_year("early", 0, 0),
      "src/lagos_heat/common.py")
E1 = ("derived", "last year of the early period, two digits", _period_year("early", -1, 1),
      "src/lagos_heat/common.py")
L0 = ("derived", "first year of the late period", _period_year("late", 0, 0),
      "src/lagos_heat/common.py")
L1 = ("derived", "last year of the late period, two digits", _period_year("late", -1, 1),
      "src/lagos_heat/common.py")


def _end_year():
    from lagos_heat.common import PERIODS
    return int(PERIODS["late"][-1].split("/")[1][:4])


END = ("derived", "last year of the late period (data end)", _end_year, "src/lagos_heat/common.py")
H1 = run("H1_primary")
NRUNS = ("json", "H1/n_runs", "x")
NMAIN = run("H1_primary", "n")
H2NEG = run("H2_primary", tr="neg")
SHARE = h3("p80", "share")
BUILT, VEG = h4("main", "built", "per_10pp"), h4("main", "vegetated", "per_10pp")
HOT = ("json", "exposure/hotspot_cells", "x")
WP = ("json", "exposure/worldpop2023_hotspot", "x")
WPS = ("json", "exposure/worldpop2023_hotspot_share", "pct")
WPL = ("label", "data-set name (WorldPop 2023)")
MMIN = ("derived", "lowest Moran's I, 2 decimals", tf._moran(np.min), "data/processed/phase2_moran.csv")
MMAX = ("derived", "highest Moran's I, 2 decimals", tf._moran(np.max), "data/processed/phase2_moran.csv")
CUT = ("json", "exposure/quintile_cut_bshare_change", "x")
P3 = ("label", "phase number")
LIC = ("label", "licence or product version")
S2 = ("label", "sensor name (Sentinel-2)")
LIST = ("label", "list item number")


def j4(path, tr="x"):
    return ("jsonf", "phase4.json", path, tr)


def j3(path, tr="x"):
    return ("jsonf", "phase3.json", path, tr)


def human_ok():
    import pandas as pd
    d = pd.read_csv(tf.PROCESSED / "phase3_spotcheck_ai.csv")
    return int((d["human_agrees_with_ai"] == "yes").sum())


HUMAN = ("derived", "rows with human_agrees_with_ai = yes", human_ok,
         "data/processed/phase3_spotcheck_ai.csv")
N40 = j3("spotcheck/n_cells")


def a1(lga):
    w = {"list": "A1_cool_now_lga", "lga": lga}
    return [("csv", LGA4, w, c, "x") for c in
            ("priority_rank", "n_hotspot_cells", "residents_hotspot", "mean_lst_rel_change")]


def bcell(cid):
    w = {"list": "B_protect_cell", "unit_id": cid}
    return [("csv", CELLS4, w, "priority_rank", "x"), ("csv", CELLS4, w, "bshare_change", "pp"),
            ("csv", CELLS4, w, "residents_nearby", "x")]


def cite(url, *phrases):
    return [("cited", url, p) for p in phrases]


NOAA = "https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/"
NAIRA = ("https://nairametrics.com/2026/05/22/"
         "lagos-records-over-27000-building-violations-in-three-years/")
AFRICA = "https://allafrica.com/stories/202607150132.html"
C40 = "https://www.c40.org/wp-content/uploads/2026/07/RfP_Lagos-Heat-Risk-Assessment.pdf"
THISDAY = ("https://www.thisdaylive.com/2026/07/23/"
           "cambridge-unilag-identify-10-locations-for-urban-heat-research-in-lagos/")
LAWNG = "https://lawnigeria.com/2019/10/15/urban-and-regional-planning-law-of-lagos-state/"
EPP = "https://epp.lagosstate.gov.ng/regulations/REVISED_LASPPPA_REGULATION_2019_1.pdf"
LASBCA = "https://lasbca.lagosstate.gov.ng/wp-content/uploads/2023/11/LASBCA-GUIDE.pdf"

LEVELS_SRC = [("cited", "https://www.c40.org/wp-content/uploads/2026/07/RfP_Lagos-Heat-Risk-Assessment.pdf",
               "The CCRA identifies Lagos Island, Lagos Mainland, Agege, Kosofe, and Oshodi-Isolo "
               "as the LGAs most severely affected (Lagos CCRA 2021)")]


def _fringe_in_top4():
    t = table(LGA4)
    top = t[(t["list"] == "A1_cool_now_lga") & (t["priority_rank"] <= 4)]["lga"].tolist()
    return {"Eti Osa", "Epe", "Alimosho"} <= set(top)


FRINGE_CLAIM = ("concentrated on the expanding fringe (for example Eti-Osa, Epe and Alimosho)",
                "Eti Osa, Epe and Alimosho are among the top 4 Cool-now LGAs "
                "(priority_areas_lga.csv)", _fringe_in_top4)

# ------------------------------------------------------------------ docs/findings.md
_unchanged = [e for e in tf.SPEC if not e[0].startswith(("# Findings draft",
                                                          "not causes. Intervals",
                                                          "- **", "The before/after maps"))]
FINDINGS = _unchanged + [
    ("# Findings: building density", [E0, END]),
    ("Across 9,157 cells of 500 m in Lagos State", [
        NMAIN, CELL, OB16, OB23, E0, E1, L0, L1, H1, STEP, P95, *ci_only("H1_primary"), NRUNS,
        H2NEG, STEP, SHARE, OB16, BUILT, STEP, VEG, HOT, WP, WPL, WPS, MMIN, MMAX]),
    ("This document reports the analysis fixed in", [P3, P3, CELL, P95,
                                                     ("const", "phase2.BLOCK", "x"),
                                                     ("const", "phase2.BLOCK", "x")]),
    ("- Without the sand-filling candidates, the H1 coefficient", [
        run("H1_excl_sandfill"), run("H2_excl_sandfill"), STEP]),
    ("- Without cells within 1 km of the coast", [("const", "phase2.COAST_M", "km"),
                                                   *ci("H1_excl_coast_lagoon_1km")]),
    ("- Without industrial cells, H1 is", [*ci("H1_excl_industrial"),
                                           ("const", "phase2.ROAD_M", "x"),
                                           *ci("H1_excl_road_250m")]),
    ("- With Landsat 8-only LST, which removes", [L8, *ci("H1_landsat8_only"),
                                                  run("H1_landsat8_only", "n")]),
    ("Haze in the Sentinel-2 composites", [S2]),
    ("- With LST observation counts as controls", [
        run("H1_obs_controls"), ("const", "phase2.LST_OBS_STRICT", "x"), run("H1_obs_min5"),
        run("H1_obs_min5", "n")]),
    ("- Without early relative LST as a control, H1 is", [
        run("H1_no_baseline_lst"), run("H1_lst_early_tercile1"), run("H1_lst_early_tercile2"),
        run("H1_lst_early_tercile3")]),
    ("- With GHSL built-up change (2015–2020)", [
        ("const", "buildings.GHSL_EPOCHS[0]", "x"), ("const", "buildings.GHSL_EPOCHS[1]", "x"),
        *ci("H1_ghsl_2015_2020")]),
    ("- By WorldPop 2023, 896,165", [WPL, WP, WPS, ("json", "exposure/worldpop2023_main", "x")]),
    ("- By GHS-POP 2020, used as a cross-check", [
        ("label", "data-set name (GHS-POP 2020)"), ("json", "exposure/ghspop2020_hotspot", "x"),
        ("json", "exposure/ghspop2020_hotspot_share", "pct"),
        ("json", "exposure/ghspop2020_main", "x")]),
    ("Area names are OpenStreetMap names near each cluster centre", [
        ("csv", CASE, {"rank": 1}, "rank", "x"), ("csv", CASE, {"rank": 2}, "rank", "x"),
        ("csv", CASE, {"rank": 3}, "rank", "x"), ("plan", "7. **Case-study maps (7)")]),
    ('- "Cool now" ranks the', [j4("hotspot_lgas"), j4("hotspot_clusters"),
                                ("csv", CASE, {"rank": 1}, "rank", "x"), *LEVELS_SRC]),
    ('- "Protect before it is lost" lists', [
        j4("protect_cells"), j4("protect_lgas"), L0, L1, ("const", "phase4.VEG_LATE_MIN", "x"),
        CUT, ("csv", LGA4, {"list": "B_protect_by_lga", "lga": "Ikorodu"}, "n_cells", "x")]),
    ("- Building-change values in cells that were already densely built", [OB16]),
    ("- The residuals of the models are strongly spatially autocorrelated", [MMIN, MMAX]),
    ("- Population comes from WorldPop, a modelled", [POP16, POP23]),
    ("- The satellite data end in February", [END, OB16, OB23, E0, E1, L0, L1, L8]),
    ("- The first late dry season overlaps", cite(
        NOAA, "2023 (RONI above 0.5 from JAS 2023)", "24 (RONI above 0.5 through JFM 2024)",
        "Warm value: RONI greater than 0.5 degrees Celsius (El Niño)",
        "JAS 2023: 0.8, first season above 0.5", "JFM 2024: 0.8, last season above 0.5",
        "OND 2023: 1.4, NDJ 2023: 1.4")),
    ("- An AI coding agent (Claude Code) wrote", [N40, HUMAN, N40]),
]


def _a1_rank1():
    t = table(LGA4)
    return t[(t["list"] == "A1_cool_now_lga") & (t["priority_rank"] == 1)]["lga"].iloc[0]


def _a2_rank1_is_case1():
    t = table(CELLS4)
    cl = t[(t["list"] == "A2_cool_now_cluster") & (t["priority_rank"] == 1)]["unit_id"].iloc[0]
    return int(cl) == int(csv_value(CASE, {"rank": 1}, "cluster"))


def _b_top_lga():
    t = table(LGA4)
    return t[(t["list"] == "B_protect_by_lga") & (t["priority_rank"] == 1)]["lga"].iloc[0]


FINDINGS_CLAIMS = [c for c in tf.CLAIMS if not c[0].startswith("- **")] + [
    ('- "Cool now" ranks the', *FRINGE_CLAIM),
    ("The first-ranked LGA is Eti Osa", "first-ranked LGA is Eti Osa",
     "A1 rank 1 is Eti Osa (priority_areas_lga.csv)", lambda: _a1_rank1() == "Eti Osa"),
    ("The first-ranked LGA is Eti Osa", "first-ranked cluster is case 1",
     "A2 rank 1 is the case 1 cluster (priority_cells.csv, phase2_case_studies.csv)",
     _a2_rank1_is_case1),
    ('- "Protect before it is lost" lists', "Ikorodu has the most such cells",
     "B by LGA rank 1 is Ikorodu (priority_areas_lga.csv)", lambda: _b_top_lga() == "Ikorodu"),
    ("- The satellite data end in February", "about half the primary estimate",
     "ratio H1_landsat8_only / H1_primary between 0.4 and 0.6 (phase2_models.csv)",
     tf.CLAIMS[-1][3]),
]
# the Phase 3 claim on the Landsat ratio pointed at the old "- **Sensor mix:**" line;
# it is carried over above with the new line key.

# ------------------------------------------------------------------ docs/policy_brief.md
BRIEF = [
    ("Policy brief for the Lagos State Parks", [("label", "date of the brief")]),
    ("Lagos has no dedicated heat action plan yet", cite(
        C40, "Urban Heat Action Plan planned by the city in 2027", "July 2026 (RfP cover)")),
    ("The study compares 500 m cells", [CELL, OB16, OB23, E0, E1, L0, L1]),
    ("1. Where buildings were added", [LIST, NMAIN, STEP, H1, P95, *ci_only("H1_primary"), NRUNS]),
    ("2. Building growth went with vegetation loss", [LIST, SHARE]),
    ("3. The warming is strongest where green land", [LIST, STEP, OB16, VEG, BUILT]),
    ("4. About 896,165 people live in heat hotspots", [LIST, WP, HOT, WPL, WPS]),
    ("| 1 | Eti Osa |", a1("Eti Osa")),
    ("| 2 | Epe |", a1("Epe")),
    ("| 3 | Oshodi/Isolo |", a1("Oshodi/Isolo")),
    ("| 4 | Alimosho |", a1("Alimosho")),
    ("| 5 | Shomolu |", a1("Shomolu")),
    ('"Protect before it is lost" lists cells that are still mostly green', [
        L0, L1, ("json", "exposure/quintile_cut_bshare_change", "pp"), j4("protect_cells"),
        j4("protect_lgas"),
        ("csv", LGA4, {"list": "B_protect_by_lga", "lga": "Ikorodu"}, "n_cells", "x"),
        ("derived", "3 x 3 cells x grid cell size, in km",
         lambda: 3 * tf.const_value("grid_spec.cell") / 1000, "data/interim/grid_spec.json")]),
    ("| 1 | Trade Fair Complex (Amuwo Odofin) | 13.3", bcell(19365)),
    ("| 2 | near Ijegun (Alimosho) | 8.7 | 15,647", bcell(16780)),
    ("| 3 | near Ijegun (Alimosho) | 8.7 | 16,150", bcell(16410)),
    ("| 4 | Trade Fair Complex (Amuwo Odofin) | 9.3", bcell(18994)),
    ("| 5 | near Ikotun (Oshodi/Isolo) |", bcell(11234)),
    ("The 2021 Climate Change Risk Assessment, as cited", LEVELS_SRC),
    ("- Enforce setback and open-space rules", cite(
        NAIRA, "27,233 building violations", "between May 2023", "and April 2026",
        "May 22, 2026 (article date)", "May 22, 2026 (article date)")),
    ("- Target tree planting at the", cite(
        AFRICA, "at least 20,000 trees", "in 2026", "13,401 trees", "as of 30 June",
        "14 July 2026 (dateline)", "14 July 2026 (dateline)")),
    ("- Add heat to planning approvals", cite(C40, "heat action plan, being scheduled for 2027")),
    ("The study measures ground-surface temperature from satellites", [END]),
]
BRIEF_CLAIMS = [
    ("The 2021 Climate Change Risk Assessment, as cited", *FRINGE_CLAIM),
    ("2. Building growth went with vegetation loss", "falls by about half",
     "share of M1 coefficient between 0.4 and 0.6 (phase2_h3.csv)",
     lambda: table(H3)["share"].between(0.4, 0.6).all()),
    ("1. Where buildings were added", "held in all 13 model runs",
     "H1/runs_ci_excludes_zero == H1/n_runs (phase2.json)",
     lambda: tf.p2json()["H1"]["runs_ci_excludes_zero"] == tf.p2json()["H1"]["n_runs"]),
]

# ------------------------------------------------------------------ README.md
README = [
    ("# Building out the breeze", [E0, END]),
    ("This project examines how changes in building density", [CELL]),
    ("Across 9,157 cells, building share grew", [
        NMAIN, OB16, OB23, E0, E1, L0, L1, H1, STEP, P95, *ci_only("H1_primary"), NRUNS]),
    ("Building growth went with vegetation loss. Adding", [SHARE, OB16, BUILT, STEP, VEG]),
    ("329 hotspot cells are in the top fifth", [HOT, WPL, WP, WPS, MMIN, MMAX]),
    ("Case 1, Lekki Peninsula II", [("csv", CASE, {"rank": 1}, "rank", "x")]),
    ("![Case 1 before and after]", [("csv", CASE, {"rank": 1}, "rank", "x")]),
    ("Mean relative LST change by decile", [P95]),
    ("- Google Open Buildings 2.5D Temporal v1", [("label", "product name (2.5D)"), LIC, LIC]),
    ("- Sentinel-2 Level-2A", [S2, E0, END]),
    ("- Landsat Collection 2 Level-2", [LIC, LIC]),
    ("- ESA WorldCover 2021 v200", [("label", "product year (WorldCover 2021)"), LIC]),
    ("- JRC GHSL GHS-BUILT-S and GHS-POP R2023A", [LIC]),
    ("- WorldPop Global2 R2025A", [LIC]),
    ("- OpenStreetMap through Overpass and Nominatim", [LIC]),
    ("- geoBoundaries gbOpen NGA ADM2", [LIC]),
    ("The code in this repository is AI-assisted.", [P3, HUMAN, N40]),
    ("- Building-change values are less reliable", [OB16]),
    ("- The satellite data end in February", [END]),
    ("- The 2023/24 El Niño overlaps", cite(NOAA, "2023 (JAS 0.8 ...)", "2024 (DJF 1.1, JFM 0.8)")),
    ("Badmus, M. (2026).", [("label", "citation year"), E0, END]),
]

# ------------------------------------------------------------------ docs/policy_context.md
CONTEXT = [
    ("Nairametrics, \"Lagos records over 27,000", cite(
        NAIRA, "Lagos records over 27,000 building violations (title)", "May 22, 2026",
        "May 22, 2026")),
    ("The article reports figures disclosed at the Lagos State Government's 2026", cite(
        NAIRA, "2026 Annual Ministerial Press Briefing")),
    ("- 27,233 building violations were recorded", cite(
        NAIRA, "27,233 building violations", "between May 2023", "and April 2026")),
    ("- 4,711 developers obtained", cite(NAIRA, "4,711 developers obtained")),
    ("- The state has set a target of up to 45,000", cite(
        NAIRA, "up to 45,000 building permits annually", "May 2025 ministerial briefing",
        "37,266 building planning permits")),
    ("AllAfrica, republishing Premium Times", cite(
        AFRICA, "Targets 20,000 This Year (title)", "14 July 2026", "14 July 2026")),
    ("The article reports statements made at the 2026 Tree Planting Day", cite(
        AFRICA, "2026 Tree Planting Day celebration")),
    ("- The state aims to plant at least 20,000", cite(
        AFRICA, "at least 20,000 trees", "in 2026", "as of 30 June", "13,401 trees",
        "about 67 per cent of the target")),
    ("- A tree-planting campaign covers the 57", cite(AFRICA, "57 Local Government Areas and "
                                                              "Local Council Development Areas",
                                                      "for at least 12 months")),
    ("C40 Cities Climate Leadership Group, \"Request for Proposal", cite(C40, "July 2026")),
    ("- A consultant is sought to carry out", cite(C40, "Lagos Climate Action Plan (CAP 2026)")),
    ("- The assessment \"needs to establish", cite(
        C40, "planned by the city in 2027", "page 3 (footer page number)",
        "being scheduled for 2027", "5 – 10 years' timeline", "5 – 10 years' timeline")),
    ("- According to the document, \"no dedicated", cite(C40, "page 3 (footer page number)")),
    ("- The 2021 Climate Change Risk Assessment identified", cite(
        C40, "The CCRA (2021)", "across all 20 LGAs", "page 3 (footer page number)")),
    ("- Land-use change, densification and loss of green cover", cite(
        C40, "page 3 (footer page number)", "page 4 (footer page number)")),
    ("- Final report submission and project completion", cite(
        C40, "[29] [January] 2027", "[29] [January] 2027")),
    ("ThisDay, \"Cambridge, UNILAG Identify 10", cite(
        THISDAY, "Identify 10 Locations (title)", "2026/07/23 (URL date)", "2026/07/23 (URL date)")),
    ("- The project is a three-year study", cite(THISDAY, "across 10 selected settlements")),
    ("- The article does not name the 10", cite(THISDAY, "10 selected settlements")),
    ("- the Lagos State Urban and Regional Planning and Development Law (2010)", cite(
        LAWNG, "commencement 5 July 2010 (legal publisher's page)")),
    ("- the Lagos State Physical Planning Permit Regulations, 2019", cite(
        EPP, "REVISED_LASPPPA_REGULATION_2019 (file name; document could not be opened)")),
    ("LASBCA publishes a guide that restates setback rules", cite(
        LASBCA, "The Building Control Agency Regulations, 2019")),
]

# ------------------------------------------------------------------ docs/phase4_priority_rules.md
P2L, RANK1 = ("label", "phase number"), ("label", "rank definition (rank 1 = most severe)")
KM15 = ("derived", "3 x 3 cells x grid cell size, in km",
        lambda: 3 * tf.const_value("grid_spec.cell") / 1000, "data/interim/grid_spec.json")
RULES = [
    ("# Phase 4: rules for the priority-area lists", [P2L]),
    ("Population is WorldPop 2023", [POP23, P2L]),
    ("For each list, every unit", [RANK1]),
    ("The units are the hotspot cells defined in Phase 2", [P2L]),
    ("- A2, hotspot clusters:", [P2L]),
    ("| Residents in hotspot cells |", [RANK1]),
    ("| Mean relative LST change |", [RANK1]),
    ("| Mean relative NDVI loss |", [RANK1]),
    ("- Still mostly vegetated in the late period:", [
        ("label", "Sentinel-2 composite pixel size (20 m, docs/phase1_notes.md)"),
        ("const", "vegetation.NDVI_VEGETATED", "x")]),
    ("- Building share grew above the hotspot cut:", [P2L]),
    ("| Building change |", [RANK1]),
    ("| Residents nearby |", [KM15, KM15, RANK1]),
    ("Area names are added for the top 10", [("const", "phase4.TOP", "x"),
                                             ("const", "phase4.TOP", "x")]),
    ("- Each cluster centre or cell centre is reverse-geocoded",
     [("label", "Nominatim zoom level (16, in phase4.names_for)")]),
    ("- For the three case-study clusters", [P2L]),
    ("- Where the zoom-16 lookup returned no local name", [
        ("label", "Nominatim zoom level"), ("const", "phase4.ZOOMS[1]", "x"),
        ("const", "phase4.ZOOMS[2]", "x")]),
    ("1. Area names (labels only; ranks unchanged).", [LIST]),
    ("- Reason: three of the top 10 protect cells", [("const", "phase4.TOP", "x"),
                                                    ("const", "phase4.ZOOMS[0]", "x")]),
]

# ------------------------------------------------------------------ docs/data_sources.md
OBCAT = "https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_Research_open-buildings-temporal_v1"
WC = "https://esa-worldcover.org/en/data-access"
GHSL = "https://human-settlement.emergency.copernicus.eu/GHSLhowToCite.php"
WPOP = "https://hub.worldpop.org/doi/10.5258/SOTON/WP00839"
GEOB = "https://www.geoboundaries.org/"
LS = "https://planetarycomputer.microsoft.com/api/stac/v1/collections/landsat-c2-l2"
SOURCES = [
    ("## Google Open Buildings 2.5D Temporal", [("label", "product name (2.5D)")]),
    ("Use: building presence and height for mid-2016", [OB16, OB23]),
    ("Licence: \"The data is shared under the Creative Commons", cite(OBCAT, "CC-BY 4.0", "4.0")),
    ("Citation: Wojciech Sirko", cite(OBCAT, "Sentinel-2", "2023")),
    ("## Sentinel-2 Level-2A", [S2]),
    ("Use: NDVI and BSI composites, dry seasons", [E0, E1, L0, L1]),
    ("Attribution, as required by the legal notice", [E0, END]),
    ("## Landsat Collection 2 Level-2", [LIC, LIC]),
    ("Use: land surface temperature (ST_B10)", [L8, ("label", "sensor name (Landsat 9)")]),
    ("Citation (DOI given by the collection)", cite(LS, "Landsat 8-9 OLI/TIRS C2 L2", "8-9",
                                                    "Collection 2", "Level-2")),
    ("## ESA WorldCover 2021 (v200)", [("label", "product year (WorldCover 2021)")]),
    ("Licence: CC BY 4.0. \"The ESA WorldCover product", cite(WC, "CC BY 4.0")),
    ("Attribution: \"© ESA WorldCover project 2021", cite(WC, "project 2021", "data (2021)")),
    ("Citation: Zanaga, D.", cite(WC, "2022", "10 m", "2021 v200")),
    ("- GHS-BUILT-S R2023A, 2015 and 2020", [("const", "buildings.GHSL_EPOCHS[0]", "x"),
                                            ("const", "buildings.GHSL_EPOCHS[1]", "x")]),
    ("- GHS-POP R2023A, epoch 2020", [("label", "GHS-POP epoch (phase2.GHSPOP_STEM E2020)")]),
    ("Licence: \"The GHSL data are provided under", cite(GHSL, "4.0", "4.0")),
    ("Citation: Pesaresi, M.", cite(GHSL, "(2024)", "17(1)", "17(1)")),
    ("## WorldPop Global2 R2025A, constrained, 100 m", [("label", "product resolution (100 m)")]),
    ("Use: residents per cell, 2016 and 2023", [POP16, POP23]),
    ("Licence: Creative Commons Attribution 4.0 International. Source", cite(WPOP, "4.0")),
    ("Citation: Bondarenko M.", cite(WPOP, "2025", "2015-2030", "2015-2030", "3 arc",
                                     "INV-045237", "DOI:10.5258/SOTON/WP00839")),
    ("- place names, through Nominatim reverse geocoding", [P2L, P2L]),
    ("Licence: Open Data Commons Open Database License (ODbL) 1.0.", [LIC]),
    ("Licence: CC BY 4.0. \"Our license requires", cite(GEOB, "CC BY 4.0")),
    ("Citation: Runfola, D. et al. (2020)", cite(GEOB, "(2020)", "15(4)", "15(4)")),
    ("Use: visual spot-check of 40 cells in Phase 3", [N40, P2L]),
]

DOCS = {
    "docs/phase4_priority_rules.md": (RULES, []),
    "docs/data_sources.md": (SOURCES, []),
    "docs/findings.md": (FINDINGS, FINDINGS_CLAIMS),
    "docs/policy_brief.md": (BRIEF, BRIEF_CLAIMS),
    "README.md": (README, []),
    "docs/policy_context.md": (CONTEXT, []),
}
