# Building out the breeze: where Lagos is heating up as it builds

Policy brief for the Lagos State Parks and Gardens Agency (LASPARK), the Lagos State Building Control Agency (LASBCA), the Ministry of Physical Planning and Urban Development, the Lagos Resilience Office and C40 Cities. November 2026.

## The issue

Lagos is building fast on its fringe, often on land that was bush, farmland or wetland. Buildings and paving replace vegetation, which cools the ground through shade and evaporation. As a result, the surfaces of newly built areas can heat up compared with the rest of the city.

Lagos has no dedicated heat action plan yet. A heat assessment commissioned through C40 is to provide the evidence for an Urban Heat Action Plan that the city plans for 2027 (C40 request for proposals, July 2026). This brief summarises satellite evidence on where building growth, vegetation loss and surface warming coincide, and where action could cool areas now or keep green areas from being lost.

## What the evidence shows

The study compares 500 m cells across Lagos State between 2016 and 2023 (buildings) and between the dry seasons of 2016–18 and 2023–25 (vegetation and land surface temperature). Temperatures are measured relative to the Lagos-wide median on each date, so a warmer cell is one that warmed more than the city as a whole.

1. Where buildings were added, the ground warmed more than the city as a whole. Across 9,157 cells, each 10 percentage points of added building cover went with 0.579 °C more relative surface warming (95% interval 0.447 to 0.711). The result held in all 13 model runs. Besides the main run, these dropped sand-filled sites, coastal cells, industrial cells or cells near major roads, among other checks.
2. Building growth went with vegetation loss, and vegetation loss accounts for about half of the warming association. When vegetation change is included in the model, the building coefficient falls by about half (share 0.522).
3. The warming is strongest where green land is built over. Per 10 percentage points of added building cover, cells that were mostly vegetated in 2016 warmed by 0.731 °C. Cells that were already built up warmed by 0.190 °C.
4. About 896,165 people live in heat hotspots. 329 cells are in the top fifth of both building growth and relative surface warming. By WorldPop 2023 estimates, they hold 8.5% of residents in the cells studied.

## Where to act

"Cool now" ranks LGAs that contain hotspot cells. The ranks of three measures get equal weight: residents in hotspot cells, average relative warming, and average vegetation loss.

| Rank | LGA | Hotspot cells | Residents in hotspot cells | Mean relative warming (°C) |
|---|---|---|---|---|
| 1 | Eti Osa | 56 | 56,761 | 2.12 |
| 2 | Epe | 30 | 35,412 | 2.19 |
| 3 | Oshodi/Isolo | 24 | 77,077 | 1.64 |
| 4 | Alimosho | 74 | 353,967 | 1.62 |
| 5 | Shomolu | 9 | 47,569 | 1.86 |

The highest-ranked hotspot cluster lies in Lekki Peninsula II, around the Chevron area of Eti-Osa.

The 2021 Climate Change Risk Assessment, as cited in the C40 request for proposals, named Lagos Island, Lagos Mainland, Agege, Kosofe and Oshodi-Isolo as the LGAs most affected by heat, which describes where heat levels are highest, whereas this study measures where surface heat grew fastest relative to the rest of the city, which is concentrated on the expanding fringe (for example Eti-Osa, Epe and Alimosho); the two views are complementary, not contradictory.

"Protect before it is lost" lists cells that are still mostly green, with at least half their area vegetated in 2023–25, but where building cover grew by at least 5.8 percentage points. There are 104 such cells in 10 LGAs; Ikorodu has the most, 33. The cells are ranked by building growth and by residents within 1.5 km.

| Rank | Area (LGA) | Building growth (percentage points) | Residents nearby |
|---|---|---|---|
| 1 | Trade Fair Complex (Amuwo Odofin) | 13.3 | 19,852 |
| 2 | near Ijegun (Alimosho) | 8.7 | 15,647 |
| 3 | near Ijegun (Alimosho) | 8.7 | 16,150 |
| 4 | Trade Fair Complex (Amuwo Odofin) | 9.3 | 11,877 |
| 5 | near Ikotun (Oshodi/Isolo) | 7.5 | 24,232 |

Names come from OpenStreetMap. "near" marks either the nearest named place from a wider lookup, or a name that OpenStreetMap places in another LGA. Ikotun, for example, is placed in Alimosho, while this cell lies in Oshodi/Isolo in the boundary data. LGAs follow the boundary data.

The full lists and maps are in `data/processed/priority_areas_lga.csv`, `data/processed/priority_cells.csv` and `figures/phase4/`.

## Recommendations

- Protect green space on the expanding fringe, starting with the "protect" cells and the fringe LGAs with the most of them (Ikorodu, Eti Osa and Ibeju Lekki). Responsible: Ministry of Physical Planning and Urban Development, LASPARK.
- Enforce setback and open-space rules in these fringe areas, where new building is replacing vegetation. LASBCA recorded 27,233 building violations between May 2023 and April 2026 (Nairametrics, 22 May 2026). Responsible: LASBCA, Lagos State Physical Planning Permit Authority.
- Target tree planting at the "cool now" LGAs, within the state's goal of at least 20,000 trees in 2026, of which 13,401 had been planted by 30 June (Premium Times via AllAfrica, 14 July 2026). Responsible: LASPARK, Ministry of Environment and Water Resources, LGA "Tree Custodians".
- Publish building-violation data by LGA, so that enforcement can be compared with where building growth and warming occur. Responsible: LASBCA.
- Add heat to planning approvals: consider vegetation loss and surface heat when approving development in vegetated areas. Use these lists as one input to the 2027 heat action plan and the C40 heat assessment. Responsible: Ministry of Physical Planning and Urban Development, Lagos Resilience Office, C40.

## Limits

The study measures ground-surface temperature from satellites, not the air temperature people feel. It shows associations, not proven causes. Building-change figures are less reliable in areas that were already dense. Nearby cells are similar to one another, so the uncertainty ranges may be too narrow. Population figures are modelled estimates, and the data end in February 2025. The analysis was coded and partly checked with an AI coding assistant (Claude Code). A person confirmed a sample of the AI's imagery checks. The full findings and limits are in `docs/findings.md`, and the sources for the policy facts are in `docs/policy_context.md`.
