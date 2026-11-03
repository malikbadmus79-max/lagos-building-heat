import json

import numpy as np
import pandas as pd

from lagos_heat import phase5_map as p5


def test_heat_class_bins_and_edges():
    v = [-3.0, -2.0, -1.5, -0.25, 0.0, 0.25, 1.5, 2.0, 5.0]
    assert p5.heat_class(v).tolist() == [0, 1, 1, 3, 3, 4, 5, 6, 6]
    assert len(p5.COLOURS) == len(p5.BIN_LABELS) == len(p5.BINS) + 1


def test_encode_nodes_roundtrip_to_1e5_degrees():
    lat = np.array([[6.712345, 6.71], [6.40001, 6.4]])
    lon = np.array([[2.7, 2.80456], [2.7, 4.399999]])
    enc = p5.encode_nodes(lat, lon)
    for i, (la, lo) in enumerate(zip(lat.ravel(), lon.ravel())):
        dla, dlo = p5.decode_node(enc, i)
        assert abs(dla - la) <= 1e-5 and abs(dlo - lo) <= 1e-5


def test_lga_stats_figures():
    cells = pd.DataFrame({"lga": ["A", "A", "A", "B"], "hotspot": [True, False, True, False],
                          "pop_2023": [10.0, 99.0, 5.0, 7.0],
                          "lst_rel_change": [2.0, 0.0, 1.0, -1.0]})
    cool = pd.DataFrame({"lga": ["A"], "priority_rank": [1], "mean_lst_rel_change": [1.5]})
    prot = pd.DataFrame({"lga": ["B"], "n_cells": [3]})
    s = p5.lga_stats(cells, cool, prot).set_index("lga")
    assert s.loc["A", "residents_in_hotspots"] == 15 and s.loc["A", "hotspot_squares"] == 2
    assert s.loc["A", "mean_heat_change_all_squares"] == 1.0
    assert s.loc["A", "cool_now_rank"] == 1 and np.isnan(s.loc["B", "cool_now_rank"])
    assert s.loc["A", "mean_heat_change_hotspot_squares"] == 1.5
    assert np.isnan(s.loc["B", "mean_heat_change_hotspot_squares"])
    assert s.loc["B", "protect_squares"] == 3 and s.loc["A", "protect_squares"] == 0


def test_page_fills_numbers_from_data_and_inlines_it():
    data = {"text": {"per10": 0.6, "hotspots": 329, "residents_k": 896000, "protect": 104,
                     "squares": 9157}, "cells": {"v": [12]}}
    html = p5.page(data)
    assert "about 0.6 &deg;C" in html and "about 896,000 people" in html
    assert "Across 9,157 squares" in html and "(104)" in html
    assert json.dumps(data, separators=(",", ":")) in html
    assert f"leaflet@{p5.LEAFLET}/dist/leaflet.js" in html and p5.LEAFLET_JS_SRI in html
    assert "relative LST" not in html and "NDVI" not in html  # plain language on the page
