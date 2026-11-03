import numpy as np
import pandas as pd

from lagos_heat import phase4 as p4


def test_rank_table_equal_weights_and_tiebreak():
    t = pd.DataFrame({"id": ["a", "b", "c"], "pop": [10.0, 30.0, 20.0],
                      "heat": [3.0, 1.0, 2.0]})
    r = p4.rank_table(t, ["pop", "heat"], "pop", "id")
    # a: ranks 3,1 -> 2; b: 1,3 -> 2; c: 2,2 -> 2; all tied, so ordered by pop
    assert r["id"].tolist() == ["b", "c", "a"]
    assert r["priority_rank"].tolist() == [1, 2, 3]
    assert (r["mean_rank"] == 2).all()


def test_rank_table_ties_share_lowest_rank():
    t = pd.DataFrame({"id": [1, 2, 3], "x": [5.0, 5.0, 1.0]})
    r = p4.rank_table(t, ["x"], "x", "id")
    assert r["rank_x"].tolist() == [1, 1, 3] and r["id"].tolist() == [1, 2, 3]


def test_neighbourhood_sum_3x3_with_edges_and_nan():
    v = np.arange(1.0, 10.0).reshape(3, 3)
    v[0, 0] = np.nan
    s = p4.neighbourhood_sum([1, 0], [1, 2], v, 3, 3)
    assert s[0] == np.nansum(v)            # centre sees all 9 (NaN as 0)
    assert s[1] == 2 + 3 + 5 + 6           # corner (0, 2) sees 4 cells


def test_local_name_prefers_most_local():
    assert p4.local_name({"suburb": "S", "neighbourhood": "N"}) == "N"
    assert p4.local_name({"county": "Alimosho"}) == ""


def test_cool_now_metrics():
    hot = pd.DataFrame({"cell_id": [1, 2, 3], "lga": ["A", "A", "B"], "cluster": [7, 7, 8],
                        "pop_2023": [10.0, 20.0, 5.0], "lst_rel_change": [1.0, 3.0, 4.0],
                        "ndvi_rel_change_p80": [-0.1, -0.3, 0.0]})
    lga, cl = p4.cool_now(hot)
    a = lga.set_index("lga").loc["A"]
    assert a["residents_hotspot"] == 30 and a["mean_lst_rel_change"] == 2.0
    assert np.isclose(a["mean_ndvi_loss"], 0.2)
    assert set(cl["cluster"]) == {7, 8}


def test_protect_selection_rule():
    cells = pd.DataFrame({"cell_id": [1, 2, 3, 4], "lga": ["A", "A", "B", "B"],
                          "veg_share_late": [0.6, 0.4, 0.5, 0.9],
                          "bshare_change": [0.10, 0.20, 0.058, 0.01],
                          "residents_nearby": [5.0, 9.0, 7.0, 1.0],
                          "pop_2023": [1.0, 1.0, 1.0, 1.0]})
    ranked, by_lga = p4.protect(cells, 0.058)
    assert sorted(ranked["cell_id"]) == [1, 3]   # 2 fails vegetation, 4 fails change
    assert by_lga["n_cells"].sum() == 2


def test_label_name_near_rules():
    assert p4.label_name("Ikotun", "Alimosho", "Alimosho", False) == ("Ikotun", "")
    name, note = p4.label_name("Ikotun", "Alimosho", "Oshodi/Isolo", False)
    assert name == "near Ikotun" and "Alimosho" in note and "Oshodi/Isolo" in note
    name, note = p4.label_name("Ijegun", "Alimosho", "Alimosho", True)
    assert name == "near Ijegun" and "wider" in note
    assert p4.label_name("Akoka", "", "Shomolu", False) == ("Akoka", "")   # no OSM LGA
    assert p4.label_name("", "Alimosho", "Alimosho", False)[0] == ""
    assert p4.label_name("Ajah", "Eti-Osa", "Eti Osa", False) == ("Ajah", "")  # spelling
