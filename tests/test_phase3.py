import numpy as np
import pandas as pd
from affine import Affine

from lagos_heat import trace_findings as tf
from lagos_heat import validate as v

SPEC = {"xmin": 1000.0, "ymax": 5000.0, "ncol": 10, "nrow": 8, "cell": 500}


def test_cell_bounds_row_major_from_top():
    assert v.cell_bounds(0, SPEC) == (1000.0, 4500.0, 1500.0, 5000.0)
    assert v.cell_bounds(12, SPEC) == (2000.0, 4000.0, 2500.0, 4500.0)  # row 1, col 2


def test_centres_in_cell_edge_rule():
    # 25 m pixels offset by 12.5 m: pixel centres fall exactly on all four cell edges
    tr = Affine(25, 0, 987.5, 0, -25, 5012.5)
    m = v.centres_in_cell(tr, (30, 30), (1000.0, 4500.0, 1500.0, 5000.0))
    xs = 987.5 + (np.arange(30) + 0.5) * 25   # 1000, 1025, ...
    ys = 5012.5 - (np.arange(30) + 0.5) * 25  # 5000, 4975, ...
    assert xs[0] == 1000 and m[:, 0].any()     # x == x0 inside
    assert xs[20] == 1500 and not m[:, 20].any()  # x == x1 outside
    assert ys[0] == 5000 and m[0, :].any()     # y == y1 inside
    assert ys[20] == 4500 and not m[20, :].any()  # y == y0 outside
    assert m.sum() == 20 * 20


def test_stratified_sample_sizes_and_strata():
    rng = np.random.default_rng(0)
    n = 400
    df = pd.DataFrame({"cell_id": np.arange(n), "hotspot": np.arange(n) < 40,
                       "bshare_change": np.r_[np.full(40, 0.2), np.full(60, 0.1),
                                              np.zeros(300)]})
    s = v.stratified_sample(df, 0.058, 0.006, rng)
    assert s.groupby("stratum").size().to_dict() == {
        "high_change_not_hotspot": 15, "hotspot": 15, "low_change": 10}
    assert s["cell_id"].is_unique
    assert s.loc[s.stratum == "hotspot", "cell_id"].lt(40).all()
    assert s.loc[s.stratum == "high_change_not_hotspot", "cell_id"].between(40, 99).all()
    again = v.stratified_sample(df, 0.058, 0.006, np.random.default_rng(0))
    assert again.equals(s)  # fixed seed, same draw


def test_compare_tolerance_and_nan():
    d, ok = v.compare(np.array([0.1, 0.1, np.nan, 0.1]), np.array([0.105, 0.2, np.nan, np.nan]),
                      0.01)
    assert ok.tolist() == [True, False, True, False]


def test_ndvi_direction():
    assert [v.ndvi_direction(x) for x in (-0.05, 0.03, 0.01, np.nan)] == [
        "decrease", "increase", "little change", "no data"]


def test_wayback_pick_closest():
    rel = [(1, "2016-01-13"), (2, "2016-06-27"), (3, "2017-01-11")]
    assert v.wayback_pick(rel, "2016-06-30")[0] == 2


def test_name_candidates_tally():
    resp = [{"case": 1, "query": "centre_z18", "response": {"address": {"suburb": "A"}}},
            {"case": 1, "query": "cell_5_z16", "response": {"address": {"suburb": "A",
                                                                         "neighbourhood": "B"}}},
            {"case": 1, "query": "cell_6_z16", "response": {"error": "timeout"}}]
    t = v.name_candidates(resp)
    a = t[t["name"] == "A"].iloc[0]
    assert a["n_points"] == 2 and a["from_centre"]
    assert not t[t["name"] == "B"].iloc[0]["from_centre"]


def test_numbers_parse_minus_commas_and_skip_code():
    toks = tf.numbers("x −0.542 (−1,008.5 to 9,157) `H1_primary 0.3` D10 25th 2016–18, 4.")
    assert [t[0] for t in toks] == ["−0.542", "−1,008.5", "9,157", "2016", "18", "4"]
    assert toks[0][1:] == (-0.542, 3) and toks[1][1] == -1008.5


def test_matches_rounding():
    assert tf.matches(0.579, 3, 0.579070060309045)
    assert not tf.matches(0.578, 3, 0.579070060309045)
    assert tf.matches(896165, 0, 896164.5627)
    assert not tf.matches(1.0, 2, float("nan"))


def test_trace_flags_wrong_number(monkeypatch):
    lines = ["Main sample: 9,157 cells and 9,999 blocks."]
    monkeypatch.setattr(tf, "SPEC", [("Main sample:", [("const", "grid_spec.cell", "x"),
                                                        ("label", "test")])])
    monkeypatch.setattr(tf, "CLAIMS", [])
    monkeypatch.setattr(tf, "const_value", lambda name: 9157.0)
    out, untraced = tf.trace(lines)
    assert out["match"].tolist() == ["yes", "n/a"] and untraced == []
    monkeypatch.setattr(tf, "const_value", lambda name: 9158.0)
    assert tf.trace(lines)[0]["match"].tolist() == ["no", "n/a"]


def test_keep_manual_columns_preserves_entries():
    new = pd.DataFrame({"cell_id": [1, 2, 3], "x": [0.1, 0.2, 0.3]})
    old = pd.DataFrame({"cell_id": ["2", "1"], "seen_2016": ["dense buildings", ""],
                        "notes": ["", "cloud"], "human_agrees_with_ai": ["yes", ""]})
    out = v.keep_manual_columns(new, old)
    assert out["seen_2016"].tolist() == ["", "dense buildings", ""]
    assert out["notes"].tolist() == ["cloud", "", ""]
    assert out["human_agrees_with_ai"].tolist() == ["", "yes", ""]
    assert out["agrees_with_data"].tolist() == ["", "", ""]
    assert v.keep_manual_columns(new, None)["seen_2023"].tolist() == ["", "", ""]
