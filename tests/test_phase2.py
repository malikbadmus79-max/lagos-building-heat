import numpy as np
import pandas as pd

from lagos_heat import phase2 as p2


def test_block_ids_5x5():
    b = p2.block_ids(np.array([0, 4, 5, 0, 9]), np.array([0, 4, 0, 5, 9]))
    assert b[0] == b[1]
    assert len({b[0], b[2], b[3], b[4]}) == 4


def test_rank_deciles_split_ties_evenly():
    d = p2.rank_deciles(np.zeros(100))
    assert np.bincount(d)[1:].tolist() == [10] * 10


def test_moran_i_hand_case():
    # 1 x 4 line, values 1, 1, -1, -1: z'Wz = 2, z'z = 4, n = S0 = 4, so I = 0.5
    r = p2.moran_i(np.array([1.0, 1.0, -1.0, -1.0]), np.zeros(4, int), np.arange(4),
                   n_perm=99)
    assert np.isclose(r["moran_i"], 0.5)
    assert r["n_islands"] == 0 and np.isclose(r["expected"], -1 / 3)


def test_moran_i_drops_islands():
    r = p2.moran_i(np.array([1.0, 1.0, -1.0, -1.0, 5.0]), np.array([0, 0, 0, 0, 9]),
                   np.array([0, 1, 2, 3, 9]), n_perm=9)
    assert r["n_islands"] == 1 and np.isclose(r["moran_i"], 0.5)


def test_group_means_ci_constant_within_groups():
    y = np.array([1.0, 1.0, 3.0, 3.0])
    g = p2.group_means_ci(y, np.array([1, 1, 2, 2]), np.array([0, 1, 0, 1]), n_boot=50)
    assert g["mean"].tolist() == [1.0, 3.0]
    assert np.allclose(g["ci_lo"], g["mean"]) and np.allclose(g["ci_hi"], g["mean"])


def test_weighted_coefs_unit_weights_equal_ols():
    rng = np.random.default_rng(0)
    x = np.column_stack([np.ones(50), rng.normal(size=50)])
    y = 2 + 3 * x[:, 1] + rng.normal(size=50)
    ols = np.linalg.lstsq(x, y, rcond=None)[0]
    assert np.allclose(p2.weighted_coefs(x, y, np.ones((1, 50)))[0], ols)


def test_bootstrap_weights_sum_to_cells_per_drawn_blocks():
    w = p2.bootstrap_weights(np.array([0, 0, 1, 2]), 20, np.random.default_rng(1))
    assert w.shape == (20, 4)
    assert (w[:, 0] == w[:, 1]).all()  # same block, same weight


def test_top_quintile_both():
    a = np.arange(10.0)
    hot, ca, cb = p2.top_quintile_both(a, a[::-1])
    assert not hot.any()
    hot, _, _ = p2.top_quintile_both(a, a)
    assert hot.sum() == 2


def test_label_clusters_rook():
    row, col = np.array([0, 0, 1, 2]), np.array([0, 1, 2, 2])
    lab = p2.label_clusters(row, col, np.ones(4, bool), 3, 3)
    assert lab[0] == lab[1] and lab[2] == lab[3] and lab[1] != lab[2]  # diagonal not joined


def test_choose_cases_third_needs_four_cells():
    t = pd.DataFrame({"cluster": [1, 2, 3], "n_cells": [5, 9, 3],
                      "sum_bshare_change": [1.0, 2.0, 0.5]})
    c = p2.choose_cases(t)
    assert c.loc[c["chosen"], "cluster"].tolist() == [2, 1]
    t.loc[2, "n_cells"] = 4
    assert p2.choose_cases(t)["chosen"].sum() == 3


def test_choose_cases_tie_broken_by_building_change():
    t = pd.DataFrame({"cluster": [1, 2], "n_cells": [5, 5], "sum_bshare_change": [1.0, 2.0]})
    assert p2.choose_cases(t)["cluster"].tolist() == [2, 1]


def test_h4_groups_exclusive():
    df = pd.DataFrame({"bshare_2016": [0.5, 0.1, 0.5, 0.1], "tree_frac_wc2021": [0.1, 0.7, 0.7, 0.1],
                       "veg_share_early": [0, 0, 0, 0], "main": [True] * 4})
    g, n_both = p2.h4_groups(df)
    assert g.tolist() == ["built", "vegetated", "neither", "neither"] and n_both == 1


def test_sign_label():
    assert p2.sign_label(0.1, 0.2) == "positive"
    assert p2.sign_label(-0.2, -0.1) == "negative"
    assert p2.sign_label(-0.1, 0.1) == "not distinguishable from zero"
