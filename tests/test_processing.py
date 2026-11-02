"""Unit tests for the Phase 1 processing functions (synthetic inputs only, no network)."""

import numpy as np
import pytest
from affine import Affine

from lagos_heat.common import GridSpec, aggregate, safe_div
from lagos_heat.temperature import lst_celsius, qa_clear, relative_lst
from lagos_heat.vegetation import (
    baseline_offset, bsi, ndvi, parse_boa_offsets, resolve_offset, s2_reflectance,
)


# ------------------------------------------------------------------ LST
def test_lst_scaling_known_values():
    # K = DN * 0.00341802 + 149.0 ; C = K - 273.15
    dn = np.array([40000, 44000, 0], dtype="uint16")
    c = lst_celsius(dn)
    assert c[0] == pytest.approx(40000 * 0.00341802 + 149.0 - 273.15, abs=1e-4)
    assert c[1] == pytest.approx(44000 * 0.00341802 + 149.0 - 273.15, abs=1e-4)
    assert np.isnan(c[2])  # DN 0 is fill


def test_lst_scaling_plausible_range():
    # DN 43000 corresponds to about 22.8 deg C, a plausible surface temperature
    assert lst_celsius(np.array([43000]))[0] == pytest.approx(22.82, abs=0.01)


def test_qa_clear_bits():
    clear = 1 << 6  # "clear" bit only
    cloud = (1 << 3) | (1 << 6)
    shadow = 1 << 4
    water = 1 << 7
    fill = 1
    dilated = 1 << 1
    qa = np.array([clear, cloud, shadow, water, fill, dilated], dtype="uint16")
    assert qa_clear(qa).tolist() == [True, False, False, False, False, False]


def test_relative_lst_subtracts_land_median():
    lst = np.array([[30.0, 32.0], [34.0, np.nan]], dtype="float32")
    ref = np.array([[True, True], [True, True]])
    rel, med = relative_lst(lst, ref)
    assert med == pytest.approx(32.0)
    assert rel[0, 0] == pytest.approx(-2.0)
    assert rel[1, 0] == pytest.approx(2.0)
    assert np.isnan(rel[1, 1])


def test_relative_lst_reference_mask_only():
    lst = np.array([10.0, 20.0, 100.0], dtype="float32")
    rel, med = relative_lst(lst, np.array([True, True, False]))
    assert med == pytest.approx(15.0)  # the masked-out 100 is not in the reference


def test_date_outliers_two_sd_rule():
    from lagos_heat.temperature import date_outliers

    medians = {f"d{i}": v for i, v in enumerate([32.0, 33.0, 31.5, 32.5, 33.5, 32.0, 24.0])}
    excluded, mean, sd = date_outliers(medians, n_sd=2.0)
    assert excluded == {"d6"}
    assert mean == pytest.approx(np.mean(list(medians.values())))
    assert date_outliers({"a": 30.0, "b": 31.0}, n_sd=2.0)[0] == set()  # too few dates


def test_hazy_dates_one_sided():
    from lagos_heat.s2_haze import hazy_dates  # AOT attempt (documented, not applied)

    base = [0.30, 0.32, 0.28, 0.31, 0.29, 0.30, 0.33, 0.27]
    medians = {f"d{i}": v for i, v in enumerate(base + [0.90])}
    excluded, mean, sd = hazy_dates(medians, n_sd=2.0)
    assert excluded == {"d8"}
    low = {f"d{i}": v for i, v in enumerate(base + [0.0])}  # unusually clear: kept
    assert hazy_dates(low, n_sd=2.0)[0] == set()
    assert hazy_dates({"a": 0.3, "b": 0.9}, n_sd=2.0)[0] == set()  # too few dates


def test_reference_forest_screen():
    from lagos_heat.s2_haze import reference_screen

    med = {"e/1": 0.70, "e/2": 0.55, "l/1": 0.72, "l/2": 0.69, "l/3": 0.62, "l/4": 0.10}
    n = {k: 1000 for k in med}
    n["l/4"] = 50  # cloud-covered reference forest: cannot be assessed, excluded
    level, thr, excluded, assessable = reference_screen(med, n, min_px=200, clearest_n=3,
                                                        max_drop=0.10)
    assert level == pytest.approx(0.70)  # median of 0.72, 0.70, 0.69
    assert thr == pytest.approx(0.60)
    assert excluded == {"e/2", "l/4"}  # 0.62 is within 0.10 of the level
    assert "l/4" not in assessable


def test_daystack_percentile(tmp_path, monkeypatch):
    import lagos_heat.common as common

    monkeypatch.setattr(common, "WORK", tmp_path)
    st = common.DayStack("p", ["a", "b", "c", "d", "e"], (1, 1), 100)
    for k, v in zip(st.keys, [0.1, 0.2, 0.3, 0.4, 0.5]):
        st.write(k, np.array([[v]]))
    assert st.median(q=80)[0][0, 0] == pytest.approx(0.42)  # linear interpolation
    assert st.median()[0][0, 0] == pytest.approx(0.3)


def test_landsat9_date_key():
    from types import SimpleNamespace
    from datetime import datetime

    from lagos_heat.temperature import date_key

    def item(platform):
        return SimpleNamespace(datetime=datetime(2024, 1, 7),
                               properties={"landsat:wrs_path": 191, "platform": platform})
    assert date_key(item("landsat-8")) == "2024-01-07_p191"
    assert date_key(item("landsat-9")) == "2024-01-07_p191_l9"


def test_daystack_range_rule(tmp_path, monkeypatch):
    import lagos_heat.common as common

    monkeypatch.setattr(common, "WORK", tmp_path)
    st = common.DayStack("lst", ["d1", "d2"], (1, 3), 100)
    st.write("d1", np.array([[5.0, 30.0, 70.0]]))
    st.write("d2", np.array([[20.0, 34.0, 40.0]]))
    m = np.array([[True, True, True]])
    assert st.count_outside(["d1", "d2"], m, 15, 65) == {"d1": (1, 1), "d2": (0, 0)}
    med, cnt = st.median(vmin=15, vmax=65)
    assert cnt.tolist() == [[1, 2, 1]]
    assert med[0, 0] == pytest.approx(20.0) and med[0, 2] == pytest.approx(40.0)
    assert st.slice_medians(["d1"], m, 15, 65)["d1"] == pytest.approx(30.0)


def test_daystack_extend_keeps_stored_days(tmp_path, monkeypatch):
    import lagos_heat.common as common

    monkeypatch.setattr(common, "WORK", tmp_path)
    st = common.DayStack("s", ["d1", "d3"], (1, 2), 100)
    st.write("d1", np.array([[1.0, 2.0]]))
    st.write("d3", np.array([[3.0, np.nan]]))
    del st
    st = common.DayStack("s", ["d1", "d2", "d3"], (1, 2), 100)  # d2 added in between
    assert st.is_done("d1") and st.is_done("d3") and not st.is_done("d2")
    med, cnt = st.median(["d1", "d3"])
    assert med[0, 0] == pytest.approx(2.0) and cnt.tolist() == [[2, 1]]
    st.write("d2", np.array([[5.0, 5.0]]))
    assert st.median()[0][0, 0] == pytest.approx(3.0)


# ---------------------------------------------------- Sentinel-2 offset
@pytest.mark.parametrize("baseline,expected", [
    ("02.12", 0), ("03.00", 0), ("03.01", 0), ("04.00", -1000), ("05.09", -1000),
    ("05.11", -1000),
])
def test_baseline_offset(baseline, expected):
    assert baseline_offset(baseline) == expected


def test_baseline_offset_missing_raises():
    with pytest.raises(ValueError):
        baseline_offset(None)


def test_reflectance_with_and_without_offset():
    dn = np.array([2000, 0], dtype="uint16")
    old = s2_reflectance(dn, "02.12")
    new = s2_reflectance(dn, "05.11")
    assert old[0] == pytest.approx(0.2)
    assert new[0] == pytest.approx(0.1)  # (2000 - 1000) / 10000
    assert np.isnan(old[1]) and np.isnan(new[1])  # DN 0 is no data


def test_offset_changes_ndvi_as_expected():
    # The same surface (red 0.05, NIR 0.30) encoded under each baseline gives the same NDVI
    # only if the offset is applied.
    red_old, nir_old = np.array([500]), np.array([3000])
    red_new, nir_new = red_old + 1000, nir_old + 1000
    n_old = ndvi(s2_reflectance(nir_old, "02.12"), s2_reflectance(red_old, "02.12"))
    n_new = ndvi(s2_reflectance(nir_new, "05.11"), s2_reflectance(red_new, "05.11"))
    assert n_new[0] == pytest.approx(n_old[0], abs=1e-6)
    n_wrong = ndvi(s2_reflectance(nir_new, "02.12"), s2_reflectance(red_new, "02.12"))
    assert n_wrong[0] < n_old[0] - 0.1  # ignoring the offset biases NDVI low


MTD_WITH_OFFSET = """<n1:Level-2A_User_Product><BOA_ADD_OFFSET_VALUES_LIST>
<BOA_ADD_OFFSET band_id="0">-1000</BOA_ADD_OFFSET>
<BOA_ADD_OFFSET band_id="1">-1000</BOA_ADD_OFFSET>
<BOA_ADD_OFFSET band_id="12">-1000</BOA_ADD_OFFSET>
</BOA_ADD_OFFSET_VALUES_LIST></n1:Level-2A_User_Product>"""


def test_parse_boa_offsets():
    assert parse_boa_offsets(MTD_WITH_OFFSET) == {0: -1000, 1: -1000, 12: -1000}
    assert parse_boa_offsets("<n1:Level-2A_User_Product></n1:Level-2A_User_Product>") == {}


def test_resolve_offset_metadata_wins_over_date_and_baseline():
    # An archive scene acquired in 2017 but reprocessed with baseline 05.00 carries the
    # offset in its metadata; the acquisition date plays no role.
    r = resolve_offset("05.00", {0: -1000, 1: -1000})
    assert r["boa_add_offset"] == -1000 and r["offset_source"].startswith("metadata")
    # Old baseline, no metadata offset -> 0
    r = resolve_offset("02.12", {})
    assert r["boa_add_offset"] == 0 and r["offset_sources_agree"]
    # Disagreement is recorded, metadata value applied
    r = resolve_offset("02.12", {0: -1000})
    assert r["boa_add_offset"] == -1000 and not r["offset_sources_agree"]


def test_reflectance_explicit_offset_overrides_baseline():
    dn = np.array([2000], dtype="uint16")
    assert s2_reflectance(dn, "02.12", offset=-1000)[0] == pytest.approx(0.1)


def test_daystack_relative_median(tmp_path, monkeypatch):
    import lagos_heat.common as common

    monkeypatch.setattr(common, "WORK", tmp_path)
    st = common.DayStack("t", ["d1", "d2", "d3"], (1, 2), 100)
    st.write("d1", np.array([[0.2, 0.4]]))
    st.write("d2", np.array([[0.3, np.nan]]))
    st.write("d3", np.array([[0.5, 0.6]]))
    med, cnt = st.median()
    assert med[0, 0] == pytest.approx(0.3) and cnt.tolist() == [[3, 2]]
    meds = st.slice_medians(["d1", "d2", "d3"], np.array([[True, True]]))
    assert meds["d1"] == pytest.approx(0.3) and meds["d2"] == pytest.approx(0.3)
    rel, _ = st.median(offsets=meds)
    # pixel 0: d1 0.2-0.3=-0.1, d2 0.3-0.3=0, d3 0.5-0.55=-0.05 -> median -0.05
    assert rel[0, 0] == pytest.approx(-0.05, abs=1e-6)


# ------------------------------------------------------- index formulas
def test_ndvi_formula():
    nir, red = np.array([0.5, 0.2, 0.0]), np.array([0.1, 0.2, 0.0])
    out = ndvi(nir, red)
    assert out[0] == pytest.approx((0.5 - 0.1) / (0.5 + 0.1))
    assert out[1] == pytest.approx(0.0)
    assert np.isnan(out[2])  # zero denominator


def test_bsi_formula():
    blue, red, nir, swir = 0.05, 0.10, 0.30, 0.20
    expected = ((swir + red) - (nir + blue)) / ((swir + red) + (nir + blue))
    assert bsi([blue], [red], [nir], [swir])[0] == pytest.approx(expected)
    # bare soil: high SWIR and red, low NIR -> positive
    assert bsi([0.10], [0.25], [0.28], [0.40])[0] > 0
    # vegetation: high NIR -> negative
    assert bsi([0.03], [0.04], [0.40], [0.15])[0] < 0


# ------------------------------------------------------ grid aggregation
def test_pixel_cells_and_aggregate_mean():
    # 2 x 2 grid of 500 m cells; 10 x 10 raster of 100 m pixels covering it exactly
    spec = GridSpec(xmin=1000.0, ymax=5000.0, ncol=2, nrow=2)
    transform = Affine(100, 0, 1000, 0, -100, 5000)
    cells = spec.pixel_cells(transform, (10, 10))
    assert cells[0, 0] == 0 and cells[0, 9] == 1 and cells[9, 0] == 2 and cells[9, 9] == 3
    values = np.zeros((10, 10))
    values[:5, :5], values[:5, 5:], values[5:, :5], values[5:, 5:] = 1, 2, 3, 4
    total, count = aggregate(cells, values=values, ncells=spec.ncells)
    assert count.tolist() == [25, 25, 25, 25]
    assert safe_div(total, count).tolist() == [1, 2, 3, 4]


def test_aggregate_with_mask_and_outside_pixels():
    spec = GridSpec(xmin=0.0, ymax=1000.0, ncol=1, nrow=1, cell=1000)
    # raster extends beyond the grid by one pixel column on the right
    transform = Affine(500, 0, 0, 0, -500, 1000)
    cells = spec.pixel_cells(transform, (2, 3))
    assert (cells[:, 2] == -1).all()
    values = np.array([[1.0, 3.0, 99.0], [np.nan, 5.0, 99.0]])
    mask = np.isfinite(values)
    total, count = aggregate(cells, values=np.where(mask, values, 0), mask=mask,
                             ncells=spec.ncells)
    assert count[0] == 3 and total[0] == pytest.approx(9.0)


def test_cell_of_xy_and_geobox_alignment():
    spec = GridSpec(xmin=466500.0, ymax=741000.0, ncol=370, nrow=74)
    assert spec.cell_of_xy(466500.0 + 1, 741000.0 - 1) == 0
    assert spec.cell_of_xy(466500.0 + 501, 741000.0 - 1) == 1
    assert spec.cell_of_xy(0.0, 0.0) == -1
    gb = spec.geobox(30, anchor=15)
    assert gb.transform.c % 30 == 15 and gb.transform.f % 30 == 15
    gb20 = spec.geobox(20)
    assert gb20.shape == (74 * 25, 370 * 25)


def test_share_above_ignores_nan_and_outside_pixels():
    from lagos_heat.vegetation import share_above

    cells = np.array([[0, 0, 1], [1, -1, 2]])
    arr = np.array([[0.5, 0.3, np.nan], [0.6, 0.9, np.nan]], dtype="float32")
    got = share_above(cells, arr, 0.4, np.array([0, 1, 2]), 3)
    assert got[0] == 0.5 and got[1] == 1.0 and np.isnan(got[2])
