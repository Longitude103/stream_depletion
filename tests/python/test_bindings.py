"""Lock the in-crate identities through the PyO3 bindings.

Numbers match src/sdf.rs, src/glover_infinite.rs, and the URF cases
called out in the request (month-gap, duplicates, day normalization).
"""

from datetime import date

import pytest

from stream_depletion import (
    combined_urf_results,
    combined_urf_results_dated,
    daily_urf_mass_by_reach,
    gpd_per_ft_to_ft2_per_day,
    lag_glover_infinite,
    lag_sdf,
    sdf_from_glover,
    urf_lagging,
    urf_lagging_daily,
    urf_mass_by_reach,
)


def test_january_100_af_sdf_265_matches_jenkins_volume():
    rows = lag_sdf({"2025-01-01": 100.0}, sdf=265, total_months=6)
    expected = [
        ("2025-01-01", 0.944405),
        ("2025-02-01", 7.169295),
        ("2025-03-01", 10.023970),
        ("2025-04-01", 7.801918),
        ("2025-05-01", 6.279747),
        ("2025-06-01", 4.833422),
    ]
    assert len(rows) == len(expected)
    for (got_date, got), (exp_date, exp) in zip(rows, expected):
        assert got_date == exp_date
        assert got == pytest.approx(exp, abs=5e-6)


def test_recharge_accretion_is_negative_mirror():
    pump = lag_sdf({"2025-01-01": 100.0}, 265, 6)
    rech = lag_sdf({"2025-01-01": -100.0}, 265, 6)
    assert all(v < 0.0 for _, v in rech)
    assert rech[0][1] == pytest.approx(-0.944405, abs=5e-6)
    for (_, p), (_, r) in zip(pump, rech):
        assert p + r == pytest.approx(0.0, abs=1e-12)


def test_farmers_pawnee_july_recharge_uses_negative_volumes():
    """July FP recharge: one monthly pulse, SDF, accretion (negative)."""
    rows = lag_sdf({"2025-07-01": -100.0}, sdf=265, total_months=6)
    assert rows[0][0] == "2025-07-01"
    assert all(v < 0.0 for _, v in rows)
    assert rows[0][1] == pytest.approx(-0.944405, abs=5e-6)


def test_glover_infinite_january_100_af_matches_volume_formula():
    t = gpd_per_ft_to_ft2_per_day(261_800.0)
    rows = lag_glover_infinite(
        {"2025-01-01": 100.0},
        distance_ft=4000.0,
        specific_yield=0.2,
        transmissivity_ft2_per_day=t,
        total_months=6,
    )
    expected = [9.231024, 20.792231, 13.118033, 7.592177, 5.342934, 3.804862]
    assert [v for _, v in rows] == pytest.approx(expected, abs=5e-6)


def test_glover_matches_sdf_when_sdf_equals_a2s_over_t():
    t = gpd_per_ft_to_ft2_per_day(261_800.0)
    sdf = sdf_from_glover(4000.0, 0.2, t)
    assert sdf == pytest.approx(91.437, abs=0.01)
    glover = lag_glover_infinite({"2025-01-01": 100.0}, 4000.0, 0.2, t, 6)
    jenkins = lag_sdf({"2025-01-01": 100.0}, round(sdf), 6)
    for (_, g), (_, j) in zip(glover, jenkins):
        assert abs(g - j) / max(j, 1e-6) < 0.02


def test_urf_month_gap_is_a_zero_lag_not_a_collapsed_lag():
    # Month 2 omitted: must not treat month 3 as “second listed row”.
    lagged = urf_lagging({"2025-01-01": 100.0}, [(1, 1, 0.5), (3, 1, 0.5)])
    combined = combined_urf_results(lagged)
    assert combined == [("2025-01-01", 50.0), ("2025-03-01", 50.0)]


def test_urf_duplicate_month_rows_sum_coefficients():
    urf = [(1, 1, 0.25), (1, 1, 0.25), (2, 1, 0.50)]
    combined = combined_urf_results(urf_lagging({"2025-01-01": 80.0}, urf))
    assert combined == [("2025-01-01", 40.0), ("2025-02-01", 40.0)]


def test_urf_mid_month_usage_date_normalizes_to_month_start():
    usage = {date(2025, 1, 15): 10.0, date(2025, 1, 1): 5.0}
    combined = combined_urf_results(urf_lagging(usage, [{"month": 1, "reach": 1, "urf_val": 1.0}]))
    assert combined == [("2025-01-01", 15.0)]


def test_urf_year_month_string_and_tuple_keys():
    a = combined_urf_results(urf_lagging({"2025-01": 10.0}, [(1, 1, 1.0)]))
    b = combined_urf_results(urf_lagging({(2025, 1): 10.0}, [(1, 1, 1.0)]))
    assert a == b == [("2025-01-01", 10.0)]


def test_daily_urf_gap_is_a_zero_lag():
    urf = [{"day": 1, "reach": 1, "urf_val": 0.4}, {"day": 3, "reach": 1, "urf_val": 0.6}]
    combined = combined_urf_results_dated(urf_lagging_daily({"2025-01-01": 50.0}, urf))
    assert combined == [("2025-01-01", 20.0), ("2025-01-03", 30.0)]


def test_daily_first_coefficient_is_stress_day():
    combined = combined_urf_results_dated(
        urf_lagging_daily({date(2025, 6, 15): 4.0}, [(1, 1, 1.0)])
    )
    assert combined == [("2025-06-15", 4.0)]


def test_daily_unit_pulse_conserves_mass():
    urf = [(1, 1, 0.50), (2, 1, 0.30), (3, 1, 0.20)]
    assert daily_urf_mass_by_reach(urf)[1] == pytest.approx(1.0)
    combined = combined_urf_results_dated(urf_lagging_daily({"2025-01-10": 10.0}, urf))
    assert combined == [
        ("2025-01-10", 5.0),
        ("2025-01-11", 3.0),
        ("2025-01-12", 2.0),
    ]


def test_urf_mass_by_reach_and_dict_rows():
    urf = [
        {"month": 1, "reach": 1, "urf_val": 0.20},
        {"month": 2, "reach": 1, "urf_val": 0.50},
        {"month": 3, "reach": 1, "urf_val": 0.30},
    ]
    assert urf_mass_by_reach(urf)[1] == pytest.approx(1.0)


def test_rust_name_aliases_are_the_same_callables():
    from stream_depletion import calculate_streamflow_depletion_sdf

    assert calculate_streamflow_depletion_sdf is lag_sdf
