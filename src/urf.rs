//! Discrete unit-response (URF) convolution for monthly **or daily** stream
//! effects.
//!
//! A URF is a tabulated unit hydrograph: the fraction of a one-step stress
//! that appears as stream depletion (or accretion) in each later step.
//! Convolution is the discrete-kernel / algebraic-technological-function
//! construction of Maddock (1972, 1974) and Morel-Seytoux and Daly (1975),
//! not an analytical Glover or SDF evaluation. See `docs/URF.md`.
//!
//! * Monthly: [`UrfValue::month`] is **1-based lag** (month 1 = stress month).
//!   Input is monthly volumes. Use [`urf_lagging`].
//! * Daily: [`DailyUrfValue::day`] is **1-based lag** (day 1 = stress day).
//!   Input is daily volumes. Use [`urf_lagging_daily`].
//!
//! A native daily kernel is not a monthly kernel interpolated onto days.
//! [`disaggregate_monthly_urf_uniform`] is an explicit, opt-in conversion
//! that splits each monthly coefficient evenly across the calendar days of
//! that response month (measured from the first of the stress month).
//! Missing periods between 1 and the last listed period are implicit zeros.
//! Duplicate `(reach, period)` rows are summed.

use crate::utils::add_months;
use chrono::{Datelike, Duration, Months, NaiveDate};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::ops::Add;

pub type LaggedUrfByDate = HashMap<NaiveDate, f64>;
pub type LaggedUrfResult = HashMap<i32, LaggedUrfByDate>;

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct UrfValue {
    /// 1-based response month (1 = same calendar month as the stress).
    pub month: i32,
    pub reach: i32,
    pub urf_val: f64,
}

impl UrfValue {
    pub fn new(month: i32, reach: i32, urf_val: f64) -> Self {
        UrfValue {
            month,
            reach,
            urf_val,
        }
    }

    /// Build a coefficient from a 0-based lag (0 = stress month), the indexing
    /// used by IDS AWAS `m_URF_MonthlyData`.
    pub fn from_lag(lag: u32, reach: i32, urf_val: f64) -> Self {
        UrfValue::new(lag as i32 + 1, reach, urf_val)
    }
}

/// One term of a **daily** unit-response series.
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct DailyUrfValue {
    /// 1-based response day (1 = same calendar day as the stress).
    pub day: i32,
    pub reach: i32,
    pub urf_val: f64,
}

impl DailyUrfValue {
    pub fn new(day: i32, reach: i32, urf_val: f64) -> Self {
        DailyUrfValue {
            day,
            reach,
            urf_val,
        }
    }

    /// 0-based lag (0 = stress day), matching IDS AWAS `m_URF_DailyData`.
    pub fn from_lag(lag: u32, reach: i32, urf_val: f64) -> Self {
        DailyUrfValue::new(lag as i32 + 1, reach, urf_val)
    }
}

/// First day of the calendar month containing `date`.
pub fn month_start(date: NaiveDate) -> NaiveDate {
    NaiveDate::from_ymd_opt(date.year(), date.month(), 1).unwrap()
}

/// Sum of URF coefficients on each reach.
///
/// For a complete unit hydrograph whose only source/sink is the stream,
/// each reach-sum (or the sum across reaches, if they partition one stream)
/// equals 1. Truncated model tails and other sinks make the sum smaller.
/// This crate does **not** rescale the series to 1.
pub fn urf_mass_by_reach(urf: &[UrfValue]) -> HashMap<i32, f64> {
    mass_by_reach(urf.iter().map(|r| (r.month, r.reach, r.urf_val)))
}

/// Sum of **daily** URF coefficients on each reach. Same unit-hydrograph
/// meaning as [`urf_mass_by_reach`].
pub fn daily_urf_mass_by_reach(urf: &[DailyUrfValue]) -> HashMap<i32, f64> {
    mass_by_reach(urf.iter().map(|r| (r.day, r.reach, r.urf_val)))
}

fn mass_by_reach<I>(rows: I) -> HashMap<i32, f64>
where
    I: IntoIterator<Item = (i32, i32, f64)>,
{
    let mut masses = HashMap::new();
    for (period, reach, val) in rows {
        if period < 1 {
            continue;
        }
        *masses.entry(reach).or_insert(0.0) += val;
    }
    masses
}

/// Dense 1-based series per reach: index 0 is period 1 (stress period).
fn dense_series_by_reach<I>(rows: I) -> BTreeMap<i32, Vec<f64>>
where
    I: IntoIterator<Item = (i32, i32, f64)>,
{
    let mut by_reach: BTreeMap<i32, BTreeMap<i32, f64>> = BTreeMap::new();
    for (period, reach, val) in rows {
        if period < 1 {
            continue;
        }
        *by_reach
            .entry(reach)
            .or_default()
            .entry(period)
            .or_insert(0.0) += val;
    }

    let mut series = BTreeMap::new();
    for (reach, periods) in by_reach {
        let max_period = *periods.keys().max().unwrap_or(&1);
        let mut dense = vec![0.0; max_period as usize];
        for (period, val) in periods {
            dense[(period - 1) as usize] = val;
        }
        series.insert(reach, dense);
    }
    series
}

/// Discrete convolution of monthly stress with a per-reach URF.
///
/// ```text
/// effect_r(t) = Σ_τ  usage(τ) · URF_r(t − τ + 1)
/// ```
///
/// with `URF_r(1)` applied in the stress month. Stress dates are normalized
/// to the first of the month. Positive usage is pumping (depletion); negative
/// usage is recharge (accretion). Linear superposition; zeros are skipped.
pub fn urf_lagging(usage: &HashMap<NaiveDate, f64>, urf: Vec<UrfValue>) -> LaggedUrfResult {
    let series = dense_series_by_reach(urf.iter().map(|r| (r.month, r.reach, r.urf_val)));

    let mut month_usage: BTreeMap<NaiveDate, f64> = BTreeMap::new();
    for (date, volume) in usage {
        *month_usage.entry(month_start(*date)).or_insert(0.0) += *volume;
    }

    convolve(&month_usage, &series, |start, lag| {
        start.add(Months::new(lag as u32))
    })
}

/// Discrete convolution of **daily** stress with a per-reach daily URF.
///
/// ```text
/// effect_r(t) = Σ_τ  usage(τ) · URF_r(t − τ + 1)
/// ```
///
/// with `URF_r(1)` applied on the stress day. Dates are used as given (not
/// rolled to month start). Positive usage is pumping; negative is recharge.
pub fn urf_lagging_daily(
    usage: &HashMap<NaiveDate, f64>,
    urf: Vec<DailyUrfValue>,
) -> LaggedUrfResult {
    let series = dense_series_by_reach(urf.iter().map(|r| (r.day, r.reach, r.urf_val)));
    let mut day_usage: BTreeMap<NaiveDate, f64> = BTreeMap::new();
    for (date, volume) in usage {
        *day_usage.entry(*date).or_insert(0.0) += *volume;
    }
    convolve(&day_usage, &series, |start, lag| {
        start + Duration::days(lag as i64)
    })
}

fn convolve(
    usage: &BTreeMap<NaiveDate, f64>,
    series: &BTreeMap<i32, Vec<f64>>,
    shift: impl Fn(NaiveDate, usize) -> NaiveDate,
) -> LaggedUrfResult {
    let mut lagged_result = HashMap::new();
    for (reach, reach_urf) in series {
        let mut reach_lagged = HashMap::new();
        for (usage_date, volume) in usage {
            if *volume == 0.0 {
                continue;
            }
            for (i, coeff) in reach_urf.iter().enumerate() {
                if *coeff == 0.0 {
                    continue;
                }
                *reach_lagged.entry(shift(*usage_date, i)).or_insert(0.0) += volume * coeff;
            }
        }
        lagged_result.insert(*reach, reach_lagged);
    }
    lagged_result
}

/// Opt-in conversion of a monthly URF into a rolling daily series.
///
/// Each monthly coefficient is split **evenly** across the calendar days of
/// that response month, counting from the first of `stress_month`. This is
/// **not** a native daily kernel: it is a uniform-in-month interpolation
/// that is only equivalent to the monthly convolution for a pulse on the
/// first of `stress_month` (see tests). Prefer a model-generated daily URF
/// when the plan has one.
pub fn disaggregate_monthly_urf_uniform(
    monthly: &[UrfValue],
    stress_month: NaiveDate,
) -> Vec<DailyUrfValue> {
    let origin = month_start(stress_month);
    let series = dense_series_by_reach(monthly.iter().map(|r| (r.month, r.reach, r.urf_val)));
    let mut daily = Vec::new();
    for (reach, coeffs) in series {
        let mut day_index = 1i32;
        for (lag, coeff) in coeffs.iter().enumerate() {
            let period = match add_months(origin, lag as i32) {
                Some(d) => d,
                None => break,
            };
            let n = period.num_days_in_month() as i32;
            let each = if n > 0 { coeff / n as f64 } else { 0.0 };
            for _ in 0..n {
                if each != 0.0 {
                    daily.push(DailyUrfValue::new(day_index, reach, each));
                }
                day_index += 1;
            }
        }
    }
    daily
}

/// Sum reach-specific lagged series onto a single monthly timeline.
///
/// Every month that received any reach contribution is returned, including
/// zeros (explicit gap months that were written) and negative accretions.
/// Months that no reach touched are omitted (they are true absences, not
/// computed zeros).
pub fn combined_urf_results(values: LaggedUrfResult) -> Vec<(NaiveDate, f64)> {
    let mut date_sums: HashMap<NaiveDate, f64> = HashMap::new();

    for inner_map in values.values() {
        for (date, value) in inner_map {
            *date_sums.entry(month_start(*date)).or_insert(0.0) += value;
        }
    }

    let mut result: Vec<(NaiveDate, f64)> = date_sums.into_iter().collect();
    result.sort_by(|a, b| a.0.cmp(&b.0));
    result
}

/// Sum reach-specific lagged series **without** rolling dates to month start.
///
/// Use this after [`urf_lagging_daily`] when daily resolution is required.
/// [`combined_urf_results`] can still be applied to a daily lag map if you
/// want calendar-month totals.
pub fn combined_urf_results_dated(values: LaggedUrfResult) -> Vec<(NaiveDate, f64)> {
    let mut date_sums: HashMap<NaiveDate, f64> = HashMap::new();
    for inner_map in values.values() {
        for (date, value) in inner_map {
            *date_sums.entry(*date).or_insert(0.0) += value;
        }
    }
    let mut result: Vec<(NaiveDate, f64)> = date_sums.into_iter().collect();
    result.sort_by(|a, b| a.0.cmp(&b.0));
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(year: i32, month: u32, day: u32) -> NaiveDate {
        NaiveDate::from_ymd_opt(year, month, day).unwrap()
    }

    fn month(year: i32, month: u32) -> NaiveDate {
        d(year, month, 1)
    }

    #[test]
    fn test_urf_lagging() {
        let urf = vec![
            UrfValue {
                month: 1,
                reach: 1,
                urf_val: 0.6,
            },
            UrfValue {
                month: 1,
                reach: 2,
                urf_val: 0.1,
            },
            UrfValue {
                month: 2,
                reach: 1,
                urf_val: 0.3,
            },
        ];
        let mut usage = HashMap::new();
        usage.insert(month(2024, 7), 100.0);
        usage.insert(month(2024, 8), 100.0);

        let mut expected_lagged = HashMap::new();
        let mut reach1 = HashMap::new();
        reach1.insert(month(2024, 7), 60.0);
        reach1.insert(month(2024, 8), 90.0);
        reach1.insert(month(2024, 9), 30.0);

        let mut reach2 = HashMap::new();
        reach2.insert(month(2024, 7), 10.0);
        reach2.insert(month(2024, 8), 10.0);
        expected_lagged.insert(1, reach1);
        expected_lagged.insert(2, reach2);

        let result = urf_lagging(&usage, urf);
        assert_eq!(result, expected_lagged);
    }

    #[test]
    fn test_urf_lagging_skip_usage_month() {
        let urf = vec![
            UrfValue {
                month: 1,
                reach: 1,
                urf_val: 0.4,
            },
            UrfValue {
                month: 1,
                reach: 2,
                urf_val: 0.2,
            },
            UrfValue {
                month: 2,
                reach: 1,
                urf_val: 0.2,
            },
            UrfValue {
                month: 2,
                reach: 2,
                urf_val: 0.1,
            },
            UrfValue {
                month: 3,
                reach: 1,
                urf_val: 0.1,
            },
        ];
        let mut usage = HashMap::new();
        usage.insert(month(2024, 5), 100.0);
        usage.insert(month(2024, 7), 100.0);
        usage.insert(month(2024, 8), 100.0);

        let mut expected_lagged = HashMap::new();
        let mut reach1 = HashMap::new();
        reach1.insert(month(2024, 5), 40.0);
        reach1.insert(month(2024, 6), 20.0);
        reach1.insert(month(2024, 7), 50.0);
        reach1.insert(month(2024, 8), 60.0);
        reach1.insert(month(2024, 9), 30.0);
        reach1.insert(month(2024, 10), 10.0);
        let mut reach2 = HashMap::new();
        reach2.insert(month(2024, 5), 20.0);
        reach2.insert(month(2024, 6), 10.0);
        reach2.insert(month(2024, 7), 20.0);
        reach2.insert(month(2024, 8), 30.0);
        reach2.insert(month(2024, 9), 10.0);
        expected_lagged.insert(1, reach1);
        expected_lagged.insert(2, reach2);

        let result = urf_lagging(&usage, urf);
        assert_eq!(result, expected_lagged);
    }

    #[test]
    fn unit_pulse_conserves_mass_when_urf_sums_to_one() {
        // Maddock / Morel-Seytoux discrete kernel: a unit stress convolved
        // with a complete unit hydrograph recovers the stressed volume.
        let urf = vec![
            UrfValue::new(1, 1, 0.20),
            UrfValue::new(2, 1, 0.50),
            UrfValue::new(3, 1, 0.30),
        ];
        assert!((urf_mass_by_reach(&urf)[&1] - 1.0).abs() < 1e-15);

        let mut usage = HashMap::new();
        usage.insert(month(2025, 1), 100.0);
        let lagged = urf_lagging(&usage, urf);
        let combined = combined_urf_results(lagged);
        assert_eq!(
            combined,
            vec![
                (month(2025, 1), 20.0),
                (month(2025, 2), 50.0),
                (month(2025, 3), 30.0),
            ]
        );
        let sum: f64 = combined.iter().map(|(_, v)| v).sum();
        assert!((sum - 100.0).abs() < 1e-12);
    }

    #[test]
    fn first_coefficient_is_stress_month_not_the_next() {
        let urf = vec![UrfValue::from_lag(0, 1, 1.0)];
        let mut usage = HashMap::new();
        usage.insert(month(2025, 6), 40.0);
        let combined = combined_urf_results(urf_lagging(&usage, urf));
        assert_eq!(combined, vec![(month(2025, 6), 40.0)]);
    }

    #[test]
    fn gap_in_month_index_is_a_zero_lag_not_a_collapsed_lag() {
        // Month 2 omitted: must not treat month 3 as “second listed row”.
        let urf = vec![UrfValue::new(1, 1, 0.5), UrfValue::new(3, 1, 0.5)];
        let mut usage = HashMap::new();
        usage.insert(month(2025, 1), 100.0);
        let combined = combined_urf_results(urf_lagging(&usage, urf));
        assert_eq!(
            combined,
            vec![(month(2025, 1), 50.0), (month(2025, 3), 50.0)]
        );
    }

    #[test]
    fn duplicate_month_rows_sum_coefficients() {
        let urf = vec![
            UrfValue::new(1, 1, 0.25),
            UrfValue::new(1, 1, 0.25),
            UrfValue::new(2, 1, 0.50),
        ];
        let mut usage = HashMap::new();
        usage.insert(month(2025, 1), 80.0);
        let combined = combined_urf_results(urf_lagging(&usage, urf));
        assert_eq!(
            combined,
            vec![(month(2025, 1), 40.0), (month(2025, 2), 40.0)]
        );
    }

    #[test]
    fn recharge_is_negative_mirror_of_pumping() {
        let urf = vec![UrfValue::new(1, 1, 0.4), UrfValue::new(2, 1, 0.6)];
        let mut pump = HashMap::new();
        pump.insert(month(2025, 1), 100.0);
        let mut rech = HashMap::new();
        rech.insert(month(2025, 1), -100.0);
        let p = combined_urf_results(urf_lagging(&pump, urf.clone()));
        let r = combined_urf_results(urf_lagging(&rech, urf));
        assert_eq!(p.len(), r.len());
        for (a, b) in p.iter().zip(r.iter()) {
            assert_eq!(a.0, b.0);
            assert!((a.1 + b.1).abs() < 1e-12);
            assert!(b.1 < 0.0);
        }
    }

    #[test]
    fn mid_month_usage_date_normalizes_to_month_start() {
        let urf = vec![UrfValue::new(1, 1, 1.0)];
        let mut usage = HashMap::new();
        usage.insert(d(2025, 1, 15), 10.0);
        usage.insert(d(2025, 1, 1), 5.0);
        let combined = combined_urf_results(urf_lagging(&usage, urf));
        assert_eq!(combined, vec![(month(2025, 1), 15.0)]);
    }

    #[test]
    fn superposition_of_two_stress_months() {
        let urf = vec![UrfValue::new(1, 1, 0.25), UrfValue::new(2, 1, 0.75)];
        let mut usage = HashMap::new();
        usage.insert(month(2025, 1), 100.0);
        usage.insert(month(2025, 2), 100.0);
        let combined = combined_urf_results(urf_lagging(&usage, urf));
        // Jan: 25; Feb: 75 (from Jan) + 25 (from Feb) = 100; Mar: 75
        assert_eq!(
            combined,
            vec![
                (month(2025, 1), 25.0),
                (month(2025, 2), 100.0),
                (month(2025, 3), 75.0),
            ]
        );
        let sum: f64 = combined.iter().map(|(_, v)| v).sum();
        assert!((sum - 200.0).abs() < 1e-12);
    }

    #[test]
    fn combined_keeps_accretion_and_does_not_invent_a_cutoff() {
        let urf = vec![UrfValue::new(1, 1, 1.0)];
        let mut usage = HashMap::new();
        usage.insert(month(2025, 1), -0.0004);
        let combined = combined_urf_results(urf_lagging(&usage, urf));
        assert_eq!(combined, vec![(month(2025, 1), -0.0004)]);
    }

    #[test]
    fn month_less_than_one_is_ignored() {
        let urf = vec![UrfValue::new(0, 1, 0.9), UrfValue::new(1, 1, 1.0)];
        let mut usage = HashMap::new();
        usage.insert(month(2025, 1), 10.0);
        let combined = combined_urf_results(urf_lagging(&usage, urf));
        assert_eq!(combined, vec![(month(2025, 1), 10.0)]);
    }

    #[test]
    fn glover_rate_fraction_is_not_a_volume_urf() {
        // Mixing bug: q/Q (Glover rate) is not the monthly volume URF.
        // Incremental volume fractions for a 31-day pulse at sdf=265 d
        // are the Jenkins v increments / 100 AF, not erfc(u) itself.
        use crate::kernel::{
            complementary_error, pulse_cumulative_volume, sdf_argument, ResponseKernel,
        };
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let q = 100.0 / 31.0;
        let v_jan = pulse_cumulative_volume(q, 31.0, 31.0, &kernel);
        let v_feb = pulse_cumulative_volume(q, 31.0, 59.0, &kernel);
        let rate_end_jan = complementary_error(sdf_argument(265.0, 31.0));
        // Rate fraction at end of January is not the January volume fraction.
        assert!((v_jan / 100.0 - 0.009444).abs() < 1e-5);
        assert!((rate_end_jan - 0.038695).abs() < 1e-5);
        assert!((v_jan / 100.0 - rate_end_jan).abs() > 0.02);
        assert!(v_feb > v_jan);
    }

    #[test]
    fn daily_unit_pulse_conserves_mass() {
        let urf = vec![
            DailyUrfValue::new(1, 1, 0.50),
            DailyUrfValue::new(2, 1, 0.30),
            DailyUrfValue::new(3, 1, 0.20),
        ];
        assert!((daily_urf_mass_by_reach(&urf)[&1] - 1.0).abs() < 1e-15);
        let mut usage = HashMap::new();
        usage.insert(d(2025, 1, 10), 10.0);
        let combined = combined_urf_results_dated(urf_lagging_daily(&usage, urf));
        assert_eq!(
            combined,
            vec![
                (d(2025, 1, 10), 5.0),
                (d(2025, 1, 11), 3.0),
                (d(2025, 1, 12), 2.0),
            ]
        );
        let sum: f64 = combined.iter().map(|(_, v)| v).sum();
        assert!((sum - 10.0).abs() < 1e-12);
    }

    #[test]
    fn daily_first_coefficient_is_stress_day_not_the_next() {
        let urf = vec![DailyUrfValue::from_lag(0, 1, 1.0)];
        let mut usage = HashMap::new();
        usage.insert(d(2025, 6, 15), 4.0);
        let combined = combined_urf_results_dated(urf_lagging_daily(&usage, urf));
        assert_eq!(combined, vec![(d(2025, 6, 15), 4.0)]);
    }

    #[test]
    fn daily_gap_is_a_zero_lag() {
        let urf = vec![DailyUrfValue::new(1, 1, 0.4), DailyUrfValue::new(3, 1, 0.6)];
        let mut usage = HashMap::new();
        usage.insert(d(2025, 1, 1), 50.0);
        let combined = combined_urf_results_dated(urf_lagging_daily(&usage, urf));
        assert_eq!(combined, vec![(d(2025, 1, 1), 20.0), (d(2025, 1, 3), 30.0)]);
    }

    #[test]
    fn daily_recharge_is_negative_mirror() {
        let urf = vec![
            DailyUrfValue::new(1, 1, 0.25),
            DailyUrfValue::new(2, 1, 0.75),
        ];
        let mut pump = HashMap::new();
        pump.insert(d(2025, 3, 1), 8.0);
        let mut rech = HashMap::new();
        rech.insert(d(2025, 3, 1), -8.0);
        let p = combined_urf_results_dated(urf_lagging_daily(&pump, urf.clone()));
        let r = combined_urf_results_dated(urf_lagging_daily(&rech, urf));
        for (a, b) in p.iter().zip(r.iter()) {
            assert_eq!(a.0, b.0);
            assert!((a.1 + b.1).abs() < 1e-12);
            assert!(b.1 < 0.0);
        }
    }

    #[test]
    fn daily_superposition_of_two_stress_days() {
        let urf = vec![DailyUrfValue::new(1, 1, 0.5), DailyUrfValue::new(2, 1, 0.5)];
        let mut usage = HashMap::new();
        usage.insert(d(2025, 1, 1), 10.0);
        usage.insert(d(2025, 1, 2), 10.0);
        let combined = combined_urf_results_dated(urf_lagging_daily(&usage, urf));
        assert_eq!(
            combined,
            vec![
                (d(2025, 1, 1), 5.0),
                (d(2025, 1, 2), 10.0),
                (d(2025, 1, 3), 5.0),
            ]
        );
    }

    #[test]
    fn jan1_daily_pulse_with_disaggregated_monthly_urf_matches_monthly() {
        // Native daily kernel is rolling from the stress *day*. A monthly URF
        // uniformly split across the calendar days of each response month,
        // applied to a pulse on the 1st, is the discrete equivalent of the
        // monthly convolution.
        let monthly = vec![
            UrfValue::new(1, 1, 0.20),
            UrfValue::new(2, 1, 0.50),
            UrfValue::new(3, 1, 0.30),
        ];
        let mut monthly_usage = HashMap::new();
        monthly_usage.insert(month(2025, 1), 100.0);
        let monthly_out = combined_urf_results(urf_lagging(&monthly_usage, monthly.clone()));

        let daily_urf = disaggregate_monthly_urf_uniform(&monthly, month(2025, 1));
        assert!((daily_urf_mass_by_reach(&daily_urf)[&1] - 1.0).abs() < 1e-12);
        let mut daily_usage = HashMap::new();
        daily_usage.insert(d(2025, 1, 1), 100.0);
        let daily_out = combined_urf_results(urf_lagging_daily(&daily_usage, daily_urf));

        assert_eq!(monthly_out.len(), daily_out.len());
        for (m, dy) in monthly_out.iter().zip(daily_out.iter()) {
            assert_eq!(m.0, dy.0);
            assert!((m.1 - dy.1).abs() < 1e-12);
        }
    }

    #[test]
    fn uniform_january_daily_stress_with_rolling_daily_urf_matches_mass_only() {
        // Spreading a monthly volume across January days and applying a
        // *rolling* daily kernel (even one disaggregated from the monthly
        // URF) is not the same monthly-shape as one January pulse. Totals
        // still match when the kernel sums to 1.
        let monthly = vec![
            UrfValue::new(1, 1, 0.20),
            UrfValue::new(2, 1, 0.50),
            UrfValue::new(3, 1, 0.30),
        ];
        let mut monthly_usage = HashMap::new();
        monthly_usage.insert(month(2025, 1), 100.0);
        let monthly_out = combined_urf_results(urf_lagging(&monthly_usage, monthly.clone()));

        let daily_urf = disaggregate_monthly_urf_uniform(&monthly, month(2025, 1));
        let mut daily_usage = HashMap::new();
        let n = month(2025, 1).num_days_in_month() as f64;
        for day in 1..=month(2025, 1).num_days_in_month() {
            daily_usage.insert(d(2025, 1, day as u32), 100.0 / n);
        }
        let daily_out = combined_urf_results(urf_lagging_daily(&daily_usage, daily_urf));
        let m_sum: f64 = monthly_out.iter().map(|(_, v)| *v).sum();
        let d_sum: f64 = daily_out.iter().map(|(_, v)| *v).sum();
        assert!((m_sum - 100.0).abs() < 1e-10);
        assert!((d_sum - 100.0).abs() < 1e-10);
        // Shape differs: late-January daily pulses leak into April.
        assert!(daily_out.len() > monthly_out.len());
        let jan_m = monthly_out
            .iter()
            .find(|(dt, _)| *dt == month(2025, 1))
            .unwrap()
            .1;
        let jan_d = daily_out
            .iter()
            .find(|(dt, _)| *dt == month(2025, 1))
            .unwrap()
            .1;
        assert!(jan_d < jan_m - 1.0);
    }
}
