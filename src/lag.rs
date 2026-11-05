//! Jenkins (1968) superposition of rectangular monthly pulses.
//!
//! Each month of stress is a rectangular pulse. Linear superposition applies
//! to both pumping (+) and recharge (−). Monthly output is the increment of
//! cumulative stream volume across that month (`v(end) − v(start)`), in the
//! same units as the input volumes.
//!
//! This is the volume form of Jenkins’ residual-effects construction (USGS
//! TWRI 4-D1), not a finite-difference of the *rate* curve.
//!
//! The length of a month is a modelling choice, set by [`MonthLength`]. Jenkins
//! only needs `t` in consistent units, so both options are the same published
//! method applied on a different time axis.

use crate::kernel::{pulse_cumulative_volume, ParameterError, ResponseKernel};
use crate::utils::add_months;
use chrono::{Datelike, NaiveDate};
use std::collections::HashMap;

/// Length of one month under [`MonthLength::Average`], in days: `365/12`
/// rounded to the value IDS AWAS uses (`SiteData::GetFact`). The literal is kept
/// rather than `365.0 / 12.0` so results stay comparable with AWAS and its
/// `q / 7.43952` gallons-per-minute to acre-feet-per-month factor.
pub const AVERAGE_DAYS_PER_MONTH: f64 = 30.41667;

/// How long a month is when building pulses and reporting windows.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MonthLength {
    /// Actual calendar lengths (January = 31, February = 28 or 29, …).
    #[default]
    Calendar,
    /// Every month lasts [`AVERAGE_DAYS_PER_MONTH`] days, so the response to a
    /// monthly pulse does not depend on when in the year it falls. Matches IDS
    /// AWAS with “use average days per month” on. Output dates are still the
    /// calendar month starts.
    Average,
}

impl MonthLength {
    fn days_in(self, month: NaiveDate) -> f64 {
        match self {
            MonthLength::Calendar => month.num_days_in_month() as f64,
            MonthLength::Average => AVERAGE_DAYS_PER_MONTH,
        }
    }

    /// Model days from the start of `origin` to the start of `month`.
    fn elapsed_days(self, origin: NaiveDate, month: NaiveDate) -> f64 {
        match self {
            MonthLength::Calendar => (month - origin).num_days() as f64,
            MonthLength::Average => {
                let months = (month.year() - origin.year()) * 12 + month.month() as i32
                    - origin.month() as i32;
                months as f64 * AVERAGE_DAYS_PER_MONTH
            }
        }
    }
}

/// Lag monthly pumping (+) or recharge (−) volumes through `kernel`, using
/// calendar month lengths.
///
/// * `stress_volumes_monthly` — acre-feet (or any consistent volume) keyed by
///   a date in the stress month. The day-of-month is ignored; the pulse covers
///   the calendar month.
/// * `total_months` — number of calendar months to report, starting at the
///   earliest stress month.
///
/// `_days_per_month` is accepted for source compatibility and has no effect.
/// Callers pass values such as `30.42` and rely on calendar results.
///
/// This function does **not** check its input, because its return type cannot
/// report an error. Prefer [`lag_monthly_volumes_with_month_length`], which
/// checks it and also offers [`MonthLength::Average`].
pub fn lag_monthly_volumes(
    stress_volumes_monthly: &HashMap<NaiveDate, f64>,
    kernel: &ResponseKernel,
    _days_per_month: f64,
    total_months: usize,
) -> Vec<(NaiveDate, f64)> {
    lag_unchecked(
        stress_volumes_monthly,
        kernel,
        MonthLength::Calendar,
        total_months,
    )
}

/// Lag monthly pumping (+) or recharge (−) volumes through `kernel`, with an
/// explicit [`MonthLength`]. Input volumes are per month under either option;
/// only the time axis changes. See [`lag_monthly_volumes`] for the other
/// arguments.
///
/// Returns an error, and computes nothing, if the kernel’s parameters are
/// invalid (see [`ResponseKernel::validate`]) or a stress volume is not finite.
pub fn lag_monthly_volumes_with_month_length(
    stress_volumes_monthly: &HashMap<NaiveDate, f64>,
    kernel: &ResponseKernel,
    month_length: MonthLength,
    total_months: usize,
) -> Result<Vec<(NaiveDate, f64)>, ParameterError> {
    kernel.validate()?;
    let bad_stress = stress_volumes_monthly
        .iter()
        .filter(|(_, volume)| !volume.is_finite())
        .min_by_key(|(month, _)| **month);
    if let Some((month, value)) = bad_stress {
        return Err(ParameterError::InvalidStress {
            month: *month,
            value: *value,
        });
    }
    Ok(lag_unchecked(
        stress_volumes_monthly,
        kernel,
        month_length,
        total_months,
    ))
}

fn lag_unchecked(
    stress_volumes_monthly: &HashMap<NaiveDate, f64>,
    kernel: &ResponseKernel,
    month_length: MonthLength,
    total_months: usize,
) -> Vec<(NaiveDate, f64)> {
    if stress_volumes_monthly.is_empty() || total_months == 0 {
        return Vec::new();
    }

    let origin = *stress_volumes_monthly.keys().min().unwrap();
    let start_month = NaiveDate::from_ymd_opt(origin.year(), origin.month(), 1).unwrap();

    // Collapse multiple entries in the same calendar month.
    let mut month_stress: HashMap<NaiveDate, f64> = HashMap::new();
    for (date, volume) in stress_volumes_monthly {
        let month = NaiveDate::from_ymd_opt(date.year(), date.month(), 1).unwrap();
        *month_stress.entry(month).or_insert(0.0) += *volume;
    }

    let mut pulses: Vec<Pulse> = month_stress
        .iter()
        .filter(|(_, v)| **v != 0.0)
        .map(|(month, volume)| {
            let days = month_length.days_in(*month);
            Pulse {
                start_day: month_length.elapsed_days(start_month, *month),
                duration_days: days,
                rate: *volume / days,
            }
        })
        .collect();
    // HashMap order varies by run; a fixed sum order keeps results reproducible.
    pulses.sort_by(|a, b| a.start_day.total_cmp(&b.start_day));

    let mut results = Vec::with_capacity(total_months);
    for month_index in 0..total_months {
        let month_start = match add_months(start_month, month_index as i32) {
            Some(d) => d,
            None => break,
        };
        let month_end = match add_months(month_start, 1) {
            Some(d) => d,
            None => break,
        };
        let window_start = month_length.elapsed_days(start_month, month_start);
        let window_end = month_length.elapsed_days(start_month, month_end);
        let mut increment = 0.0;
        for pulse in &pulses {
            increment += pulse_volume_between(pulse, window_start, window_end, kernel);
        }
        results.push((month_start, increment));
    }
    results
}

struct Pulse {
    start_day: f64,
    duration_days: f64,
    rate: f64,
}

fn pulse_volume_between(
    pulse: &Pulse,
    window_start: f64,
    window_end: f64,
    kernel: &ResponseKernel,
) -> f64 {
    let t_end = window_end - pulse.start_day;
    let t_start = window_start - pulse.start_day;
    pulse_cumulative_volume(pulse.rate, pulse.duration_days, t_end, kernel)
        - pulse_cumulative_volume(pulse.rate, pulse.duration_days, t_start, kernel)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::ResponseKernel;

    fn month(year: i32, month: u32) -> NaiveDate {
        NaiveDate::from_ymd_opt(year, month, 1).unwrap()
    }

    #[test]
    fn empty_input_is_empty() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let empty = HashMap::new();
        assert!(lag_monthly_volumes(&empty, &kernel, 30.42, 12).is_empty());
    }

    #[test]
    fn reports_every_month_including_zeros_and_negatives() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let mut stress = HashMap::new();
        stress.insert(month(2025, 1), -100.0); // recharge
        let out = lag_monthly_volumes(&stress, &kernel, 30.42, 6);
        assert_eq!(out.len(), 6);
        assert!(out.iter().all(|(_, v)| *v < 0.0));
        assert_eq!(out[0].0, month(2025, 1));
        assert_eq!(out[5].0, month(2025, 6));
    }

    #[test]
    fn two_equal_months_superpose() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 50.0 };
        let mut ten = HashMap::new();
        ten.insert(month(2025, 1), 10.0);
        let mut five = HashMap::new();
        five.insert(month(2025, 1), 5.0);

        let a = lag_monthly_volumes(&ten, &kernel, 30.0, 4);
        let half = lag_monthly_volumes(&five, &kernel, 30.0, 4);
        for (left, right) in a.iter().zip(half.iter()) {
            assert_eq!(left.0, right.0);
            assert!((left.1 - 2.0 * right.1).abs() < 1e-12);
        }

        // January + February pulses: February total = Jan residual + a February-only
        // pulse. Durations differ (31 vs 28 days), so the February pulse is computed
        // from a start in February, not by shifting January.
        let mut jan_feb = HashMap::new();
        jan_feb.insert(month(2025, 1), 10.0);
        jan_feb.insert(month(2025, 2), 10.0);
        let mut feb_only = HashMap::new();
        feb_only.insert(month(2025, 2), 10.0);
        let ab = lag_monthly_volumes(&jan_feb, &kernel, 30.0, 4);
        let feb = lag_monthly_volumes(&feb_only, &kernel, 30.0, 4);
        assert!((ab[1].1 - (a[1].1 + feb[0].1)).abs() < 1e-10);
    }

    #[test]
    fn average_month_matches_published_volume() {
        // Independent evaluation of Jenkins v = Q t · 4i²erfc(√(sdf/4t)) with
        // t = k × 30.41667 days and a one-month pulse of 100 acre-feet at
        // sdf = 265 days. Computed from the formula, not from AWAS.
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let mut stress = HashMap::new();
        stress.insert(month(2025, 1), 100.0);
        let out = lag_monthly_volumes_with_month_length(&stress, &kernel, MonthLength::Average, 6)
            .unwrap();
        let expected = [0.890056, 7.94116, 9.759836, 7.795874, 6.079456, 4.845931];
        assert_eq!(out.len(), expected.len());
        for (i, ((date, got), exp)) in out.iter().zip(expected.iter()).enumerate() {
            assert_eq!(*date, add_months(month(2025, 1), i as i32).unwrap());
            assert!((got - exp).abs() < 5e-6, "{date}: {got} vs {exp}");
        }
    }

    #[test]
    fn average_month_response_does_not_depend_on_the_calendar() {
        // A February pulse in a leap year and a July pulse give the same response.
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let run = |start: NaiveDate, mode: MonthLength| {
            let stress = HashMap::from([(start, 100.0)]);
            lag_monthly_volumes_with_month_length(&stress, &kernel, mode, 24).unwrap()
        };
        let feb = run(month(2024, 2), MonthLength::Average);
        let jul = run(month(2025, 7), MonthLength::Average);
        for (a, b) in feb.iter().zip(jul.iter()) {
            assert!((a.1 - b.1).abs() < 1e-12);
        }
        let feb = run(month(2024, 2), MonthLength::Calendar);
        let jul = run(month(2025, 7), MonthLength::Calendar);
        assert!((feb[0].1 - jul[0].1).abs() > 1e-3);
    }

    #[test]
    fn average_month_superposes_later_pulses_on_month_indexes() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 100.0 };
        let run = |stress: HashMap<NaiveDate, f64>| {
            lag_monthly_volumes_with_month_length(&stress, &kernel, MonthLength::Average, 240)
                .unwrap()
        };
        let first = run(HashMap::from([(month(1983, 11), 80.0)]));
        let later = run(HashMap::from([(month(1984, 2), 30.0)]));
        let both = run(HashMap::from([
            (month(1983, 11), 80.0),
            (month(1984, 2), 30.0),
        ]));
        assert_eq!(both[0].0, month(1983, 11));
        for i in 0..first.len() - 3 {
            let shifted = if i >= 3 { later[i - 3].1 } else { 0.0 };
            assert!((both[i].1 - first[i].1 - shifted).abs() < 1e-9, "month {i}");
        }
    }

    #[test]
    fn average_month_conserves_the_stress_volume() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 10.0 };
        let stress = HashMap::from([(month(2025, 1), -100.0)]);
        let out =
            lag_monthly_volumes_with_month_length(&stress, &kernel, MonthLength::Average, 1200)
                .unwrap();
        let total: f64 = out.iter().map(|(_, v)| v).sum();
        // After 100 years the unreturned part is 2u/√π ≈ 0.9% of the pulse.
        assert!((total + 100.0).abs() < 1.0, "total {total}");
        assert!(total > -100.0, "total {total}");
    }

    #[test]
    fn existing_entry_point_stays_calendar_whatever_days_per_month_says() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let stress = HashMap::from([(month(2025, 1), 100.0)]);
        let calendar =
            lag_monthly_volumes_with_month_length(&stress, &kernel, MonthLength::Calendar, 12)
                .unwrap();
        for days_per_month in [30.0, 30.42, AVERAGE_DAYS_PER_MONTH] {
            assert_eq!(
                lag_monthly_volumes(&stress, &kernel, days_per_month, 12),
                calendar
            );
        }
        assert_eq!(MonthLength::default(), MonthLength::Calendar);
    }

    #[test]
    fn bad_parameters_return_an_error_instead_of_numbers() {
        let stress = HashMap::from([(month(2025, 1), 100.0)]);
        for sdf_days in [0.0, -1.0, f64::NAN] {
            let kernel = ResponseKernel::EffectiveSdf { sdf_days };
            for mode in [MonthLength::Calendar, MonthLength::Average] {
                assert!(matches!(
                    lag_monthly_volumes_with_month_length(&stress, &kernel, mode, 6),
                    Err(ParameterError::Invalid {
                        name: "sdf_days",
                        ..
                    })
                ));
            }
        }
    }

    #[test]
    fn non_finite_stress_returns_an_error_naming_the_earliest_month() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let stress = HashMap::from([
            (month(2025, 1), 100.0),
            (month(2025, 3), f64::NAN),
            (month(2025, 2), f64::INFINITY),
        ]);
        assert_eq!(
            lag_monthly_volumes_with_month_length(&stress, &kernel, MonthLength::Average, 6),
            Err(ParameterError::InvalidStress {
                month: month(2025, 2),
                value: f64::INFINITY
            })
        );
    }

    #[test]
    fn bad_parameters_are_reported_even_with_nothing_to_lag() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 0.0 };
        assert!(lag_monthly_volumes_with_month_length(
            &HashMap::new(),
            &kernel,
            MonthLength::Calendar,
            6
        )
        .is_err());
    }

    #[test]
    fn average_month_length_is_the_value_files_rely_on() {
        assert_eq!(AVERAGE_DAYS_PER_MONTH, 30.41667);
    }

    #[test]
    fn calendar_mode_uses_the_length_of_a_leap_february() {
        // Independent evaluation with a 29-day pulse, then a 31-day March, of
        // 100 acre-feet at sdf = 265 days.
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let stress = HashMap::from([(month(2024, 2), 100.0)]);
        let out = lag_monthly_volumes_with_month_length(&stress, &kernel, MonthLength::Calendar, 2)
            .unwrap();
        assert!((out[0].1 - 0.76400847).abs() < 1e-8, "{}", out[0].1);
        assert!((out[1].1 - 8.02425554).abs() < 1e-8, "{}", out[1].1);
    }

    #[test]
    fn many_pulses_give_identical_results_every_time() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let run = || {
            // Each HashMap gets its own random iteration order.
            let stress: HashMap<NaiveDate, f64> = (0..48)
                .map(|i| (add_months(month(2020, 1), i).unwrap(), 1.1 + i as f64 / 7.0))
                .collect();
            lag_monthly_volumes(&stress, &kernel, 30.42, 96)
        };
        let first = run();
        for _ in 0..10 {
            assert_eq!(run(), first);
        }
    }
}
