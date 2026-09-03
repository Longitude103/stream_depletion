//! Jenkins (1968) superposition of rectangular monthly pulses.
//!
//! Each calendar month of stress is a pulse of duration equal to the number of
//! days in that month. Linear superposition applies to both pumping (+) and
//! recharge (−). Monthly output is the increment of cumulative stream volume
//! across that calendar month (`v(end) − v(start)`), in the same units as the
//! input volumes.
//!
//! This is the volume form of Jenkins’ residual-effects construction (USGS
//! TWRI 4-D1), not a finite-difference of the *rate* curve.

use crate::kernel::{pulse_cumulative_volume, ResponseKernel};
use crate::utils::add_months;
use chrono::{Datelike, NaiveDate};
use std::collections::HashMap;

/// Lag monthly pumping (+) or recharge (−) volumes through `kernel`.
///
/// * `stress_volumes_monthly` — acre-feet (or any consistent volume) keyed by
///   a date in the stress month. The day-of-month is ignored; the pulse covers
///   the calendar month.
/// * `total_months` — number of calendar months to report, starting at the
///   earliest stress month.
///
/// `days_per_month` is accepted for source compatibility with the previous
/// API. Accounting uses actual calendar month lengths, which is the
/// research-correct rectangular pulse. IDS AWAS monthly mode instead uses
/// 30.41667 days (`365/12`) unless “actual days” is selected.
pub fn lag_monthly_volumes(
    stress_volumes_monthly: &HashMap<NaiveDate, f64>,
    kernel: &ResponseKernel,
    _days_per_month: f64,
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

    let pulses: Vec<Pulse> = month_stress
        .iter()
        .filter(|(_, v)| **v != 0.0)
        .map(|(month, volume)| {
            let days = month.num_days_in_month() as f64;
            Pulse {
                start: *month,
                duration_days: days,
                rate: *volume / days,
            }
        })
        .collect();

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
        let mut increment = 0.0;
        for pulse in &pulses {
            increment += pulse_volume_between(pulse, month_start, month_end, kernel);
        }
        results.push((month_start, increment));
    }
    results
}

struct Pulse {
    start: NaiveDate,
    duration_days: f64,
    rate: f64,
}

fn pulse_volume_between(
    pulse: &Pulse,
    window_start: NaiveDate,
    window_end: NaiveDate,
    kernel: &ResponseKernel,
) -> f64 {
    let t_end = (window_end - pulse.start).num_days() as f64;
    let t_start = (window_start - pulse.start).num_days() as f64;
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
}
