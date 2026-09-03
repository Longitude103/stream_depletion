use crate::kernel::ResponseKernel;
use crate::lag::lag_monthly_volumes;
use chrono::NaiveDate;
use std::collections::HashMap;

/// Monthly stream depletion / accretion using Jenkins’ (1968) SDF method.
///
/// Jenkins, C.T., 1968, Techniques for computing rate and volume of stream
/// depletion by wells: *Ground Water*, v. 6, no. 2, p. 37–46; and USGS
/// Techniques of Water-Resources Investigations, Book 4, Chapter D1.
///
/// For an idealized aquifer `sdf = a² S / T` [days]. In South Platte practice
/// the sdf is often a mapped effective value (Hurr, Schneider, and others).
///
/// # Parameters
///
/// * `pumping_volumes_monthly` — monthly volumes, positive = pumping (depletion),
///   negative = recharge (accretion), typically acre-feet. Keys are dates in
///   the stress month (day-of-month is ignored).
/// * `sdf` — stream depletion factor in days.
/// * `days_per_month` — retained for API compatibility; pulse lengths use
///   calendar days (see [`crate::lag::lag_monthly_volumes`]).
/// * `total_months` — number of calendar months to report from the first stress.
///
/// # Returns
///
/// One entry per calendar month: `(month_start, volume)`. Sign follows the
/// stress: pumping produces positive depletion, recharge produces negative
/// values (accretion).
pub fn calculate_streamflow_depletion_sdf(
    pumping_volumes_monthly: &HashMap<NaiveDate, f64>,
    sdf: u32,
    days_per_month: f64,
    total_months: u32,
) -> Vec<(NaiveDate, f64)> {
    let kernel = ResponseKernel::EffectiveSdf {
        sdf_days: sdf as f64,
    };
    lag_monthly_volumes(
        pumping_volumes_monthly,
        &kernel,
        days_per_month,
        total_months as usize,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn month(year: i32, month: u32) -> NaiveDate {
        NaiveDate::from_ymd_opt(year, month, 1).unwrap()
    }

    #[test]
    fn january_100_af_sdf_265_matches_jenkins_volume() {
        // Independent evaluation of Jenkins v = Q t · 4i²erfc(√(sdf/4t))
        // with a 31-day rectangular pulse of 100 acre-feet, sdf = 265 days.
        // These literals were computed from the published volume formula, not
        // from AWAS and not from a previous revision of this crate.
        let mut pumping = HashMap::new();
        pumping.insert(month(2025, 1), 100.0);
        let value = calculate_streamflow_depletion_sdf(&pumping, 265, 30.42, 6);

        let expected = [
            (month(2025, 1), 0.944405),
            (month(2025, 2), 7.169295),
            (month(2025, 3), 10.023970),
            (month(2025, 4), 7.801918),
            (month(2025, 5), 6.279747),
            (month(2025, 6), 4.833422),
        ];
        assert_eq!(value.len(), expected.len());
        for (got, (date, exp)) in value.iter().zip(expected.iter()) {
            assert_eq!(got.0, *date);
            assert!(
                (got.1 - exp).abs() < 5e-6,
                "{}: got {:.8} expected {:.8}",
                got.0,
                got.1,
                exp
            );
        }
    }

    #[test]
    fn recharge_accretion_is_negative_mirror() {
        let mut recharge = HashMap::new();
        recharge.insert(month(2025, 1), -100.0);
        let value = calculate_streamflow_depletion_sdf(&recharge, 265, 30.42, 6);
        assert!(value.iter().all(|(_, v)| *v < 0.0));
        assert!((value[0].1 + 0.944405).abs() < 5e-6);
    }
}
