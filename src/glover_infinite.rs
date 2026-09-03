use crate::kernel::ResponseKernel;
use crate::lag::lag_monthly_volumes;
use chrono::NaiveDate;
use std::collections::HashMap;

/// Monthly stream depletion / accretion using Glover and Balmer (1954).
///
/// Glover, R.E., and Balmer, G.G., 1954, River depletion resulting from
/// pumping a well near a river: *Eos, Transactions American Geophysical
/// Union*, v. 35, no. 3, p. 468–470.
///
/// ```text
/// q/Q = erfc( √(a² S / (4 T t)) )
/// ```
///
/// Monthly volumes use Glover’s (1960) integral, `v = Q t · 4 i²erfc(u)`,
/// with Jenkins (1968) superposition for a rectangular monthly pulse.
/// This is identical to [`crate::sdf::calculate_streamflow_depletion_sdf`]
/// when `sdf = a² S / T`.
///
/// # Parameters
///
/// * `pumping_volumes_monthly` — monthly volumes; positive = pumping,
///   negative = recharge.
/// * `distance_to_well` — perpendicular well-to-stream distance `a`, feet.
/// * `specific_yield` — storativity / specific yield `S`, dimensionless.
/// * `transmissivity` — `T` in **ft²/day** (convert GPD/ft with
///   [`crate::kernel::gpd_per_ft_to_ft2_per_day`]).
/// * `days_per_month` — retained for API compatibility; pulse lengths use
///   calendar days.
/// * `total_months` — calendar months to report from the first stress.
pub fn calculate_streamflow_depletion_infinite(
    pumping_volumes_monthly: &HashMap<NaiveDate, f64>,
    distance_to_well: f64,
    specific_yield: f64,
    transmissivity: f64,
    days_per_month: f64,
    total_months: usize,
) -> Vec<(NaiveDate, f64)> {
    let kernel = ResponseKernel::GloverInfinite {
        distance_ft: distance_to_well,
        specific_yield,
        transmissivity_ft2_per_day: transmissivity,
    };
    lag_monthly_volumes(
        pumping_volumes_monthly,
        &kernel,
        days_per_month,
        total_months,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{gpd_per_ft_to_ft2_per_day, sdf_from_glover};
    use crate::sdf::calculate_streamflow_depletion_sdf;

    fn month(year: i32, month: u32) -> NaiveDate {
        NaiveDate::from_ymd_opt(year, month, 1).unwrap()
    }

    #[test]
    fn matches_sdf_method_for_ideal_sdf() {
        let a = 4000.0;
        let s = 0.2;
        let t = gpd_per_ft_to_ft2_per_day(261_800.0);
        let sdf = sdf_from_glover(a, s, t);
        assert!((sdf - 91.437).abs() < 0.01);

        let mut pumping = HashMap::new();
        pumping.insert(month(2025, 1), 100.0);

        let glover = calculate_streamflow_depletion_infinite(&pumping, a, s, t, 30.42, 6);
        let jenkins = calculate_streamflow_depletion_sdf(&pumping, sdf.round() as u32, 30.42, 6);
        // sdf is rounded to an integer day for the u32 API; allow that quantization.
        for (g, j) in glover.iter().zip(jenkins.iter()) {
            assert_eq!(g.0, j.0);
            assert!((g.1 - j.1).abs() / j.1.max(1e-6) < 0.02);
        }
    }

    #[test]
    fn january_100_af_matches_glover_volume_formula() {
        let a = 4000.0;
        let s = 0.2;
        let t = gpd_per_ft_to_ft2_per_day(261_800.0);
        let mut pumping = HashMap::new();
        pumping.insert(month(2025, 1), 100.0);
        let value = calculate_streamflow_depletion_infinite(&pumping, a, s, t, 30.42, 6);

        // Independent Jenkins/Glover volume evaluation of the same pulse.
        let expected = [9.231024, 20.792231, 13.118033, 7.592177, 5.342934, 3.804862];
        assert_eq!(value.len(), expected.len());
        for (got, exp) in value.iter().zip(expected.iter()) {
            assert!(
                (got.1 - exp).abs() < 5e-6,
                "{}: got {:.8} expected {:.8}",
                got.0,
                got.1,
                exp
            );
        }
    }
}
