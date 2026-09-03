use crate::kernel::ResponseKernel;
use crate::lag::lag_monthly_volumes;
use chrono::NaiveDate;
use std::collections::HashMap;

/// Monthly stream depletion / accretion for an alluvial (strip) aquifer.
///
/// Image-well extension of Glover and Balmer (1954) for a constant-head
/// stream at `x = 0` and a parallel impermeable boundary at `x = W`
/// (Glover 1960, 1977; Schroeder 1987 SEODEP; McWhorter and Sunada 1977;
/// Miller, Durnford, Halstead, Altenhofen, and Flory 2007).
///
/// `aquifer_width` is the **stream-to-wall** distance `W` (the same quantity
/// as IDS AWAS `m_w` / SEODEP `W`). The well is at distance `a` from the
/// stream and must satisfy `0 < a < W`.
///
/// # Parameters
///
/// * `pumping_volumes_monthly` — monthly volumes; positive = pumping,
///   negative = recharge.
/// * `distance_to_well` — well-to-stream distance `a`, feet.
/// * `aquifer_width` — stream-to-no-flow-boundary distance `W`, feet.
/// * `specific_yield` — `S`, dimensionless.
/// * `transmissivity` — `T` in **ft²/day**.
/// * `days_per_month` — retained for API compatibility; pulse lengths use
///   calendar days.
/// * `total_months` — calendar months to report from the first stress.
pub fn calculate_streamflow_depletion_alluvial(
    pumping_volumes_monthly: &HashMap<NaiveDate, f64>,
    distance_to_well: f64,
    aquifer_width: f64,
    specific_yield: f64,
    transmissivity: f64,
    days_per_month: f64,
    total_months: usize,
) -> Vec<(NaiveDate, f64)> {
    let kernel = ResponseKernel::GloverAlluvial {
        distance_ft: distance_to_well,
        aquifer_width_ft: aquifer_width,
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
    use crate::glover_infinite::calculate_streamflow_depletion_infinite;
    use crate::kernel::gpd_per_ft_to_ft2_per_day;

    fn month(year: i32, month: u32) -> NaiveDate {
        NaiveDate::from_ymd_opt(year, month, 1).unwrap()
    }

    #[test]
    fn january_100_af_alluvial_matches_image_series_volume() {
        let a = 4000.0;
        let w = 8000.0;
        let s = 0.2;
        let t = gpd_per_ft_to_ft2_per_day(261_800.0);
        let mut pumping = HashMap::new();
        pumping.insert(month(2025, 1), 100.0);
        let value = calculate_streamflow_depletion_alluvial(&pumping, a, w, s, t, 30.42, 6);

        let expected = [
            9.234067, 21.070681, 14.755124, 10.312423, 8.462012, 6.632475,
        ];
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

    #[test]
    fn alluvial_is_near_infinite_early_then_larger() {
        let a = 4000.0;
        let w = 8000.0;
        let s = 0.2;
        let t = gpd_per_ft_to_ft2_per_day(261_800.0);
        let mut pumping = HashMap::new();
        pumping.insert(month(2025, 1), 100.0);
        let inf = calculate_streamflow_depletion_infinite(&pumping, a, s, t, 30.42, 6);
        let al = calculate_streamflow_depletion_alluvial(&pumping, a, w, s, t, 30.42, 6);
        assert!((al[0].1 - inf[0].1).abs() < 0.01);
        assert!(al[5].1 > inf[5].1);
    }
}
