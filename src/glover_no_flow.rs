use crate::kernel::{ParameterError, ResponseKernel};
use crate::lag::{lag_monthly_volumes_with_month_length, MonthLength};
use chrono::NaiveDate;
use std::collections::HashMap;

/// Monthly stream depletion / accretion for a well near a no-flow boundary
/// perpendicular to the stream.
///
/// The boundary adds a same-sign image well (Theis 1941; Ferris, Knowles, Brown,
/// and Stallman 1962; Hantush 1964). The image supplies exactly the depletion
/// the real well would have caused on the stream beyond the boundary, so
/// depletion on the **entire stream** equals Glover and Balmer (1954) for any
/// `boundary_distance`. See [`crate::kernel`] for the derivation.
///
/// # Parameters
///
/// * `pumping_volumes_monthly` — monthly volumes; positive = pumping,
///   negative = recharge.
/// * `distance_to_well` — perpendicular well-to-stream distance `a`, feet.
/// * `boundary_distance` — well-to-boundary distance `b` along the stream, feet
///   (IDS AWAS `m_b`). This is not the parallel-wall width used by
///   [`crate::glover_alluvial`].
/// * `specific_yield` — `S`, dimensionless.
/// * `transmissivity` — `T` in **ft²/day**.
/// * `month_length` — calendar or average months.
/// * `total_months` — months to report from the first stress.
///
/// Returns an error for invalid parameters or stress volumes (see
/// [`ResponseKernel::validate`]).
pub fn calculate_streamflow_depletion_no_flow(
    pumping_volumes_monthly: &HashMap<NaiveDate, f64>,
    distance_to_well: f64,
    boundary_distance: f64,
    specific_yield: f64,
    transmissivity: f64,
    month_length: MonthLength,
    total_months: usize,
) -> Result<Vec<(NaiveDate, f64)>, ParameterError> {
    let kernel = ResponseKernel::GloverNoFlow {
        distance_ft: distance_to_well,
        boundary_distance_ft: boundary_distance,
        specific_yield,
        transmissivity_ft2_per_day: transmissivity,
    };
    lag_monthly_volumes_with_month_length(
        pumping_volumes_monthly,
        &kernel,
        month_length,
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
    fn calendar_results_equal_glover_infinite_and_ignore_the_boundary_distance() {
        let (a, s, t) = (4000.0, 0.2, gpd_per_ft_to_ft2_per_day(261_800.0));
        let stress = HashMap::from([(month(2025, 1), 100.0), (month(2025, 3), -40.0)]);
        let infinite = calculate_streamflow_depletion_infinite(&stress, a, s, t, 30.42, 24);
        for b in [500.0, 20_000.0] {
            let no_flow = calculate_streamflow_depletion_no_flow(
                &stress,
                a,
                b,
                s,
                t,
                MonthLength::Calendar,
                24,
            )
            .unwrap();
            assert_eq!(no_flow, infinite);
        }
    }

    #[test]
    fn average_months_match_published_volume() {
        // Independent evaluation of v = Q t · 4i²erfc(a / √(4 T t / S)) with
        // t = k × 30.41667 days. Computed from the formula, not from AWAS.
        let (a, s, t) = (4000.0, 0.2, gpd_per_ft_to_ft2_per_day(261_800.0));
        let stress = HashMap::from([(month(2025, 1), 100.0)]);
        let out = calculate_streamflow_depletion_no_flow(
            &stress,
            a,
            12_000.0,
            s,
            t,
            MonthLength::Average,
            6,
        )
        .unwrap();
        let expected = [8.981528, 22.246955, 12.397944, 7.507171, 5.143469, 3.80267];
        assert_eq!(out.len(), expected.len());
        for ((date, got), exp) in out.iter().zip(expected.iter()) {
            assert!((got - exp).abs() < 5e-6, "{date}: {got} vs {exp}");
        }
    }

    #[test]
    fn bad_parameters_return_an_error() {
        let stress = HashMap::from([(month(2025, 1), 100.0)]);
        let run = |a: f64, b: f64, s: f64, t: f64| {
            calculate_streamflow_depletion_no_flow(&stress, a, b, s, t, MonthLength::Calendar, 6)
        };
        assert!(run(1000.0, 500.0, 0.2, 30_000.0).is_ok());
        assert!(run(0.0, 500.0, 0.2, 30_000.0).is_err());
        assert!(run(1000.0, -500.0, 0.2, 30_000.0).is_err());
        assert!(run(1000.0, 500.0, 0.0, 30_000.0).is_err());
        assert!(run(1000.0, 500.0, 0.2, -30_000.0).is_err());
        assert!(run(1000.0, f64::NAN, 0.2, 30_000.0).is_err());
    }
}
