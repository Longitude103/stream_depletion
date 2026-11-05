//! Regression against sanitized IDS AWAS output in `tests/fixtures`.
//!
//! AWAS is not the authority for this crate (see `docs/RESEARCH_VS_AWAS.md`).
//! These lock the option settings under which the published formulas and AWAS
//! agree, so that an accidental change to either mode is caught.

use chrono::NaiveDate;
use std::collections::HashMap;
use stream_depletion::{
    calculate_streamflow_depletion_no_flow, gpd_per_ft_to_ft2_per_day,
    lag_monthly_volumes_with_month_length, utils::add_months, MonthLength, ResponseKernel,
};

const TOLERANCE_AF: f64 = 1e-5;

struct Fixture {
    stress: Vec<f64>,
    awas: Vec<f64>,
}

impl Fixture {
    fn parse(text: &str) -> Self {
        let mut stress = Vec::new();
        let mut awas = Vec::new();
        for line in text.lines().filter(|l| !l.starts_with('#')).skip(1) {
            let mut fields = line.split(',').skip(1).map(|v| v.parse::<f64>().unwrap());
            stress.push(fields.next().unwrap());
            awas.push(fields.next().unwrap());
        }
        assert_eq!(stress.len(), 96);
        Fixture { stress, awas }
    }

    fn stress_by_month(&self, start: NaiveDate) -> HashMap<NaiveDate, f64> {
        self.stress
            .iter()
            .enumerate()
            .map(|(i, v)| (add_months(start, i as i32).unwrap(), *v))
            .collect()
    }

    fn assert_matches(&self, lagged: &[(NaiveDate, f64)]) {
        assert_eq!(lagged.len(), self.awas.len());
        for (i, ((_, got), expected)) in lagged.iter().zip(&self.awas).enumerate() {
            assert!(
                (got - expected).abs() < TOLERANCE_AF,
                "month {i}: got {got}, AWAS {expected}"
            );
        }
    }
}

#[test]
fn average_month_sdf_recharge_matches_awas() {
    let fixture = Fixture::parse(include_str!("fixtures/awas_average_month_sdf.csv"));
    let start = NaiveDate::from_ymd_opt(2001, 1, 1).unwrap();
    let lagged = lag_monthly_volumes_with_month_length(
        &fixture.stress_by_month(start),
        &ResponseKernel::EffectiveSdf { sdf_days: 750.0 },
        MonthLength::Average,
        fixture.awas.len(),
    )
    .unwrap();
    fixture.assert_matches(&lagged);
}

#[test]
fn calendar_month_no_flow_pumping_matches_awas() {
    let fixture = Fixture::parse(include_str!("fixtures/awas_no_flow_calendar.csv"));
    let start = NaiveDate::from_ymd_opt(2002, 1, 1).unwrap();
    let lagged = calculate_streamflow_depletion_no_flow(
        &fixture.stress_by_month(start),
        5797.5,
        10_000.0,
        0.2,
        gpd_per_ft_to_ft2_per_day(426_562.5),
        MonthLength::Calendar,
        fixture.awas.len(),
    )
    .unwrap();
    fixture.assert_matches(&lagged);
}
