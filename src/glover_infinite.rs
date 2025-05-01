use crate::glover_alluvial::{
    create_monthly_depletion, create_results_vector, monthly_pumping_to_daily,
};
use chrono::NaiveDate;
use scirs2_special::erfc;
use std::collections::HashMap;

/// Calculates streamflow depletion for an infinite aquifer using the Glover solution.
///
/// This function computes the monthly streamflow depletion based on given pumping volumes and aquifer parameters.
/// It uses the Glover solution for an infinite aquifer to determine the depletion fractions and applies them to the pumping rates.
///
/// # Parameters
///
/// * `pumping_volumes_monthly`: A HashMap containing monthly pumping volumes in acre-ft/month, keyed by date.
/// * `distance_to_well`: The distance from the well to the stream in feet.
/// * `specific_yield`: The specific yield of the aquifer (dimensionless).
/// * `transmissivity`: The transmissivity of the aquifer in ft²/day.
/// * `days_per_month`: The average number of days per month used in calculations.
/// * `total_months`: The total number of months to calculate depletion for.
///
/// # Returns
///
/// A Vec of tuples, where each tuple contains a date and the corresponding monthly streamflow depletion in acre-ft/month.
pub fn calculate_streamflow_depletion_infinite(
    pumping_volumes_monthly: &HashMap<NaiveDate, f64>, // Monthly pumping volumes in acre-ft / month
    distance_to_well: f64,
    specific_yield: f64,
    transmissivity: f64,
    days_per_month: f64,
    total_months: usize,
) -> Vec<(NaiveDate, f64)> {
    // get total days
    let total_days = (total_months as f64 * days_per_month).ceil() as usize;

    // 1. calculate the depletion fraction for each time step
    let mut base_depletion_fraction = vec![0.0; total_days];
    for m in 0..total_days {
        base_depletion_fraction[m] = calculate_depletion_fraction(
            distance_to_well,
            specific_yield,
            transmissivity,
            m as f64,
        );
    }

    let pumping_rates_daily = monthly_pumping_to_daily(pumping_volumes_monthly);

    // 3. Create a daily results Hashmap with daily time steps to hold the daily depletion amounts
    let mut daily_depletion_amount = HashMap::new();
    for (date, pumping_rate) in pumping_rates_daily {
        if pumping_rate <= 0.0 {
            continue;
        }
        let mut day_depletion = vec![0.0; total_days];
        for base_depletion_index in 0..base_depletion_fraction.len() {
            day_depletion[base_depletion_index] =
                pumping_rate * base_depletion_fraction[base_depletion_index];
        }

        // add the day depletion to the daily depletion amount for the corresponding date and forward
        for depletion_index in 0..day_depletion.len() {
            let depletion_date = date + chrono::Duration::days(depletion_index as i64 + 1i64); // depletion is always the day after the pumping occurs
            if depletion_index == 0 {
                *daily_depletion_amount.entry(depletion_date).or_insert(0.0) +=
                    day_depletion[depletion_index];
                continue;
            }

            *daily_depletion_amount.entry(depletion_date).or_insert(0.0) +=
                day_depletion[depletion_index] - day_depletion[depletion_index - 1];
        }
    }

    let monthly_depletion_amount = create_monthly_depletion(&daily_depletion_amount);
    let results = create_results_vector(
        pumping_volumes_monthly,
        total_months,
        &monthly_depletion_amount,
    );

    results
}

/// Calculates the depletion fraction for streamflow depletion using the Glover solution.
///
/// This function computes the fraction of pumping that has been captured from the stream
/// at a given time, based on aquifer properties and the distance to the stream.
///
/// # Parameters
///
/// * `d`: Distance from the well to the stream (in length units, typically feet).
/// * `s`: Storativity of the aquifer (dimensionless).
/// * `t`: Transmissivity of the aquifer (in length²/time units, typically ft²/day).
/// * `time`: Time since pumping began (in time units, typically days).
///
/// # Returns
///
/// Returns the depletion fraction as a `f64`, representing the proportion of pumping
/// that has been captured from the stream at the given time.
fn calculate_depletion_fraction(d: f64, s: f64, t: f64, time: f64) -> f64 {
    // Calculate the argument of the complementary error function
    let z = ((s * d.powi(2)) / (4.0 * t * time)).sqrt();
    // Calculate erfc(z)
    erfc(z)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to round a float to 5 decimal places
    fn round_to_5_decimals(value: f64) -> f64 {
        (value * 100_000.0).round() / 100_000.0
    }

    #[test]
    fn test_with_infinite_aquifer() {
        // Aquifer parameters (in feet-based units)
        let d: f64 = 4000.0; // Distance to stream (ft)
        let s: f64 = 0.2; // Storativity (dimensionless)
        let t: f64 = 261_800.0; // Transmissivity (GPD/ft)

        // Pumping rates in acre-feet/month for month 1
        let mut pumping_volumes = HashMap::new();
        pumping_volumes.insert(NaiveDate::from_ymd_opt(2025, 1, 1).unwrap(), 100.0); // acre-feet for month 1
        let days_per_month = 30.42; // Average days per month
        let total_months = 120; // 10 years

        let converted_t = t / 7.481; // Convert GPD to ft2/day
        let value = calculate_streamflow_depletion_infinite(
            &pumping_volumes,
            d,
            s,
            converted_t,
            days_per_month,
            total_months,
        );
        // println!("Monthly depletion amounts");
        // for month in 0..value.len() {
        //     println!("{}: {}", value[month].0, value[month].1);
        // }

        assert!(value.len() <= total_months); // Test if results vector has correct length

        let tolerance = 0.00001; // 10^-5 for 5 decimal places
        
        // values that should be checked are:
        // 2025-01-01: 8.169915278703847
        // 2025-02-01: 20.979264088137487
        // 2025-03-01: 13.514164851251204
        // 2025-04-01: 7.75855587035028
        // 2025-05-01: 5.433551969020377
        // 2025-06-01: 3.857354439468754
        assert!((round_to_5_decimals(value[0].1) - 8.16991).abs() < tolerance);
        assert!((round_to_5_decimals(value[1].1) - 20.97926).abs() < tolerance);
        assert!((round_to_5_decimals(value[2].1) - 13.51416).abs() < tolerance);
        assert!((round_to_5_decimals(value[3].1) - 7.75856).abs() < tolerance);
        assert!((round_to_5_decimals(value[4].1) - 5.43355).abs() < tolerance);
        assert!((round_to_5_decimals(value[5].1) - 3.85735).abs() < tolerance);
    }
}
