use crate::glover_alluvial::{
    create_monthly_depletion, create_results_vector, monthly_pumping_to_daily,
};
use chrono::NaiveDate;
use scirs2_special::erfc;
use std::collections::HashMap;

/// Calculates streamflow depletion using the Stream Depletion Factor (SDF) method.
///
/// This function computes the monthly streamflow depletion based on given monthly pumping volumes,
/// SDF value, and time parameters. It uses the SDF method to estimate the impact of groundwater
/// pumping on streamflow over time.
///
/// # Parameters
///
/// * `pumping_volumes_monthly`: A HashMap containing monthly pumping volumes in acre-feet per month,
///   with NaiveDate keys representing the start of each month.
/// * `sdf`: The Stream Depletion Factor in days, representing the time it takes for stream depletion
///   to reach about 28% of the pumping rate.
/// * `days_per_month`: The average number of days per month used in calculations.
/// * `total_months`: The total number of months for which to calculate depletion.
///
/// # Returns
///
/// A Vec of tuples, where each tuple contains:
/// * A NaiveDate representing the start of a month
/// * A f64 value representing the calculated streamflow depletion for that month in acre-feet
pub fn calculate_streamflow_depletion_sdf(
    pumping_volumes_monthly: &HashMap<NaiveDate, f64>, // Monthly pumping volumes in acre-ft / month
    sdf: u32,
    days_per_month: f64,
    total_months: u32,
) -> Vec<(NaiveDate, f64)> {
    let total_days = (total_months as f64 * days_per_month).ceil() as usize;

    // 1. calculate the depletion fraction for each time step
    let mut base_depletion_fraction = vec![0.0; total_days as usize];
    for m in 0..total_days {
        base_depletion_fraction[m as usize] = calculate_depletion_fraction_sdf(sdf, m);
    }

    // println!("Base Depletion Fractions: {:?}", base_depletion_fraction);

    let pumping_rates_daily = monthly_pumping_to_daily(pumping_volumes_monthly);

    // println!("Pumping Rates Daily: {:?}", pumping_rates_daily);

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

            // println!(
            //     "Daily Depletion Amount for Date: {}, Day {} - Previous {}",
            //     depletion_date,
            //     day_depletion[depletion_index],
            //     day_depletion[depletion_index - 1]
            // );
            *daily_depletion_amount.entry(depletion_date).or_insert(0.0) +=
                day_depletion[depletion_index] - day_depletion[depletion_index - 1];
        }
    }

    // println!("Daily Depletion Amounts: {:?}", daily_depletion_amount);

    let monthly_depletion_amount = create_monthly_depletion(&daily_depletion_amount);
    let results = create_results_vector(
        pumping_volumes_monthly,
        total_months as usize,
        &monthly_depletion_amount,
    );

    results
}

/// Calculates the depletion fraction using the Stream Depletion Factor (SDF) method.
///
/// This function computes the fraction of pumping that has been depleted from the stream
/// at a given time step, based on the SDF value.
///
/// # Parameters
///
/// * `sdf`: The Stream Depletion Factor in days, representing the time it takes for stream
///   depletion to reach about 28% of the pumping rate.
/// * `time_step`: The current time step (in days) for which the depletion fraction is being calculated.
///
/// # Returns
///
/// A `f64` value representing the calculated depletion fraction at the given time step.
/// This value ranges from 0 to 1, where 0 means no depletion and 1 means complete depletion.
fn calculate_depletion_fraction_sdf(sdf: u32, time_step: usize) -> f64 {
    let u = (sdf as f64 / (4.0 * time_step as f64)).sqrt(); // u factor
    erfc(u)
}

#[cfg(test)]
mod tests {
    use super::*;


    // Helper function to round a float to 5 decimal places
    fn round_to_5_decimals(value: f64) -> f64 {
        (value * 100_000.0).round() / 100_000.0
    }

    #[test]
    fn test_with_sdf() {
        // Aquifer parameters (in feet-based units)
        let sdf: u32 = 265; // SDF Value in days

        // Pumping rates in acre-feet/month for month 1
        let mut pumping_volumes = HashMap::new();
        pumping_volumes.insert(NaiveDate::from_ymd_opt(2025, 1, 1).unwrap(), 100.0); // acre-feet for month 1
        let days_per_month = 30.42; // Average days per month
        let total_months = 120; // 10 years

        let value =
            calculate_streamflow_depletion_sdf(&pumping_volumes, sdf, days_per_month, total_months);
        // println!("Monthly depletion amounts");
        // for month in 0..value.len() {
        //     println!("{}: {}", value[month].0, value[month].1);
        // }
        
        assert!(value.len() <= total_months as usize); // Test if results vector has correct length

        let tolerance = 0.00001; // 10^-5 for 5 decimal places

        // values that should be checked are:
        // 2025-01-01: 0.7680351810235876
        // 2025-02-01: 6.842148015576688
        // 2025-03-01: 10.08458908541661
        // 2025-04-01: 7.8994824450947645
        // 2025-05-01: 6.35488902954355
        // 2025-06-01: 4.885147733308235
        assert!((round_to_5_decimals(value[0].1) - 0.76803).abs() < tolerance);
        assert!((round_to_5_decimals(value[1].1) - 6.84215).abs() < tolerance);
        assert!((round_to_5_decimals(value[2].1) - 10.08459).abs() < tolerance);
        assert!((round_to_5_decimals(value[3].1) - 7.89948).abs() < tolerance);
        assert!((round_to_5_decimals(value[4].1) - 6.35489).abs() < tolerance);
        assert!((round_to_5_decimals(value[5].1) - 4.88515).abs() < tolerance);
    }
}
