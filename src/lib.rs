mod glover_alluvial;
mod glover_infinite;
mod sdf;
mod urf;

use chrono::{NaiveDate, Datelike};

/// Adds a specified number of months to a given date.
///
/// This function calculates a new date by adding (or subtracting) a number of months
/// to the input date. It handles year rollovers and adjusts for varying month lengths.
///
/// # Parameters
///
/// * `date`: The starting `NaiveDate` to which months will be added.
/// * `months`: The number of months to add. Can be positive (to add months) or
///             negative (to subtract months).
///
/// # Returns
///
/// Returns an `Option<NaiveDate>`:
/// * `Some(NaiveDate)` if the resulting date is valid.
/// * `None` if the resulting date is invalid (e.g., February 30th).
fn add_months(date: NaiveDate, months: i32) -> Option<NaiveDate> {
    let mut year = date.year() + (date.month() as i32 + months - 1) / 12;
    let mut month = (date.month() as i32 + months - 1) % 12 + 1;
    if month <= 0 {
        month += 12;
        year -= 1;
    }
    NaiveDate::from_ymd_opt(year, month as u32, date.day())
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use super::*;
    use crate::glover_infinite::calculate_streamflow_depletion_infinate;
    use crate::glover_alluvial::calculate_streamflow_depletion_alluvial;

    #[test]
    fn test_with_infinite_aquifer() {
        // Aquifer parameters (in feet-based units)
        let d: f64 = 4000.0;  // Distance to stream (ft)
        let s: f64 = 0.2;      // Storativity (dimensionless)
        let t: f64 = 67000.0;  // Transmissivity (ft²/day)

        // Pumping rates in acre-feet/month for month 1
        let pumping_volumes = vec![100.0]; // acre-feet for month 1
        let days_per_month = 30.42; // Average days per month
        let total_months = 120; // 10 years

        // Convert pumping rates to ft³/day and calculate rate changes
        let mut pumping_rates: Vec<f64> = vec![0.0]; // Start with Q_0 = 0
        for &vol in &pumping_volumes {
            let rate = vol * 43_560.0 / days_per_month; // Convert acre-feet to ft³/day
            pumping_rates.push(rate);
        }
        // Set pumping to 0 after month 1
        for _ in pumping_volumes.len()..total_months {
            pumping_rates.push(0.0);
        }

        // Calculate depletion for each month
        let mut results = Vec::new();
        for month in 1..=total_months {
            let time = month as f64 * days_per_month; // Time at end of month (days)
            let qs_ft3_per_day = calculate_streamflow_depletion_infinate(&pumping_rates, d, s, t, time, days_per_month);
            // Convert to acre-feet/month
            let qs_af_per_month = qs_ft3_per_day * days_per_month / 43_560.0;
            results.push((month, qs_af_per_month));
        }

        // Print results
        println!("Infinite Aquifer Results:");
        println!("Month | Streamflow Depletion (acre-feet/month)");
        println!("------|-------------------------------------");
        for (month, qs) in results {
            println!("{:>5} | {:>35.2}", month, qs);
        }
    }
    
    #[test]
    fn test_with_alluvial_aquifer() {
        // Aquifer parameters (in feet-based units)
        let d: f64 = 4000.0;  // Distance to stream (ft)
        let b: f64 = 8000.0;     // Distance from well to boundary (ft)
        let s: f64 = 0.2;      // Storativity (dimensionless)
        let t: f64 = 261_800.0;  // Transmissivity (GPD/ft)

        // Pumping rates in acre-feet/month for month 1
        let mut pumping_volumes = HashMap::new();
        pumping_volumes.insert(NaiveDate::from_ymd_opt(2025, 1,1).unwrap(), 100.0); // acre-feet for month 1
        let days_per_month = 30.42; // Average days per month
        let total_months = 120; // 10 years

        let converted_t = t / 7.481; // Convert GPD to ft2/day
        let value = calculate_streamflow_depletion_alluvial(&pumping_volumes, d, b, s, converted_t, days_per_month, 120);
        println!("Monthly depletion amounts");
        for month in 0..value.len() {
            println!("{}: {}", value[month].0, value[month].1);
        }

    }
}
