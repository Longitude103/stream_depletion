use std::collections::HashMap;
use scirs2_special::erfc; // Crate for error function calculations
use chrono::{NaiveDate, Datelike};

fn calculate_streamflow_depletion_infinate(
    pumping_rates: &[f64],
    d: f64,
    s: f64,
    t: f64,
    time: f64,
    days_per_month: f64,
) -> f64 {
    let mut qs_total = 0.0;
    // Sum depletion contributions from each rate change
    for i in 1..pumping_rates.len() {
        let delta_q = pumping_rates[i] - pumping_rates[i - 1]; // Rate change
        let start_time = (i - 1) as f64 * days_per_month; // Start of month i
        if time > start_time {
            let qs = delta_q * calculate_depletion_fraction(d, s, t, time - start_time);
            qs_total += qs;
        }
    }
    qs_total.max(0.0) // Ensure depletion is non-negative
}


/// Calculates streamflow depletion for an alluvial aquifer based on monthly pumping volumes.
///
/// This function computes the streamflow depletion over time in an alluvial aquifer setting,
/// taking into account various hydrogeological parameters and pumping volumes.
///
/// # Parameters
///
/// * `pumping_volumes_monthly`: A HashMap containing monthly pumping volumes (in acre-ft/month)
///   indexed by their corresponding dates.
/// * `distance_to_well`: The distance from the well to the stream (in feet).
/// * `distance_to_boundary`: The distance to the aquifer boundary (in feet).
/// * `specific_yield`: The specific yield of the aquifer (dimensionless).
/// * `transmissivity`: The transmissivity of the aquifer (in ft²/day).
/// * `days_per_month`: The average number of days per month used in calculations.
/// * `total_months`: The total number of months to calculate depletion for.
///
/// # Returns
///
/// A vector of tuples, where each tuple contains:
/// * A `NaiveDate` representing the start of a month.
/// * A `f64` value representing the streamflow depletion for that month (in acre-ft/month).
///
/// The vector only includes months when the depletion is greater than 0.001 acre-ft/month.
/// The calculation stops if a negative depletion value is encountered, indicating complete aquifer depletion.
fn calculate_streamflow_depletion_alluvial(
    pumping_volumes_monthly: &HashMap<NaiveDate, f64>,  // Monthly pumping volumes in acre-ft / month
    distance_to_well: f64,
    distance_to_boundary: f64,
    specific_yield: f64,
    transmissivity: f64,
    days_per_month: f64,
    total_months: usize,
) -> Vec<(NaiveDate, f64)> {
    let total_days = (total_months as f64 * days_per_month).ceil() as usize;
    
    // 1. calculate the depletion fraction for each time step
    let mut base_depletion_fraction = vec![0.0; total_days];
    for m in 0..total_days {
        base_depletion_fraction[m] = calculate_depletion_fraction_alluvial_aquifer(distance_to_well, distance_to_boundary, specific_yield, transmissivity, m as f64);
    }

    // println!("Base Depletion Fraction");
    // for step in 0..120 {
    //     println!("{}: {}", step, base_depletion_fraction[step]);
    // }

    // total up base_depletion_fraction
    // let mut total_base_depletion_fraction = 0.0;
    // for step in 0..total_days {
    //     total_base_depletion_fraction += base_depletion_fraction[step];
    // }
    //
    // println!("Total Base Depletion Fraction: {}", total_base_depletion_fraction);
    
    // 2. convert pumping_volumes_monthly to pumping_rates_daily using the number of days in the month of the NaiveDate
    let mut pumping_rates_daily = HashMap::new();
    for (date, pumping_volume) in pumping_volumes_monthly {
        let m = date.month() as usize - 1; // Month index (0-11)
        let days_in_month = date.num_days_in_month();

        // for each day in the month, calculate the daily pumping rate, and store it in pumping_rates_daily by NaiveDate and amount
        for d in 0..days_in_month {
            let date_daily = NaiveDate::from_ymd_opt(date.year(), date.month(), (d + 1u8) as u32).unwrap();
            let daily_pumping_rate = pumping_volume * 43_560f64 / (days_in_month as f64);
            *pumping_rates_daily.entry(date_daily).or_insert(0.0) += daily_pumping_rate;
        }
    }

    // println!("{:?}", pumping_rates_daily);  // order is not sorted
    // println!("Daily pumping rates");
    // let start_date = NaiveDate::from_ymd_opt(2025, 1, 1).unwrap();
    // for day in 0..35 {
    //     let pump_date = start_date + chrono::Duration::days(day as i64);  // depletion is always the day after the pumping occurs
    //     println!("{}: {}", pump_date, pumping_rates_daily.get(&pump_date).unwrap_or(&0.0));
    // }

    // 3. Create a daily results Hashmap with daily time steps to hold the daily depletion amounts
    let mut daily_depletion_amount = HashMap::new();

    for (date, pumping_rate) in pumping_rates_daily {
        if pumping_rate <= 0.0 {
            continue;
        }
        let mut day_depletion = vec![0.0; total_days];
        for base_depletion_index in 0..base_depletion_fraction.len() {
            day_depletion[base_depletion_index] = pumping_rate * base_depletion_fraction[base_depletion_index];
        }

        // add the day depletion to the daily depletion amount for the corresponding date and forward
        for depletion_index in 0..day_depletion.len() {
            let depletion_date = date + chrono::Duration::days(depletion_index as i64 + 1i64);  // depletion is always the day after the pumping occurs
            if depletion_index == 0 {
                *daily_depletion_amount.entry(depletion_date).or_insert(0.0) += day_depletion[depletion_index];
                continue;
            }

            *daily_depletion_amount.entry(depletion_date).or_insert(0.0) += day_depletion[depletion_index] - day_depletion[depletion_index - 1];
        }
    }

    // println!("Daily depletion amounts");
    // let start_date = NaiveDate::from_ymd_opt(2025, 1, 1).unwrap();
    // for day in 0..35 {
    //     let depletion_date = start_date + chrono::Duration::days(day as i64);  // depletion is always the day after the pumping occurs
    //     println!("{}: {}", depletion_date, daily_depletion_amount.get(&depletion_date).unwrap_or(&0.0));
    // }

    // println!("{:?}", daily_depletion_amount);  // order is not sorted, this is ft³/day

    // 4. sum the daily depletion amounts to monthly depletion totals and convert to acre-ft / month from ft³/month
    let mut monthly_depletion_amount = HashMap::new();
    for (date, depletion_amount) in daily_depletion_amount {
        let monthly_date = NaiveDate::from_ymd_opt(date.year(), date.month(), 1).unwrap();  // Monthly date
        *monthly_depletion_amount.entry(monthly_date).or_insert(0.0) += depletion_amount / 43560f64;  // Convert ft³ to acre-ft
    }

    // println!("{:?}", monthly_depletion_amount);  // order is sorted

    // println!("Monthly depletion amounts");
    // let start_date = NaiveDate::from_ymd_opt(2025, 1, 1).unwrap();
    // for month in 0..60 {
    //     let result_date = add_months(start_date, month).unwrap();  // depletion is always the day after the pumping occurs
    //     println!("{}: {}", result_date, monthly_depletion_amount.get(&result_date).unwrap_or(&0.0));
    // }
    //
    // // 5. sum the monthly depletion amounts to get the total depletion
    // let total_depletion = monthly_depletion_amount.values().sum::<f64>();
    //
    // println!("Total depletion: {}", total_depletion);

    let mut results: Vec<(NaiveDate, f64)> = vec![];
    // start date should be the oldest date key in the pumping_volumes_monthly HashMap
    let start_date = pumping_volumes_monthly.keys().min().unwrap().clone();
    results.reserve(total_months);  // Reserve space for results to avoid reallocating
    // let start_date = NaiveDate::from_ymd_opt(2025, 1, 1).unwrap();  // should get from the input parameters
    for month in 0..total_months {
        let result_date = add_months(start_date, month as i32).unwrap();  // depletion is always the day after the pumping occurs
        let monthly_depletion = *monthly_depletion_amount.get(&result_date).unwrap_or(&0.0);

        if monthly_depletion < 0.0 {
            // The depletion is negative, which means the aquifer has been depleted completely,
            // so we should stop the simulation and return the results.
            break;
        }

        if monthly_depletion > 0.001 {
            results.push((result_date, monthly_depletion));
        }
    }

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

/// Calculates the depletion fraction for streamflow depletion in an alluvial aquifer.
///
/// This function computes the fraction of pumping that has been captured from the stream
/// at a given time, based on aquifer properties and the distance to the stream, specifically
/// for an alluvial aquifer scenario. It also includes the image wells that are determined by
/// the factor that is created.
///
/// # Parameters
///
/// * `distance_to_well`: Distance from the well to the stream (in length units, typically feet).
/// * `distance_to_boundary`: Distance from the well to boundary (in length units, typically feet).
/// * `specific_yield`: Storativity of the aquifer (dimensionless).
/// * `transmissivity`: Transmissivity of the aquifer (in length²/time units, typically ft²/day).
/// * `time`: Time since pumping began (in time units, typically days).
///
/// # Returns
///
/// Returns the depletion fraction as a `f64`, representing the proportion of pumping
/// that has been captured from the stream at the given time in an alluvial aquifer setting.
fn calculate_depletion_fraction_alluvial_aquifer(distance_to_well: f64, distance_to_boundary: f64,
    specific_yield: f64, transmissivity: f64, time: f64) -> f64 {
    let mut total_depletion_fraction = 0.0;
    let mut image_factor = 1.0;
    let mut well_distance = -distance_to_well;  // distance is negative to account for first loop
    
    loop {
        // Real well or positive image well
        well_distance += 2.0 * distance_to_well;
        let u = well_distance / (4.0 * transmissivity * time / (specific_yield)).sqrt();
        let depletion_fraction = if u > 2.9 { 0.0 } else { erfc(u) };
        total_depletion_fraction += depletion_fraction * image_factor;

        if depletion_fraction == 0.0 {
            break;
        }
        
        // Negative image well
        well_distance = well_distance - 2.0 * distance_to_well + 2.0 * distance_to_boundary;
        let u = well_distance / (4.0 * transmissivity * time / (specific_yield)).sqrt();
        let depletion_fraction = if u > 2.9 { 0.0 } else { erfc(u) };
        total_depletion_fraction += depletion_fraction * image_factor;

        if depletion_fraction == 0.0 {
            break;
        }

        image_factor *= -1.0; // Alternate sign for next pair of image wells
    }

    total_depletion_fraction
}

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
    use super::*;

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
