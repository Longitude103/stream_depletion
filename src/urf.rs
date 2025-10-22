use chrono::{Months, NaiveDate};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ops::Add;

#[derive(Serialize, Deserialize, Debug)]
pub struct UrfValue {
    pub month: i32,
    pub reach: i32,
    pub urf_val: f64,
}

impl UrfValue {
    pub fn new(month: i32, reach: i32, urf_val: f64) -> Self {
        UrfValue {
            month,
            reach,
            urf_val,
        }
    }
}

/// Computes the lagged usage rate factor (URF) for different reaches over time.
///
/// This function takes a vector of `UrfValue` and a usage map, and calculates the lagged URF
/// for each reach over the given usage dates. The result is a nested `HashMap` where the outer
/// key is the reach identifier, and the inner key-value pairs represent the date and the corresponding
/// lagged URF value.
///
/// # Parameters
///
/// - `urf`: A vector of `UrfValue` structs, each containing the month, reach, and URF value.
/// - `usage`: A reference to a `HashMap` where the keys are `NaiveDate` objects representing the usage dates,
///   and the values are the usage amounts for those dates.
///
/// # Returns
///
/// A `HashMap` where the keys are reach identifiers (`i32`), and the values are another `HashMap` with
/// `NaiveDate` keys and `f64` values representing the lagged URF for each date.
pub fn urf_lagging(
    usage: &HashMap<NaiveDate, f64>,
    urf: Vec<UrfValue>,
) -> HashMap<i32, HashMap<NaiveDate, f64>> {
    let reaches = urf.iter().map(|u| u.reach).unique().collect::<Vec<_>>();
    let usage_dates: Vec<&NaiveDate> = usage.keys().into_iter().sorted().collect();

    let mut lagged_result = HashMap::new();
    for reach in reaches {
        let mut reach_lagged = HashMap::new();
        let reach_urf = urf
            .iter()
            .filter(|u| u.reach == reach)
            .sorted_by_key(|u| u.month)
            .map(|u| u.urf_val)
            .collect::<Vec<f64>>();

        for usage_date in &usage_dates {
            let month_usage = usage.get(*usage_date).unwrap_or(&0.0);
            for (i, urf) in reach_urf.iter().enumerate() {
                let urf_date = usage_date.add(Months::new(i as u32));
                let urf_dep = month_usage * urf;
                *reach_lagged.entry(urf_date).or_insert(0.0) += urf_dep;
            }
        }

        lagged_result.insert(reach, reach_lagged);
    }

    lagged_result
}

/// Creates a combined result of depletion from the lagged_result provided by the urf_lagging function
///
/// This function takes a hashmap of the lagged URF results and combines all the values for each reach
/// into a single value. If there are multiple values across the reaches for the same month, they are summed up.
/// The result is then sorted by date and returned as a vector of tuples.
///
/// # Parameters
///
/// - `values`: A `HashMap` where the keys are reach identifiers (`i32`), and the values are another `HashMap` with
/// `NaiveDate` keys and `f64` values representing the lagged URF for each date.
///
/// # Returns
///
/// A vector of tuples, where each tuple contains:
/// * A `NaiveDate` representing the start of a month.
/// * A `f64` value representing the streamflow depletion for that month (in acre-ft/month).
///
/// The vector only includes months when the depletion is greater than 0.001 acre-ft/month.
/// The calculation stops if a negative depletion value is encountered, indicating complete aquifer depletion.
pub fn combined_urf_results(
    values: HashMap<i32, HashMap<NaiveDate, f64>>,
) -> Vec<(NaiveDate, f64)> {
    // Aggregate f64 values by NaiveDate
    let mut date_sums: HashMap<NaiveDate, f64> = HashMap::new();

    // Iterate over all inner HashMaps
    for inner_map in values.values() {
        for (date, value) in inner_map {
            // Sum values for each date
            *date_sums.entry(*date).or_insert(0.0) += value;
        }
    }

    // Convert to Vec and sort by date
    let mut result: Vec<(NaiveDate, f64)> = date_sums.into_iter().collect();
    result.sort_by(|a, b| a.0.cmp(&b.0));

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_urf_lagging() {
        let urf = vec![
            UrfValue {
                month: 1,
                reach: 1,
                urf_val: 0.6,
            },
            UrfValue {
                month: 1,
                reach: 2,
                urf_val: 0.1,
            },
            UrfValue {
                month: 2,
                reach: 1,
                urf_val: 0.3,
            },
        ];
        let mut usage = HashMap::new();
        usage.insert(NaiveDate::from_ymd_opt(2024, 7, 1).unwrap(), 100.0);
        usage.insert(NaiveDate::from_ymd_opt(2024, 8, 1).unwrap(), 100.0);

        let mut expected_lagged = HashMap::new();
        let mut reach1 = HashMap::new();
        reach1.insert(NaiveDate::from_ymd_opt(2024, 7, 1).unwrap(), 60.0);
        reach1.insert(NaiveDate::from_ymd_opt(2024, 8, 1).unwrap(), 90.0);
        reach1.insert(NaiveDate::from_ymd_opt(2024, 9, 1).unwrap(), 30.0);

        let mut reach2 = HashMap::new();
        reach2.insert(NaiveDate::from_ymd_opt(2024, 7, 1).unwrap(), 10.0);
        reach2.insert(NaiveDate::from_ymd_opt(2024, 8, 1).unwrap(), 10.0);
        expected_lagged.insert(1, reach1);
        expected_lagged.insert(2, reach2);

        let result = urf_lagging(&usage, urf);
        assert_eq!(result, expected_lagged);
    }

    #[test]
    fn test_urf_lagging_skip_usage_month() {
        let urf = vec![
            UrfValue {
                month: 1,
                reach: 1,
                urf_val: 0.4,
            },
            UrfValue {
                month: 1,
                reach: 2,
                urf_val: 0.2,
            },
            UrfValue {
                month: 2,
                reach: 1,
                urf_val: 0.2,
            },
            UrfValue {
                month: 2,
                reach: 2,
                urf_val: 0.1,
            },
            UrfValue {
                month: 3,
                reach: 1,
                urf_val: 0.1,
            },
        ];
        let mut usage = HashMap::new();
        usage.insert(NaiveDate::from_ymd_opt(2024, 5, 1).unwrap(), 100.0);
        usage.insert(NaiveDate::from_ymd_opt(2024, 7, 1).unwrap(), 100.0);
        usage.insert(NaiveDate::from_ymd_opt(2024, 8, 1).unwrap(), 100.0);

        let mut expected_lagged = HashMap::new();
        let mut reach1 = HashMap::new();
        reach1.insert(NaiveDate::from_ymd_opt(2024, 5, 1).unwrap(), 40.0);
        reach1.insert(NaiveDate::from_ymd_opt(2024, 6, 1).unwrap(), 20.0);
        reach1.insert(NaiveDate::from_ymd_opt(2024, 7, 1).unwrap(), 50.0);
        reach1.insert(NaiveDate::from_ymd_opt(2024, 8, 1).unwrap(), 60.0);
        reach1.insert(NaiveDate::from_ymd_opt(2024, 9, 1).unwrap(), 30.0);
        reach1.insert(NaiveDate::from_ymd_opt(2024, 10, 1).unwrap(), 10.0);
        let mut reach2 = HashMap::new();
        reach2.insert(NaiveDate::from_ymd_opt(2024, 5, 1).unwrap(), 20.0);
        reach2.insert(NaiveDate::from_ymd_opt(2024, 6, 1).unwrap(), 10.0);
        reach2.insert(NaiveDate::from_ymd_opt(2024, 7, 1).unwrap(), 20.0);
        reach2.insert(NaiveDate::from_ymd_opt(2024, 8, 1).unwrap(), 30.0);
        reach2.insert(NaiveDate::from_ymd_opt(2024, 9, 1).unwrap(), 10.0);
        expected_lagged.insert(1, reach1);
        expected_lagged.insert(2, reach2);

        let result = urf_lagging(&usage, urf);
        assert_eq!(result, expected_lagged);
    }
}
