use chrono::{Datelike, NaiveDate};

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
pub(crate) fn add_months(date: NaiveDate, months: i32) -> Option<NaiveDate> {
    let mut year = date.year() + (date.month() as i32 + months - 1) / 12;
    let mut month = (date.month() as i32 + months - 1) % 12 + 1;
    if month <= 0 {
        month += 12;
        year -= 1;
    }
    NaiveDate::from_ymd_opt(year, month as u32, date.day())
}
