use scirs2_special::erfc;

pub fn calculate_streamflow_depletion_infinate(
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