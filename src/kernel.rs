//! Glover–Balmer / Jenkins stream-depletion kernels.
//!
//! Rate (Glover and Balmer 1954; Jenkins 1968, USGS TWRI 4-D1, curve A):
//!
//! ```text
//! q/Q = erfc(u),    u = sqrt(sdf / (4 t)) = a / sqrt(4 T t / S)
//! ```
//!
//! Volume (Glover 1960; Hantush 1964; Jenkins 1968, curve B):
//!
//! ```text
//! v / (Q t) = 4 i²erfc(u) = (1 + 2 u²) erfc(u) − (2 u / √π) exp(−u²)
//! ```
//!
//! `sdf` is Jenkins’ stream-depletion factor [T]. In an idealized semi-infinite
//! aquifer, `sdf = a² S / T`. Jenkins also *defines* sdf as the time at which
//! `v = 0.28 Qt` (because `4 i²erfc(1/2) ≈ 0.27986`).
//!
//! Alluvial (strip) aquifers add Glover’s image-well series for a constant-head
//! stream at `x = 0` and a parallel no-flow boundary at `x = W` (Glover 1960/1977;
//! Schroeder 1987 SEODEP; McWhorter and Sunada 1977, p. 127; Miller et al. 2007):
//!
//! ```text
//! q/Q = Σ_{n=0}^∞ (−1)^n [ erfc(u(2nW + a)) + erfc(u(2(n+1)W − a)) ]
//! ```
//!
//! and the same coefficients on `4 i²erfc` for `v/(Q t)`. `W` is the aquifer
//! width (stream to impermeable valley wall), **not** the well-to-wall distance.
//!
//! A no-flow boundary **perpendicular** to the stream (distance `b` along the
//! stream from the well) adds a same-sign image well at the same distance `a`
//! from the stream, `2b` away along it (Theis 1941; Ferris et al. 1962; Hantush
//! 1964). The stream’s flux density from a well at distance `a`, at distance
//! `x` along the stream, is
//!
//! ```text
//! f(x) = a / (π (a² + x²)) · exp(−(a² + x²) S / (4 T t)),   ∫ f dx = erfc(u)
//! ```
//!
//! Depletion on the real stream (`x < b`) is therefore
//! `∫_{x<b} f(x) dx + ∫_{x<b} f(x − 2b) dx`. Mirroring the second integral
//! across `x = b` makes it `∫_{x>b} f(x) dx`, so the two sum to the whole-line
//! integral: `q/Q = erfc(u)` for every `b`. The boundary moves depletion along
//! the stream but does not change the total for the entire stream.

use chrono::NaiveDate;
use libm::erfc;
use std::fmt;

/// US gallon per cubic foot used by IDS AWAS C++ (`7.48051945`).
/// Jenkins (1968) rounded this to 7.48; SEODEP BASIC used 7.481.
pub const GALLONS_PER_CUBIC_FOOT: f64 = 7.48051945;

/// Convert transmissivity from gallons/day/ft to ft²/day.
pub fn gpd_per_ft_to_ft2_per_day(transmissivity_gpd_per_ft: f64) -> f64 {
    transmissivity_gpd_per_ft / GALLONS_PER_CUBIC_FOOT
}

/// Complementary error function.
#[inline]
pub fn complementary_error(u: f64) -> f64 {
    if !u.is_finite() {
        return if u.is_sign_positive() { 0.0 } else { 2.0 };
    }
    erfc(u)
}

/// Jenkins / Glover volume fraction `v/(Q t) = 4 i²erfc(u)`.
///
/// Equivalent to Abramowitz and Stegun (1964) §7.2 via
/// `4 i²erfc(u) = (1 + 2u²) erfc(u) − 2u exp(−u²) / √π`.
pub fn four_i2erfc(u: f64) -> f64 {
    if !u.is_finite() {
        return if u.is_sign_positive() { 0.0 } else { 1.0 };
    }
    if u <= 0.0 {
        return 1.0;
    }
    // erfc(u) and exp(−u²) underflow to 0 well before u = 26.
    if u > 26.0 {
        return 0.0;
    }
    let erfc_u = complementary_error(u);
    let exp_term = (-u * u).exp();
    (1.0 + 2.0 * u * u) * erfc_u - (2.0 * u * exp_term) / std::f64::consts::PI.sqrt()
}

/// `u = √(sdf / (4 t))`. `t = 0` → `+∞` so both fractions are 0.
pub fn sdf_argument(sdf_days: f64, time_days: f64) -> f64 {
    if time_days <= 0.0 || sdf_days < 0.0 || !time_days.is_finite() {
        return f64::INFINITY;
    }
    if sdf_days == 0.0 {
        return 0.0;
    }
    (sdf_days / (4.0 * time_days)).sqrt()
}

/// `u = a / √(4 T t / S) = √(a² S / (4 T t))`.
pub fn glover_argument(
    distance_ft: f64,
    specific_yield: f64,
    transmissivity_ft2_per_day: f64,
    time_days: f64,
) -> f64 {
    if time_days <= 0.0 || transmissivity_ft2_per_day <= 0.0 || specific_yield < 0.0 {
        return f64::INFINITY;
    }
    if distance_ft == 0.0 || specific_yield == 0.0 {
        return 0.0;
    }
    distance_ft / (4.0 * transmissivity_ft2_per_day * time_days / specific_yield).sqrt()
}

/// Idealized Jenkins sdf [days] = `a² S / T`.
pub fn sdf_from_glover(
    distance_ft: f64,
    specific_yield: f64,
    transmissivity_ft2_per_day: f64,
) -> f64 {
    if transmissivity_ft2_per_day <= 0.0 {
        return f64::INFINITY;
    }
    distance_ft * distance_ft * specific_yield / transmissivity_ft2_per_day
}

/// A kernel parameter or stress value that cannot produce a meaningful result.
#[derive(Clone, Debug, PartialEq)]
pub enum ParameterError {
    /// `name` is non-finite or outside its allowed range.
    Invalid {
        name: &'static str,
        value: f64,
        requirement: &'static str,
    },
    /// The well is farther from the stream than the no-flow wall is.
    DistanceBeyondAquiferWidth { distance_ft: f64, width_ft: f64 },
    /// A monthly volume is NaN or infinite.
    InvalidStress { month: NaiveDate, value: f64 },
}

impl fmt::Display for ParameterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParameterError::Invalid {
                name,
                value,
                requirement,
            } => write!(f, "{name} must be {requirement}, got {value}"),
            ParameterError::DistanceBeyondAquiferWidth {
                distance_ft,
                width_ft,
            } => write!(
                f,
                "well-to-stream distance ({distance_ft} ft) is beyond the aquifer width ({width_ft} ft)"
            ),
            ParameterError::InvalidStress { month, value } => {
                write!(f, "stress volume for {month} must be finite, got {value}")
            }
        }
    }
}

impl std::error::Error for ParameterError {}

fn require(
    name: &'static str,
    value: f64,
    requirement: &'static str,
    valid: bool,
) -> Result<(), ParameterError> {
    if value.is_finite() && valid {
        Ok(())
    } else {
        Err(ParameterError::Invalid {
            name,
            value,
            requirement,
        })
    }
}

fn require_glover(a: f64, s: f64, t: f64) -> Result<(), ParameterError> {
    require("distance_ft", a, "finite and greater than 0", a > 0.0)?;
    require(
        "specific_yield",
        s,
        "finite, greater than 0, and at most 1",
        s > 0.0 && s <= 1.0,
    )?;
    require(
        "transmissivity_ft2_per_day",
        t,
        "finite and greater than 0",
        t > 0.0,
    )
}

/// Response kernels used for rate (`q/Q`) and volume (`v/Qt`) fractions.
#[derive(Clone, Debug)]
pub enum ResponseKernel {
    /// Jenkins effective SDF (semi-infinite Glover with `sdf` replacing `a²S/T`).
    EffectiveSdf { sdf_days: f64 },
    /// Glover–Balmer (1954) infinite / semi-infinite aquifer.
    GloverInfinite {
        distance_ft: f64,
        specific_yield: f64,
        transmissivity_ft2_per_day: f64,
    },
    /// Glover / Schroeder alluvial strip: stream at 0, no-flow at `aquifer_width_ft`.
    GloverAlluvial {
        distance_ft: f64,
        aquifer_width_ft: f64,
        specific_yield: f64,
        transmissivity_ft2_per_day: f64,
    },
    /// Glover–Balmer with a no-flow boundary perpendicular to the stream,
    /// `boundary_distance_ft` along the stream from the well.
    ///
    /// Whole-stream depletion does not depend on that distance (see the module
    /// docs), so this gives the same fractions as [`Self::GloverInfinite`]. The
    /// distance is kept so a site’s geometry stays explicit. Depletion on a
    /// *segment* of stream would depend on it, and is not implemented.
    GloverNoFlow {
        distance_ft: f64,
        boundary_distance_ft: f64,
        specific_yield: f64,
        transmissivity_ft2_per_day: f64,
    },
}

impl ResponseKernel {
    /// Check that the parameters describe a real site.
    ///
    /// [`Self::rate_fraction`] and [`Self::volume_fraction`] do not check, so
    /// they can return a meaningless number for bad input. Call this first, or
    /// use [`crate::lag::lag_monthly_volumes_with_month_length`], which does.
    ///
    /// A boundary distance of 0 (well on the boundary) is valid for
    /// [`Self::GloverNoFlow`]. An alluvial well may sit at the wall, not past it.
    pub fn validate(&self) -> Result<(), ParameterError> {
        match *self {
            ResponseKernel::EffectiveSdf { sdf_days } => require(
                "sdf_days",
                sdf_days,
                "finite and greater than 0",
                sdf_days > 0.0,
            ),
            ResponseKernel::GloverInfinite {
                distance_ft,
                specific_yield,
                transmissivity_ft2_per_day,
            } => require_glover(distance_ft, specific_yield, transmissivity_ft2_per_day),
            ResponseKernel::GloverAlluvial {
                distance_ft,
                aquifer_width_ft,
                specific_yield,
                transmissivity_ft2_per_day,
            } => {
                require_glover(distance_ft, specific_yield, transmissivity_ft2_per_day)?;
                require(
                    "aquifer_width_ft",
                    aquifer_width_ft,
                    "finite and greater than 0",
                    aquifer_width_ft > 0.0,
                )?;
                if distance_ft > aquifer_width_ft {
                    return Err(ParameterError::DistanceBeyondAquiferWidth {
                        distance_ft,
                        width_ft: aquifer_width_ft,
                    });
                }
                Ok(())
            }
            ResponseKernel::GloverNoFlow {
                distance_ft,
                boundary_distance_ft,
                specific_yield,
                transmissivity_ft2_per_day,
            } => {
                require_glover(distance_ft, specific_yield, transmissivity_ft2_per_day)?;
                require(
                    "boundary_distance_ft",
                    boundary_distance_ft,
                    "finite and not negative",
                    boundary_distance_ft >= 0.0,
                )
            }
        }
    }

    pub fn rate_fraction(&self, time_days: f64) -> f64 {
        self.fold_images(time_days, complementary_error)
    }

    pub fn volume_fraction(&self, time_days: f64) -> f64 {
        self.fold_images(time_days, four_i2erfc)
    }

    fn fold_images(&self, time_days: f64, term: fn(f64) -> f64) -> f64 {
        match *self {
            ResponseKernel::EffectiveSdf { sdf_days } => term(sdf_argument(sdf_days, time_days)),
            ResponseKernel::GloverInfinite {
                distance_ft,
                specific_yield,
                transmissivity_ft2_per_day,
            } => term(glover_argument(
                distance_ft,
                specific_yield,
                transmissivity_ft2_per_day,
                time_days,
            )),
            ResponseKernel::GloverAlluvial {
                distance_ft,
                aquifer_width_ft,
                specific_yield,
                transmissivity_ft2_per_day,
            } => alluvial_series(
                distance_ft,
                aquifer_width_ft,
                specific_yield,
                transmissivity_ft2_per_day,
                time_days,
                term,
            ),
            ResponseKernel::GloverNoFlow {
                distance_ft,
                specific_yield,
                transmissivity_ft2_per_day,
                ..
            } => term(glover_argument(
                distance_ft,
                specific_yield,
                transmissivity_ft2_per_day,
                time_days,
            )),
        }
    }
}

/// Schroeder / SEODEP image-well series (entire stream, option 2).
///
/// Positions and signs match `SEODEP2.BAS` lines 3360–3560 and AWAS
/// `SDFdata.cpp` (non-segment alluvial path):
/// `+f(a) + f(2W−a) − f(2W+a) − f(4W−a) + f(4W+a) + f(6W−a) − …`
///
/// Terms are taken in pairs and the loop stops when a pair’s contribution is
/// negligible. That avoids the Cesàro oscillation of the infinite `+2, 0, +2, …`
/// partial sums as `t → ∞` (each `erfc → 1`).
fn alluvial_series(
    distance_ft: f64,
    aquifer_width_ft: f64,
    specific_yield: f64,
    transmissivity_ft2_per_day: f64,
    time_days: f64,
    term: fn(f64) -> f64,
) -> f64 {
    if time_days <= 0.0 || transmissivity_ft2_per_day <= 0.0 || aquifer_width_ft <= 0.0 {
        return 0.0;
    }

    let mut total = 0.0;
    let mut image_factor = 1.0;
    let mut well_distance = -distance_ft;
    // Pair-wise tail: last pair at u ≈ 8 is already < 1e-15 for both erfc and 4i²erfc.
    const PAIR_TOL: f64 = 1e-15;
    const MAX_PAIRS: usize = 10_000;

    for _ in 0..MAX_PAIRS {
        well_distance += 2.0 * distance_ft;
        let u1 = glover_argument(
            well_distance,
            specific_yield,
            transmissivity_ft2_per_day,
            time_days,
        );
        let t1 = term(u1);

        well_distance = well_distance - 2.0 * distance_ft + 2.0 * aquifer_width_ft;
        let u2 = glover_argument(
            well_distance,
            specific_yield,
            transmissivity_ft2_per_day,
            time_days,
        );
        let t2 = term(u2);

        let pair = (t1 + t2) * image_factor;
        total += pair;
        if pair.abs() < PAIR_TOL && u2 > 1.0 {
            break;
        }
        image_factor *= -1.0;
    }
    total
}

/// Cumulative stream-volume change at elapsed time `t_eval` caused by a
/// rectangular pulse of rate `q` from `t = 0` to `t = duration` (Jenkins
/// superposition, TWRI 4-D1 “residual effects”).
///
/// `v(t) = q t · 4i²erfc(u(t))` while the pulse is on;
/// after shutoff, subtract an opposite pulse that started at `duration`.
pub fn pulse_cumulative_volume(
    rate: f64,
    duration_days: f64,
    t_eval_days: f64,
    kernel: &ResponseKernel,
) -> f64 {
    if rate == 0.0 || t_eval_days <= 0.0 || duration_days <= 0.0 {
        return 0.0;
    }
    let v_on = rate * t_eval_days * kernel.volume_fraction(t_eval_days);
    if t_eval_days <= duration_days {
        return v_on;
    }
    let t_off = t_eval_days - duration_days;
    v_on - rate * t_off * kernel.volume_fraction(t_off)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    #[test]
    fn jenkins_sdf_definition_is_28_percent_at_t_equals_sdf() {
        // Jenkins (1968): sdf is the time where v = 28% of Qt.
        // Analytically that is 4 i²erfc(1/2).
        let v_over_qt = four_i2erfc(0.5);
        assert!((v_over_qt - 0.2798588938127078).abs() < 1e-12);
        assert!((complementary_error(0.5) - 0.4795001221869535).abs() < 1e-12);

        let sdf = 100.0;
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: sdf };
        assert!((kernel.volume_fraction(sdf) - v_over_qt).abs() < TOL);
        assert!((kernel.rate_fraction(sdf) - complementary_error(0.5)).abs() < TOL);
    }

    #[test]
    fn glover_infinite_matches_sdf_when_sdf_equals_a2s_over_t() {
        // Jenkins TWRI 4-D1 rounds a=3660 ft and T/S=134000 ft²/day to sdf=100 d;
        // the exact identity is sdf = a² / (T/S) = 3660²/134000 ≈ 99.967 d.
        // Use T/S = a²/100 so the kernels are algebraically identical.
        let a = 3660.0;
        let s = 0.2;
        let t = a * a * s / 100.0;
        let sdf = sdf_from_glover(a, s, t);
        assert!((sdf - 100.0).abs() < 1e-9);

        let glover = ResponseKernel::GloverInfinite {
            distance_ft: a,
            specific_yield: s,
            transmissivity_ft2_per_day: t,
        };
        let jenkins = ResponseKernel::EffectiveSdf { sdf_days: sdf };
        for days in [1.0, 10.0, 35.0, 100.0, 365.0, 3650.0] {
            assert!((glover.rate_fraction(days) - jenkins.rate_fraction(days)).abs() < 1e-12);
            assert!((glover.volume_fraction(days) - jenkins.volume_fraction(days)).abs() < 1e-12);
        }
    }

    #[test]
    fn jenkins_twri_sample_35_days_at_sdf_100() {
        // Jenkins TWRI 4-D1: a = 3660 ft, T/S = 134000 ft²/day, sdf = 100 days,
        // Q = 10 acre-ft/day, tp = 35 days.
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 100.0 };
        let u = sdf_argument(100.0, 35.0);
        assert!((u - (100.0_f64 / 140.0).sqrt()).abs() < 1e-14);
        let q_over_q = complementary_error(u);
        let v_over_qt = four_i2erfc(u);
        assert!((q_over_q - 0.2319977236287341).abs() < 1e-12);
        assert!((v_over_qt - 0.09656945903460429).abs() < 1e-12);

        let q = 10.0; // acre-ft/day
        let v = pulse_cumulative_volume(q, 35.0, 35.0, &kernel);
        assert!((v - q * 35.0 * v_over_qt).abs() < 1e-12);
        assert!((v - 33.7993106621115).abs() < 1e-9);
    }

    #[test]
    fn residual_volume_after_shutoff_uses_superposition() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 100.0 };
        let q = 10.0;
        // 10 days after a 35-day pulse: v = v_on(45) − v_on(10)
        let v = pulse_cumulative_volume(q, 35.0, 45.0, &kernel);
        let v_on_45 = q * 45.0 * kernel.volume_fraction(45.0);
        let v_on_10 = q * 10.0 * kernel.volume_fraction(10.0);
        assert!((v - (v_on_45 - v_on_10)).abs() < 1e-12);
        assert!(v > 33.7); // residual volume continues to grow after shutoff
        assert!(v < q * 35.0);
    }

    #[test]
    fn recharge_is_the_negative_of_pumping() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        let pump = pulse_cumulative_volume(5.0, 31.0, 90.0, &kernel);
        let rech = pulse_cumulative_volume(-5.0, 31.0, 90.0, &kernel);
        assert!((pump + rech).abs() < 1e-12);
        assert!(rech < 0.0);
    }

    #[test]
    fn t_zero_is_zero_not_nan() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: 265.0 };
        assert_eq!(kernel.rate_fraction(0.0), 0.0);
        assert_eq!(kernel.volume_fraction(0.0), 0.0);
        assert_eq!(pulse_cumulative_volume(10.0, 31.0, 0.0, &kernel), 0.0);
    }

    #[test]
    fn alluvial_exceeds_infinite_and_approaches_one() {
        let a = 4000.0;
        let w = 8000.0;
        let s = 0.2;
        let t = gpd_per_ft_to_ft2_per_day(261_800.0);
        let inf = ResponseKernel::GloverInfinite {
            distance_ft: a,
            specific_yield: s,
            transmissivity_ft2_per_day: t,
        };
        let al = ResponseKernel::GloverAlluvial {
            distance_ft: a,
            aquifer_width_ft: w,
            specific_yield: s,
            transmissivity_ft2_per_day: t,
        };
        // Early time: wall not felt.
        assert!((inf.rate_fraction(30.0) - al.rate_fraction(30.0)).abs() < 1e-3);
        // Later: no-flow wall increases depletion rate toward 1.
        assert!(al.rate_fraction(365.0) > inf.rate_fraction(365.0) + 0.1);
        assert!((al.rate_fraction(3650.0) - 1.0).abs() < 1e-4);
        assert!(al.volume_fraction(3650.0) > inf.volume_fraction(3650.0));
        assert!(al.volume_fraction(3650.0) < 1.0);
    }

    #[test]
    fn gallons_conversion_matches_awas_cpp_factor() {
        let t = gpd_per_ft_to_ft2_per_day(261_800.0);
        assert!((t - 261_800.0 / 7.48051945).abs() < 1e-12);
    }

    /// `∫ f(x) dx` from `x = −∞` to `x = upper` for the stream flux density in
    /// the module docs. With `x = a tan θ` the integrand becomes
    /// `exp(−a² sec²θ · S / (4 T t)) / π`, which Simpson’s rule handles well.
    fn stream_flux_up_to(upper: f64, a: f64, s: f64, t: f64, days: f64) -> f64 {
        const STEPS: usize = 20_000;
        let lo = -std::f64::consts::FRAC_PI_2;
        let hi = (upper / a).atan();
        let h = (hi - lo) / STEPS as f64;
        let f = |theta: f64| {
            let cos = theta.cos();
            if cos <= 1e-12 {
                return 0.0;
            }
            (-(a * a) * s / (4.0 * t * days) / (cos * cos)).exp() / std::f64::consts::PI
        };
        let mut sum = f(lo) + f(hi);
        for i in 1..STEPS {
            let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
            sum += weight * f(lo + i as f64 * h);
        }
        sum * h / 3.0
    }

    #[test]
    fn no_flow_equals_real_plus_image_flux_on_the_real_stream() {
        // Integrates the image-well flux directly, without using erfc.
        let (a, s, t) = (5000.0, 0.2, gpd_per_ft_to_ft2_per_day(250_000.0));
        for b in [500.0, 15_000.0, 100_000.0] {
            let kernel = ResponseKernel::GloverNoFlow {
                distance_ft: a,
                boundary_distance_ft: b,
                specific_yield: s,
                transmissivity_ft2_per_day: t,
            };
            for days in [5.0, 30.0, 365.0, 3650.0] {
                let real = stream_flux_up_to(b, a, s, t, days);
                // Image at x = 2b: real-stream x < b is image-relative x' < −b.
                let image = stream_flux_up_to(-b, a, s, t, days);
                let rate = kernel.rate_fraction(days);
                assert!(
                    (real + image - rate).abs() < 1e-9,
                    "b={b} t={days}: {} vs {rate}",
                    real + image
                );
            }
        }
    }

    #[test]
    fn no_flow_boundary_changes_where_depletion_lands_but_not_the_total() {
        // Between the well and a nearby boundary, the image adds real depletion,
        // so a segment answer would depend on b even though the total does not.
        let (a, s, t, days) = (1000.0, 0.2, gpd_per_ft_to_ft2_per_day(250_000.0), 365.0);
        let b = 500.0;
        let kernel = ResponseKernel::GloverNoFlow {
            distance_ft: a,
            boundary_distance_ft: b,
            specific_yield: s,
            transmissivity_ft2_per_day: t,
        };
        let real_only = stream_flux_up_to(b, a, s, t, days);
        let image_share = stream_flux_up_to(-b, a, s, t, days);
        assert!(image_share > 0.1);
        assert!((kernel.rate_fraction(days) - real_only - image_share).abs() < 1e-9);
    }

    #[test]
    fn no_flow_matches_infinite_for_any_boundary_distance() {
        let (a, s, t) = (7000.0, 0.2, gpd_per_ft_to_ft2_per_day(200_000.0));
        let infinite = ResponseKernel::GloverInfinite {
            distance_ft: a,
            specific_yield: s,
            transmissivity_ft2_per_day: t,
        };
        for b in [0.0, 12_000.0, 1e7] {
            let no_flow = ResponseKernel::GloverNoFlow {
                distance_ft: a,
                boundary_distance_ft: b,
                specific_yield: s,
                transmissivity_ft2_per_day: t,
            };
            for days in [0.0, 1.0, 31.0, 365.0, 36_500.0] {
                assert_eq!(no_flow.rate_fraction(days), infinite.rate_fraction(days));
                assert_eq!(
                    no_flow.volume_fraction(days),
                    infinite.volume_fraction(days)
                );
            }
        }
    }

    fn glover(a: f64, s: f64, t: f64) -> [ResponseKernel; 3] {
        [
            ResponseKernel::GloverInfinite {
                distance_ft: a,
                specific_yield: s,
                transmissivity_ft2_per_day: t,
            },
            ResponseKernel::GloverAlluvial {
                distance_ft: a,
                aquifer_width_ft: 1e6,
                specific_yield: s,
                transmissivity_ft2_per_day: t,
            },
            ResponseKernel::GloverNoFlow {
                distance_ft: a,
                boundary_distance_ft: 500.0,
                specific_yield: s,
                transmissivity_ft2_per_day: t,
            },
        ]
    }

    fn invalid_name(kernel: &ResponseKernel) -> &'static str {
        match kernel.validate() {
            Err(ParameterError::Invalid { name, .. }) => name,
            other => panic!("{kernel:?}: expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn validate_accepts_real_sites() {
        for kernel in glover(1000.0, 0.2, 30_000.0) {
            assert_eq!(kernel.validate(), Ok(()));
        }
        let sdf = ResponseKernel::EffectiveSdf { sdf_days: 750.0 };
        assert_eq!(sdf.validate(), Ok(()));
    }

    #[test]
    fn validate_rejects_bad_glover_parameters() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            for kernel in glover(bad, 0.2, 30_000.0) {
                assert_eq!(invalid_name(&kernel), "distance_ft");
            }
            for kernel in glover(1000.0, bad, 30_000.0) {
                assert_eq!(invalid_name(&kernel), "specific_yield");
            }
            for kernel in glover(1000.0, 0.2, bad) {
                assert_eq!(invalid_name(&kernel), "transmissivity_ft2_per_day");
            }
        }
        for kernel in glover(1000.0, 1.5, 30_000.0) {
            assert_eq!(invalid_name(&kernel), "specific_yield");
        }
    }

    #[test]
    fn validate_rejects_bad_sdf() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let kernel = ResponseKernel::EffectiveSdf { sdf_days: bad };
            assert_eq!(invalid_name(&kernel), "sdf_days");
        }
    }

    #[test]
    fn validate_rejects_a_well_beyond_the_aquifer_width() {
        let alluvial = |w: f64| ResponseKernel::GloverAlluvial {
            distance_ft: 1000.0,
            aquifer_width_ft: w,
            specific_yield: 0.2,
            transmissivity_ft2_per_day: 30_000.0,
        };
        assert_eq!(
            alluvial(999.0).validate(),
            Err(ParameterError::DistanceBeyondAquiferWidth {
                distance_ft: 1000.0,
                width_ft: 999.0
            })
        );
        assert_eq!(alluvial(1000.0).validate(), Ok(()));
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(invalid_name(&alluvial(bad)), "aquifer_width_ft");
        }
    }

    #[test]
    fn validate_rejects_a_bad_no_flow_boundary_distance() {
        for (bad, ok) in [
            (-1.0, false),
            (f64::NAN, false),
            (f64::INFINITY, false),
            (0.0, true),
        ] {
            let kernel = ResponseKernel::GloverNoFlow {
                distance_ft: 1000.0,
                boundary_distance_ft: bad,
                specific_yield: 0.2,
                transmissivity_ft2_per_day: 30_000.0,
            };
            if ok {
                assert_eq!(kernel.validate(), Ok(()));
            } else {
                assert_eq!(invalid_name(&kernel), "boundary_distance_ft");
            }
        }
    }

    #[test]
    fn parameter_errors_say_what_is_wrong() {
        let kernel = ResponseKernel::EffectiveSdf { sdf_days: -5.0 };
        assert_eq!(
            kernel.validate().unwrap_err().to_string(),
            "sdf_days must be finite and greater than 0, got -5"
        );
    }
}
