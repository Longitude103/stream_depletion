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

use libm::erfc;

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
}

impl ResponseKernel {
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
}
