//! Analytical stream depletion / lagged accretion for South Platte-style
//! monthly accounting.
//!
//! Kernels follow the published Glover–Balmer (1954) rate equation and the
//! Glover (1960) / Jenkins (1968) volume equation `v = Q t · 4 i²erfc(u)`,
//! not IDS AWAS output. See `docs/RESEARCH_VS_AWAS.md`.

pub mod glover_alluvial;
pub mod glover_infinite;
pub mod kernel;
pub mod lag;
pub mod sdf;
pub mod urf;
pub mod utils;

#[cfg(feature = "python")]
mod python;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Maturin loads `stream_depletion._native` (see `pyproject.toml`).
/// The `python` feature is off by default so Rust dependents do not link PyO3.
#[cfg(feature = "python")]
#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    python::register(module)
}

pub use glover_alluvial::calculate_streamflow_depletion_alluvial;
pub use glover_infinite::calculate_streamflow_depletion_infinite;
pub use kernel::{
    four_i2erfc, gpd_per_ft_to_ft2_per_day, sdf_from_glover, ResponseKernel, GALLONS_PER_CUBIC_FOOT,
};
pub use sdf::calculate_streamflow_depletion_sdf;
pub use urf::{
    combined_urf_results, combined_urf_results_dated, daily_urf_mass_by_reach,
    disaggregate_monthly_urf_uniform, month_start, urf_lagging, urf_lagging_daily,
    urf_mass_by_reach, DailyUrfValue, LaggedUrfByDate, LaggedUrfResult, UrfValue,
};
pub use utils::add_months;
