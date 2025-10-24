pub mod glover_alluvial;
pub mod glover_infinite;
pub mod sdf;
pub mod urf;
pub mod utils;

pub use glover_alluvial::calculate_streamflow_depletion_alluvial;
pub use glover_infinite::calculate_streamflow_depletion_infinite;
pub use sdf::calculate_streamflow_depletion_sdf;
pub use urf::{LaggedUrfByDate, LaggedUrfResult, UrfValue, combined_urf_results, urf_lagging};
pub use utils::add_months;
