mod glover_alluvial;
mod glover_infinite;
mod sdf;
mod urf;
mod utils;

pub use glover_alluvial::calculate_streamflow_depletion_alluvial;
pub use glover_infinite::calculate_streamflow_depletion_infinite;
pub use sdf::calculate_streamflow_depletion_sdf;
pub use urf::{UrfValue, urf_lagging};
