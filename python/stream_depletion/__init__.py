"""Python bindings for the stream_depletion Rust crate.

Monthly (and daily URF) lagged stream effects for South Platte-style
accounting. Farmers Pawnee recharge uses :func:`lag_sdf` with **negative**
monthly volumes (accretion). Other ARI Glover plans can call
:func:`lag_glover_infinite` or :func:`lag_glover_alluvial`.

The native module is built with maturin / PyO3. Rust consumers of this
repository (lagging-api and any ``use stream_depletion`` CLI) are
unchanged: the ``python`` Cargo feature is off by default.
"""

from stream_depletion._native import (
    __version__,
    calculate_streamflow_depletion_alluvial,
    calculate_streamflow_depletion_infinite,
    calculate_streamflow_depletion_sdf,
    combined_urf_results,
    combined_urf_results_dated,
    daily_urf_mass_by_reach,
    disaggregate_monthly_urf_uniform,
    gpd_per_ft_to_ft2_per_day,
    lag_glover_alluvial,
    lag_glover_infinite,
    lag_sdf,
    sdf_from_glover,
    urf_lagging,
    urf_lagging_daily,
    urf_mass_by_reach,
)

__all__ = [
    "__version__",
    "calculate_streamflow_depletion_alluvial",
    "calculate_streamflow_depletion_infinite",
    "calculate_streamflow_depletion_sdf",
    "combined_urf_results",
    "combined_urf_results_dated",
    "daily_urf_mass_by_reach",
    "disaggregate_monthly_urf_uniform",
    "gpd_per_ft_to_ft2_per_day",
    "lag_glover_alluvial",
    "lag_glover_infinite",
    "lag_sdf",
    "sdf_from_glover",
    "urf_lagging",
    "urf_lagging_daily",
    "urf_mass_by_reach",
]
