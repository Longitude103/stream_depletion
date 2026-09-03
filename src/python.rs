//! Optional PyO3 surface. Compiled only with `--features python` (maturin).
//!
//! Function signatures of the Rust library are unchanged. These wrappers
//! accept Python dates / dicts / sequences and return ISO date strings.

use crate::{
    calculate_streamflow_depletion_alluvial, calculate_streamflow_depletion_infinite,
    calculate_streamflow_depletion_sdf, combined_urf_results, combined_urf_results_dated,
    daily_urf_mass_by_reach, disaggregate_monthly_urf_uniform, gpd_per_ft_to_ft2_per_day,
    sdf_from_glover, urf_lagging, urf_lagging_daily, urf_mass_by_reach, DailyUrfValue,
    LaggedUrfResult, UrfValue,
};
use chrono::NaiveDate;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyModule};
use std::collections::HashMap;

fn parse_date(obj: &Bound<'_, PyAny>) -> PyResult<NaiveDate> {
    if let Ok(s) = obj.extract::<String>() {
        return parse_date_str(&s);
    }
    if obj.hasattr("year")? && obj.hasattr("month")? {
        let year: i32 = obj.getattr("year")?.extract()?;
        let month: u32 = obj.getattr("month")?.extract()?;
        let day: u32 = if obj.hasattr("day")? {
            obj.getattr("day")?.extract()?
        } else {
            1
        };
        return NaiveDate::from_ymd_opt(year, month, day).ok_or_else(|| {
            PyValueError::new_err(format!("invalid calendar date {year}-{month:02}-{day:02}"))
        });
    }
    if let Ok((year, month, day)) = obj.extract::<(i32, u32, u32)>() {
        return NaiveDate::from_ymd_opt(year, month, day).ok_or_else(|| {
            PyValueError::new_err(format!("invalid calendar date {year}-{month:02}-{day:02}"))
        });
    }
    if let Ok((year, month)) = obj.extract::<(i32, u32)>() {
        return NaiveDate::from_ymd_opt(year, month, 1)
            .ok_or_else(|| PyValueError::new_err(format!("invalid year-month ({year}, {month})")));
    }
    Err(PyValueError::new_err(
        "date must be 'YYYY-MM-DD', 'YYYY-MM', datetime.date, or (year, month[, day])",
    ))
}

fn parse_date_str(raw: &str) -> PyResult<NaiveDate> {
    let s = raw.trim();
    if let Ok(d) = NaiveDate::parse_from_str(s, "%Y-%m-%d") {
        return Ok(d);
    }
    if let Ok(d) = NaiveDate::parse_from_str(s, "%Y-%m") {
        return Ok(d);
    }
    if s.len() == 7 && s.as_bytes().get(4) == Some(&b'-') {
        if let Ok(d) = NaiveDate::parse_from_str(&format!("{s}-01"), "%Y-%m-%d") {
            return Ok(d);
        }
    }
    Err(PyValueError::new_err(format!(
        "could not parse date '{s}' (use YYYY-MM-DD or YYYY-MM)"
    )))
}

fn iso(date: NaiveDate) -> String {
    date.format("%Y-%m-%d").to_string()
}

fn parse_inflows(obj: &Bound<'_, PyAny>) -> PyResult<HashMap<NaiveDate, f64>> {
    let mut map = HashMap::new();
    if let Ok(dict) = obj.downcast::<PyDict>() {
        for (key, value) in dict.iter() {
            let volume: f64 = value.extract()?;
            *map.entry(parse_date(&key)?).or_insert(0.0) += volume;
        }
        return Ok(map);
    }
    let iter = obj.try_iter().map_err(|_| {
        PyValueError::new_err(
            "inflows must be a dict of date → volume or a sequence of (date, volume)",
        )
    })?;
    for item in iter {
        let item = item?;
        if let Ok((date_obj, volume)) = item.extract::<(Bound<'_, PyAny>, f64)>() {
            *map.entry(parse_date(&date_obj)?).or_insert(0.0) += volume;
            continue;
        }
        if let Ok(mapping) = item.downcast::<PyDict>() {
            let date_obj = mapping
                .get_item("date")?
                .or_else(|| mapping.get_item("month").ok().flatten())
                .ok_or_else(|| PyValueError::new_err("inflow row missing 'date'"))?;
            let volume: f64 = mapping
                .get_item("volume")?
                .or_else(|| mapping.get_item("value").ok().flatten())
                .ok_or_else(|| PyValueError::new_err("inflow row missing 'volume'"))?
                .extract()?;
            *map.entry(parse_date(&date_obj)?).or_insert(0.0) += volume;
            continue;
        }
        return Err(PyValueError::new_err(
            "each inflow row must be (date, volume) or {date, volume}",
        ));
    }
    Ok(map)
}

fn monthly_rows(rows: Vec<(NaiveDate, f64)>) -> Vec<(String, f64)> {
    rows.into_iter().map(|(d, v)| (iso(d), v)).collect()
}

fn field_i32(item: &Bound<'_, PyAny>, names: &[&str]) -> PyResult<i32> {
    if let Ok(mapping) = item.downcast::<PyDict>() {
        for name in names {
            if let Ok(Some(v)) = mapping.get_item(*name) {
                return v.extract();
            }
        }
    } else {
        for name in names {
            if item.hasattr(*name)? {
                return item.getattr(*name)?.extract();
            }
        }
    }
    Err(PyValueError::new_err(format!(
        "URF row missing integer field ({})",
        names.join("/")
    )))
}

fn field_f64(item: &Bound<'_, PyAny>, names: &[&str]) -> PyResult<f64> {
    if let Ok(mapping) = item.downcast::<PyDict>() {
        for name in names {
            if let Ok(Some(v)) = mapping.get_item(*name) {
                return v.extract();
            }
        }
    } else {
        for name in names {
            if item.hasattr(*name)? {
                return item.getattr(*name)?.extract();
            }
        }
    }
    Err(PyValueError::new_err(format!(
        "URF row missing numeric field ({})",
        names.join("/")
    )))
}

fn parse_monthly_urf(obj: &Bound<'_, PyAny>) -> PyResult<Vec<UrfValue>> {
    let mut rows = Vec::new();
    let iter = obj
        .try_iter()
        .map_err(|_| PyValueError::new_err("urf must be a sequence of (month, reach, urf_val)"))?;
    for item in iter {
        let item = item?;
        if let Ok((month, reach, urf_val)) = item.extract::<(i32, i32, f64)>() {
            rows.push(UrfValue::new(month, reach, urf_val));
            continue;
        }
        rows.push(UrfValue::new(
            field_i32(&item, &["month"])?,
            field_i32(&item, &["reach"])?,
            field_f64(&item, &["urf_val", "value", "coefficient"])?,
        ));
    }
    Ok(rows)
}

fn parse_daily_urf(obj: &Bound<'_, PyAny>) -> PyResult<Vec<DailyUrfValue>> {
    let mut rows = Vec::new();
    let iter = obj
        .try_iter()
        .map_err(|_| PyValueError::new_err("urf must be a sequence of (day, reach, urf_val)"))?;
    for item in iter {
        let item = item?;
        if let Ok((day, reach, urf_val)) = item.extract::<(i32, i32, f64)>() {
            rows.push(DailyUrfValue::new(day, reach, urf_val));
            continue;
        }
        rows.push(DailyUrfValue::new(
            field_i32(&item, &["day"])?,
            field_i32(&item, &["reach"])?,
            field_f64(&item, &["urf_val", "value", "coefficient"])?,
        ));
    }
    Ok(rows)
}

fn lagged_to_py(py: Python<'_>, result: LaggedUrfResult) -> PyResult<Bound<'_, PyDict>> {
    let outer = PyDict::new(py);
    for (reach, dates) in result {
        let inner = PyDict::new(py);
        for (date, volume) in dates {
            inner.set_item(iso(date), volume)?;
        }
        outer.set_item(reach, inner)?;
    }
    Ok(outer)
}

fn lagged_from_py(obj: &Bound<'_, PyAny>) -> PyResult<LaggedUrfResult> {
    let dict = obj
        .downcast::<PyDict>()
        .map_err(|_| PyValueError::new_err("lagged URF result must be {reach: {date: volume}}"))?;
    let mut result = HashMap::new();
    for (reach_key, inner_obj) in dict.iter() {
        let reach: i32 = reach_key.extract()?;
        let inner_dict = inner_obj
            .downcast::<PyDict>()
            .map_err(|_| PyValueError::new_err("each reach map must be {date: volume}"))?;
        let mut inner = HashMap::new();
        for (date_key, volume) in inner_dict.iter() {
            inner.insert(parse_date(&date_key)?, volume.extract()?);
        }
        result.insert(reach, inner);
    }
    Ok(result)
}

fn mass_to_py(py: Python<'_>, masses: HashMap<i32, f64>) -> PyResult<Bound<'_, PyDict>> {
    let out = PyDict::new(py);
    for (reach, mass) in masses {
        out.set_item(reach, mass)?;
    }
    Ok(out)
}

/// Monthly Jenkins SDF lagging. Positive volumes = pumping (depletion);
/// negative = recharge (accretion). Farmers Pawnee recharge uses the
/// negative-volume path with a mapped sdf (days).
///
/// `days_per_month` is accepted for Rust-API compatibility and is unused;
/// pulse lengths are calendar month lengths.
#[pyfunction]
#[pyo3(name = "lag_sdf", signature = (inflows, sdf, total_months, days_per_month=30.42))]
fn py_lag_sdf(
    inflows: Bound<'_, PyAny>,
    sdf: u32,
    total_months: u32,
    days_per_month: f64,
) -> PyResult<Vec<(String, f64)>> {
    let map = parse_inflows(&inflows)?;
    Ok(monthly_rows(calculate_streamflow_depletion_sdf(
        &map,
        sdf,
        days_per_month,
        total_months,
    )))
}

#[pyfunction]
#[pyo3(
    name = "lag_glover_infinite",
    signature = (inflows, distance_ft, specific_yield, transmissivity_ft2_per_day, total_months, days_per_month=30.42)
)]
fn py_lag_glover_infinite(
    inflows: Bound<'_, PyAny>,
    distance_ft: f64,
    specific_yield: f64,
    transmissivity_ft2_per_day: f64,
    total_months: usize,
    days_per_month: f64,
) -> PyResult<Vec<(String, f64)>> {
    let map = parse_inflows(&inflows)?;
    Ok(monthly_rows(calculate_streamflow_depletion_infinite(
        &map,
        distance_ft,
        specific_yield,
        transmissivity_ft2_per_day,
        days_per_month,
        total_months,
    )))
}

#[pyfunction]
#[pyo3(
    name = "lag_glover_alluvial",
    signature = (inflows, distance_ft, aquifer_width_ft, specific_yield, transmissivity_ft2_per_day, total_months, days_per_month=30.42)
)]
fn py_lag_glover_alluvial(
    inflows: Bound<'_, PyAny>,
    distance_ft: f64,
    aquifer_width_ft: f64,
    specific_yield: f64,
    transmissivity_ft2_per_day: f64,
    total_months: usize,
    days_per_month: f64,
) -> PyResult<Vec<(String, f64)>> {
    let map = parse_inflows(&inflows)?;
    Ok(monthly_rows(calculate_streamflow_depletion_alluvial(
        &map,
        distance_ft,
        aquifer_width_ft,
        specific_yield,
        transmissivity_ft2_per_day,
        days_per_month,
        total_months,
    )))
}

#[pyfunction]
#[pyo3(name = "urf_lagging")]
fn py_urf_lagging<'py>(
    py: Python<'py>,
    usage: Bound<'py, PyAny>,
    urf: Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyDict>> {
    let usage_map = parse_inflows(&usage)?;
    let rows = parse_monthly_urf(&urf)?;
    lagged_to_py(py, urf_lagging(&usage_map, rows))
}

#[pyfunction]
#[pyo3(name = "urf_lagging_daily")]
fn py_urf_lagging_daily<'py>(
    py: Python<'py>,
    usage: Bound<'py, PyAny>,
    urf: Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyDict>> {
    let usage_map = parse_inflows(&usage)?;
    let rows = parse_daily_urf(&urf)?;
    lagged_to_py(py, urf_lagging_daily(&usage_map, rows))
}

#[pyfunction]
#[pyo3(name = "combined_urf_results")]
fn py_combined_urf_results(values: Bound<'_, PyAny>) -> PyResult<Vec<(String, f64)>> {
    Ok(monthly_rows(combined_urf_results(lagged_from_py(&values)?)))
}

#[pyfunction]
#[pyo3(name = "combined_urf_results_dated")]
fn py_combined_urf_results_dated(values: Bound<'_, PyAny>) -> PyResult<Vec<(String, f64)>> {
    Ok(monthly_rows(combined_urf_results_dated(lagged_from_py(
        &values,
    )?)))
}

#[pyfunction]
#[pyo3(name = "urf_mass_by_reach")]
fn py_urf_mass_by_reach<'py>(
    py: Python<'py>,
    urf: Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyDict>> {
    mass_to_py(py, urf_mass_by_reach(&parse_monthly_urf(&urf)?))
}

#[pyfunction]
#[pyo3(name = "daily_urf_mass_by_reach")]
fn py_daily_urf_mass_by_reach<'py>(
    py: Python<'py>,
    urf: Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyDict>> {
    mass_to_py(py, daily_urf_mass_by_reach(&parse_daily_urf(&urf)?))
}

#[pyfunction]
#[pyo3(name = "disaggregate_monthly_urf_uniform")]
fn py_disaggregate_monthly_urf_uniform(
    monthly: Bound<'_, PyAny>,
    stress_month: Bound<'_, PyAny>,
) -> PyResult<Vec<(i32, i32, f64)>> {
    let daily =
        disaggregate_monthly_urf_uniform(&parse_monthly_urf(&monthly)?, parse_date(&stress_month)?);
    Ok(daily
        .into_iter()
        .map(|row| (row.day, row.reach, row.urf_val))
        .collect())
}

#[pyfunction]
#[pyo3(name = "sdf_from_glover")]
fn py_sdf_from_glover(
    distance_ft: f64,
    specific_yield: f64,
    transmissivity_ft2_per_day: f64,
) -> f64 {
    sdf_from_glover(distance_ft, specific_yield, transmissivity_ft2_per_day)
}

#[pyfunction]
#[pyo3(name = "gpd_per_ft_to_ft2_per_day")]
fn py_gpd_per_ft_to_ft2_per_day(transmissivity_gpd_per_ft: f64) -> f64 {
    gpd_per_ft_to_ft2_per_day(transmissivity_gpd_per_ft)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(py_lag_sdf, module)?)?;
    module.add_function(wrap_pyfunction!(py_lag_glover_infinite, module)?)?;
    module.add_function(wrap_pyfunction!(py_lag_glover_alluvial, module)?)?;
    module.add_function(wrap_pyfunction!(py_urf_lagging, module)?)?;
    module.add_function(wrap_pyfunction!(py_urf_lagging_daily, module)?)?;
    module.add_function(wrap_pyfunction!(py_combined_urf_results, module)?)?;
    module.add_function(wrap_pyfunction!(py_combined_urf_results_dated, module)?)?;
    module.add_function(wrap_pyfunction!(py_urf_mass_by_reach, module)?)?;
    module.add_function(wrap_pyfunction!(py_daily_urf_mass_by_reach, module)?)?;
    module.add_function(wrap_pyfunction!(
        py_disaggregate_monthly_urf_uniform,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(py_sdf_from_glover, module)?)?;
    module.add_function(wrap_pyfunction!(py_gpd_per_ft_to_ft2_per_day, module)?)?;
    // Rust-name aliases so pyOWW / lagging-api style call sites map 1:1.
    module.add(
        "calculate_streamflow_depletion_sdf",
        module.getattr("lag_sdf")?,
    )?;
    module.add(
        "calculate_streamflow_depletion_infinite",
        module.getattr("lag_glover_infinite")?,
    )?;
    module.add(
        "calculate_streamflow_depletion_alluvial",
        module.getattr("lag_glover_alluvial")?,
    )?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
