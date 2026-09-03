# Unit response functions (URF)

This crate’s URF path is a **tabulated discrete convolution**, not Glover or
SDF. The series is an input (from a numerical model, a mapped table, or a
pre-computed kernel). Tests lock the convolution identities below, not IDS
AWAS output.

## Named methods

### Maddock algebraic technological function

Maddock, T., 1972, Algebraic technological function from a simulation
model: *Water Resources Research*, 8(1), 129–134.

Maddock, T., 1974, The operation of a stream-aquifer system under
stochastic demands: *Water Resources Research*, 10(1), 1–10.

Drawdown or stream exchange is the convolution of pumping with a
pre-computed response. Maddock (1974) treats the stream as a constant-head
boundary and recovers stream volume from storage-vs-pumpage mass balance.

### Morel-Seytoux discrete kernel

Morel-Seytoux, H.J., and Daly, C.J., 1975, A discrete kernel generator
for stream-aquifer studies: *Water Resources Research*, 11(2), 253–260.

Illangasekare, T., and Morel-Seytoux, H.J., 1982, Stream-aquifer
influence coefficients as tools for simulation and management: *Water
Resources Research*, 18(1), 168–176.

The Green’s function is sampled at discrete points in time and space.
Once the kernels exist, any linear response (reach return flow in month
*n*, etc.) is an explicit function of the pumping/recharge sequence:

```text
effect(n) = Σ_{ν}  δ(n − ν + 1) · stress(ν)
```

`δ(1)` is the response during the stress period itself.

### Unit hydrograph / fraction function

Dreizin and Haimes (1977) call the stream-exchange kernels **fraction
functions**: the fraction of a unit stress that appears as
stream–aquifer exchange in each later period. If the stream is the sole
source/sink and the series is carried to infinity, the fractions sum to 1.
A truncated model tail or other sinks (ET, bounding aquifers) make the
sum smaller. **Rescaling a truncated URF to 1 invents mass** and is not
done here.

Jenkins (1968) SDF/Glover curves are one way to *generate* a URF (the
incremental `Δv / V` of a unit monthly pulse). They are not themselves a
URF. Using `q/Q = erfc(u)` as monthly URF coefficients is a mixing error
(rate vs incremental volume).

## What AWAS / SEODEP implement

SEODEP BASIC has no URF option. IDS AWAS C++ (`SDFdata.cpp`,
`m_bi == URF`) does:

```text
vd[j+x] += URF[x] * Q[j] * pumpingSign * GetFact(j) * dela
```

then accumulates `vd` for a cumulative printout. Incremental monthly
volume is `URF[x] * (rate × days)` — a volume times a fraction.

- `URF[0]` is the stress timestep (0-based array).
- Site type sets sign (`WELL = +1`, recharge = −1).
- Horizon is truncated at `min(npa, URF.size())`.
- Monthly vs daily arrays; optional “use monthly URF for daily.”
- No reach split; one series per site.
- No interpolation of the tabulated values.
- No check that the series sums to 1.
- `GetFact * dela` is required because AWAS `Q` is a **rate**.

## What this crate implements

| Topic | Research | AWAS | This crate |
| --- | --- | --- | --- |
| Operation | Discrete convolution of stress with `δ` | Same, `URF[x] * volume` | `usage * urf_val` |
| Time origin | `δ(1)` in the stress period | `URF[0]` in period `j` | `month = 1` in the stress month |
| Indexing | 1-based period in the papers | 0-based `vector` | `UrfValue.month` is **1-based**; `from_lag(0, …)` for AWAS arrays |
| Gaps | Missing period ⇒ `δ = 0` | Array is dense | Missing months 1…N are zeros (not collapsed) |
| Superposition | Linear | Linear | Linear; same-month dates are summed |
| Recharge | Opposite sign | `m_pumpingSign` | Negative `usage` |
| Units | Consistent volume × fraction | Rate × days × fraction | Caller passes **volume** (acre-feet) |
| Normalization | Sum = 1 only if complete unit hydrograph | Not enforced | `urf_mass_by_reach`; no rescale |
| Interpolation | None (discrete kernel) | None (optional monthly-as-daily copy) | None |
| Reaches | Influence coeffs can be per reach | One series / site | `reach` key; `combined_urf_results` sums them |
| SDF / Glover mix | Generate URF from `Δv/V`, never from `q/Q` | Separate `EFFECTIVE_SDF` / Glover paths | Separate modules; test locks `q/Q ≠ Δv/V` |

## Bugs found in the previous URF code

1. **`month` was only a sort key.** Lag was `enumerate()` after sort. A
   table with month 1 and month 3 (no month 2) applied month 3 at a
   one-month lag. Research and AWAS treat that as a two-month lag with a
   zero in between.
2. **Duplicate `(reach, month)` rows** became two consecutive lags
   instead of one summed coefficient.
3. **Usage day-of-month was kept.** 15 January + 1 month is 15 February,
   so January-1 and January-15 stresses did not combine and could miss a
   monthly join with Glover/SDF output (which keys the 1st).
4. **`combined_urf_results` docs were wrong.** They claimed a 0.001
   cutoff and a stop-on-negative (copied from the old Glover helper).
   The code never did that; docs now match the code. Accretion and small
   volumes are kept.
5. **No mass check.** Callers could not see whether a table summed to 1.

Convolution itself (first bin = stress month, linear superposition,
signed usage) was already the Maddock / AWAS construction. Existing
dense-table tests still pass.

## Remaining risks

- **This crate does not generate URFs.** A mapped or MODFLOW-derived
  table is only as good as the model and the truncation. A sum of 0.92
  means 8% of the stress never appears in the lagged series.
- **Do not stack URF on Glover or SDF** for the same stress. Do not load
  `erfc(u)` or `4 i²erfc(u)` into `UrfValue` unless you have converted
  them to *incremental volume fractions of a rectangular monthly pulse*.
- **No daily URF and no interpolation.** AWAS can store a daily series
  or reuse monthly values on a daily step (`m_useMonthlyURFForDaily`).
  Spreading a monthly URF across days here would be an unstated
  assumption.
- **AWAS rate units.** If you compare to Windows AWAS, convert AWAS `Q`
  to a monthly volume (`Q × GetFact × dela`) before comparing to this
  crate’s acre-foot usage.
- **Reach partition.** If several reaches are meant to be a partition of
  one stream, check `Σ_reach mass ≈ 1`, not each reach separately.
- **`month < 1` is ignored.** Import AWAS arrays with
  `UrfValue::from_lag(i, reach, val)`.
- **Horizon.** AWAS drops URF tail past the pumping-record length. This
  crate applies the full series, so a long URF will write months after
  the last stress (correct discrete-kernel behavior).

## Worked example (locked in tests)

Unit hydrograph `(0.20, 0.50, 0.30)` on one reach, 100 acre-feet in
January 2025:

| Month | Effect (ac-ft) |
| --- | ---: |
| 2025-01 | 20 |
| 2025-02 | 50 |
| 2025-03 | 30 |
| **sum** | **100** |

Two successive 100 acre-foot months with `(0.25, 0.75)`:

| Month | Effect |
| --- | ---: |
| Jan | 25 |
| Feb | 75 + 25 = 100 |
| Mar | 75 |
| **sum** | **200** |
