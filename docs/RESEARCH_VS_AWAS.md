# Research vs IDS AWAS vs this crate

South Platte augmentation accounting needs a lagged stream response that
can be defended from the published analytical methods, not from a
particular Windows build of IDS AWAS. This note records:

1. The formulas actually used (authors, year, equation).
2. What IDS AWAS / SEODEP implement.
3. What this crate implements after the review.
4. Where Windows AWAS will disagree with a research-correct result.

AWAS output is **not** treated as an oracle. Tests lock the published
identities (Jenkins’ 28% definition, Glover–Balmer `erfc`, Glover/Hantush
`4 i²erfc`, Jenkins superposition).

## Named methods

### Glover–Balmer rate (infinite / semi-infinite aquifer)

Glover, R.E., and Balmer, G.G., 1954, River depletion resulting from
pumping a well near a river: *Eos, Transactions, American Geophysical
Union*, 35(3), 468–470.

Derived from Theis (1941). For steady pumping at rate `Q` beginning at
`t = 0`:

```text
q(t) / Q = erfc(u),     u = √(a² S / (4 T t)) = a / √(4 T t / S)
```

- `q` — stream-depletion **rate** [L³/T]
- `Q` — net pumping (or recharge) rate [L³/T]
- `a` — perpendicular well-to-stream distance [L]
- `S` — specific yield / storativity [–]
- `T` — transmissivity [L²/T]
- `t` — time since the stress began [T]

Assumptions (Glover and Balmer 1954; Jenkins 1968): homogeneous isotropic
semi-infinite aquifer; straight fully penetrating stream with no streambed
resistance; constant `T`; instantaneous release from storage; Dupuit flow.

### Glover / Hantush volume

The rate equation is integrable in closed form (Glover 1960, *Ground
water–surface water relationships*, CSU CER60REG45; Hantush 1964,
*Hydraulics of wells*; Jenkins 1968, USGS TWRI 4-D1, curve B / eqs. 7–10):

```text
v(t) / (Q t) = 4 i²erfc(u) = (1 + 2u²) erfc(u) − (2u / √π) exp(−u²)
```

`v` is the **cumulative** stream-volume change. This is the quantity
monthly acre-foot accounting needs. Differencing `erfc(u)` at integer
days and treating the difference as a volume is only an Euler
approximation of `∫ q dt`.

### Jenkins SDF

Jenkins, C.T., 1968, Techniques for computing rate and volume of stream
depletion by wells: *Ground Water*, 6(2), 37–46.

Also: Jenkins, C.T., 1968 (rev. 1977), Computation of rate and volume of
stream depletion by wells: USGS Techniques of Water-Resources
Investigations, Book 4, Chapter D1.

```text
sdf ≡ time at which v = 0.28 Q t          [T]
```

In an idealized aquifer `sdf = a² S / T`, because `4 i²erfc(1/2) ≈ 0.27986`.
In a real valley the sdf is often a **mapped effective** value (Hurr and
others, South Platte / Arkansas SDF maps) that already folds in some
boundary and heterogeneity effects.

Substituting `sdf` into Glover:

```text
q/Q = erfc(√(sdf / 4t))
v/(Q t) = 4 i²erfc(√(sdf / 4t))
```

Jenkins is explicit that recharge is the same problem with opposite sign
(“recharging” / “accretion” for “pumping” / “depletion”).

### Residual effects and intermittent stress

Jenkins TWRI 4-D1: after a pulse of duration `t_p`, assume pumping
continues and superimpose an equal opposite stress starting at `t_p`.

```text
q_pulse(t) = Q [ erfc(u(t)) − erfc(u(t − t_p)) ]     (t > t_p)
v_pulse(t) = Q t · 4i²erfc(u(t)) − Q (t − t_p) · 4i²erfc(u(t − t_p))
```

Monthly acre-feet in `[t_a, t_b]` = `v_pulse(t_b) − v_pulse(t_a)`.
Several months superpose linearly.

### Alluvial (strip) aquifer — image wells

Stream at `x = 0` (constant head), impermeable valley wall at `x = W`,
well at `x = a`, `0 < a < W`.

Glover (1960, 1966, 1977 *Transient Ground Water Hydraulics*); Schroeder,
D.R., 1987, Analytical stream depletion model, Colorado DWR (SEODEP);
McWhorter and Sunada (1977, p. 127); Miller, Durnford, Halstead,
Altenhofen, and Flory (2007), *Ground Water* 45(4); Contor (2011), IWRRI
201101.

Glover already contains the image across the stream. The remaining series
(the series SEODEP/AWAS code actually loops) is:

```text
q/Q = Σ_{n=0}^∞ (−1)^n [ erfc(u(2nW + a)) + erfc(u(2(n+1)W − a)) ]
```

and the same coefficients on `4 i²erfc` for `v/(Q t)`.

`W` is **stream to wall** (aquifer width). AWAS field `m_w` / SEODEP `W`
use this meaning. A single extra well at `2W − a` is only the first pair;
the infinite series is required for a no-flow wall that remains a no-flow
wall after the stream images are added.

A **perpendicular** no-flow boundary (AWAS `NO_FLOW` / SEODEP option 3,
field `m_b`) is a different geometry (stream-segment images). This crate
does not implement that option.

### URF

Unit response functions are tabulated monthly fractions, not an
analytical kernel. This crate’s `urf_lagging` is a discrete convolution
and was not part of the Glover/SDF review.

## Comparison matrix

| Topic | Research (source of truth) | IDS AWAS / SEODEP | This crate |
| --- | --- | --- | --- |
| Rate kernel | `q/Q = erfc(u)` — Glover and Balmer 1954; Jenkins 1968 eq. (5) / curve A | `CalcErrorFunc` / BASIC 5490: series for `erfc(u)` stored as `BQQ` | `libm::erfc` |
| Volume kernel | `v/(Q t) = 4 i²erfc(u)` — Glover 1960; Hantush 1964; Jenkins 1968 eqs. (7)–(10) / curve B | Same closed form as `SVV` in BASIC 5630 / `SDFdata.cpp` | Same closed form (`four_i2erfc`) |
| SDF argument | `u = √(sdf / 4t)` | `u = sqrt(m_sdf / (4 * GetFact * t))` | `sqrt(sdf / (4 t))` with `t` in days |
| Glover argument | `u = a / √(4 T t / S)` | `u = m_dxx / sqrt(4 * m_tr * fact * t / (m_s * 7.48x))` — `m_tr` is GPD/ft | Caller passes `T` in **ft²/day**; helper `gpd_per_ft_to_ft2_per_day` uses 7.48051945 |
| Alluvial `W` | Stream-to-wall aquifer width | `m_w` / `W`: “distance from the stream to the parallel impermeable boundary” | `aquifer_width` — same meaning. Previous comments that called this well-to-wall were wrong. |
| Alluvial series | Infinite image pairs, signs `(+ + − − + + …)` | SEODEP BASIC 3360–3560 (entire stream). AWAS C++ has a later “Kenny” rewrite that only builds the kernel from the first non-zero `Q` month and uses `n` in place of `t` in `GetFact(x)` | Same pair sequence as SEODEP BASIC; pair-wise convergence (no `u > 2.9` kill) |
| No-flow (perp.) | Contor 2011 / stream-segment images | AWAS `NO_FLOW` + `CalcStreamQ` Simpson integration | **Not implemented** |
| Superposition | Jenkins residual construction on **volume** | `vd += SVV * Q * fact * N * dela` then subtract an opposite pulse at shutoff (BASIC 3600–3760) | `v_on(t) − v_on(t − t_p)` per monthly pulse, then `v(end) − v(start)` |
| Time origin | `t = 0` at the start of the pulse; `erfc(∞) = 0` | First step is `t = dela` (never evaluates `t = 0`) | `t ≤ 0` → fraction 0; first month uses the calendar length of that month |
| Time step | Any consistent `T`. Monthly filing: rectangular pulse over the calendar month | Monthly mode: `GetFact` = 30.41667 (`365/12`) or actual days. Daily mode: `dela` days | Calendar days in the stress month (January = 31, February = 28/29, …) |
| Units in / out | Any consistent set. Jenkins example: acre-ft/day and acre-ft | Pumping often GPM; volumes converted with `1440/325851` ac-ft (C++) or `1440/325900` (BASIC) | Monthly **volumes** in, monthly **volumes** out (acre-feet if that is what you put in) |
| `T` conversion | `1 ft³ = 7.48 gal` (Jenkins 1968) | BASIC `7.481`; AWAS C++ `7.48051945` | `7.48051945` in `GALLONS_PER_CUBIC_FOOT` |
| Pumping vs recharge | Identical except sign (Jenkins TWRI 4-D1 p. 2) | Site type sets `m_pumpingSign` (`WELL = +1`, recharge = `−1`) | Signed monthly volumes: `+` pumping / depletion, `−` recharge / accretion |
| `erfc` cutoff | None (tail is well defined) | BASIC: `U > 2.9` → 0. AWAS `CalcErrorFunc`: `U > 3.9` → 0. `CalcErrorFuncNew` (Kenny): continued-fraction for `U ≥ 2` | No cutoff; `libm::erfc` |
| Failed kernel | — | `m_qq` / `m_vv` initialized to **1.0** (“instantaneous effect if the calculation fails”) | Fractions are 0 at `t = 0` and computed otherwise |
| Output filter | Report the signed monthly volume, including small and negative | GUI / print path varies | Every month in the window is returned, including zeros and accretions |
| URF | Tabulated unit response | AWAS `URF` path convolves `m_URF_DataCalc` | Unchanged discrete convolution |

## What was wrong in this crate before the review

1. **Rate-differencing used as volume.** The code formed
   `Q · [erfc(u(t)) − erfc(u(t−1))]` and summed it as acre-feet. That is
   Jenkins’ *rate* superposition with a 1-day Euler step, not
   `v = Q t · 4 i²erfc(u)`. For a January 100 acre-foot pulse at
   `sdf = 265 d` the published volume in January is **0.944 acre-feet**;
   the old kernel produced **0.768 acre-feet** (~19% low that month). Later
   months were closer (a few percent).
2. **Recharge dropped.** `if pumping_rate <= 0.0 { continue; }` plus
   `create_results_vector` aborting on negative totals. Farmers Pawnee is
   ditch recharge with no wells — the previous code returned nothing
   useful.
3. **`+1` day shift and `t = 0` slot.** The daily loop started at
   `time_step = 0` (`erfc(∞) = 0`) and assigned the increment to
   `date + index + 1`, sliding response by an extra day.
4. **Alluvial `W` documented as well-to-wall.** The loop was already the
   SEODEP stream-to-wall series (`2W − a`); the README formula
   `erfc(d) + erfc(2W − d)` described only the first pair and used the
   wrong name for `W`.
5. **`u > 2.9` kill on the alluvial path only.** SEODEP does this for
   `erfc`; it is an AWAS/SEODEP numerical shortcut, not part of Glover.
6. **Self-referential tests** locked the old rate-differenced numbers.

## Remaining AWAS deviations (Windows AWAS will disagree)

These are places where a research-correct result **should** differ from
IDS AWAS. They are not crate bugs.

| AWAS / SEODEP behavior | Research | Typical effect |
| --- | --- | --- |
| Monthly `t` in units of 30.41667-day “months” when average-days is on | Calendar month lengths (31, 28, …) | First-month volume and the lag of February/March in leap vs common years |
| `erfc` / `4i²erfc` forced to 0 for `u > 2.9` (BASIC) or `u > 3.9` (C++ `CalcErrorFunc`) | Full tail | Very early time or very large `a` / sdf: AWAS is slightly later |
| Series `erfc` (`1e-8` BASIC, `1e-14` C++) vs `libm` | Either is fine at filing precision | Sub-0.001 acre-foot |
| `7.481` (BASIC) vs `7.48051945` (C++) vs Jenkins’ `7.48` | State the factor | ~0.01% in `T`, smaller in `u` |
| `325900` gal/ac-ft (BASIC) vs `325851` (C++) | Not used here (we never convert GPM) | Only if you compare AWAS GPM prints to acre-feet |
| `m_qq`/`m_vv` default **1.0** | 0 at `t = 0` | If transmissivity or sdf is missing, AWAS can report 100% instantaneous depletion |
| Comment in AWAS C++: “error in line 4810” on the **segment** alluvial path | Entire-stream series is the SEODEP loop we match | Segment-of-stream runs in AWAS are a different, known-fragile path |
| “Kenny” alluvial rewrite (`SDFdata.cpp` ~3680–3724): kernel rebuilt from the first non-zero `Q`, `GetFact(x)` with `n` instead of elapsed `t` | One kernel vs elapsed time since each pulse | Alluvial monthly AWAS can drift from SEODEP BASIC and from this crate |
| `CalcErrorFuncNew` continued fraction for `u ≥ 2` | Same `erfc` everywhere | Depends which function pointer the GUI passes |
| Stream-segment Simpson (`CalcStreamQ`) | Not used for whole-stream Glover/SDF | Irrelevant to Farmers Pawnee / plan-wide Glover |
| URF and “effective sdf” sharing the infinite `erfc` path | SDF is just Glover with mapped sdf | Agrees when both use the volume formula |

### Worked numbers (research, not AWAS)

January 2025 pulse of **100 acre-feet**, calendar days, Jenkins volume.

| Month | SDF = 265 d | Glover infinite `a=4000 ft`, `S=0.2`, `T=261800 GPD/ft` | Same + alluvial `W=8000 ft` |
| --- | ---: | ---: | ---: |
| 2025-01 | 0.944405 | 9.231024 | 9.234067 |
| 2025-02 | 7.169295 | 20.792231 | 21.070681 |
| 2025-03 | 10.023970 | 13.118033 | 14.755124 |
| 2025-04 | 7.801918 | 7.592177 | 10.312423 |
| 2025-05 | 6.279747 | 5.342934 | 8.462012 |
| 2025-06 | 4.833422 | 3.804862 | 6.632475 |

Jenkins TWRI 4-D1 sample: `a = 3660 ft`, `T/S = 134000 ft²/day` →
`sdf = 100 d`. After 35 days of `Q = 10` acre-ft/day,
`q/Q = erfc(√(100/140)) ≈ 0.2320`,
`v = 10 × 35 × 4i²erfc(…) ≈ 33.80` acre-feet.

Identity: `4 i²erfc(1/2) ≈ 0.2798588938` (the 28% in the sdf definition).

## Filing notes for ARI

- **Farmers Pawnee** (ditch recharge, SDF): call
  `calculate_streamflow_depletion_sdf` with **negative** monthly recharge
  volumes. Mapped sdf values already include a great deal of valley
  geometry; do **not** also apply the alluvial image series on top of a
  mapped sdf unless the sdf has been “boundary-adjusted” as in Miller et
  al. (2007).
- **Other ARI Glover plans:** `calculate_streamflow_depletion_infinite`
  or `_alluvial`. Pass `T` in ft²/day. For alluvial, `aquifer_width` is
  stream-to-wall.
- This crate is a **library**. There is still no command-line front end
  in this repository.
- Do not treat a difference versus Windows AWAS as a crate regression
  until it has been traced to one of the rows above.

## References

- Theis, C.V., 1941, The effect of a well on the flow of a nearby stream:
  *Eos*, 22(3), 734–738.
- Glover, R.E., and Balmer, G.G., 1954, River depletion resulting from
  pumping a well near a river: *Eos*, 35(3), 468–470.
- Glover, R.E., 1960, Ground water–surface water relationships: Colorado
  State University paper CER60REG45.
- Glover, R.E., 1964, Ground-water movement: USBR Engineering Monograph 31.
- Hantush, M.S., 1964, Hydraulics of wells, in Chow, V.T., ed.,
  *Advances in Hydroscience*, v. 1.
- Jenkins, C.T., 1968a, Techniques for computing rate and volume of
  stream depletion by wells: *Ground Water*, 6(2), 37–46.
- Jenkins, C.T., 1968 (TWRI 4-D1), Computation of rate and volume of
  stream depletion by wells: USGS.
- Schroeder, D.R., 1987, Analytical stream depletion model: Colorado
  Division of Water Resources (SEODEP; source in the AWAS tree as
  `SEODEP2.BAS`).
- McWhorter, D.B., and Sunada, D.K., 1977, *Ground-water hydrology and
  hydraulics*, p. 127 (image pattern).
- Miller, C.D., Durnford, D., Halstead, M.R., Altenhofen, J., and Flory,
  V., 2007, Stream depletion in alluvial valleys using the SDF
  semianalytical model: *Ground Water*, 45(4), 506–514.
- Contor, B.A., 2011, Adaptation of the Glover/Balmer/Jenkins analytical
  stream-depletion methods for no-flow and recharge boundaries: IWRRI
  Technical Completion Report 201101.
- Gautschi, W., 1964, Error function and Fresnel integrals, in
  Abramowitz and Stegun, *Handbook of mathematical functions*, NBS AMS 55.
- IDS AWAS: Integrated Decision Support Group, Colorado State
  University; C++ port in Longitude103/AWAS (`SDFdata.cpp`).
