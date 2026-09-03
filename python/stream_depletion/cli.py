"""Thin CLI around the PyO3 bindings.

Farmers Pawnee monthly recharge (replaces IDS AWAS SDF on Windows)::

    stream-depletion \\
        --inflows examples/farmers_pawnee_july_recharge.csv \\
        --method sdf --sdf 265 --as-recharge \\
        --total-months 12 --output july_fp_accretions.csv

Also: ``python -m stream_depletion ...``
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence, TextIO

from stream_depletion import (
    combined_urf_results,
    combined_urf_results_dated,
    gpd_per_ft_to_ft2_per_day,
    lag_glover_alluvial,
    lag_glover_infinite,
    lag_sdf,
    urf_lagging,
    urf_lagging_daily,
)

VOLUME_ALIASES = ("volume", "inflow", "inflows", "af", "acre_feet", "acre-feet", "value")
RECHARGE_ALIASES = ("recharge", "accretion")
PUMPING_ALIASES = ("pumping", "depletion")


def _open_text(path: str) -> TextIO:
    return Path(path).open(newline="", encoding="utf-8")


def _strip_comments(text: str) -> str:
    kept = [
        line
        for line in text.splitlines(keepends=True)
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return "".join(kept)


def _norm_fields(names: Iterable[str] | None) -> list[str]:
    return [name.strip().lower().replace(" ", "_") for name in (names or [])]


def _get(row: Mapping[str, str], *names: str) -> str | None:
    lower = {key.strip().lower().replace(" ", "_"): value for key, value in row.items() if key}
    for name in names:
        if name in lower and lower[name] != "":
            return lower[name]
    return None


def read_inflows(path: str) -> dict[str, float]:
    """Read ``date,volume`` or ``year,month[,day],volume`` CSV.

    Duplicate dates are summed (same as the Rust monthly collapse).
    """
    raw = Path(path).read_text(encoding="utf-8")
    reader = csv.DictReader(io.StringIO(_strip_comments(raw)))
    fields = _norm_fields(reader.fieldnames)
    if not fields:
        raise SystemExit(f"{path}: missing header row")

    inflows: dict[str, float] = {}
    for row in reader:
        date = _get(row, "date", "month_start", "month")
        if date is None:
            year = _get(row, "year")
            month = _get(row, "month")
            if year is None or month is None:
                raise SystemExit(
                    f"{path}: each row needs 'date' or 'year'+'month' (got {list(row)})"
                )
            day = _get(row, "day") or "1"
            date = f"{int(year):04d}-{int(month):02d}-{int(day):02d}"
        volume_s = None
        for alias in VOLUME_ALIASES + RECHARGE_ALIASES + PUMPING_ALIASES:
            volume_s = _get(row, alias)
            if volume_s is not None:
                break
        if volume_s is None:
            raise SystemExit(f"{path}: no volume column in {list(row)}")
        inflows[date] = inflows.get(date, 0.0) + float(volume_s)
    if not inflows:
        raise SystemExit(f"{path}: no inflow rows")
    return inflows


def read_urf(path: str, timestep: str) -> list[tuple[int, int, float]]:
    raw = Path(path).read_text(encoding="utf-8")
    reader = csv.DictReader(io.StringIO(_strip_comments(raw)))
    period_names = ("day",) if timestep == "daily" else ("month",)
    rows: list[tuple[int, int, float]] = []
    for row in reader:
        period_s = _get(row, *period_names)
        reach_s = _get(row, "reach")
        value_s = _get(row, "urf_val", "value", "coefficient", "urf")
        if period_s is None or reach_s is None or value_s is None:
            raise SystemExit(
                f"{path}: URF rows need {period_names[0]},reach,urf_val (got {list(row)})"
            )
        rows.append((int(period_s), int(reach_s), float(value_s)))
    if not rows:
        raise SystemExit(f"{path}: no URF rows")
    return rows


def _write_pairs(out: TextIO, rows: Sequence[tuple[str, float]]) -> None:
    writer = csv.writer(out)
    writer.writerow(["date", "volume"])
    for date, volume in rows:
        writer.writerow([date, f"{volume:.10g}"])


def _write_reach_map(out: TextIO, lagged: Mapping[int, Mapping[str, float]]) -> None:
    writer = csv.writer(out)
    writer.writerow(["date", "reach", "volume"])
    triples: list[tuple[str, int, float]] = []
    for reach, series in lagged.items():
        for date, volume in series.items():
            triples.append((date, int(reach), volume))
    triples.sort(key=lambda item: (item[0], item[1]))
    for date, reach, volume in triples:
        writer.writerow([date, reach, f"{volume:.10g}"])


def _glover_transmissivity(args: argparse.Namespace) -> float:
    if args.transmissivity is not None:
        return args.transmissivity
    if args.gpd_per_ft is not None:
        return gpd_per_ft_to_ft2_per_day(args.gpd_per_ft)
    raise SystemExit("--transmissivity (ft²/day) or --gpd-per-ft is required for Glover")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="stream-depletion",
        description=(
            "Lag pumping (+) or recharge (−) volumes to monthly stream "
            "depletion / accretion. Farmers Pawnee: --method sdf with "
            "negative recharge (or --as-recharge). Replaces Windows IDS AWAS."
        ),
    )
    parser.add_argument(
        "--inflows",
        "-i",
        required=True,
        help="CSV of monthly/daily volumes (date,volume). Positive=pumping, negative=recharge.",
    )
    parser.add_argument(
        "--method",
        "-m",
        required=True,
        choices=("sdf", "glover", "glover-alluvial", "urf"),
        help="Lagging kernel. Farmers Pawnee uses sdf; other ARI plans may use glover.",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output CSV (date,volume). Default: stdout.",
    )
    parser.add_argument(
        "--sdf",
        type=int,
        help="Stream depletion factor in days (required for --method sdf).",
    )
    parser.add_argument(
        "--total-months",
        type=int,
        default=24,
        help="Calendar months to report from the first stress (SDF/Glover). Default: 24.",
    )
    parser.add_argument(
        "--days-per-month",
        type=float,
        default=30.42,
        help="Accepted for API compatibility; pulse lengths use calendar days.",
    )
    parser.add_argument("--distance", type=float, help="Well-to-stream distance a, feet (Glover).")
    parser.add_argument("--specific-yield", type=float, help="Storativity / specific yield S (Glover).")
    parser.add_argument(
        "--transmissivity",
        type=float,
        help="Transmissivity T in ft²/day (Glover).",
    )
    parser.add_argument(
        "--gpd-per-ft",
        type=float,
        help="Transmissivity in GPD/ft; converted with the AWAS C++ factor 7.48051945.",
    )
    parser.add_argument(
        "--aquifer-width",
        type=float,
        help="Stream-to-wall aquifer width W, feet (--method glover-alluvial).",
    )
    parser.add_argument(
        "--urf",
        help="URF CSV: month,reach,urf_val (monthly) or day,reach,urf_val (daily).",
    )
    parser.add_argument(
        "--timestep",
        choices=("monthly", "daily"),
        default="monthly",
        help="URF stress step. Daily SDF is not in the crate; use --method urf --timestep daily.",
    )
    parser.add_argument(
        "--as-recharge",
        action="store_true",
        help="Negate inflow volumes (use when the CSV stores recharge as positive acre-feet).",
    )
    parser.add_argument(
        "--by-reach",
        action="store_true",
        help="With --method urf, write date,reach,volume instead of combining reaches.",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    inflows = read_inflows(args.inflows)
    if args.as_recharge:
        inflows = {date: -volume for date, volume in inflows.items()}

    out = Path(args.output).open("w", newline="", encoding="utf-8") if args.output else sys.stdout
    try:
        if args.method == "sdf":
            if args.sdf is None:
                raise SystemExit("--sdf DAYS is required for --method sdf")
            rows = lag_sdf(inflows, args.sdf, args.total_months, args.days_per_month)
            _write_pairs(out, rows)
        elif args.method == "glover":
            if args.distance is None or args.specific_yield is None:
                raise SystemExit("--distance and --specific-yield are required for --method glover")
            rows = lag_glover_infinite(
                inflows,
                args.distance,
                args.specific_yield,
                _glover_transmissivity(args),
                args.total_months,
                args.days_per_month,
            )
            _write_pairs(out, rows)
        elif args.method == "glover-alluvial":
            if args.distance is None or args.specific_yield is None or args.aquifer_width is None:
                raise SystemExit(
                    "--distance, --specific-yield, and --aquifer-width are required "
                    "for --method glover-alluvial"
                )
            rows = lag_glover_alluvial(
                inflows,
                args.distance,
                args.aquifer_width,
                args.specific_yield,
                _glover_transmissivity(args),
                args.total_months,
                args.days_per_month,
            )
            _write_pairs(out, rows)
        else:
            if not args.urf:
                raise SystemExit("--urf FILE is required for --method urf")
            kernel = read_urf(args.urf, args.timestep)
            if args.timestep == "daily":
                lagged = urf_lagging_daily(inflows, kernel)
                combined = combined_urf_results_dated(lagged)
            else:
                lagged = urf_lagging(inflows, kernel)
                combined = combined_urf_results(lagged)
            if args.by_reach:
                _write_reach_map(out, lagged)
            else:
                _write_pairs(out, combined)
    finally:
        if args.output:
            out.close()
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
