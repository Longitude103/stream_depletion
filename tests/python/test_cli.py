"""CLI: Farmers Pawnee monthly recharge and the URF month-gap case."""

from pathlib import Path

import pytest

from stream_depletion.cli import main

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "farmers_pawnee_july_recharge.csv"
URF_GAP = ROOT / "examples" / "urf_month_gap.csv"


def test_cli_farmers_pawnee_july_sdf(tmp_path):
    out = tmp_path / "july_fp_accretions.csv"
    rc = main(
        [
            "--inflows",
            str(EXAMPLE),
            "--method",
            "sdf",
            "--sdf",
            "265",
            "--as-recharge",
            "--total-months",
            "6",
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "date,volume"
    first = lines[1].split(",")
    assert first[0] == "2025-07-01"
    assert float(first[1]) < 0.0
    assert abs(float(first[1]) + 0.944405) < 5e-6
    assert len(lines) == 7  # header + 6 months


def test_cli_urf_month_gap(tmp_path):
    usage = tmp_path / "usage.csv"
    usage.write_text("date,volume\n2025-01-01,100\n", encoding="utf-8")
    out = tmp_path / "lagged.csv"
    rc = main(
        [
            "--inflows",
            str(usage),
            "--method",
            "urf",
            "--urf",
            str(URF_GAP),
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    body = [line for line in out.read_text(encoding="utf-8").splitlines() if line]
    assert body[1].startswith("2025-01-01,")
    assert "2025-02-01" not in out.read_text(encoding="utf-8")
    assert body[2].startswith("2025-03-01,")


def test_cli_module_help(capsys):
    try:
        main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0
    captured = capsys.readouterr().out
    assert "--inflows" in captured
    assert "--method" in captured
    assert "sdf" in captured


def test_cli_average_month_sdf(tmp_path):
    out = tmp_path / "average.csv"
    rc = main(
        [
            "--inflows",
            str(EXAMPLE),
            "--method",
            "sdf",
            "--sdf",
            "265",
            "--as-recharge",
            "--total-months",
            "2",
            "--month-length",
            "average",
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert abs(float(lines[1].split(",")[1]) + 0.890056) < 5e-6
    assert abs(float(lines[2].split(",")[1]) + 7.94116) < 5e-6


def test_cli_no_flow_matches_glover(tmp_path):
    usage = tmp_path / "usage.csv"
    usage.write_text("date,volume\n2025-01-01,100\n", encoding="utf-8")
    outputs = {}
    for method, extra in (
        ("glover", []),
        ("glover-no-flow", ["--boundary-distance", "12000"]),
    ):
        out = tmp_path / f"{method}.csv"
        rc = main(
            [
                "--inflows",
                str(usage),
                "--method",
                method,
                "--distance",
                "4000",
                "--specific-yield",
                "0.2",
                "--gpd-per-ft",
                "261800",
                "--total-months",
                "6",
                "--output",
                str(out),
                *extra,
            ]
        )
        assert rc == 0
        outputs[method] = out.read_text(encoding="utf-8")
    assert outputs["glover"] == outputs["glover-no-flow"]


def test_cli_bad_parameter_exits_cleanly_and_keeps_existing_output(tmp_path):
    out = tmp_path / "keep.csv"
    out.write_text("previous results\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="sdf_days"):
        main(
            [
                "--inflows",
                str(EXAMPLE),
                "--method",
                "sdf",
                "--sdf",
                "0",
                "--output",
                str(out),
            ]
        )
    assert out.read_text(encoding="utf-8") == "previous results\n"


@pytest.mark.parametrize(
    ("volume", "extra"),
    [
        ("abc", ["--sdf", "265"]),
        ("100", ["--sdf", "-1"]),
        ("100", ["--sdf", "265", "--total-months", "-1"]),
    ],
)
def test_cli_other_bad_input_exits_cleanly(tmp_path, volume, extra):
    usage = tmp_path / "usage.csv"
    usage.write_text(f"date,volume\n2025-01-01,{volume}\n", encoding="utf-8")
    out = tmp_path / "keep.csv"
    out.write_text("previous results\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="error:"):
        main(["--inflows", str(usage), "--method", "sdf", "--output", str(out), *extra])
    assert out.read_text(encoding="utf-8") == "previous results\n"
