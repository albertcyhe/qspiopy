from pathlib import Path

import pandas as pd
import pytest

from scripts.summarize_equivalence import main as summarize_main


def test_summarize_equivalence_generates_tables(tmp_path):
    validation_dir = tmp_path / "validation"
    validation_dir.mkdir()

    metrics_df = pd.DataFrame(
        [
            {
                "scenario": "example1_control",
                "observable": "tumour_volume_l",
                "max_fractional_error": 1e-7,
                "relative_rmse": 5e-8,
                "rmse": 1e-9,
                "r2": 0.99,
                "delta_auc": 0.0,
            },
            {
                "scenario": "example1_control",
                "observable": "pd1_occupancy",
                "max_fractional_error": 2e-7,
                "relative_rmse": 1e-7,
                "rmse": 1e-9,
                "r2": 0.98,
                "delta_auc": 0.0,
            },
        ]
    )
    metrics_df.to_csv(validation_dir / "metrics.csv", index=False)

    event_df = pd.DataFrame(
        [
            {
                "scenario": "example1_control",
                "event_index": 1,
                "time_fire": 1.5,
                "time_trigger": 1.0,
                "delay": 0.5,
                "type": "delayed",
                "assignments": "A=0",
            },
            {
                "scenario": "example1_control",
                "event_index": 2,
                "time_fire": 1.0,
                "time_trigger": 1.0,
                "delay": 0.0,
                "type": "immediate",
                "assignments": "B=1",
            },
        ]
    )
    event_df.to_csv(validation_dir / "events_example1_control_python.csv", index=False)

    output_dir = tmp_path / "analysis"
    args = [
        "--validation-dir",
        str(validation_dir),
        "--snapshot-root",
        "artifacts/matlab_frozen_model",
        "--output-dir",
        str(output_dir),
    ]
    summarize_main(args)

    alignment_metrics = pd.read_csv(output_dir / "alignment_metrics.csv")
    assert not alignment_metrics.empty

    event_logs = pd.read_csv(output_dir / "event_logs.csv")
    assert "time_residual" in event_logs.columns

    event_summary = pd.read_csv(output_dir / "event_summary.csv")
    assert "max_time_residual" in event_summary.columns
    row = event_summary.iloc[0]
    assert pytest.approx(row["max_time_residual"], rel=1e-9) == 0.5
    assert pytest.approx(row["median_time_residual"], rel=1e-9) == 0.25
    assert pytest.approx(row["max_delay_residual"], rel=1e-9) == 0.0
