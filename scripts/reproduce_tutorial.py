"""Run the reduced-order tutorial simulations from the command line."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd

from src.offline import tutorial_sim


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reproduce tutorial scenarios using the Python surrogate")
    parser.add_argument(
        "--parameters",
        type=Path,
        action="append",
        required=True,
        help="Path to a JSON parameter catalogue (can be provided multiple times)",
    )
    parser.add_argument(
        "--days",
        type=float,
        default=400.0,
        help="Duration of the simulation horizon in days (default: 400)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional destination CSV.  When omitted the result is printed to stdout",
    )
    parser.add_argument(
        "--therapy",
        choices=["none", "anti_pd1"],
        default="anti_pd1",
        help="Select which checkpoint scenario to run (default: anti_pd1)",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    result = tutorial_sim.simulate_tutorial(
        args.parameters,
        days=args.days,
        therapy=args.therapy,
    )
    frame = result.to_frame()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.output, index=False)
    else:
        pd.set_option("display.max_rows", 20)
        print(frame)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

