"""Command line interface for the deterministic QSP-I/O surrogate.

The tool mirrors the MATLAB-orientated workflow by accepting a tidy
``EPOCH_FEATURES`` table, validating it against the offline schema, and
executing :func:`src.offline.qsp_io.run_qsp` to generate the synthetic
trajectories used in contract tests.

The CLI is intentionally strict:
* Missing columns, inconsistent mechanism assignments, empty inputs, or unit
  mismatches are treated as fatal errors.
* A compact mechanism registry keeps the plumbing honest; the registry is
  loaded from ``parameters/MECHANISM_REGISTRY.csv`` when available and falls
  back to the feature schema otherwise.
* Every successful run records reproducibility metadata (seed, git hash,
  source version override, input SHA-256) alongside the simulated values so the
  contract tests can diff the CSV artefacts end-to-end.

Typical usage::

    python -m scripts.qsp_io_cli --input artifacts/EPOCH_FEATURES.csv \
        --output artifacts/QSP_OUT.csv --seed 42

The ``--dry-run`` flag performs the full validation/simulation cycle without
persisting the output artefact; this is useful when integrating the CLI in
continuous-integration pipelines where the calling harness is responsible for
capturing the generated frame.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, Mapping, MutableMapping

import pandas as pd

from src.offline import qsp_io

LOGGER = logging.getLogger("qsp_io_cli")


_DEFAULT_INPUT = Path("artifacts") / "EPOCH_FEATURES.csv"
_DEFAULT_OUTPUT = Path("artifacts") / "QSP_OUT.csv"

_REQUIRED_COLUMNS = {"mechanism_id", "feature", "time_h", "value", "units"}


def _configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(levelname)s:%(name)s:%(message)s",
    )


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the QSP-I/O surrogate over a CSV input")
    parser.add_argument(
        "--input",
        type=Path,
        default=_DEFAULT_INPUT,
        help=f"Input CSV containing EPOCH feature rows (default: {_DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_DEFAULT_OUTPUT,
        help=f"Destination CSV for the simulated trajectories (default: {_DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Deterministic seed forwarded to src.offline.qsp_io.run_qsp (default: 0)",
    )
    parser.add_argument(
        "--source-version",
        type=str,
        default=None,
        help="Optional override for the source_version metadata attached to the output",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and simulate without writing the output CSV",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Emit debug logging (useful when auditing validation failures)",
    )
    return parser.parse_args(argv)


def _read_input_frame(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input CSV {path} does not exist")
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError("Input CSV does not contain any rows")
    missing = _REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)!r}")
    return frame


def _load_mechanism_registry(schema: Mapping[str, qsp_io.FeatureSchema]) -> set[str]:
    """Load the mechanism registry from disk or fall back to the feature schema."""

    candidate_paths = []
    env_path = os.environ.get("QSPI_MECHANISM_REGISTRY")
    if env_path:
        candidate_paths.append(Path(env_path))
    repo_root = Path(__file__).resolve().parents[1]
    candidate_paths.append(repo_root / "parameters" / "MECHANISM_REGISTRY.csv")

    for path in candidate_paths:
        if not path:
            continue
        if path.is_file():
            registry = pd.read_csv(path)
            if "mechanism_id" not in registry.columns:
                raise ValueError(
                    f"Mechanism registry {path} must contain a 'mechanism_id' column"
                )
            mechanisms = {
                str(value).strip()
                for value in registry["mechanism_id"].dropna().tolist()
                if str(value).strip()
            }
            if not mechanisms:
                raise ValueError(f"Mechanism registry {path} is empty")
            LOGGER.info(
                "Loaded %d mechanism identifiers from %s", len(mechanisms), path
            )
            return mechanisms
    LOGGER.info(
        "Mechanism registry not found; defaulting to feature identifiers from the schema"
    )
    return {str(key) for key in schema.keys()}


def _normalise_input(
    frame: pd.DataFrame,
    mechanism_registry: set[str],
) -> tuple[pd.DataFrame, Dict[str, str]]:
    canonical = frame.copy()
    canonical["mechanism_id"] = canonical["mechanism_id"].astype(str).str.strip()
    canonical["feature"] = canonical["feature"].astype(str).str.strip()

    unknown = sorted(set(canonical["mechanism_id"]) - mechanism_registry)
    if unknown:
        raise ValueError(
            "Unknown mechanism_id values encountered: "
            + ", ".join(unknown)
        )

    mapping_series = canonical.groupby("feature")["mechanism_id"].nunique()
    inconsistent = mapping_series[mapping_series != 1]
    if not inconsistent.empty:
        problems = ", ".join(
            f"{feature} -> {count} ids" for feature, count in inconsistent.items()
        )
        raise ValueError(
            "Each feature must map to exactly one mechanism_id; violations: " + problems
        )
    feature_to_mechanism = (
        canonical.groupby("feature")["mechanism_id"].first().to_dict()
    )
    return canonical, feature_to_mechanism


def _log_unit_dictionary(schema: Mapping[str, qsp_io.FeatureSchema], features: Iterable[str]) -> None:
    units_dict: Dict[str, str] = {}
    for feature in sorted(set(features)):
        if feature in schema:
            units_dict[feature] = schema[feature].units
    if units_dict:
        LOGGER.info("Unit dictionary: %s", json.dumps(units_dict, sort_keys=True))


def _compute_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _build_output_frame(
    result: MutableMapping[str, object],
    feature_to_mechanism: Mapping[str, str],
    *,
    seed: int,
    input_sha256: str,
    source_version_override: str | None = None,
) -> pd.DataFrame:
    try:
        summary: Mapping[str, object] = result["summary"]  # type: ignore[index]
        summary_features: Mapping[str, Mapping[str, object]] = summary[  # type: ignore[index]
            "features"
        ]
    except KeyError as exc:  # pragma: no cover - defensive, should never trigger
        raise ValueError("run_qsp did not return the expected summary payload") from exc

    trajectory_records = []
    trajectories = result.get("trajectories", {})
    if not isinstance(trajectories, Mapping):
        raise TypeError("'trajectories' entry from run_qsp must be a mapping")

    for time_key, feature_map in trajectories.items():
        try:
            time_value = float(time_key)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid time index {time_key!r} in trajectories") from None
        if not isinstance(feature_map, Mapping):
            raise TypeError("Each trajectory entry must be a mapping of feature->value")
        for feature, value in feature_map.items():
            feature_name = str(feature)
            if feature_name not in feature_to_mechanism:
                raise ValueError(
                    f"Feature {feature_name!r} missing mechanism mapping in output construction"
                )
            if feature_name not in summary_features:
                raise ValueError(
                    f"Feature {feature_name!r} missing from summary metadata"
                )
            feature_summary = summary_features[feature_name]
            trajectory_records.append(
                {
                    "mechanism_id": feature_to_mechanism[feature_name],
                    "feature": feature_name,
                    "time_h": float(time_value),
                    "value": float(value),
                    "units": str(feature_summary.get("units", "")),
                    "transform": str(feature_summary.get("transform", "")),
                    "long_tail": bool(feature_summary.get("long_tail", False)),
                    "mean": float(feature_summary.get("mean", 0.0)),
                    "std": float(feature_summary.get("std", 0.0)),
                    "lo95": float(feature_summary.get("lo95", 0.0)),
                    "hi95": float(feature_summary.get("hi95", 0.0)),
                    "n_eff": float(feature_summary.get("n_eff", 0.0)),
                    "posterior_weight": float(
                        feature_summary.get("posterior_weight", summary.get("posterior_weight", 0.0))
                    ),
                }
            )

    if not trajectory_records:
        raise ValueError("run_qsp produced an empty trajectory payload")

    output_frame = pd.DataFrame.from_records(trajectory_records)
    output_frame.sort_values(["mechanism_id", "feature", "time_h"], inplace=True)
    output_frame.reset_index(drop=True, inplace=True)

    source_model = str(summary.get("source_model", "offline.qsp_io"))
    source_version = source_version_override or str(summary.get("source_version", "unknown"))

    output_frame["source_model"] = source_model
    output_frame["source_version"] = source_version
    output_frame["input_sha256"] = input_sha256
    output_frame["seed"] = int(seed)

    column_order = [
        "mechanism_id",
        "feature",
        "time_h",
        "value",
        "units",
        "transform",
        "long_tail",
        "mean",
        "std",
        "lo95",
        "hi95",
        "n_eff",
        "posterior_weight",
        "source_model",
        "source_version",
        "input_sha256",
        "seed",
    ]
    output_frame = output_frame[column_order]

    # Enforce dtypes explicitly for deterministic contract testing.
    float_columns = ["time_h", "value", "mean", "std", "lo95", "hi95", "n_eff", "posterior_weight"]
    for column in float_columns:
        output_frame[column] = output_frame[column].astype(float)
    output_frame["seed"] = output_frame["seed"].astype(int)
    output_frame["long_tail"] = output_frame["long_tail"].astype(bool)
    string_columns = [
        "mechanism_id",
        "feature",
        "units",
        "transform",
        "source_model",
        "source_version",
        "input_sha256",
    ]
    for column in string_columns:
        output_frame[column] = output_frame[column].astype(str)

    return output_frame


def _write_output(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    _configure_logging(args.verbose)

    LOGGER.info("Starting QSP-I/O CLI")
    LOGGER.info("Input CSV: %s", args.input)
    LOGGER.info("Output CSV: %s", args.output)
    LOGGER.info("Seed: %d", args.seed)

    try:
        schema = qsp_io._load_epoch_schema()  # pylint: disable=protected-access
        mechanism_registry = _load_mechanism_registry(schema)
        frame = _read_input_frame(args.input)
        canonical_frame, feature_to_mechanism = _normalise_input(frame, mechanism_registry)
        _log_unit_dictionary(schema, canonical_frame["feature"].tolist())
        input_sha = _compute_sha256(args.input)
        result = qsp_io.run_qsp(canonical_frame, seed=args.seed)
        output_frame = _build_output_frame(
            result,
            feature_to_mechanism,
            seed=args.seed,
            input_sha256=input_sha,
            source_version_override=args.source_version,
        )
        if args.dry_run:
            LOGGER.info("Dry run requested; skipping write step after successful validation")
            return 0
        _write_output(args.output, output_frame)
        LOGGER.info("Wrote %d rows to %s", len(output_frame), args.output)
        return 0
    except Exception as exc:  # pragma: no cover - exercised via CLI
        LOGGER.error("%s", exc)
        LOGGER.debug("Full exception", exc_info=True)
        return 1


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    sys.exit(main())
