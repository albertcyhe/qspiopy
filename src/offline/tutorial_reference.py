"""High-fidelity reference simulations orchestrated via MATLAB SimBiology."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

from .frozen_model import simulate_frozen_model
from .tutorial_sim import ScenarioResult


@dataclass(frozen=True)
class ReferenceSettings:
    """Compatibility shim for the legacy API."""

    rtol: float = 1e-7
    atol: float = 1e-10
    method: str = "BDF"
    dense_output: bool = True


def _infer_script(parameter_paths: Iterable[Path | str]) -> str:
    for raw in parameter_paths:
        name = Path(raw).name
        if "example2" in name.lower():
            return "example2"
    return "example1"


def simulate_reference(
    parameter_paths: Iterable[Path | str],
    *,
    days: float = 400.0,
    therapy: Literal["none", "anti_pd1"] = "anti_pd1",
    settings: ReferenceSettings | None = None,  # noqa: ARG001 - preserved for API compatibility
    seed: int | None = None,
    emit_diagnostics: bool = False,
    run_label: str | None = None,
    event_log: list[dict[str, object]] | None = None,
    rtol_override: float | None = None,
    atol_override: float | None = None,
) -> ScenarioResult:
    """Simulate the SimBiology tutorial model using MATLAB."""

    script = _infer_script(parameter_paths)
    return simulate_frozen_model(
        script,
        days=days,
        therapy=therapy,
        seed=seed,
        emit_diagnostics=emit_diagnostics,
        run_label=run_label,
        event_log=event_log,
        rtol_override=rtol_override,
        atol_override=atol_override,
    )
