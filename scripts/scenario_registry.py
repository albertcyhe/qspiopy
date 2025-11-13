"""Scenario registry for alignment tuning (A-series / B-series)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Sequence

from src.offline import DoseEntry
from src.offline.units import convert_amount

HOUR = 1.0
DAY = 24.0

PD1_TARGET_DEFAULT = "V_C.nivolumab"
CTLA4_TARGET_DEFAULT = "V_C.ipililumab"

DRUG_PROPERTIES = {
    PD1_TARGET_DEFAULT: ("nivolumab", 1.436e5),  # grams per mole
    CTLA4_TARGET_DEFAULT: ("ipilimumab", 1.486349e5),
}


@dataclass(frozen=True)
class Dose:
    """High-level dose description used by the scenario registry."""

    time_hours: float
    amount_mg: float
    drug: str
    target: str
    molecular_weight_g_per_mol: float
    interval_hours: float = 0.0
    repeat: int = 0


@dataclass(frozen=True)
class ScenarioSpec:
    """Scenario definition for extended alignment sweeps."""

    name: str
    label: str
    snapshot: str
    matlab_script: str
    days: float
    sample_interval_hours: float
    therapy: str
    doses: Sequence[Dose]
    context_outputs: Mapping[str, str]


def _build_repeat_series(interval_hours: float, occurrences: int, dose_mg: float, target: str) -> List[Dose]:
    drug, mw = _resolve_drug(target)
    return [
        Dose(time_hours=i * interval_hours, amount_mg=dose_mg, drug=drug, target=target, molecular_weight_g_per_mol=mw)
        for i in range(occurrences)
    ]


def _resolve_drug(target: str) -> tuple[str, float]:
    if target not in DRUG_PROPERTIES:
        raise ValueError(f"No drug metadata for target '{target}'")
    return DRUG_PROPERTIES[target]


def a_series(snapshot: str = "example1") -> List[ScenarioSpec]:
    """Return A1–A6 scenario specs."""
    return [
        ScenarioSpec(
            name="A1",
            label="200 mg Q3W ×4",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(21 * DAY, 4, 200.0, PD1_TARGET_DEFAULT),
            context_outputs=_pd1_outputs(),
        ),
        ScenarioSpec(
            name="A2",
            label="240 mg Q2W ×6",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(14 * DAY, 6, 240.0, PD1_TARGET_DEFAULT),
            context_outputs=_pd1_outputs(),
        ),
        ScenarioSpec(
            name="A3",
            label="400 mg Q6W ×2",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(42 * DAY, 2, 400.0, PD1_TARGET_DEFAULT),
            context_outputs=_pd1_outputs(),
        ),
        ScenarioSpec(
            name="A4",
            label="3 mg/kg (70 kg) Q2W ×6",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(14 * DAY, 6, 3.0 * 70.0, PD1_TARGET_DEFAULT),
            context_outputs=_pd1_outputs(),
        ),
        ScenarioSpec(
            name="A5",
            label="800 mg load + 200 mg Q3W",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=[
                Dose(
                    time_hours=0.0,
                    amount_mg=800.0,
                    drug=_resolve_drug(PD1_TARGET_DEFAULT)[0],
                    target=PD1_TARGET_DEFAULT,
                    molecular_weight_g_per_mol=_resolve_drug(PD1_TARGET_DEFAULT)[1],
                )
            ]
            + _build_repeat_series(21 * DAY, 3, 200.0, PD1_TARGET_DEFAULT),
            context_outputs=_pd1_outputs(),
        ),
        ScenarioSpec(
            name="A6",
            label="50 mg Q6W ×2",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(42 * DAY, 2, 50.0, PD1_TARGET_DEFAULT),
            context_outputs=_pd1_outputs(),
        ),
    ]


def microdose(snapshot: str = "example1") -> ScenarioSpec:
    """Return a 1 mg single-dose, 48 h scenario for diagnostics."""
    drug, mw = _resolve_drug(PD1_TARGET_DEFAULT)
    return ScenarioSpec(
        name="microdose_1mg_48h",
        label="1 mg single dose (48 h)",
        snapshot=snapshot,
        matlab_script="example1",
        days=2.0,
        sample_interval_hours=0.5,
        therapy="anti_pd1",
        doses=[
            Dose(
                time_hours=0.0,
                amount_mg=1.0,
                drug=drug,
                target=PD1_TARGET_DEFAULT,
                molecular_weight_g_per_mol=mw,
            )
        ],
        context_outputs=_pd1_outputs(),
    )


def b_series(snapshot: str = "example1") -> List[ScenarioSpec]:
    """Return combination therapy scenario specs (PD-1 + CTLA-4)."""
    return [
        ScenarioSpec(
            name="B1",
            label="PD-1 200 mg Q3W ×4 + CTLA-4 1 mg/kg Q6W ×2",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(21 * DAY, 4, 200.0, PD1_TARGET_DEFAULT)
            + _build_repeat_series(42 * DAY, 2, 70.0, CTLA4_TARGET_DEFAULT),
            context_outputs=_pd1_outputs({"ctla4_plasma_molar": "V_C.ipililumab"}),
        ),
        ScenarioSpec(
            name="B2",
            label="PD-1 200 mg Q3W ×4 + CTLA-4 delay 1 week",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(21 * DAY, 4, 200.0, PD1_TARGET_DEFAULT)
            + [
                Dose(
                    time_hours=7 * DAY,
                    amount_mg=70.0,
                    drug=_resolve_drug(CTLA4_TARGET_DEFAULT)[0],
                    target=CTLA4_TARGET_DEFAULT,
                    molecular_weight_g_per_mol=_resolve_drug(CTLA4_TARGET_DEFAULT)[1],
                ),
                Dose(
                    time_hours=49 * DAY,
                    amount_mg=70.0,
                    drug=_resolve_drug(CTLA4_TARGET_DEFAULT)[0],
                    target=CTLA4_TARGET_DEFAULT,
                    molecular_weight_g_per_mol=_resolve_drug(CTLA4_TARGET_DEFAULT)[1],
                ),
            ],
            context_outputs=_pd1_outputs({"ctla4_plasma_molar": "V_C.ipililumab"}),
        ),
        ScenarioSpec(
            name="B3",
            label="PD-1 240 mg Q2W ×6 + CTLA-4 3 mg/kg Q3W ×2",
            snapshot=snapshot,
            matlab_script="example1",
            days=84.0,
            sample_interval_hours=12.0,
            therapy="anti_pd1",
            doses=_build_repeat_series(14 * DAY, 6, 240.0, PD1_TARGET_DEFAULT)
            + _build_repeat_series(21 * DAY, 2, 3.0 * 70.0, CTLA4_TARGET_DEFAULT),
            context_outputs=_pd1_outputs({"ctla4_plasma_molar": "V_C.ipililumab"}),
        ),
    ]


def doses_to_entries(
    doses: Iterable[Dose],
    *,
    index_offset: int = 50_000,
) -> List[DoseEntry]:
    """Convert registry doses to DoseEntry objects (mg + MW → mol conversion)."""
    entries: List[DoseEntry] = []
    for idx, dose in enumerate(doses, start=1):
        amount_moles = convert_amount(dose.amount_mg, "milligram", dose.molecular_weight_g_per_mol)
        entries.append(
            DoseEntry(
                index=index_offset + idx,
                name=f"registry_{dose.drug}_{idx}",
                dose_type="RepeatDose",
                target=dose.target,
                amount=amount_moles,
                amount_units="mole",
                start_time=dose.time_hours / 24.0,
                interval=dose.interval_hours / 24.0 if dose.interval_hours else 0.0,
                repeat_count=dose.repeat,
                rate=None,
                rate_units="",
                duration=None,
                amount_mg=dose.amount_mg,
            )
        )
    return entries


__all__ = [
    "Dose",
    "ScenarioSpec",
    "a_series",
    "b_series",
    "microdose",
    "doses_to_entries",
    "HOUR",
    "DAY",
]
PD1_CONTEXT_OUTPUTS: Mapping[str, str] = {
    "drug_plasma_molar": "V_C.nivolumab",
    "drug_tumor_molar": "V_T.nivolumab",
    "pdl1_occupancy": "H_PD1_C1",
    "syn_pd1_total": "syn_T1_C1.PD1",
    "syn_pdl1_total": "syn_T1_C1.PDL1",
    "syn_pdl2_total": "syn_T1_C1.PDL2",
    "syn_pd1_pdl1": "syn_T1_C1.PD1_PDL1",
    "syn_pd1_pdl2": "syn_T1_C1.PD1_PDL2",
    "syn_pd1_apd1": "syn_T1_C1.PD1_aPD1",
    "syn_pd1_apd1_pd1": "syn_T1_C1.PD1_aPD1_PD1",
    "tcell_tumour": "V_T.T1",
    "tcell_treg": "V_T.T0",
    "tcell_ln": "V_LN.T1",
    "tcell_peripheral": "V_P.T1",
    "tcell_central": "V_C.T1",
    "tcell_total": "T_total",
    "tcell_total_ln": "T_total_LN",
    "tcell_kill_rate": "R_Tcell",
}


def _merge_context_outputs(*layers: Mapping[str, str]) -> Dict[str, str]:
    merged: Dict[str, str] = {}
    for layer in layers:
        if not layer:
            continue
        merged.update(layer)
    return merged


def _pd1_outputs(extra: Mapping[str, str] | None = None) -> Dict[str, str]:
    return _merge_context_outputs(extra or {}, PD1_CONTEXT_OUTPUTS)
