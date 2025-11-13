"""Canonical context alias helpers."""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

AliasSpec = Tuple[str, Tuple[str, ...]]

ALIAS_GROUPS: Tuple[AliasSpec, ...] = (
    ("H_PD1_C1", ("H_PD1", "PD1_OCCUPANCY", "pd1_occupancy", "H_PD1_TOTAL")),
    ("V_T.T1", ("T", "T_total", "T_tumour", "CD8_T_cells")),
    ("V_T", ("tumor_volume_l", "tumour_volume_l", "Vtumour")),
    ("C_x", ("C_dead", "Cdead", "C_d", "C1_dead", "dead_cells")),
    ("syn_pd1_total", ("syn_T1_C1.PD1",)),
    ("syn_pdl1_total", ("syn_T1_C1.PDL1",)),
    ("syn_pdl2_total", ("syn_T1_C1.PDL2",)),
    ("syn_pd1_pdl1", ("syn_T1_C1.PD1_PDL1",)),
    ("syn_pd1_pdl2", ("syn_T1_C1.PD1_PDL2",)),
    ("syn_pd1_apd1", ("syn_T1_C1.PD1_aPD1",)),
    ("syn_pd1_apd1_pd1", ("syn_T1_C1.PD1_aPD1_PD1",)),
)


def inject_output_aliases(context: Dict[str, float]) -> None:
    """Populate canonical observable keys using known aliases."""

    for canonical, candidates in ALIAS_GROUPS:
        if canonical in context:
            continue
        for candidate in candidates:
            if candidate in context:
                context[canonical] = context[candidate]
                break


__all__ = ["inject_output_aliases"]
