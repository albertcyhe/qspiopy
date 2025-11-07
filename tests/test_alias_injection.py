from __future__ import annotations

from src.offline.aliases import inject_output_aliases


def test_alias_injection_populates_missing_keys():
    ctx = {
        "H_PD1": 0.42,
        "T_total": 1500.0,
        "C_dead": 12.0,
        "tumor_volume_l": 1e-3,
    }
    inject_output_aliases(ctx)
    assert ctx["H_PD1_C1"] == 0.42
    assert ctx["V_T.T1"] == 1500.0
    assert ctx["C_x"] == 12.0
    assert ctx["V_T"] == 1e-3
