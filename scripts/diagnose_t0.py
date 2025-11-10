from __future__ import annotations

import argparse
import numpy as np

from src.offline.snapshot import load_frozen_model


def _inject_time(ctx: dict[str, float], t: float) -> None:
    ctx["time"] = t
    ctx["time_days"] = t
    ctx["t"] = t


def main(snapshot_name: str) -> None:
    model = load_frozen_model(snapshot_name)
    state = model.initial_state().astype(float)
    model.apply_initial_assignments_to_state(state)
    context = model.build_context_from_state(state.copy())
    _inject_time(context, 0.0)
    model.evaluate_repeated_assignments(context)
    model.apply_algebraic_rules(context, state, mutate=False)
    model.sync_state_from_context(context, state)

    print("=== t=0 trigger samples ===")
    for entry in model.events:
        try:
            raw = float(entry.trigger_compiled.evaluate(context))
        except Exception as exc:  # pragma: no cover - diagnostic
            print(f"[{entry.index:03d}] {entry.name:25s} expr=ERROR {exc}")
            continue
        bool_val = None
        if entry.trigger_boolean_compiled is not None:
            try:
                bool_val = bool(entry.trigger_boolean_compiled.evaluate(context))
            except Exception as exc:  # pragma: no cover - diagnostic
                bool_val = f"ERR({exc})"
        print(
            f"[{entry.index:03d}] {entry.name:25s} expr={raw:+.3e} bool={bool_val}"
        )

    print("\n=== RHS finite check ===")
    rhs = model.rhs(0.0, state.copy())
    finite = np.all(np.isfinite(rhs))
    print(f"all_finite={finite}")
    if not finite:
        bad_idx = np.where(~np.isfinite(rhs))[0]
        print(f"non-finite derivative indices: {bad_idx}")


if __name__ == "__main__":  # pragma: no cover - CLI utility
    parser = argparse.ArgumentParser(description="Diagnose t=0 trigger/RHS values.")
    parser.add_argument("snapshot", help="Snapshot name (e.g. example2)")
    args = parser.parse_args()
    main(args.snapshot)
