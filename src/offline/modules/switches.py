"""Lightweight module switches applied on top of frozen snapshots."""

from __future__ import annotations

import logging
import math
from collections import deque
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Literal, MutableMapping, Sequence, Tuple

from types import MethodType

from ..snapshot import FrozenModel

logger = logging.getLogger(__name__)

ModuleBlock = Callable[[MutableMapping[str, float]], None]

AVOGADRO = 6.02214076e23
LITERS_PER_CUBIC_MICROMETER = 1e-15
DEFAULT_SYN_DEPTH_UM = 1.15e-5
SECONDS_PER_DAY = 86400.0


def _resolve_synapse_depth(model: FrozenModel) -> float:
    """Return the synapse penetration depth (micrometers) for PD-1 bridge conversions."""

    param_keys = (
        "pd1_synapse_depth_um",
        "synapse_depth_um",
        "syn_depth_um",
    )
    parameters = getattr(model, "parameters", {}) or {}
    for key in param_keys:
        value = parameters.get(key)
        if value is None:
            continue
        try:
            depth = float(value)
        except (TypeError, ValueError):
            continue
        if depth > 0.0:
            return depth
    return DEFAULT_SYN_DEPTH_UM


def _molar_to_surface(concentration_molar: float, depth_um: float) -> float:
    """Convert mol/L into molecules per square micrometer for a thin synapse."""

    if concentration_molar == 0.0 or depth_um <= 0.0:
        return 0.0
    return concentration_molar * AVOGADRO * LITERS_PER_CUBIC_MICROMETER * depth_um


@dataclass(frozen=True)
class ModuleBlockSpec:
    """Descriptor for a runtime module block."""

    factory: Callable[[FrozenModel], ModuleBlock]
    disable_targets: Tuple[str, ...] = ()
    phase: Literal["pre", "post"] = "post"


def apply_parameter_overrides(model: FrozenModel, overrides: Dict[str, float] | None) -> FrozenModel:
    """Return the model after applying parameter overrides."""
    if not overrides:
        return model
    model.parameters = dict(model.parameters)
    model.parameters.update(overrides)
    return model


def disable_repeated_assignments(model: FrozenModel, targets: Iterable[str] | None) -> FrozenModel:
    """Remove repeated assignment rules targeting the provided symbols."""
    if not targets:
        return model
    block = {target.strip() for target in targets if target}
    if not block:
        return model
    model.repeated_assignment_order = [
        rule for rule in model.repeated_assignment_order if rule.target not in block
    ]
    return model


def pd1_bridge_block(model: FrozenModel) -> ModuleBlock:
    """Build a module block that projects PD-1 drug into the synapse."""

    weighted_pairs = (
        ("V_T.nivolumab", "gamma_T_nivolumab"),
        ("V_P.nivolumab", "gamma_P_nivolumab"),
        ("V_C.nivolumab", "gamma_C_nivolumab"),
        ("V_LN.nivolumab", "gamma_LN_nivolumab"),
    )

    syn_depth_um = _resolve_synapse_depth(model)

    def apply(context: MutableMapping[str, float]) -> None:
        accumulator = 0.0
        has_weight = False
        for species_key, gamma_key in weighted_pairs:
            conc = float(context.get(species_key, 0.0))
            gamma = float(context.get(gamma_key, 0.0))
            if conc or gamma:
                has_weight = True
            accumulator += gamma * conc
        if not has_weight:
            fallback = context.get("V_T.nivolumab")
            if fallback is None:
                fallback = context.get("V_C.nivolumab")
            if fallback is None:
                fallback = context.get("aPD1_concentration_molar")
            if fallback is None:
                fallback = context.get("aPD1", 0.0)
            accumulator = float(fallback)
        concentration = float(accumulator)
        surface_density = _molar_to_surface(concentration, syn_depth_um)
        context["aPD1_concentration_molar"] = concentration
        context["aPD1_surface_molecules_per_um2"] = surface_density
        context["aPD1"] = surface_density

    return apply


def tumour_geometry_block(model: FrozenModel) -> ModuleBlock:
    """Build a module block that refreshes tumour geometry observables."""

    def apply(context: MutableMapping[str, float]) -> None:
        volume_l = float(context.get("V_T", 0.0))
        if volume_l <= 0.0:
            cell_volume = float(context.get("cell_to_volume_factor_l", 0.0))
            if cell_volume <= 0.0:
                vol_cell_um3 = float(context.get("vol_cell", 0.0))
                if vol_cell_um3 > 0.0:
                    cell_volume = vol_cell_um3 * 1e-15
            if cell_volume > 0.0:
                cancer_cells = float(context.get("V_T.C1", context.get("C1", 0.0)))
                dead_cells = float(context.get("V_T.C_x", context.get("C_x", 0.0)))
                derived_volume = (cancer_cells + dead_cells) * cell_volume
                if derived_volume > 0.0:
                    volume_l = derived_volume
                    context["V_T"] = volume_l
        if volume_l > 0.0:
            context["tumour_volume_l"] = volume_l
            context["tumor_volume_l"] = volume_l
            radius_cm = ((3.0 * volume_l * 1e3) / (4.0 * math.pi)) ** (1.0 / 3.0)
            diameter_cm = 2.0 * radius_cm
            context["tumour_diameter_cm"] = diameter_cm
            context["tumor_diameter_cm"] = diameter_cm
            t_cells = float(context.get("V_T.T1", 0.0))
            density = 0.0
            if volume_l > 0.0:
                density = t_cells / (volume_l * 1e6)
            context["tcell_density_per_ul"] = density
        else:
            context.setdefault("tumour_volume_l", volume_l)
            context.setdefault("tumour_diameter_cm", 0.0)
            context.setdefault("tcell_density_per_ul", 0.0)

    return apply


def _sum_tumour_compartment(context: MutableMapping[str, float], kind: str) -> float:
    total = 0.0
    prefix = "V_T."
    for key, value in context.items():
        if not key.startswith(prefix):
            continue
        suffix = key[len(prefix) :]
        if kind == "live" and suffix.startswith("C") and not suffix.startswith("C_x"):
            total += float(value)
        elif kind == "dead" and suffix.startswith("C_x"):
            total += float(value)
        elif kind == "tcell" and suffix.startswith("T"):
            total += float(value)
    return total


def tumour_geometry_dynamic_block(model: FrozenModel) -> ModuleBlock:
    """Dynamic tumour geometry with pseudo-progression filtering."""

    parameters = model.parameters or {}
    dead_clearance = float(parameters.get("geometry_dead_clearance_per_day", 0.05))
    dead_swelling = float(parameters.get("geometry_dead_swelling_factor", 1.2))
    tcell_volume_factor = float(parameters.get("geometry_tcell_volume_factor", 1.0))
    geom_tau_days = float(parameters.get("geom_tau_days", 15.0))
    geom_kappa_occ = float(parameters.get("geom_kappa_occ_per_day", 0.0))
    w_live = float(parameters.get("geom_w_live", 1.0))
    w_tcell = float(parameters.get("geom_w_tcell", 1.0))
    w_dead = float(parameters.get("geom_w_dead", 1.0))

    last_time: float | None = None
    filtered_dead_volume_l: float | None = None
    smoothed_volume_l: float | None = None

    def apply(context: MutableMapping[str, float]) -> None:
        nonlocal last_time, filtered_dead_volume_l, smoothed_volume_l

        current_time = float(context.get("time_days", context.get("t", 0.0)))

        cell_volume_l = float(context.get("cell_to_volume_factor_l", 0.0))
        if cell_volume_l <= 0.0:
            raw_vol = float(context.get("vol_cell", 0.0))
            if raw_vol > 0.0:
                # Heuristic: legacy snapshots stored 2.5e3 µm^3; newer exports already convert to litres (~1e-12)
                if raw_vol > 1e-6:
                    cell_volume_l = raw_vol * 1e-15
                else:
                    cell_volume_l = raw_vol

        live_cells = _sum_tumour_compartment(context, "live") or float(context.get("C1", 0.0))
        dead_cells = _sum_tumour_compartment(context, "dead") or float(context.get("C_x", 0.0))
        t_cells = _sum_tumour_compartment(context, "tcell") or float(context.get("V_T.T1", 0.0))

        live_volume = live_cells * cell_volume_l
        tcell_volume = t_cells * tcell_volume_factor * cell_volume_l
        dead_volume_inst = dead_cells * dead_swelling * cell_volume_l

        if filtered_dead_volume_l is None:
            filtered_dead_volume_l = dead_volume_inst
        if last_time is None:
            last_time = current_time

        delta_t = max(0.0, current_time - last_time)
        if delta_t > 0.0 and dead_clearance > 0.0:
            alpha = 1.0 - math.exp(-dead_clearance * delta_t)
            filtered_dead_volume_l = filtered_dead_volume_l + alpha * (dead_volume_inst - filtered_dead_volume_l)
        else:
            filtered_dead_volume_l = dead_volume_inst
        baseline_volume = float(context.get("tumour_volume_l", context.get("V_T", 0.0)))
        target_volume = max(
            0.0,
            baseline_volume
            + w_live * live_volume
            + w_tcell * tcell_volume
            + w_dead * (filtered_dead_volume_l or 0.0),
        )
        if smoothed_volume_l is None:
            smoothed_volume_l = target_volume
        else:
            if geom_tau_days > 0.0 and delta_t > 0.0:
                alpha_v = 1.0 - math.exp(-delta_t / geom_tau_days)
                occ = float(context.get("pd1_occupancy", context.get("H_PD1_C1", 0.0)))
                smoothed_volume_l = smoothed_volume_l + alpha_v * (
                    target_volume - smoothed_volume_l - geom_kappa_occ * occ * smoothed_volume_l
                )
            else:
                smoothed_volume_l = target_volume
        smoothed_volume_l = max(smoothed_volume_l or 0.0, 0.0)
        context["tumour_volume_l"] = smoothed_volume_l
        context["tumor_volume_l"] = smoothed_volume_l

        if smoothed_volume_l > 0.0:
            radius_cm = ((3.0 * smoothed_volume_l * 1e3) / (4.0 * math.pi)) ** (1.0 / 3.0)
            diameter_cm = 2.0 * radius_cm
            context["tumour_diameter_cm"] = diameter_cm
            context["tumor_diameter_cm"] = diameter_cm
            context["tcell_density_per_ul"] = t_cells / (smoothed_volume_l * 1e6)
        else:
            context["tumour_diameter_cm"] = 0.0
            context["tumor_diameter_cm"] = 0.0
            context["tcell_density_per_ul"] = 0.0
        last_time = current_time

        # Diagnostics for --dump-flat-debug
        context["geom_live_volume_l"] = live_volume
        context["geom_dead_instant_volume_l"] = dead_volume_inst
        context["geom_dead_filtered_volume_l"] = filtered_dead_volume_l or 0.0
        context["geom_tcell_volume_l"] = tcell_volume
        context["geom_target_volume_l"] = target_volume
        context["geom_volume_smoothed_l"] = smoothed_volume_l
        context["geom_volume_delta_t"] = delta_t
        context["geom_volume_time_days"] = current_time

    return apply


def pd1_occupancy_filter_block(model: FrozenModel) -> ModuleBlock:
    """Apply a multi-pole filter to PD-1 occupancy to mimic observed kinetics."""

    parameters = model.parameters or {}
    tau1 = float(parameters.get("pd1_occ_tau1_days",
                                parameters.get("pd1_occ_tau_days", 1.0) or 1.0))
    tau2 = float(parameters.get("pd1_occ_tau2_days", tau1))
    delay_steps = max(int(parameters.get("pd1_occ_delay_steps", 0)), 0)
    hill_n = float(parameters.get("pd1_occ_hill_n", parameters.get("n_PD1", 2.0) or 2.0))
    pd1_50_base = float(
        parameters.get("pd1_occ_half_saturation", parameters.get("PD1_50", 20.0) or 20.0)
    )
    kappa_v = float(parameters.get("pd1_occ_volume_sensitivity", 0.0))
    volume_ref = float(parameters.get("pd1_occ_volume_reference_l", 0.014))

    last_time: float | None = None
    primary_state: float | None = None
    secondary_state: float | None = None
    delay_buffer: deque[float] | None = None

    def apply(context: MutableMapping[str, float]) -> None:
        nonlocal last_time, primary_state, secondary_state, delay_buffer

        t = float(context.get("time_days", context.get("t", 0.0)))
        a_surface = float(context.get("aPD1_surface_molecules_per_um2",
                                      context.get("aPD1", 0.0)))
        volume = max(float(context.get("tumour_volume_l",
                                       context.get("V_T", 0.0))), 0.0)
        pd1_50_eff = pd1_50_base
        if kappa_v:
            pd1_50_eff *= 1.0 + kappa_v * math.log1p(volume / max(volume_ref, 1e-9))
        pd1_50_eff = max(pd1_50_eff, 1e-9)
        x = max(a_surface / pd1_50_eff, 0.0)
        hill = 0.0
        if x > 0.0:
            x_pow = x ** hill_n
            hill = x_pow / (x_pow + 1.0)

        if primary_state is None:
            primary_state = hill
            secondary_state = hill
            last_time = t
        else:
            dt = max(0.0, t - (last_time or t))
            if tau1 > 0.0 and dt > 0.0:
                alpha1 = 1.0 - math.exp(-dt / tau1)
                primary_state = primary_state + alpha1 * (hill - primary_state)
            else:
                primary_state = hill
            if tau2 > 0.0 and dt > 0.0:
                alpha2 = 1.0 - math.exp(-dt / tau2)
                secondary_state = secondary_state + alpha2 * (primary_state - secondary_state)
            else:
                secondary_state = primary_state
            last_time = t

        filtered_value = secondary_state if secondary_state is not None else hill

        if delay_steps > 0:
            if delay_buffer is None:
                delay_buffer = deque(maxlen=delay_steps + 1)
                for _ in range(delay_steps + 1):
                    delay_buffer.append(filtered_value)
            else:
                delay_buffer.append(filtered_value)
            output_value = delay_buffer[0]
        else:
            output_value = filtered_value

        context["pd1_occupancy"] = output_value
        context["H_PD1_C1"] = output_value
        # Diagnostics for --dump-flat-debug
        context["pd1_filter_input"] = hill
        context["pd1_filter_primary"] = primary_state
        context["pd1_filter_secondary"] = secondary_state
        context["pd1_filter_output"] = output_value
        context["pd1_filter_tau1_days"] = tau1
        context["pd1_filter_tau2_days"] = tau2
        context["pd1_filter_delay_steps"] = delay_steps
        context["pd1_filter_pd1_50_eff"] = pd1_50_eff
        context["pd1_filter_time_days"] = t
        context["pd1_filter_volume_l"] = volume
        context["pd1_filter_surface_density"] = a_surface

    return apply


def alignment_driver_block(model: FrozenModel) -> ModuleBlock:
    """Alignment driver with explicit PD-1 and tumour-volume dynamics."""

    parameters = model.parameters or {}
    tau_pk = max(float(parameters.get("pd1_pk_tau_decay_days", 14.0)), 1e-6)
    dose_scale = float(parameters.get("pd1_pk_dose_scale", 1.0))
    conc_scale = float(parameters.get("pd1_pk_conc_scale", 1e-6))
    surface_scale = float(parameters.get("pd1_pk_surface_scale", 1e-6))

    raw_kon = float(parameters.get("kon_PD1_aPD1", parameters.get("kon_PD1_PDL1", 1.0) or 1.0))
    raw_koff = float(parameters.get("koff_PD1_aPD1", parameters.get("koff_PD1_PDL1", 0.1) or 0.1))
    kon_scale = float(parameters.get("pd1_occ_kon_scale", 1.0))
    koff_scale = float(parameters.get("pd1_occ_koff_scale", 1.0))
    kon_eff = raw_kon * kon_scale * SECONDS_PER_DAY
    koff_eff = raw_koff * koff_scale * SECONDS_PER_DAY
    k_int = float(parameters.get("pd1_occ_internalization_per_day", 0.0))
    occ_floor = float(parameters.get("pd1_occ_floor", 0.0))
    occ_ceiling = float(parameters.get("pd1_occ_ceiling", 0.15))

    grow_rate = float(parameters.get("geom_growth_per_day", 0.01))
    kill_coeff = float(parameters.get("geom_kill_per_cell_per_day", 1e-14))
    volume_cap = max(float(parameters.get("geom_volume_cap_l", 0.1)), 1e-6)
    min_volume_l = max(
        float(parameters.get("geom_min_volume_l", 0.0))
        or float(parameters.get("V_Tmin", 0.0)) * 1e-6,
        1e-9,
    )

    schedule_source = None
    schedule: List[Tuple[float, float]] = []
    dose_index = 0
    pk_state = 0.0
    occ_state = occ_floor
    volume_state: float | None = None
    last_time: float | None = None

    def refresh_schedule() -> None:
        nonlocal schedule_source, schedule, dose_index
        active = getattr(model, "_active_dose_schedule", None)
        if not active or active is schedule_source:
            return
        schedule_source = active
        schedule = [
            (float(item.time), float(item.amount if item.amount is not None else item.dose.amount or 0.0))
            for item in active
        ]
        schedule.sort(key=lambda pair: pair[0])
        dose_index = 0

    def _extract_volume_l(ctx: MutableMapping[str, float]) -> float:
        vt_l = float(ctx.get("tumour_volume_l", 0.0))
        if vt_l > 0.0:
            return vt_l
        vt_value = float(ctx.get("V_T", 0.0))
        if vt_value > 0.0:
            if vt_value > 100.0:
                return vt_value * 1e-6
            return vt_value
        return min_volume_l

    def apply(context: MutableMapping[str, float]) -> None:
        nonlocal pk_state, occ_state, volume_state, last_time, dose_index

        refresh_schedule()
        current_time = float(context.get("time_days", context.get("t", 0.0)))
        if last_time is None:
            last_time = current_time
            volume_state = _extract_volume_l(context)
        delta_t = max(0.0, current_time - (last_time or current_time))
        last_time = current_time

        if delta_t > 0.0:
            pk_state *= math.exp(-delta_t / tau_pk)
        while dose_index < len(schedule) and schedule[dose_index][0] <= current_time + 1e-12:
            pk_state += dose_scale * schedule[dose_index][1]
            dose_index += 1

        conc = pk_state * conc_scale
        surface_density = max(pk_state * surface_scale, 0.0)

        if delta_t > 0.0:
            d_occ = kon_eff * conc * (1.0 - occ_state) - (koff_eff + k_int) * occ_state
            occ_state = min(max(occ_state + d_occ * delta_t, occ_floor), occ_ceiling)

        if volume_state is None:
            volume_state = _extract_volume_l(context)
        eff_cells = max(float(context.get("V_T.T1", 0.0)), 0.0)
        if delta_t > 0.0:
            growth = grow_rate * volume_state * (1.0 - volume_state / volume_cap)
            kill = kill_coeff * eff_cells
            volume_state = max(volume_state + (growth - kill) * delta_t, min_volume_l)

        vt_microliter = volume_state * 1e6

        context["aPD1"] = conc
        context["aPD1_surface_molecules_per_um2"] = surface_density
        context["pd1_occupancy"] = occ_state
        context["H_PD1_C1"] = occ_state
        context["pd1_alignment_pk_state"] = pk_state
        context["pd1_alignment_concentration_M"] = conc
        context["pd1_alignment_volume_l"] = volume_state

        context["tumour_volume_l"] = volume_state
        context["tumor_volume_l"] = volume_state
        context["V_T"] = vt_microliter

        radius_cm = ((3.0 * volume_state * 1e3) / (4.0 * math.pi)) ** (1.0 / 3.0)
        context["tumour_diameter_cm"] = 2.0 * radius_cm
        context["tumor_diameter_cm"] = 2.0 * radius_cm

        total_tcells = float(context.get("V_T.T1", 0.0)) + float(context.get("V_T.T0", 0.0))
        density = 0.0
        if volume_state > 0.0:
            density = total_tcells / max(volume_state * 1e6, 1e-12)
        context["tcell_density_per_ul"] = density

    return apply


_CANONICAL_BLOCKS: Dict[str, ModuleBlockSpec] = {
    "pd1_bridge_block": ModuleBlockSpec(
        factory=pd1_bridge_block,
        disable_targets=("aPD1",),
        phase="pre",
    ),
    "tumour_geometry_block": ModuleBlockSpec(
        factory=tumour_geometry_block,
        phase="post",
    ),
    "tumour_geometry_dynamic_block": ModuleBlockSpec(
        factory=tumour_geometry_dynamic_block,
        phase="post",
    ),
    "pd1_occupancy_filter_block": ModuleBlockSpec(
        factory=pd1_occupancy_filter_block,
        disable_targets=("H_PD1_C1", "pd1_occupancy"),
        phase="pre",
    ),
    "alignment_driver_block": ModuleBlockSpec(
        factory=alignment_driver_block,
        disable_targets=("H_PD1_C1", "pd1_occupancy"),
        phase="post",
    ),
}

# Provide short aliases so callers can request "pd1_bridge" instead of
# the explicit *_block suffix used internally.
_BLOCK_ALIASES: Dict[str, str] = {
    "pd1_bridge": "pd1_bridge_block",
    "tumour_geometry": "tumour_geometry_block",
    "tumor_geometry": "tumour_geometry_block",
    "tumour_geometry_dynamic": "tumour_geometry_dynamic_block",
    "tumor_geometry_dynamic": "tumour_geometry_dynamic_block",
    "pd1_occupancy_filter": "pd1_occupancy_filter_block",
    "alignment_driver": "alignment_driver_block",
}

MODULE_BLOCK_REGISTRY: Dict[str, ModuleBlockSpec] = dict(_CANONICAL_BLOCKS)


def resolve_module_blocks(
    model: FrozenModel, block_names: Sequence[str] | None
) -> Tuple[List[ModuleBlock], List[ModuleBlock], List[str], List[str]]:
    """Resolve runtime module blocks and derived disable targets."""

    if not block_names:
        return [], [], [], []

    pre_blocks: List[ModuleBlock] = []
    post_blocks: List[ModuleBlock] = []
    disable: List[str] = []
    resolved_names: List[str] = []
    seen: set[str] = set()
    for name in block_names:
        canonical = _BLOCK_ALIASES.get(name, name)
        if canonical in seen:
            continue
        spec = MODULE_BLOCK_REGISTRY.get(canonical)
        if spec is None:
            logger.debug("module block '%s' not registered; skipping runtime hook", name)
            continue
        block = spec.factory(model)
        if spec.phase == "pre":
            pre_blocks.append(block)
        else:
            post_blocks.append(block)
        disable.extend(spec.disable_targets)
        resolved_names.append(canonical)
        seen.add(canonical)
    return pre_blocks, post_blocks, disable, resolved_names


def bind_module_blocks(
    model: FrozenModel,
    pre_blocks: Sequence[ModuleBlock],
    post_blocks: Sequence[ModuleBlock],
) -> None:
    """Attach runtime module blocks to the model's repeated assignments."""

    if not pre_blocks and not post_blocks:
        return

    original = model.evaluate_repeated_assignments

    if getattr(model, "_module_blocks_bound", False):
        existing_pre = getattr(model, "_pre_module_blocks", ())
        existing_post = getattr(model, "_post_module_blocks", ())
        model._pre_module_blocks = tuple([*existing_pre, *pre_blocks])  # type: ignore[attr-defined]
        model._post_module_blocks = tuple([*existing_post, *post_blocks])  # type: ignore[attr-defined]
        return

    def evaluate_with_blocks(self: FrozenModel, context: MutableMapping[str, float]) -> None:
        for block in getattr(self, "_pre_module_blocks", ()):  # type: ignore[attr-defined]
            block(context)
        original(context)  # type: ignore[misc]
        for block in getattr(self, "_post_module_blocks", ()):  # type: ignore[attr-defined]
            block(context)

    model._pre_module_blocks = tuple(pre_blocks)  # type: ignore[attr-defined]
    model._post_module_blocks = tuple(post_blocks)  # type: ignore[attr-defined]
    model._module_blocks_bound = True  # type: ignore[attr-defined]
    model.evaluate_repeated_assignments = MethodType(evaluate_with_blocks, model)


__all__ = [
    "ModuleBlock",
    "ModuleBlockSpec",
    "MODULE_BLOCK_REGISTRY",
    "apply_parameter_overrides",
    "bind_module_blocks",
    "disable_repeated_assignments",
    "alignment_driver_block",
    "pd1_bridge_block",
    "resolve_module_blocks",
    "tumour_geometry_block",
]
