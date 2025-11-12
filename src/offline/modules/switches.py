"""Lightweight module switches applied on top of frozen snapshots."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Literal, MutableMapping, Sequence, Tuple

from types import MethodType

from ..snapshot import FrozenModel

logger = logging.getLogger(__name__)

ModuleBlock = Callable[[MutableMapping[str, float]], None]

AVOGADRO = 6.02214076e23
LITERS_PER_CUBIC_MICROMETER = 1e-15
DEFAULT_SYN_DEPTH_UM = 1.15e-5


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
}

# Provide short aliases so callers can request "pd1_bridge" instead of
# the explicit *_block suffix used internally.
_BLOCK_ALIASES: Dict[str, str] = {
    "pd1_bridge": "pd1_bridge_block",
    "tumour_geometry": "tumour_geometry_block",
    "tumor_geometry": "tumour_geometry_block",
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
    "pd1_bridge_block",
    "resolve_module_blocks",
    "tumour_geometry_block",
]
