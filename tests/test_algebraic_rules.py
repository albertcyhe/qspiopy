import math
from pathlib import Path

import numpy as np
import sympy as sp

from src.offline.frozen_model import CompiledExpression, FrozenModel, RuleEntry


def _compiled(expression: sp.Expr, symbols: tuple[sp.Symbol, ...]) -> CompiledExpression:
    func = sp.lambdify(symbols, expression, modules=["math"])
    tokens = tuple(str(symbol) for symbol in symbols)
    return CompiledExpression(tokens=tokens, func=func, sympy_expr=expression)


def test_algebraic_rule_solves_non_linear_equation():
    # Build an algebraic rule A^3 - B = 0 with a symbolic solver A = B**(1/3).
    a_sym, b_sym = sp.symbols("A B")
    residual_expr = a_sym**3 - b_sym
    solver_expr = b_sym ** sp.Rational(1, 3)

    algebraic_rule = RuleEntry(
        index=1,
        name="alg_rule",
        rule_type="algebraic",
        target="A",
        expression="A^3 - B",
        compiled=_compiled(residual_expr, (a_sym, b_sym)),
        dependencies=("B",),
        safe_target="SYM_A",
        algebraic_solver=_compiled(solver_expr, (b_sym,)),
    )

    model = FrozenModel(
        name="test",
        source_dir=Path("."),
        species=[],
        species_lookup={},
        species_name_lookup={},
        parameters={},
        compartments={},
        config={},
        rules=[],
        rate_rules=[],
        algebraic_rules=[algebraic_rule],
        reactions=[],
        odes=[],
        events=[],
        doses=[],
        variants=[],
        stoichiometry=[],
        constants={},
        dynamic_indices={"A": 0},
        repeated_assignment_order=[],
        initial_assignment_rules=[],
        time_unit="day",
        solver_type="BDF",
        provenance={},
    )

    context = {"B": 8.0}
    state = np.array([0.0], dtype=float)
    model.apply_algebraic_rules(context, state, mutate=True)

    assert math.isclose(context["A"], 2.0, rel_tol=1e-12)
    assert math.isclose(state[0], 2.0, rel_tol=1e-12)
