from __future__ import annotations

import pytest

from src.offline.param_graph import ParameterGraph
from src.offline.errors import SnapshotError


def test_parameter_graph_evaluates_dependencies() -> None:
    graph = ParameterGraph(lambda value, unit: value)
    graph.add_base("CL", 2.0, "liter/day")
    graph.add_base("V_C", 1.0, "liter")
    graph.add_spec("k_el", unit="1/day", expression="p(1)/p(2)", dependencies=["CL", "V_C"])
    values = graph.evaluate()
    assert values["k_el"] == pytest.approx(2.0)
    meta = graph.metadata["k_el"]
    assert meta["canonical_value"] == pytest.approx(2.0)
    assert meta["dependencies"] == ["CL", "V_C"]


def test_parameter_graph_detects_unresolved_cycle() -> None:
    graph = ParameterGraph(lambda value, unit: value)
    graph.add_spec("A", unit="", expression="p(1)", dependencies=["B"])
    graph.add_spec("B", unit="", expression="p(1)", dependencies=["A"])
    with pytest.raises(SnapshotError):
        graph.evaluate()
