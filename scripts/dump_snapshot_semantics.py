#!/usr/bin/env python3
"""Dump snapshot reactions/rules/events/parameters for white-box diff."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List, Sequence


def _read_csv(path: Path) -> List[dict]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [row for row in reader]


def _print_section(title: str, rows: Iterable[str]) -> None:
    print(f"=== {title} ===")
    any_row = False
    for row in rows:
        any_row = True
        print(row)
    if not any_row:
        print("<none>")
    print()


def _matches(text: str, keywords: Sequence[str]) -> bool:
    if not keywords:
        return True
    lower = text.lower()
    return any(keyword.lower() in lower for keyword in keywords)


def dump_snapshot(snapshot_path: Path, keywords: Sequence[str]) -> None:
    reactions = _read_csv(snapshot_path / "reactions.csv")
    rules = _read_csv(snapshot_path / "rules.csv")
    events = _read_csv(snapshot_path / "events.csv")
    parameters = _read_csv(snapshot_path / "parameters.csv")

    reaction_rows = []
    for row in reactions:
        name = row.get("name") or row.get("Reaction") or row.get("identifier") or ""
        reactants = row.get("reactants") or row.get("Reactants") or ""
        products = row.get("products") or row.get("Products") or ""
        rate = row.get("reaction_rate") or row.get("KineticLaw") or row.get("kinetic_expression") or ""
        line = f"[{name}] {reactants} -> {products} ; rate={rate}"
        if _matches(line, keywords):
            reaction_rows.append(line)
    _print_section("Reactions", reaction_rows)

    rule_rows = []
    for row in rules:
        rule_type = row.get("type") or row.get("RuleType") or ""
        expression = row.get("expression") or row.get("Rule") or ""
        target = row.get("target") or ""
        line = f"({rule_type}) {target} = {expression}" if target else f"({rule_type}) {expression}"
        if _matches(line, keywords):
            rule_rows.append(line)
    _print_section("Rules", rule_rows)

    event_rows = []
    for row in events:
        trigger = row.get("trigger") or row.get("Trigger") or ""
        delay = row.get("delay") or row.get("Delay") or "0"
        assignments = row.get("assignments") or row.get("EventFcns") or ""
        line = f"Trigger: {trigger} | delay={delay} | assignments={assignments}"
        if _matches(line, keywords):
            event_rows.append(line)
    _print_section("Events", event_rows)

    param_rows = []
    for row in parameters:
        name = row.get("name") or row.get("Name") or ""
        value = row.get("value") or row.get("Value") or ""
        units = row.get("units") or row.get("Units") or ""
        line = f"{name} = {value} {units}"
        if _matches(line, keywords):
            param_rows.append(line)
    _print_section("Parameters", param_rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Dump snapshot semantics for comparison")
    parser.add_argument("snapshot", help="Snapshot directory (e.g. artifacts/matlab_frozen_model/example1)")
    parser.add_argument(
        "--keyword",
        action="append",
        default=[],
        help="Substring filter (case-insensitive). Repeat to add multiple keywords",
    )
    args = parser.parse_args()
    dump_snapshot(Path(args.snapshot), args.keyword)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
