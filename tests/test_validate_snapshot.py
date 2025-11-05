import shutil
from pathlib import Path

import pandas as pd
import pytest

from scripts.validate_snapshot import validate_snapshot
from src.offline.snapshot import load_frozen_model


def test_validate_snapshot_example_pass():
    validate_snapshot(Path("artifacts/matlab_frozen_model/example1"))


def test_validate_snapshot_detects_invalid_species(tmp_path):
    src = Path("artifacts/matlab_frozen_model/example1")
    dst = tmp_path / "snapshot"
    shutil.copytree(src, dst)
    species_path = dst / "species.csv"
    frame = pd.read_csv(species_path)
    frame.loc[0, "interpreted_dimension"] = "invalid_dimension"
    frame.to_csv(species_path, index=False)

    with pytest.raises(SystemExit) as excinfo:
        validate_snapshot(dst)
    assert "interpreted_dimension" in str(excinfo.value)


def test_load_frozen_model_accepts_directory_path():
    path = Path("artifacts/matlab_frozen_model/example1")
    model = load_frozen_model(path)
    assert model.name == path.name
