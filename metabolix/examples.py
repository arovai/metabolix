"""Packaged configuration examples."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from metabolix.config import DEFAULT_CONFIG


def example_config() -> dict[str, Any]:
    config = deepcopy(DEFAULT_CONFIG)
    config["processing"].update({"coil_combine": "auto", "water_removal": "off", "water_ppm": [4.5, 4.9]})
    config["basis"]["path"] = "/data/basis/LCModel_GE_UnEdited_PRESS_144_.BASIS"
    config["fit"].update(
        {
            "ppmlim": [1.8, 4.0],
            "ignore": ["H2O", "CrCH2", "Ala", "Lac"],
            "combine": [["NAA", "NAAG"], ["Cr", "PCr"], ["GPC", "PCh"], ["Glu", "Gln"]],
            "internal_reference": ["Cr", "PCr"],
            "baseline": "polynomial, 2",
            "algorithm": "Newton",
        }
    )
    config["roi"].update(
        {
            "freesurfer_dir": "/data/derivatives/freesurfer_7.3.2",
            "labels": {"left_thalamus": 10, "right_thalamus": 49},
            "sampling_step_mm": 0.5,
            "minimum_overlap_percent": None,
        }
    )
    config["execution"]["n_jobs"] = 1
    config["report"]["metabolites"] = ["NAA+NAAG", "Cr+PCr", "GPC+PCh", "Glu+Gln"]
    return config