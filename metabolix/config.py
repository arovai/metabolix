"""Configuration defaults, YAML loading, validation, and serialization."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml


DEFAULT_CONFIG: dict[str, Any] = {
    "processing": {
        "coil_combine": "auto",
        "coil_reference": None,
        "water_removal": "off",
        "water_ppm": [4.5, 4.9],
    },
    "basis": {"path": None, "mapping": [], "binary_dir": None},
    "fit": {
        "ppmlim": None,
        "ignore": [],
        "combine": [],
        "internal_reference": [],
        "baseline": None,
        "algorithm": None,
        "mask": None,
    },
    "roi": {
        "freesurfer_dir": None,
        "map": None,
        "labels": {},
        "sampling_step_mm": 0.5,
        "minimum_overlap_percent": None,
    },
    "execution": {"n_jobs": 1, "resume": False, "overwrite": False},
    "qc": {"max_crlb_percent": None, "min_snr": None, "max_linewidth_hz": None},
    "report": {"enabled": True, "metabolites": []},
    "selection": {
        "participant_label": [],
        "session": [],
        "acquisition": [],
        "run": [],
    },
}

_ALLOWED_KEYS = {key: set(value) for key, value in DEFAULT_CONFIG.items()}


class ConfigError(ValueError):
    """Raised when configuration is malformed or contains invalid values."""


def _merge(target: dict[str, Any], source: dict[str, Any]) -> dict[str, Any]:
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _merge(target[key], value)
        else:
            target[key] = copy.deepcopy(value)
    return target


def _resolve_paths(config: dict[str, Any], base: Path) -> None:
    for section, key in (
        ("processing", "coil_reference"),
        ("basis", "path"),
        ("basis", "binary_dir"),
        ("fit", "mask"),
        ("roi", "freesurfer_dir"),
        ("roi", "map"),
    ):
        value = config[section].get(key)
        if value and not Path(value).is_absolute():
            config[section][key] = str((base / value).resolve())
    for mapping in config["basis"].get("mapping", []):
        if "path" in mapping and mapping["path"] and not Path(mapping["path"]).is_absolute():
            mapping["path"] = str((base / mapping["path"]).resolve())


def validate_config(config: dict[str, Any]) -> None:
    if set(config) != set(DEFAULT_CONFIG):
        raise ConfigError("Configuration sections are incomplete or invalid.")
    if config["processing"]["water_removal"] not in {"off", "on", "compare"}:
        raise ConfigError("processing.water_removal must be off, on, or compare.")
    ppm = config["processing"]["water_ppm"]
    if not isinstance(ppm, list) or len(ppm) != 2 or not 0 < float(ppm[0]) < float(ppm[1]):
        raise ConfigError("processing.water_ppm must be two increasing positive ppm values.")
    jobs = config["execution"]["n_jobs"]
    if not isinstance(jobs, int) or jobs < 1 or jobs > 256:
        raise ConfigError("execution.n_jobs must be an integer from 1 to 256.")
    step = float(config["roi"]["sampling_step_mm"])
    if not 0 < step <= 0.5:
        raise ConfigError("roi.sampling_step_mm must be > 0 and <= 0.5 mm.")
    threshold = config["roi"]["minimum_overlap_percent"]
    if threshold is not None and not 0 <= float(threshold) <= 100:
        raise ConfigError("roi.minimum_overlap_percent must be between 0 and 100.")
    if config["processing"]["coil_combine"] not in {"auto", "on", "off"}:
        raise ConfigError("processing.coil_combine must be auto, on, or off.")
    if config["fit"]["algorithm"] not in {None, "Newton", "MH"}:
        raise ConfigError("fit.algorithm must be Newton or MH.")
    if config["fit"]["baseline"] is not None and not isinstance(config["fit"]["baseline"], str):
        raise ConfigError("fit.baseline must be a string understood by FSL-MRS.")
    for key in ("ignore", "combine", "internal_reference"):
        if not isinstance(config["fit"][key], list):
            raise ConfigError(f"fit.{key} must be a list.")
    if any(not isinstance(pair, list) or len(pair) < 2 for pair in config["fit"]["combine"]):
        raise ConfigError("Each fit.combine entry must contain at least two component names.")
    if not isinstance(config["basis"]["mapping"], list):
        raise ConfigError("basis.mapping must be a list.")
    for key, value in config["qc"].items():
        if value is not None and float(value) <= 0:
            raise ConfigError(f"qc.{key} must be positive when set.")
    if config["fit"]["ppmlim"] is not None:
        limits = config["fit"]["ppmlim"]
        if not isinstance(limits, list) or len(limits) != 2 or float(limits[0]) >= float(limits[1]):
            raise ConfigError("fit.ppmlim must contain two increasing ppm values.")
    for mapping in config["basis"]["mapping"]:
        if not isinstance(mapping, dict) or "path" not in mapping:
            raise ConfigError("Each basis.mapping entry must be a mapping with a path.")


def load_config(path: Path | None, overrides: dict[str, Any]) -> dict[str, Any]:
    config = copy.deepcopy(DEFAULT_CONFIG)
    if path:
        try:
            loaded = yaml.safe_load(path.read_text()) or {}
        except (OSError, yaml.YAMLError) as exc:
            raise ConfigError(f"Could not read YAML config {path}: {exc}") from exc
        if not isinstance(loaded, dict):
            raise ConfigError("The YAML root must be a mapping.")
        unknown_sections = set(loaded) - set(DEFAULT_CONFIG)
        if unknown_sections:
            raise ConfigError(f"Unknown config section(s): {', '.join(sorted(unknown_sections))}")
        for section, values in loaded.items():
            if not isinstance(values, dict):
                raise ConfigError(f"Config section {section!r} must be a mapping.")
            unknown = set(values) - _ALLOWED_KEYS[section]
            if unknown:
                raise ConfigError(f"Unknown key(s) in {section}: {', '.join(sorted(unknown))}")
        _merge(config, loaded)
        _resolve_paths(config, path.resolve().parent)
    _merge(config, overrides)
    validate_config(config)
    return config


def dump_yaml(config: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=False))