"""BIDS-style MRS input discovery and entity handling."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


_ENTITY = re.compile(r"(?:^|_)(sub|ses|task|acq|run|rec|echo|part)-([^_]+)")


@dataclass(frozen=True)
class MRSInput:
    path: Path
    sidecar: Path | None
    entities: dict[str, str]
    standard_bids: bool

    @property
    def source_name(self) -> str:
        name = self.path.name
        for suffix in (".nii.gz", ".nii"):
            if name.endswith(suffix):
                return name[: -len(suffix)]
        return self.path.stem


def parse_entities(path: Path, bids_root: Path) -> dict[str, str]:
    entities = dict(_ENTITY.findall(path.name))
    if "sub" not in entities:
        for parent in path.parents:
            match = re.fullmatch(r"sub-(.+)", parent.name)
            if match:
                entities["sub"] = match.group(1)
                break
    if "ses" not in entities:
        for parent in path.parents:
            match = re.fullmatch(r"ses-(.+)", parent.name)
            if match:
                entities["ses"] = match.group(1)
                break
    if not entities.get("sub"):
        try:
            relative = path.relative_to(bids_root)
        except ValueError:
            relative = path
        raise ValueError(f"Could not determine a BIDS subject entity for {relative}.")
    return entities


def discover_mrs(bids_root: Path, filters: dict[str, list[str]]) -> list[MRSInput]:
    if not bids_root.is_dir():
        raise ValueError(f"BIDS directory does not exist: {bids_root}")
    paths: set[Path] = set()
    for pattern in ("*_mrs.nii", "*_mrs.nii.gz", "*_mrsi.nii", "*_mrsi.nii.gz"):
        paths.update(bids_root.rglob(pattern))
    candidates: list[MRSInput] = []
    for path in sorted(paths):
        if "derivatives" in path.relative_to(bids_root).parts:
            continue
        entities = parse_entities(path, bids_root)
        accepted = True
        for entity, requested in filters.items():
            if requested and entities.get(entity) not in requested:
                accepted = False
                break
        if not accepted:
            continue
        sidecar = path.with_name(path.name[:-7] + ".json") if path.name.endswith(".nii.gz") else path.with_suffix(".json")
        standard = path.name.endswith(("_mrs.nii", "_mrs.nii.gz")) and "mrs" in path.parent.parts
        candidates.append(MRSInput(path, sidecar if sidecar.is_file() else None, entities, standard))
    if not candidates:
        raise ValueError(
            f"No NIfTI-MRS inputs found under {bids_root}. Expected BIDS mrs/*_mrs.nii[.gz] "
            "or the documented compatibility pattern *_mrsi.nii[.gz]."
        )
    return candidates


def output_identity(item: MRSInput) -> Path:
    subject = f"sub-{item.entities['sub']}"
    parts = [subject]
    if item.entities.get("ses"):
        parts.append(f"ses-{item.entities['ses']}")
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "-", item.source_name)
    return Path(*parts, safe_name)