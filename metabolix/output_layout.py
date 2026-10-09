"""BIDS-entity filenames and tidy per-acquisition output paths."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from metabolix.discovery import MRSInput


_ENTITY_ORDER = (
    "sub", "ses", "task", "acq", "ce", "rec", "dir", "run", "mod", "echo", "part",
    "recording", "proc", "space", "split", "voi", "nuc", "mt", "inv", "flip", "tr", "te",
)


def entity_stem(item: MRSInput) -> str:
    """Return the entity prefix in a stable BIDS-oriented order."""
    parts = []
    for entity in _ENTITY_ORDER:
        value = item.entities.get(entity)
        if value:
            safe = re.sub(r"[^A-Za-z0-9]+", "", str(value))
            if safe:
                parts.append(f"{entity}-{safe}")
    return "_".join(parts)


def bids_filename(item: MRSInput, description: str, suffix: str, extension: str) -> str:
    description = re.sub(r"[^A-Za-z0-9-]+", "", description)
    if not description:
        raise ValueError("BIDS derivative description must contain at least one alphanumeric character.")
    extension = extension if extension.startswith(".") else f".{extension}"
    return f"{entity_stem(item)}_desc-{description}_{suffix}{extension}"


@dataclass(frozen=True)
class OutputLayout:
    """Locations for one acquisition bundle in a BIDS derivatives dataset."""

    output_root: Path
    item: MRSInput

    @property
    def stem(self) -> str:
        return entity_stem(self.item)

    @property
    def subject_dir(self) -> Path:
        return self.output_root / f"sub-{self.item.entities['sub']}"

    @property
    def session_dir(self) -> Path:
        session = self.item.entities.get("ses")
        return self.subject_dir / f"ses-{session}" if session else self.subject_dir

    @property
    def mrs_dir(self) -> Path:
        return self.session_dir / "mrs"

    @property
    def bundle_dir(self) -> Path:
        return self.mrs_dir / f"{self.stem}_desc-metabolix"

    @property
    def processing_dir(self) -> Path:
        return self.bundle_dir / "processing" / "fsl-mrs"

    @property
    def work_dir(self) -> Path:
        return self.bundle_dir / "work"

    @property
    def qc_dir(self) -> Path:
        return self.bundle_dir / "qc"

    @property
    def map_dir(self) -> Path:
        return self.bundle_dir / "maps"

    @property
    def figure_dir(self) -> Path:
        return self.bundle_dir / "figures"

    @property
    def roi_dir(self) -> Path:
        return self.bundle_dir / "roi"

    @property
    def log_dir(self) -> Path:
        return self.bundle_dir / "logs"

    @property
    def provenance_dir(self) -> Path:
        return self.bundle_dir / "provenance"

    def branch_dir(self, branch: str) -> Path:
        name = "original" if branch == "original" else "waterremoved"
        return self.processing_dir / f"fit-{name}"

    def branch_qc_table(self, branch: str) -> Path:
        name = "original" if branch == "original" else "waterremoved"
        return self.qc_dir / bids_filename(self.item, f"metabolix-{name}", "voxelqc", ".tsv")

    def stage_dir(self, stage: str) -> Path:
        mapping = {
            "coilcombine": self.processing_dir / "coilcombine",
            "waterremove": self.processing_dir / "water-removal",
            "fit_original": self.branch_dir("original"),
            "fit_waterremoved": self.branch_dir("waterremoved"),
        }
        return mapping[stage]

    def report_path(self, extension: str = ".html") -> Path:
        return self.bundle_dir / bids_filename(self.item, "metabolix", "report", extension)

    def manifest_path(self) -> Path:
        return self.provenance_dir / bids_filename(self.item, "metabolix", "provenance", ".json")

    def config_path(self) -> Path:
        return self.provenance_dir / bids_filename(self.item, "metabolix", "config", ".yaml")

    def log_path(self) -> Path:
        return self.log_dir / bids_filename(self.item, "metabolix", "log", ".txt")