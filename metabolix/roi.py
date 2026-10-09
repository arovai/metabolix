"""Nominal MRS voxel/ROI overlap using image-world affine sampling."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from metabolix.discovery import MRSInput
from metabolix.output_layout import bids_filename
from metabolix.validation import ValidatedInput


def _find_segmentation(config: dict[str, Any], item: MRSInput) -> tuple[Path, dict[str, int]]:
    roi = config["roi"]
    if roi.get("map"):
        path = Path(roi["map"])
        labels = roi.get("labels") or {"roi": 1}
        return path, {str(name): int(value) for name, value in labels.items()}
    root = roi.get("freesurfer_dir")
    if not root:
        raise ValueError("ROI calculation requires roi.map or roi.freesurfer_dir.")
    root_path = Path(root)
    subject = f"sub-{item.entities['sub']}"
    candidates = [root_path / subject / "mri" / "aseg.mgz", root_path / "mri" / "aseg.mgz"]
    segmentation = next((candidate for candidate in candidates if candidate.is_file()), None)
    if segmentation is None:
        raise ValueError(f"Could not find FreeSurfer aseg.mgz for {subject} under {root_path}.")
    labels = roi.get("labels") or {"left_thalamus": 10, "right_thalamus": 49}
    return segmentation, {str(name): int(value) for name, value in labels.items()}


def _sample_fractions(
    seg: np.ndarray,
    mrsi_affine: np.ndarray,
    seg_affine: np.ndarray,
    shape: tuple[int, int, int],
    mrsi_spacing: np.ndarray,
    labels: dict[str, int],
    step_mm: float,
) -> tuple[dict[str, np.ndarray], np.ndarray, list[int]]:
    counts = np.maximum(1, np.ceil(mrsi_spacing / step_mm).astype(int))
    axes = [(np.arange(count) + 0.5) / count - 0.5 for count in counts]
    offsets = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    transform = np.linalg.inv(seg_affine) @ mrsi_affine
    result = {name: np.zeros(shape, dtype=np.float32) for name in labels}
    fov = np.zeros(shape, dtype=np.float32)
    for voxel in np.ndindex(*shape):
        seg_vox = nib.affines.apply_affine(transform, offsets + np.asarray(voxel))
        indices = np.floor(seg_vox + 0.5).astype(int)
        inside = np.all((indices >= 0) & (indices < np.asarray(seg.shape)), axis=1)
        fov[voxel] = inside.mean()
        values = np.zeros(len(indices), dtype=seg.dtype)
        values[inside] = seg[tuple(indices[inside].T)]
        for name, label in labels.items():
            result[name][voxel] = np.count_nonzero(values == label) / len(values)
    return result, fov, counts.tolist()


def calculate_roi_overlap(
    item: MRSInput,
    validated: ValidatedInput,
    output_dir: Path,
    config: dict[str, Any],
    logger: logging.Logger,
) -> dict[str, Any]:
    seg_path, labels = _find_segmentation(config, item)
    seg_img = nib.load(str(seg_path))
    if len(seg_img.shape) != 3:
        raise ValueError(f"ROI image must be 3D: {seg_path}")
    seg = np.asanyarray(seg_img.dataobj)
    if not np.isfinite(seg).all():
        raise ValueError(f"ROI image contains nonfinite values: {seg_path}")
    if not np.allclose(seg, np.rint(seg)):
        raise ValueError("ROI label map must contain integer-valued labels.")
    seg = np.rint(seg).astype(np.int32)
    if config["roi"].get("map") and not config["roi"].get("labels"):
        seg = (seg != 0).astype(np.int32)
        labels = {"roi": 1}
    missing = [name for name, value in labels.items() if not np.any(seg == value)]
    if missing:
        raise ValueError(f"ROI labels absent from segmentation {seg_path}: {', '.join(missing)}")
    if len(validated.image.shape) > 3:
        spatial_shape = tuple(validated.image.shape[:3])
    else:
        spatial_shape = tuple(validated.image.shape[:3])
    affine = np.asarray(validated.image.affine)
    spacing = nib.affines.voxel_sizes(affine)[:3]
    directions = affine[:3, :3] / spacing
    if not np.allclose(directions.T @ directions, np.eye(3), atol=1e-4):
        raise ValueError("Sheared MRS affine; nominal voxel sampling assumptions need review.")
    step = float(config["roi"]["sampling_step_mm"])
    fine, fov, sample_counts = _sample_fractions(seg, affine, seg_img.affine, spatial_shape, spacing, labels, step)
    coarse, _, _ = _sample_fractions(seg, affine, seg_img.affine, spatial_shape, spacing, labels, 1.0)
    roi_dir = output_dir / "roi"
    roi_dir.mkdir(parents=True, exist_ok=True)
    qform, qcode = validated.image.get_qform(coded=True)
    sform, scode = validated.image.get_sform(coded=True)
    seg_qform, seg_qcode = seg_img.get_qform(coded=True)
    seg_sform, seg_scode = seg_img.get_sform(coded=True)
    geometry_warning = bool(
        (qcode and scode and not np.allclose(qform, sform, atol=0.01))
        or (seg_qcode and seg_scode and not np.allclose(seg_qform, seg_sform, atol=0.01))
    )
    field_coverage = float(fov.mean())
    rows = []
    mrsi_voxel_volume = abs(float(np.linalg.det(affine[:3, :3])))
    seg_voxel_volume = abs(float(np.linalg.det(seg_img.affine[:3, :3])))
    for voxel in np.ndindex(*spatial_shape):
        row: dict[str, Any] = {"i": voxel[0], "j": voxel[1], "k": voxel[2], "segmentation_fov_pct": float(fov[voxel] * 100)}
        for name in labels:
            fraction = float(fine[name][voxel])
            row[f"{name}_overlap_pct"] = fraction * 100
            row[f"{name}_overlap_mm3"] = fraction * mrsi_voxel_volume
            row[f"{name}_sampling_difference_pp"] = abs(float(fine[name][voxel] - coarse[name][voxel])) * 100
        rows.append(row)
    columns = list(rows[0])
    overlap_table = roi_dir / bids_filename(item, "metabolix-roi-overlap", "voxelqc", ".tsv")
    with overlap_table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    overlap_table.with_suffix(".json").write_text(
        json.dumps(
            {
                "Description": "Nominal per-MRS-voxel ROI overlap fractions and segmentation field-of-view coverage.",
                "Denominator": "Entire nominal MRS voxel volume.",
                "SamplingStepMM": step,
                "SensitivityComparisonStepMM": 1.0,
                "Coordinates": "Zero-based MRS spatial voxel indices.",
                "Registration": "No registration performed; existing image-world affines are assumed aligned.",
                "Labels": labels,
                "MRSIResponse": "Finite spatial response may extend beyond nominal voxel boundaries.",
            },
            indent=2,
        )
        + "\n"
    )
    summaries = []
    for name, label in labels.items():
        volume = int(np.count_nonzero(seg == label)) * seg_voxel_volume
        covered_volume = float(fine[name].sum()) * mrsi_voxel_volume
        summaries.append(
            {
                "roi": name,
                "label": label,
                "segmentation_volume_mm3": volume,
                "estimated_volume_within_mrsi_grid_mm3": covered_volume,
                "roi_fraction_covered_by_mrsi_grid_pct": 100 * covered_volume / volume if volume else None,
                "maximum_fine_vs_1mm_difference_pp": float(np.max(np.abs(fine[name] - coarse[name])) * 100),
                "voxels_at_or_above_reporting_threshold": (
                    int(np.count_nonzero(fine[name] * 100 >= float(threshold)))
                    if (threshold := config["roi"].get("minimum_overlap_percent")) is not None
                    else None
                ),
            }
        )
        map_image = nib.Nifti1Image((fine[name] * 100).astype(np.float32), affine)
        map_image.header.set_xyzt_units("mm")
        map_path = roi_dir / bids_filename(item, f"metabolix-roi-{name}-overlapPercent", "statmap", ".nii.gz")
        nib.save(map_image, str(map_path))
        map_path.with_name(map_path.name[:-7] + ".json").write_text(
            json.dumps({"Description": f"Nominal overlap percentage for ROI {name}.", "Units": "percent", "ROIlabel": label, "Denominator": "Entire nominal MRS voxel volume."}, indent=2) + "\n"
        )
    fov_image = nib.Nifti1Image((fov * 100).astype(np.float32), affine)
    fov_path = roi_dir / bids_filename(item, "metabolix-roi-segmentationFOVCoverage", "statmap", ".nii.gz")
    nib.save(fov_image, str(fov_path))
    fov_path.with_name(fov_path.name[:-7] + ".json").write_text(
        json.dumps({"Description": "Fraction of each nominal MRS voxel covered by the segmentation image field of view.", "Units": "percent"}, indent=2) + "\n"
    )
    report = {
        "segmentation": str(seg_path),
        "labels": labels,
        "sampling_step_mm": step,
        "samples_per_voxel_axis": sample_counts,
        "coarse_comparison_step_mm": 1.0,
        "mean_segmentation_fov_coverage_pct": field_coverage * 100,
        "mri_grid_affine": affine.tolist(),
        "segmentation_affine": np.asarray(seg_img.affine).tolist(),
        "registration": "No registration performed; existing image-world affines are assumed to be aligned.",
        "qform_sform_disagreement": geometry_warning,
        "interpretation": "Nominal geometric overlap only; MRSI spatial response may extend beyond nominal voxel boundaries.",
        "roi_summaries": summaries,
    }
    report["overlap_table"] = str(overlap_table)
    (roi_dir / bids_filename(item, "metabolix-roi", "summary", ".json")).write_text(json.dumps(report, indent=2) + "\n")
    logger.info("Wrote nominal ROI overlap for %s to %s", item.source_name, roi_dir)
    return report