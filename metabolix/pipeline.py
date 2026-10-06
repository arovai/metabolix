"""Participant and group workflows."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import logging
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from metabolix import __version__
from metabolix.discovery import MRSInput, discover_mrs, output_identity
from metabolix.execution import dependency_report, executable_version, find_executable, run_command, save_manifest
from metabolix.roi import calculate_roi_overlap
from metabolix.validation import (
    checksum,
    inspect_basis,
    make_voxel_diagnostic_pdf,
    prepare_finite_copy,
    validate_mrs,
    write_validity_mask,
)
from metabolix.config import dump_yaml


def _entity_filters(config: dict[str, Any]) -> dict[str, list[str]]:
    result = {}
    for key in ("participant_label", "session", "acquisition", "run"):
        values = config["selection"].get(key, [])
        if key == "participant_label":
            values = [value.removeprefix("sub-") for value in values]
            key = "sub"
        elif key == "session":
            values = [value.removeprefix("ses-") for value in values]
            key = "ses"
        elif key == "acquisition":
            key = "acq"
        result[key] = [str(value) for value in values]
    return result


def _select_basis(item: MRSInput, metadata: dict[str, Any], config: dict[str, Any]) -> Path | None:
    basis_cfg = config["basis"]
    if basis_cfg.get("path"):
        path = Path(basis_cfg["path"])
        if not path.exists():
            raise ValueError(f"Configured basis does not exist: {path}")
        return path.resolve()
    matches = []
    for entry in basis_cfg.get("mapping", []):
        criteria = entry.get("match", {})
        candidate_values = {**item.entities, **metadata}
        if all(candidate_values.get(key) == value for key, value in criteria.items()):
            path = Path(entry["path"])
            if path.exists():
                matches.append(path.resolve())
    unique = list(dict.fromkeys(matches))
    if len(unique) > 1:
        choices = "\n".join(f"  - {path}" for path in unique)
        raise ValueError(
            f"Basis mapping is ambiguous for {item.path}; multiple entries matched:\n{choices}\n"
            "Make the `match` fields more specific so exactly one mapping entry applies, or set `basis.path` / pass `--basis PATH` for a single-basis run. "
            "See README.md (Configuration) for the basis.mapping format."
        )
    if not unique:
        return None
    return unique[0]


def _run_directory(output_root: Path, item: MRSInput) -> Path:
    return output_root / output_identity(item)


def _write_dataset_description(output_root: Path) -> None:
    path = output_root / "dataset_description.json"
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "Name": "Metabolix MRSI derivatives",
                "BIDSVersion": "1.9.0",
                "DatasetType": "derivative",
                "GeneratedBy": [{"Name": "metabolix", "Version": __version__}],
                "SourceDatasets": [],
            },
            indent=2,
        )
        + "\n"
    )


def _fingerprint(item: MRSInput, basis: Path, config: dict[str, Any], fsl_version: str) -> str:
    stable_config = json.loads(json.dumps(config, default=str))
    stable_config["execution"].pop("resume", None)
    stable_config["execution"].pop("overwrite", None)
    payload = {
        "input_sha256": checksum(item.path),
        "sidecar_sha256": checksum(item.sidecar) if item.sidecar else None,
        "basis_sha256": checksum(basis),
        "coil_reference_sha256": checksum(Path(config["processing"]["coil_reference"])) if config["processing"].get("coil_reference") else None,
        "config": stable_config,
        "metabolix_version": __version__,
        "fsl_mrs_version": fsl_version,
    }
    canonical = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def _resolve_fsl(name: str, config: dict[str, Any]) -> str:
    binary = find_executable(name, config["basis"].get("binary_dir"))
    if not binary:
        raise RuntimeError(
            f"Required executable {name!r} was not found. Install FSL-MRS in its supported FSL/Conda environment, "
            "activate that environment, then run `metabolix --check-deps`. You can also set `basis.binary_dir` or "
            "pass `--binary-dir PATH` if the executables are not on PATH. `fslpy` alone is not sufficient. "
            "See README.md (Installation) and https://pages.fmrib.ox.ac.uk/fsl/fsl_mrs/install.html."
        )
    return binary


def _expected_nifti(stage_dir: Path, stem: str) -> Path:
    expected = [stage_dir / f"{stem}.nii.gz", stage_dir / f"{stem}.nii"]
    for path in expected:
        if path.is_file():
            return path
    found = sorted(stage_dir.glob(f"{stem}*.nii*"))
    if len(found) == 1:
        return found[0]
    raise RuntimeError(f"Expected one NIfTI-MRS output for {stem!r} in {stage_dir}; found {found}.")


def _validate_stage_output(
    output_path: Path,
    reference_path: Path,
    entities: dict[str, str],
    coil_removed: bool,
    logger: logging.Logger,
) -> None:
    reference = nib.load(str(reference_path))
    output_item = MRSInput(output_path, None, entities, False)
    output = validate_mrs(output_item, logger)
    expected_shape = tuple(reference.shape[:4]) if coil_removed else tuple(reference.shape)
    if tuple(output.image.shape) != expected_shape:
        raise ValueError(f"Preprocessing changed unexpected dimensions: {reference.shape} -> {output.image.shape}; expected {expected_shape}.")
    if not np.allclose(output.image.affine, reference.affine, atol=1e-5):
        raise ValueError(f"Preprocessing changed spatial affine: {output_path}")
    if not np.allclose(output.image.header.get_zooms()[:4], reference.header.get_zooms()[:4], rtol=1e-6, atol=1e-8):
        raise ValueError(f"Preprocessing changed spatial or spectral sampling: {output_path}")


def _fit_command(
    executable: str,
    data_path: Path,
    basis: Path,
    mask: Path,
    output: Path,
    config: dict[str, Any],
) -> list[str]:
    fit = config["fit"]
    argv = [executable, "--data", str(data_path), "--basis", str(basis), "--mask", str(mask), "--output", str(output)]
    if fit.get("ppmlim"):
        argv.extend(["--ppmlim", *map(str, fit["ppmlim"])])
    if fit.get("ignore"):
        argv.extend(["--ignore", *map(str, fit["ignore"])])
    for pair in fit.get("combine", []):
        argv.extend(["--combine", *map(str, pair)])
    if fit.get("internal_reference"):
        argv.extend(["--internal_ref", *map(str, fit["internal_reference"])])
    if fit.get("baseline"):
        argv.extend(["--baseline", str(fit["baseline"])])
    if fit.get("algorithm"):
        argv.extend(["--algo", str(fit["algorithm"])])
    workers = int(config["execution"]["n_jobs"])
    argv.extend(["--parallel", "local" if workers > 1 else "off"])
    if workers > 1:
        argv.extend(["--parallel-workers", str(workers)])
    if config["report"]["enabled"]:
        argv.append("--report")
    return argv


def _load_metric_maps(
    fit_dir: Path,
    spatial_shape: tuple[int, ...],
    data_path: Path,
    config: dict[str, Any],
    logger: logging.Logger,
) -> dict[str, np.ndarray]:
    metric_maps: dict[str, np.ndarray] = {}
    map_groups = (
        ("concs/raw", "raw_fit_", None),
        ("concs/internal", "internal_reference_ratio_", None),
        ("uncertainties", "crlb_percent_", "_sd"),
        ("qc", "snr_", "_snr"),
        ("qc", "linewidth_fwhm_hz_", "_fwhm"),
        ("nuisance", "nuisance_", None),
    )
    for folder, prefix, suffix in map_groups:
        directory = fit_dir / folder
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*.nii*")):
            name = path.name.removesuffix(".nii.gz").removesuffix(".nii")
            if suffix and not name.endswith(suffix):
                continue
            if suffix:
                name = name[: -len(suffix)]
            try:
                image = nib.load(str(path))
                values = np.asanyarray(image.dataobj)
            except Exception as exc:
                logger.warning("Could not read FSL-MRS QC map %s: %s", path, exc)
                continue
            if values.shape != spatial_shape:
                logger.warning("Ignoring map with unexpected spatial shape: %s (%s)", path, values.shape)
                continue
            key = f"{prefix}{name}"
            metric_maps[key] = np.asarray(values, dtype=float)

    reference = config["fit"].get("internal_reference") or ["Cr", "PCr"]
    reference_maps = [metric_maps.get(f"raw_fit_{name}") for name in reference]
    if all(values is not None for values in reference_maps):
        reference_sum = np.sum(reference_maps, axis=0)
        zero_reference = ~np.isfinite(reference_sum) | (reference_sum == 0)
        for key, values in list(metric_maps.items()):
            if key.startswith("internal_reference_ratio_"):
                values = values.copy()
                values[zero_reference] = np.nan
                metric_maps[key] = values

    try:
        data = np.asanyarray(nib.load(str(data_path)).dataobj)
        fit = np.asanyarray(nib.load(str(fit_dir / "fit" / "fit.nii.gz")).dataobj)
        residual = np.asanyarray(nib.load(str(fit_dir / "fit" / "residual.nii.gz")).dataobj)
        reduce_axes = tuple(range(3, residual.ndim))
        input_norm = np.linalg.norm(data, axis=reduce_axes)
        residual_norm = np.linalg.norm(residual, axis=reduce_axes)
        with np.errstate(divide="ignore", invalid="ignore"):
            metric_maps["residual_l2_over_input_l2"] = residual_norm / input_norm
    except (OSError, ValueError) as exc:
        logger.warning("Could not derive the documented residual norm summary: %s", exc)
    return metric_maps


def _write_qc_table(
    path: Path,
    item: MRSInput,
    valid: np.ndarray,
    reasons: np.ndarray,
    branch: str,
    fit_status: str,
    metric_maps: dict[str, np.ndarray],
    config: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    entities = {f"entity_{key}": value for key, value in item.entities.items()}
    fields = [*entities, "branch", "i", "j", "k", "status", "reason", *sorted(metric_maps), "qc_flags"]
    flags_by_voxel: dict[tuple[int, ...], list[str]] = {}
    thresholds = config["qc"]
    for index in np.ndindex(*valid.shape):
        flags = []
        if not valid[index]:
            flags_by_voxel[index] = flags
            continue
        for prefix, threshold_key, comparator, flag in (
            ("crlb_percent_", "max_crlb_percent", lambda values, limit: values > limit, "crlb_above_operational_threshold"),
            ("snr_", "min_snr", lambda values, limit: values < limit, "snr_below_operational_threshold"),
            ("linewidth_fwhm_hz_", "max_linewidth_hz", lambda values, limit: values > limit, "linewidth_above_operational_threshold"),
        ):
            limit = thresholds.get(threshold_key)
            selected = [values[index] for key, values in metric_maps.items() if key.startswith(prefix) and limit is not None and np.isfinite(values[index])]
            if selected and any(comparator(value, float(limit)) for value in selected):
                flags.append(flag)
        flags_by_voxel[index] = flags
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index in np.ndindex(*valid.shape):
            okay = bool(valid[index])
            reason = "" if okay else reasons[index]
            row = {
                **entities,
                "branch": branch,
                "i": index[0],
                "j": index[1],
                "k": index[2],
                "status": fit_status if okay else (str(reason) if reason else "outside_fit_mask"),
                "reason": reason or ("outside_fit_mask" if not okay else ""),
                **{key: ("" if not np.isfinite(values[index]) else float(values[index])) for key, values in metric_maps.items()},
                "qc_flags": ";".join(flags_by_voxel[index]),
            }
            writer.writerow(row)
    definitions = {
        "raw_fit_": "FSL-MRS raw scaling values in arbitrary fit units; not molar concentrations.",
        "internal_reference_ratio_": "FSL-MRS internal scaling relative to configured/default metabolite reference; ratio uncertainty is not automatically propagated.",
        "crlb_percent_": "FSL-MRS Newton CRLB-derived percentage uncertainty; may be capped and is not absolute SD.",
        "snr_": "FSL-MRS component-specific QC SNR map.",
        "linewidth_fwhm_hz_": "FSL-MRS component-specific FWHM map in Hz.",
        "residual_l2_over_input_l2": "Euclidean norm of complex time-domain residual divided by complex input-FID norm; descriptive fit residual summary, not metabolite loss.",
        "qc_flags": "Operational threshold flags only; estimates are retained and never automatically rejected.",
    }
    (path.parent / "voxel_qc_definitions.json").write_text(json.dumps(definitions, indent=2) + "\n")


def _write_branch_comparison(original: Path, waterremoved: Path, output_dir: Path) -> None:
    with original.open(newline="") as stream:
        original_rows = {(int(row["i"]), int(row["j"]), int(row["k"])): row for row in csv.DictReader(stream, delimiter="\t")}
    with waterremoved.open(newline="") as stream:
        water_rows = {(int(row["i"]), int(row["j"]), int(row["k"])): row for row in csv.DictReader(stream, delimiter="\t")}
    numeric_fields = [
        key for key in next(iter(original_rows.values()), {})
        if key.startswith(("raw_fit_", "internal_reference_ratio_", "crlb_percent_", "snr_", "linewidth_fwhm_hz_", "nuisance_"))
        or key == "residual_l2_over_input_l2"
    ]
    rows = []
    means: dict[str, list[float]] = {key: [] for key in numeric_fields}
    paired = 0
    invalid_counts = {"original": 0, "waterremoved": 0}
    for index in sorted(set(original_rows) & set(water_rows)):
        before, after = original_rows[index], water_rows[index]
        invalid_counts["original"] += before["status"] != "fit_completed"
        invalid_counts["waterremoved"] += after["status"] != "fit_completed"
        row: dict[str, Any] = {"i": index[0], "j": index[1], "k": index[2]}
        both_valid = before["status"] == after["status"] == "fit_completed"
        if both_valid:
            paired += 1
        for field in numeric_fields:
            try:
                delta = float(after[field]) - float(before[field]) if both_valid and before[field] and after[field] else float("nan")
            except (KeyError, ValueError):
                delta = float("nan")
            row[f"waterremoved_minus_original_{field}"] = "" if not np.isfinite(delta) else delta
            if np.isfinite(delta):
                means[field].append(delta)
        rows.append(row)
    comparison_dir = original.parent.parent / "qc"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    table = comparison_dir / "water_removal_paired_differences.tsv"
    fields = list(rows[0]) if rows else ["i", "j", "k"]
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "branches": ["original", "waterremoved"],
        "paired_valid_voxel_count": paired,
        "invalid_or_excluded_voxel_count": invalid_counts,
        "paired_difference_definition": "waterremoved - original; values are available only where both branch voxels have numeric maps and fit status is completed.",
        "mean_paired_difference": {key: float(np.mean(values)) if values else None for key, values in means.items()},
        "selection": "No branch is selected automatically from fit quality or uncertainty.",
    }
    (comparison_dir / "water_removal_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")


def _write_report(
    path: Path,
    item: MRSInput,
    config: dict[str, Any],
    validated: Any,
    basis: Path,
    branches: list[str],
    roi_summary: dict[str, Any] | None,
    warnings: list[str],
) -> None:
    total = int(validated.valid_voxels.size)
    valid_count = int(validated.valid_voxels.sum())
    lines = [
        f"# Metabolix report: {item.source_name}",
        "",
        "Processing completed. This report documents software execution and QC; successful fitting does not establish validated metabolite measurements.",
        "",
        "## Input and model",
        "",
        f"- Source: `{item.path}`",
        f"- Discovery convention: {'BIDS MRS' if item.standard_bids else 'legacy *_mrsi compatibility pattern'}",
        f"- Shape: `{validated.image.shape}`; valid spatial voxels: {valid_count}/{total}",
        f"- Basis: `{basis}` (sequence/model compatibility remains a scientific assumption)",
        f"- Fitting branches: {', '.join(branches)}",
        f"- Water removal mode: `{config['processing']['water_removal']}`; no branch is selected automatically.",
        "- Estimates are not labelled as molar concentration. Internal-reference ratios require a valid nonzero reference.",
        "- FSL-MRS's standard HTML report summarizes the average FID in the mask; it is not a summary of independent voxel fits.",
        "",
        "## QC and outputs",
        "",
        "The voxel QC TSV preserves every spatial coordinate. Input-invalid voxels are flagged and excluded from the fit mask. Numeric component metrics remain blank unless they can be read unambiguously from native FSL-MRS outputs; native output maps are retained in each fit directory.",
        "",
        "Uncertainty products named `*_sd` by the Newton fit are CRLB-derived percentage estimates, not absolute standard deviations, and may be capped. They are not automatically propagated uncertainties for metabolite/reference ratios.",
        "",
        "## Warnings",
        "",
    ]
    lines.extend([f"- {warning}" for warning in warnings] or ["- None recorded by the application."])
    lines.extend(
        [
            "",
            "## Voxel-wise distributions",
            "",
            f"Reporting focus configured: {', '.join(config['report']['metabolites']) or 'all available fitted components'}.",
            "The model fit is not restricted to this reporting list. Summary statistics below use finite values from voxels marked fit-completed.",
            "",
            "| Branch | Metric | N | Mean | SD | Median | P05 | P95 |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for branch in branches:
        fit_folder = "03_fit_waterremoved" if branch == "waterremoved" else "03_fit_original"
        qc_path = path.parent / fit_folder / "voxel_qc.tsv"
        if not qc_path.is_file():
            continue
        with qc_path.open(newline="") as stream:
            qc_rows = list(csv.DictReader(stream, delimiter="\t"))
        metric_fields = [
            key for key in (qc_rows[0] if qc_rows else {})
            if key.startswith(("raw_fit_", "internal_reference_ratio_", "crlb_percent_", "snr_", "linewidth_fwhm_hz_", "nuisance_"))
            or key == "residual_l2_over_input_l2"
        ]
        report_focus = config["report"]["metabolites"]
        if report_focus:
            metric_fields = [
                key for key in metric_fields
                if key.startswith(("nuisance_", "residual_l2_over_input_l2"))
                or any(key.endswith(f"_{name}") for name in report_focus)
            ]
        completed_rows = [row for row in qc_rows if row.get("status") == "fit_completed"]
        for metric in metric_fields:
            values = []
            for row in completed_rows:
                try:
                    value = float(row[metric])
                    if np.isfinite(value):
                        values.append(value)
                except (KeyError, ValueError):
                    continue
            if values:
                array = np.asarray(values, dtype=float)
                lines.append(
                    f"| {branch} | `{metric}` | {len(values)} | {array.mean():.5g} | {array.std(ddof=0):.5g} | "
                    f"{np.median(array):.5g} | {np.percentile(array, 5):.5g} | {np.percentile(array, 95):.5g} |"
                )
        lines.extend(
            [
                "",
                f"- Native voxel maps and individual fit, baseline, and residual FIDs: `{fit_folder}/` (open with FSLeyes or another NIfTI-MRS-aware viewer).",
                f"- Paginated all-fitted-voxel spectral diagnostics: `{fit_folder}/voxel_diagnostics.pdf`.",
                f"- Voxel table and metric definitions: `{fit_folder}/voxel_qc.tsv`, `{fit_folder}/voxel_qc_definitions.json`.",
            ]
        )
    if len(branches) == 2 and (path.parent / "qc" / "water_removal_paired_differences.tsv").is_file():
        lines.extend(
            [
                "",
                "## Water-removal comparison",
                "",
                "Paired differences are water-removed minus original, voxel-matched; see `qc/water_removal_paired_differences.tsv` and `qc/water_removal_comparison.json`. No branch is selected automatically.",
            ]
        )
    if roi_summary:
        lines.extend(
            [
                "",
                "## Anatomical overlap",
                "",
                "ROI overlap is nominal geometric overlap from image-world affines. No registration was performed or validated. MRSI spatial response can spread signal beyond nominal voxel boundaries. Overlap did not restrict fitting.",
                "",
                f"Mean segmentation field-of-view coverage: {roi_summary['mean_segmentation_fov_coverage_pct']:.2f}%.",
                "Overlap-weighted estimates, when available, are saved in `roi/roi_weighted_estimates.tsv`; they are not pure-ROI concentrations or fits of an averaged ROI spectrum.",
            ]
        )
        for summary in roi_summary["roi_summaries"]:
            lines.append(
                f"- {summary['roi']}: {summary['estimated_volume_within_mrsi_grid_mm3']:.1f} mm^3 nominal overlap; "
                f"{summary['roi_fraction_covered_by_mrsi_grid_pct']:.2f}% of label-map ROI volume lies within the MRS grid."
            )
    lines.append("")
    path.write_text("\n".join(lines))


def _add_roi_columns(qc_path: Path, roi_path: Path) -> None:
    with roi_path.open(newline="") as stream:
        overlap_rows = {
            (int(row["i"]), int(row["j"]), int(row["k"])): row
            for row in csv.DictReader(stream, delimiter="\t")
        }
    with qc_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or not overlap_rows:
        return
    roi_fields = list(next(iter(overlap_rows.values())))
    fields = list(rows[0]) + [field for field in roi_fields if field not in {"i", "j", "k"} and field not in rows[0]]
    for row in rows:
        overlap = overlap_rows.get((int(row["i"]), int(row["j"]), int(row["k"])), {})
        row.update({key: value for key, value in overlap.items() if key not in {"i", "j", "k"}})
    with qc_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_roi_weighted_estimates(run_dir: Path, branches: list[str], config: dict[str, Any]) -> None:
    overlap_path = run_dir / "roi" / "voxel_overlap.tsv"
    with overlap_path.open(newline="") as stream:
        overlap = {
            (int(row["i"]), int(row["j"]), int(row["k"])): row
            for row in csv.DictReader(stream, delimiter="\t")
        }
    threshold = config["roi"].get("minimum_overlap_percent")
    roi_columns = [key for key in next(iter(overlap.values()), {}) if key.endswith("_overlap_pct")]
    rows = []
    for branch in branches:
        fit_folder = "03_fit_waterremoved" if branch == "waterremoved" else "03_fit_original"
        with (run_dir / fit_folder / "voxel_qc.tsv").open(newline="") as stream:
            voxel_rows = list(csv.DictReader(stream, delimiter="\t"))
        if not voxel_rows:
            continue
        metrics = [key for key in voxel_rows[0] if key.startswith(("raw_fit_", "internal_reference_ratio_"))]
        for roi_column in roi_columns:
            roi_name = roi_column.removesuffix("_overlap_pct")
            for metric in metrics:
                weighted_sum = 0.0
                weight_sum = 0.0
                count = 0
                for voxel in voxel_rows:
                    if voxel["status"] != "fit_completed" or not voxel.get(metric):
                        continue
                    index = (int(voxel["i"]), int(voxel["j"]), int(voxel["k"]))
                    weight = float(overlap[index][roi_column]) / 100.0
                    if weight <= 0 or (threshold is not None and weight * 100 < float(threshold)):
                        continue
                    try:
                        estimate = float(voxel[metric])
                    except ValueError:
                        continue
                    if np.isfinite(estimate):
                        weighted_sum += estimate * weight
                        weight_sum += weight
                        count += 1
                rows.append(
                    {
                        "branch": branch,
                        "roi": roi_name,
                        "metric": metric,
                        "overlap_weighted_voxel_estimate": weighted_sum / weight_sum if weight_sum else "",
                        "weighted_voxel_count": count,
                        "minimum_overlap_threshold_percent": threshold if threshold is not None else "",
                        "interpretation": "Nominal overlap-weighted voxel estimate; not a pure-ROI concentration or an averaged-spectrum fit.",
                    }
                )
    fields = ["branch", "roi", "metric", "overlap_weighted_voxel_estimate", "weighted_voxel_count", "minimum_overlap_threshold_percent", "interpretation"]
    with (run_dir / "roi" / "roi_weighted_estimates.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _process_one(
    bids_dir: Path,
    output_root: Path,
    item: MRSInput,
    config: dict[str, Any],
    dry_run: bool,
    logger: logging.Logger,
) -> None:
    validated = validate_mrs(item, logger)
    basis = _select_basis(item, validated.metadata, config)
    if basis is None and not dry_run:
        raise ValueError(
            f"No compatible FSL-MRS basis was selected for {item.path}. Fitting cannot start without an explicit basis.\n"
            "Add a compatible basis path to your YAML configuration:\n\n"
            "basis:\n"
            "  path: /absolute/path/to/compatible_basis.BASIS\n\n"
            "Or provide it on the command line:\n\n"
            "  metabolix BIDS_DIR OUTPUT_DIR participant --basis /absolute/path/to/compatible_basis.BASIS\n\n"
            "For datasets with different acquisitions, add an unambiguous basis.mapping entry instead. "
            "Check that the basis matches the sequence, echo time, field strength, and acquisition; successful resampling alone is not proof of compatibility.\n"
            "See README.md (Configuration) and docs/LIMITATIONS.md (No basis or ambiguous basis mapping)."
        )
    if basis is not None and not basis.exists():
        raise ValueError(f"Selected basis does not exist: {basis}")
    if validated.coil_dimension is not None and config["processing"]["coil_combine"] == "off":
        raise ValueError("Input contains a tagged coil dimension but coil combination was explicitly disabled.")
    coil_reference_path = config["processing"].get("coil_reference")
    if coil_reference_path:
        reference_path = Path(coil_reference_path)
        if not reference_path.is_file():
            raise ValueError(f"Configured coil reference does not exist: {reference_path}")
        if validated.coil_dimension is None:
            raise ValueError("A coil reference was provided, but the input has no tagged coil dimension to combine.")
        if reference_path.name.endswith(".nii.gz"):
            reference_sidecar = reference_path.with_name(reference_path.name[:-7] + ".json")
        else:
            reference_sidecar = reference_path.with_suffix(".json")
        reference_item = MRSInput(reference_path, reference_sidecar if reference_sidecar.is_file() else None, item.entities, False)
        reference_data = validate_mrs(reference_item, logger)
        coil_axis = validated.coil_dimension
        if reference_data.coil_dimension != coil_axis:
            raise ValueError("Coil reference and MRS input must use the same tagged coil dimension.")
        if reference_data.data.shape[coil_axis] != validated.data.shape[coil_axis]:
            raise ValueError("Coil reference and MRS input have different coil counts.")
        if reference_data.image.shape[:3] != validated.image.shape[:3] or not np.allclose(reference_data.image.affine, validated.image.affine, atol=1e-5):
            raise ValueError("Coil reference spatial grid/affine does not match the MRS input.")
        data_nucleus = validated.metadata.get("ResonantNucleus", validated.metadata.get("Nucleus"))
        reference_nucleus = reference_data.metadata.get("ResonantNucleus", reference_data.metadata.get("Nucleus"))
        if data_nucleus != reference_nucleus:
            raise ValueError("Coil reference and MRS input have different nucleus metadata.")
        data_frequency = validated.metadata["SpectrometerFrequency"]
        reference_frequency = reference_data.metadata["SpectrometerFrequency"]
        data_frequency = float(data_frequency[0] if isinstance(data_frequency, (list, tuple)) else data_frequency)
        reference_frequency = float(reference_frequency[0] if isinstance(reference_frequency, (list, tuple)) else reference_frequency)
        if not np.isclose(data_frequency, reference_frequency, rtol=1e-6):
            raise ValueError("Coil reference and MRS input have different spectrometer frequencies.")
    explicit_mask = None
    if config["fit"].get("mask"):
        mask_path = Path(config["fit"]["mask"])
        if not mask_path.is_file():
            raise ValueError(f"Explicit fit mask does not exist: {mask_path}")
        explicit_mask = nib.load(str(mask_path))
        if explicit_mask.shape != validated.image.shape[:3] or not np.allclose(explicit_mask.affine, validated.image.affine):
            raise ValueError("Explicit fit mask shape/affine does not match the MRS spatial grid.")
        if not np.any((np.asanyarray(explicit_mask.dataobj) > 0) & validated.valid_voxels):
            raise ValueError("Explicit fit mask contains no valid nonzero MRS voxels.")
    if config["roi"].get("map") and not Path(config["roi"]["map"]).is_file():
        raise ValueError(f"ROI map does not exist: {config['roi']['map']}")
    if config["roi"].get("freesurfer_dir"):
        fs_root = Path(config["roi"]["freesurfer_dir"])
        subject = f"sub-{item.entities['sub']}"
        if not any(path.is_file() for path in (fs_root / subject / "mri" / "aseg.mgz", fs_root / "mri" / "aseg.mgz")):
            raise ValueError(f"Could not find FreeSurfer aseg.mgz for {subject} under {fs_root}.")
    basis_metadata = None
    fsl_mrsi = find_executable("fsl_mrsi", config["basis"].get("binary_dir"))
    if basis is not None and fsl_mrsi:
        basis_metadata = inspect_basis(
            basis,
            1.0 / float(validated.image.header.get_zooms()[3]),
            int(validated.image.shape[3]),
            fsl_mrsi,
        )
        logger.info(
            "Basis readable through FSL-MRS: %d components, %d points, %.6g MHz; formatted to %.6g Hz x %d points.",
            basis_metadata["component_count"],
            basis_metadata["original_points"],
            basis_metadata["central_frequency_mhz"],
            basis_metadata["target_bandwidth_hz"],
            basis_metadata["target_points"],
        )
        if isinstance(validated.metadata.get("SpectrometerFrequency"), (list, tuple)):
            data_frequency = float(validated.metadata["SpectrometerFrequency"][0])
        else:
            data_frequency = float(validated.metadata["SpectrometerFrequency"])
        basis_metadata["data_central_frequency_mhz"] = data_frequency
        basis_metadata["relative_frequency_difference_percent"] = 100.0 * (data_frequency / basis_metadata["central_frequency_mhz"] - 1.0)
        basis_metadata["sequence_compatibility_verified"] = False
    if dry_run:
        logger.info("DRY RUN input %s entities=%s shape=%s basis=%s", item.path, item.entities, validated.image.shape, basis or "UNRESOLVED")
        logger.info("Requested spatial coverage: all %d valid voxels; water removal=%s", int(validated.valid_voxels.sum()), config["processing"]["water_removal"])
        run_dir = _run_directory(output_root, item)
        clean_input = run_dir / "00_work" / f"{item.source_name}_finite.nii.gz"
        planned_input = clean_input
        coil_mode = config["processing"]["coil_combine"]
        if validated.coil_dimension is not None and coil_mode in {"auto", "on"}:
            output_name = f"{item.source_name}_desc-coilcombined_mrsi"
            planned_input = run_dir / "01_coilcombine" / f"{output_name}.nii.gz"
            logger.info(
                "DRY RUN argv: %s",
                json.dumps([find_executable("fsl_mrs_proc", config["basis"].get("binary_dir")) or "fsl_mrs_proc", "coilcombine", "--file", str(clean_input), *(["--reference", str(coil_reference_path)] if coil_reference_path else []), "--output", str(run_dir / "01_coilcombine"), "--filename", output_name, "--generateReports"]),
            )
        water_mode = config["processing"]["water_removal"]
        fit_inputs = {"original": planned_input}
        if water_mode in {"on", "compare"}:
            output_name = f"{item.source_name}_desc-waterremoved_mrsi"
            water_input = run_dir / "02_waterremove" / f"{output_name}.nii.gz"
            logger.info(
                "DRY RUN argv: %s",
                json.dumps([find_executable("fsl_mrs_proc", config["basis"].get("binary_dir")) or "fsl_mrs_proc", "remove", "--file", str(planned_input), "--ppm", *map(str, config["processing"]["water_ppm"]), "--output", str(run_dir / "02_waterremove"), "--filename", output_name, "--generateReports"]),
            )
            fit_inputs = {"waterremoved": water_input} if water_mode == "on" else {"original": planned_input, "waterremoved": water_input}
        if basis is not None and fsl_mrsi:
            planned_mask = run_dir / "00_work" / f"{item.source_name}_valid_mask.nii.gz"
            for branch, data_path in fit_inputs.items():
                fit_dir = run_dir / ("03_fit_waterremoved" if branch == "waterremoved" else "03_fit_original")
                logger.info("DRY RUN argv: %s", json.dumps(_fit_command(fsl_mrsi, data_path, basis, planned_mask, fit_dir, config)))
        if basis_metadata:
            logger.info("DRY RUN basis metadata: %s", json.dumps(basis_metadata, sort_keys=True))
        return
    assert basis is not None
    fsl_version = executable_version("fsl_mrsi", config["basis"].get("binary_dir"))
    fingerprint = _fingerprint(item, basis, config, fsl_version)
    run_dir = _run_directory(output_root, item)
    manifest_path = run_dir / "manifest.json"
    resume_manifest = None
    if run_dir.exists():
        if config["execution"]["resume"] and manifest_path.exists():
            previous = json.loads(manifest_path.read_text())
            if previous.get("fingerprint") != fingerprint:
                raise ValueError(f"Resume refused for {run_dir}: input/configuration fingerprint changed.")
            if previous.get("status") == "completed":
                logger.info("Verified completed run; resume skips %s", run_dir)
                return
            if previous.get("status") not in {"running", "failed"}:
                raise ValueError(f"Cannot resume run with manifest status {previous.get('status')!r}.")
            resume_manifest = previous
        elif config["execution"]["resume"]:
            raise ValueError(f"Resume refused for {run_dir}: no manifest exists to verify inputs/configuration.")
        elif not config["execution"]["overwrite"]:
            raise ValueError(f"Output already exists: {run_dir}. Use --resume for an identical completed run or --overwrite explicitly.")
        else:
            shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "processing.log"
    handler = logging.FileHandler(log_path, mode="a" if resume_manifest else "w")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    dump_yaml(config, run_dir / "config_resolved.yaml")
    _write_dataset_description(output_root)
    manifest: dict[str, Any] = resume_manifest or {
        "status": "running",
        "created": datetime.now(timezone.utc).isoformat(),
        "source": str(item.path),
        "source_sha256": checksum(item.path),
        "source_sidecar": str(item.sidecar) if item.sidecar else None,
        "source_sidecar_sha256": checksum(item.sidecar) if item.sidecar else None,
        "basis": str(basis),
        "basis_sha256": checksum(basis),
        "basis_metadata": basis_metadata,
        "coil_reference": coil_reference_path,
        "coil_reference_sha256": checksum(Path(coil_reference_path)) if coil_reference_path else None,
        "entities": item.entities,
        "discovery_convention": "BIDS MRS" if item.standard_bids else "legacy_mrsi_compatibility",
        "metabolix_version": __version__,
        "fsl_mrs_version": fsl_version,
        "python_version": __import__("platform").python_version(),
        "fingerprint": fingerprint,
        "resolved_config": config,
        "commands": [],
        "completed_stages": [],
        "warnings": list(validated.warnings),
    }
    manifest["status"] = "running"
    commands: list[dict[str, Any]] = list(manifest.get("commands", []))
    completed_stages = set(manifest.get("completed_stages", []))

    def record_stage(name: str, command_record: dict[str, Any]) -> None:
        commands.append(command_record)
        if command_record.get("exit_status") == 0:
            completed_stages.add(name)
        manifest["commands"] = commands
        manifest["completed_stages"] = sorted(completed_stages)
        save_manifest(manifest_path, manifest)

    def require_success(stage_name: str, command_record: dict[str, Any]) -> None:
        if command_record.get("exit_status") != 0:
            raise RuntimeError(
                f"{stage_name} failed with exit status {command_record.get('exit_status')}; "
                f"see {command_record.get('stderr_file')} and {command_record.get('stdout_file')}."
            )

    save_manifest(manifest_path, manifest)
    if basis_metadata:
        (run_dir / "basis_metadata.json").write_text(json.dumps(basis_metadata, indent=2) + "\n")
    try:
        work = run_dir / "00_work"
        clean_input = work / f"{item.source_name}_finite.nii.gz"
        prepare_finite_copy(validated, clean_input)
        mask = work / f"{item.source_name}_valid_mask.nii.gz"
        fit_mask = validated.valid_voxels.copy()
        if config["fit"].get("mask"):
            assert explicit_mask is not None
            fit_mask &= np.asanyarray(explicit_mask.dataobj) > 0
        write_validity_mask(fit_mask, validated.image.affine, mask)
        if not fit_mask.any():
            raise ValueError("No valid spatial voxels remain in the fitting mask.")

        processed_input = clean_input
        has_coils = validated.coil_dimension is not None
        coil_mode = config["processing"]["coil_combine"]
        if has_coils and coil_mode == "off":
            raise ValueError("Input contains a tagged coil dimension but coil combination was explicitly disabled.")
        do_coil = has_coils and coil_mode in {"auto", "on"}
        if coil_mode == "on" and not has_coils:
            logger.warning("coil_combine is on but input has no tagged coil dimension; using the input unchanged.")
        if do_coil:
            stage = run_dir / "01_coilcombine"
            stem = f"{item.source_name}_desc-coilcombined_mrsi"
            if "coilcombine" not in completed_stages:
                if stage.exists():
                    shutil.rmtree(stage)
                argv = [_resolve_fsl("fsl_mrs_proc", config), "coilcombine", "--file", str(clean_input)]
                if coil_reference_path:
                    argv.extend(["--reference", str(coil_reference_path)])
                argv.extend(["--output", str(stage), "--filename", stem, "--generateReports"])
                command_record = run_command(argv, stage, "coilcombine", logger)
                record_stage("coilcombine", command_record)
                require_success("coilcombine", command_record)
                for warning in command_record.get("warnings", []):
                    manifest["warnings"].append(f"FSL-MRS coilcombine: {warning}")
                processed_input = _expected_nifti(stage, stem)
                _validate_stage_output(processed_input, clean_input, item.entities, True, logger)
            else:
                processed_input = _expected_nifti(stage, stem)
                _validate_stage_output(processed_input, clean_input, item.entities, True, logger)
            processed_input = _expected_nifti(stage, stem)
            manifest["coilcombine_output"] = str(processed_input)
            save_manifest(manifest_path, manifest)

        branches: dict[str, Path] = {}
        water_mode = config["processing"]["water_removal"]
        if water_mode in {"on", "compare"}:
            stage = run_dir / "02_waterremove"
            stem = f"{item.source_name}_desc-waterremoved_mrsi"
            if "waterremove" not in completed_stages:
                if stage.exists():
                    shutil.rmtree(stage)
                argv = [_resolve_fsl("fsl_mrs_proc", config), "remove", "--file", str(processed_input), "--ppm", *map(str, config["processing"]["water_ppm"]), "--output", str(stage), "--filename", stem, "--generateReports"]
                command_record = run_command(argv, stage, "waterremove", logger)
                record_stage("waterremove", command_record)
                require_success("waterremove", command_record)
            else:
                command_record = None
            water_input = _expected_nifti(stage, stem)
            _validate_stage_output(water_input, processed_input, item.entities, False, logger)
            if water_mode == "on":
                branches["waterremoved"] = water_input
            else:
                branches["original"] = processed_input
                branches["waterremoved"] = water_input
        else:
            branches["original"] = processed_input

        fit_executable = _resolve_fsl("fsl_mrsi", config)
        for branch, data_path in branches.items():
            stage = run_dir / ("03_fit_waterremoved" if branch == "waterremoved" else "03_fit_original")
            stage_name = f"fit_{branch}"
            if stage_name not in completed_stages:
                if stage.exists():
                    shutil.rmtree(stage)
                argv = _fit_command(fit_executable, data_path, basis, mask, stage, config)
                command_record = run_command(argv, stage, stage_name, logger)
                record_stage(stage_name, command_record)
                require_success(stage_name, command_record)
                fit_reasons = validated.invalid_reasons.copy()
                fit_reasons[validated.valid_voxels & ~fit_mask] = "outside_fit_mask"
                metric_maps = _load_metric_maps(stage, tuple(validated.image.shape[:3]), data_path, config, logger)
                _write_qc_table(stage / "voxel_qc.tsv", item, fit_mask, fit_reasons, branch, "fit_completed", metric_maps, config)
            elif not (stage / "voxel_qc.tsv").is_file():
                raise ValueError(f"Manifest marks {stage_name} complete but its voxel QC table is missing.")
            fsl_mrsi_runtime = find_executable("fsl_mrsi", config["basis"].get("binary_dir"))
            if fsl_mrsi_runtime and not (stage / "voxel_diagnostics.pdf").is_file():
                make_voxel_diagnostic_pdf(
                    data_path,
                    stage,
                    stage / "voxel_diagnostics.pdf",
                    mask,
                    validated.metadata,
                    fsl_mrsi_runtime,
                )

        roi_report = None
        if config["roi"].get("freesurfer_dir") or config["roi"].get("map"):
            roi_report = calculate_roi_overlap(item, validated, run_dir, config, logger)
            for branch in branches:
                fit_dir = run_dir / ("03_fit_waterremoved" if branch == "waterremoved" else "03_fit_original")
                _add_roi_columns(fit_dir / "voxel_qc.tsv", run_dir / "roi" / "voxel_overlap.tsv")
            _write_roi_weighted_estimates(run_dir, list(branches), config)
            if config["roi"].get("minimum_overlap_percent") is not None:
                threshold = float(config["roi"]["minimum_overlap_percent"])
                manifest["warnings"].append(f"ROI summary threshold {threshold:g}% is reporting-only and did not restrict fitting.")

        if "original" in branches and "waterremoved" in branches:
            _write_branch_comparison(
                run_dir / "03_fit_original" / "voxel_qc.tsv",
                run_dir / "03_fit_waterremoved" / "voxel_qc.tsv",
                run_dir,
            )

        report_path = run_dir / "report.md"
        _write_report(report_path, item, config, validated, basis, list(branches), roi_report, manifest["warnings"])
        links = []
        for branch in branches:
            fit_folder = "03_fit_waterremoved" if branch == "waterremoved" else "03_fit_original"
            links.append(f'<li><a href="{fit_folder}/voxel_diagnostics.pdf">{html.escape(branch)} voxel diagnostics (PDF)</a></li>')
        report_html = (
            "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            "<title>Metabolix report</title><style>body{font:16px/1.5 sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#222}"
            "pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f5f2;padding:1rem;border-left:4px solid #41616a}</style><body>"
            f"<h1>Metabolix report: {html.escape(item.source_name)}</h1><h2>Voxel diagnostics</h2><ul>{''.join(links)}</ul>"
            f"<pre>{html.escape(report_path.read_text())}</pre></body></html>\n"
        )
        (run_dir / "report.html").write_text(report_html)
        manifest["commands"] = commands
        manifest["status"] = "completed"
        manifest["finished"] = datetime.now(timezone.utc).isoformat()
        save_manifest(manifest_path, manifest)
    except Exception as exc:
        manifest["commands"] = commands
        manifest["status"] = "failed"
        manifest["failure"] = str(exc)
        manifest["finished"] = datetime.now(timezone.utc).isoformat()
        save_manifest(manifest_path, manifest)
        raise
    finally:
        logger.removeHandler(handler)
        handler.close()


def run_participant(
    bids_dir: Path,
    output_dir: Path,
    config: dict[str, Any],
    dry_run: bool,
    logger: logging.Logger,
) -> int:
    inputs = discover_mrs(bids_dir, _entity_filters(config))
    identities = [output_identity(item) for item in inputs]
    if len(set(identities)) != len(identities):
        raise ValueError("Multiple MRS inputs map to the same output identity; resolve duplicate names/entities before processing.")
    if not dry_run:
        dependencies = dependency_report(config)
        missing = [name for name, path in dependencies.items() if not path]
        if missing:
            raise RuntimeError(f"Missing required FSL-MRS executable(s): {', '.join(missing)}. Use --check-deps and see installation docs.")
    failures = []
    for item in inputs:
        try:
            _process_one(bids_dir, output_dir, item, config, dry_run, logger)
        except Exception as exc:
            logger.error("Failed input %s: %s", item.path, exc)
            failures.append((item.path, str(exc)))
    if failures:
        logger.error("%d/%d input(s) failed.", len(failures), len(inputs))
        return 1
    logger.info("%s %d input(s).", "Validated" if dry_run else "Processed", len(inputs))
    return 0


def run_group(output_dir: Path, config: dict[str, Any], dry_run: bool, logger: logging.Logger) -> int:
    derivative_root = output_dir
    tables = sorted(derivative_root.rglob("voxel_qc.tsv"))
    if not tables:
        raise ValueError(f"No participant voxel_qc.tsv tables found under {derivative_root}.")
    if dry_run:
        logger.info("DRY RUN group aggregation: %d participant/run table(s) found.", len(tables))
        return 0
    summary_rows = []
    group_values: dict[tuple[str, str], list[float]] = {}
    group_subjects: dict[tuple[str, str], set[str]] = {}
    group_voxel_counts: dict[tuple[str, str], int] = {}
    for table in tables:
        with table.open(newline="") as stream:
            voxel_rows = list(csv.DictReader(stream, delimiter="\t"))
        if not voxel_rows:
            continue
        subjects = sorted({row.get("entity_sub", "") for row in voxel_rows if row.get("entity_sub")})
        branches = sorted({row.get("branch", "unknown") for row in voxel_rows})
        metrics = [
            key for key in voxel_rows[0]
            if key.startswith(("raw_fit_", "internal_reference_ratio_", "crlb_percent_", "snr_", "linewidth_fwhm_hz_"))
            or key == "residual_l2_over_input_l2"
        ]
        for branch in branches:
            branch_rows = [row for row in voxel_rows if row.get("branch") == branch]
            completed = [row for row in branch_rows if row.get("status") == "fit_completed"]
            common = {
                "source_table": str(table.relative_to(derivative_root)),
                "participant": ";".join(subjects),
                "branch": branch,
                "voxel_rows": len(branch_rows),
                "fit_completed_rows": len(completed),
                "invalid_or_excluded_rows": len(branch_rows) - len(completed),
            }
            if not metrics:
                summary_rows.append({**common, "metric": "fit_status_only", "n": len(completed), "mean": "", "sd": "", "median": ""})
            for metric in metrics:
                values = []
                for row in completed:
                    try:
                        value = float(row[metric])
                        if np.isfinite(value):
                            values.append(value)
                    except (KeyError, ValueError):
                        continue
                if values:
                    array = np.asarray(values, dtype=float)
                    summary_rows.append(
                        {
                            **common,
                            "metric": metric,
                            "n": len(values),
                            "mean": float(array.mean()),
                            "sd": float(array.std(ddof=0)),
                            "median": float(np.median(array)),
                        }
                    )
                    key = (branch, metric)
                    group_values.setdefault(key, []).extend(values)
                    group_subjects.setdefault(key, set()).update(subjects)
                    group_voxel_counts[key] = group_voxel_counts.get(key, 0) + len(branch_rows)
    for (branch, metric), values in sorted(group_values.items()):
        array = np.asarray(values, dtype=float)
        summary_rows.append(
            {
                "source_table": "ALL_PARTICIPANTS",
                "participant": "GROUP",
                "branch": branch,
                "voxel_rows": group_voxel_counts[(branch, metric)],
                "fit_completed_rows": len(values),
                "invalid_or_excluded_rows": "",
                "metric": metric,
                "n": len(values),
                "mean": float(array.mean()),
                "sd": float(array.std(ddof=0)),
                "median": float(np.median(array)),
                "participant_count": len(group_subjects[(branch, metric)]),
            }
        )
    if not summary_rows:
        raise ValueError("Participant QC tables contain no rows to aggregate.")
    group_dir = derivative_root / "group"
    group_dir.mkdir(parents=True, exist_ok=True)
    summary_path = group_dir / "group_summary.tsv"
    if summary_path.exists() and not config["execution"]["overwrite"]:
        raise ValueError(f"Group summary already exists: {summary_path}. Set execution.overwrite=true to replace it.")
    fields = ["source_table", "participant", "participant_count", "branch", "voxel_rows", "fit_completed_rows", "invalid_or_excluded_rows", "metric", "n", "mean", "sd", "median"]
    with summary_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    logger.info("Aggregated %d participant/run QC table(s) into %s; no participant refitting performed.", len(tables), summary_path)
    return 0