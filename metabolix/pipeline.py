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
import yaml

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
from metabolix.output_layout import OutputLayout, bids_filename
from metabolix.reporting import generate_report_figures


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
    return OutputLayout(output_root, item).bundle_dir


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


def _publish_fsl_maps(layout: OutputLayout, branch: str, fit_dir: Path) -> list[dict[str, str]]:
    """Expose native FSL scalar maps as BIDS-named maps and report heatmap inputs."""
    definitions = (
        ("concs/raw", f"metabolix-{branch}-raw", "Raw fit scaling", "a.u."),
        ("concs/internal", f"metabolix-{branch}-internal", "Internal-reference ratio", "relative units"),
        ("uncertainties", f"metabolix-{branch}-crlb-percent", "CRLB-derived uncertainty", "percent"),
        ("qc", f"metabolix-{branch}-qc", "Fit QC metric", "metric-specific units"),
    )
    layout.map_dir.mkdir(parents=True, exist_ok=True)
    layout.figure_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for source_folder, description, label, units in definitions:
        source_dir = fit_dir / source_folder
        if not source_dir.is_dir():
            continue
        for source in sorted(source_dir.glob("*.nii*")):
            name = source.name.removesuffix(".nii.gz").removesuffix(".nii")
            if source_folder == "uncertainties" and not name.endswith("_sd"):
                continue
            if source_folder == "qc" and not name.endswith(("_snr", "_fwhm")):
                continue
            if source_folder == "uncertainties":
                name = name.removesuffix("_sd")
                kind = "crlb"
                map_description = description
                title = f"{branch}: {label} {name}"
                map_units = units
            elif source_folder == "qc":
                kind = "snr" if name.endswith("_snr") else "linewidth"
                name = name.removesuffix("_snr").removesuffix("_fwhm")
                map_description = f"{description}-{kind}"
                title = f"{branch}: component {kind} {name}"
                map_units = "component SNR metric" if kind == "snr" else "Hz"
            else:
                kind = "metabolite-map"
                map_description = description
                title = f"{branch}: {label} {name}"
                map_units = units
            metabolite = name.replace("+", "-").replace("_", "-")
            desc = f"{map_description}-{metabolite}"
            output = layout.map_dir / bids_filename(layout.item, desc, "statmap", ".nii.gz")
            image = nib.load(str(source))
            nib.save(image, str(output))
            map_metadata = {
                "Description": label,
                "Metabolite": name,
                "ProcessingBranch": branch,
                "Scaling": source_folder.removeprefix("concs/"),
                "Units": map_units,
                "SourceMap": str(source),
            }
            output.with_name(output.name[:-7] + ".json").write_text(json.dumps(map_metadata, indent=2) + "\n")
            # Concentration maps are heatmapped for immediate visual inspection; QC maps remain downloadable.
            figure_path = layout.figure_dir / bids_filename(layout.item, f"{description}-{metabolite}-heatmap", "map", ".png")
            results.append({
                "path": str(output),
                "figure_path": str(figure_path),
                "title": title,
                "branch": branch,
                "scale": source_folder,
                "metabolite": name,
                "kind": kind,
                "units": map_units,
            })
    return results


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
    definitions_path = path.with_name(path.name.replace("_voxelqc.tsv", "_voxelqc_definitions.json"))
    definitions_path.write_text(json.dumps(definitions, indent=2) + "\n")


def _write_branch_comparison(original: Path, waterremoved: Path, output_dir: Path, item: MRSInput) -> tuple[Path, Path]:
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
    comparison_dir = output_dir / "qc"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    table = comparison_dir / bids_filename(item, "metabolix-waterRemovalComparison", "voxelqc", ".tsv")
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
    summary_path = comparison_dir / bids_filename(item, "metabolix-waterRemovalComparison", "summary", ".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return table, summary_path


def _write_report(
    layout: OutputLayout,
    item: MRSInput,
    config: dict[str, Any],
    validated: Any,
    basis: Path,
    branches: list[str],
    branch_qc_paths: dict[str, Path],
    map_outputs: list[dict[str, str]],
    figures: list[dict[str, str]],
    roi_summary: dict[str, Any] | None,
    warnings: list[str],
    manifest: dict[str, Any],
) -> tuple[Path, Path]:
    report_md = layout.report_path(".md")
    report_html = layout.report_path(".html")
    valid_count = int(validated.valid_voxels.sum())
    total_count = int(validated.valid_voxels.size)
    qc_rows: dict[str, list[dict[str, str]]] = {}
    statistics: list[dict[str, Any]] = []
    for branch, qc_path in branch_qc_paths.items():
        with qc_path.open(newline="") as stream:
            qc_rows[branch] = list(csv.DictReader(stream, delimiter="\t"))
        rows = qc_rows[branch]
        columns = list(rows[0]) if rows else []
        selected = config["report"]["metabolites"]
        metric_columns = [
            column for column in columns
            if column.startswith(("raw_fit_", "internal_reference_ratio_", "crlb_percent_", "snr_", "linewidth_fwhm_hz_"))
            and (not selected or any(column.endswith(f"_{name}") for name in selected))
        ]
        for metric in metric_columns:
            values = []
            for row in rows:
                if row.get("status") != "fit_completed":
                    continue
                try:
                    value = float(row[metric])
                    if np.isfinite(value):
                        values.append(value)
                except (KeyError, ValueError):
                    continue
            if values:
                array = np.asarray(values)
                statistics.append({
                    "branch": branch,
                    "metric": metric,
                    "n": len(values),
                    "mean": float(array.mean()),
                    "sd": float(array.std(ddof=0)),
                    "median": float(np.median(array)),
                    "p05": float(np.percentile(array, 5)),
                    "p95": float(np.percentile(array, 95)),
                })

    def rel(path: str | Path) -> str:
        return Path(path).resolve().relative_to(layout.bundle_dir.resolve()).as_posix()

    repro = {
        "created": manifest.get("finished", manifest.get("created")),
        "metabolix_version": manifest.get("metabolix_version", __version__),
        "python_version": manifest.get("python_version"),
        "fsl_mrs_version": manifest.get("fsl_mrs_version"),
        "source": manifest.get("source"),
        "source_sha256": manifest.get("source_sha256"),
        "basis": str(basis),
        "basis_sha256": manifest.get("basis_sha256"),
        "basis_metadata": manifest.get("basis_metadata"),
        "fingerprint": manifest.get("fingerprint"),
        "commands": manifest.get("commands", []),
    }
    basis_metadata = manifest.get("basis_metadata") or {}
    basis_facts = []
    if basis_metadata:
        basis_facts = [
            f"Basis sampling: {basis_metadata.get('original_points', 'unknown')} points at "
            f"{basis_metadata.get('original_bandwidth_hz', 'unknown')} Hz; fitted data formatted to "
            f"{basis_metadata.get('target_points', 'unknown')} points at {basis_metadata.get('target_bandwidth_hz', 'unknown')} Hz.",
            f"Basis frequency: {basis_metadata.get('central_frequency_mhz', 'unknown')} MHz; data frequency: "
            f"{basis_metadata.get('data_central_frequency_mhz', 'unknown')} MHz "
            f"(difference {basis_metadata.get('relative_frequency_difference_percent', 'unknown')}%).",
            "Basis component names: " + ", ".join(basis_metadata.get("component_names", [])),
            "Sequence/pulse-sequence equivalence is not verified by successful parsing or resampling.",
        ]
    config_text = yaml.safe_dump(config, sort_keys=False, allow_unicode=False)
    qc_figures = [entry for entry in figures if entry["kind"] == "preprocessing-qc"]
    basis_figures = [entry for entry in figures if entry["kind"] == "basis"]
    map_figures = [entry for entry in figures if entry["kind"] == "metabolite-heatmap"]

    md_lines = [
        f"# Metabolix report: {layout.stem}",
        "",
        "> Processing/QC report. Successful fitting does not establish validated metabolite measurements.",
        "",
        "## Quick links",
        "",
        f"- [All voxel maps](maps/) and [voxel QC tables](qc/)",
        f"- [Figures](figures/) and [native FSL-MRS products](processing/fsl-mrs/)",
        f"- [Reproducibility/configuration records](provenance/)",
        "",
        "## Acquisition and fitting summary",
        "",
        f"- Source: `{item.path}` ({'standard BIDS MRS' if item.standard_bids else 'legacy *_mrsi compatibility layout'})",
        f"- Spatial/spectral shape: `{validated.image.shape}`; valid voxels: {valid_count}/{total_count}",
        f"- Basis: `{basis}` (sequence compatibility must be assessed independently)",
        *[f"- {fact}" for fact in basis_facts],
        f"- Fit branches: {', '.join(branches)}; water-removal mode: `{config['processing']['water_removal']}`",
        f"- Report focus: {', '.join(config['report']['metabolites']) or 'all available metabolites'}; this does not restrict the fitted model.",
        "- Internal-reference outputs are relative estimates, not molar concentrations. Raw fit values are arbitrary units.",
        "",
        "## Preprocessing QC figures",
        "",
    ]
    for entry in qc_figures:
        caption = entry["title"]
        md_lines.append(f"### {caption}")
        md_lines.append("")
        md_lines.append(f"![{caption}]({rel(entry['path'])})")
        md_lines.append("")
    md_lines.extend(["## Basis figures", ""])
    for entry in basis_figures:
        md_lines.append(f"![{entry['title']}]({rel(entry['path'])})")
        md_lines.append("")
    md_lines.extend(["## Metabolite map heatmaps", ""])
    for entry in map_figures:
        md_lines.append(f"### {entry['title']}")
        md_lines.append("")
        md_lines.append(f"![{entry['title']}]({rel(entry['path'])})")
        md_lines.append(f"[Download NIfTI map]({rel(entry['map_path'])})")
        md_lines.append("")
    md_lines.extend(["## Map files", ""])
    for entry in map_outputs:
        md_lines.append(f"- [{entry['title']}]({rel(entry['path'])})")
    md_lines.append("")
    md_lines.extend(["## All-voxel diagnostics", ""])
    for entry in figures:
        if entry["kind"] == "voxel-diagnostics":
            md_lines.append(f"- [{entry['title']}]({rel(entry['path'])})")
    md_lines.append("")
    md_lines.extend(["## Voxel distributions", "", "| Branch | Metric | N | Mean | SD | Median | P05 | P95 |", "|---|---|---:|---:|---:|---:|---:|---:|"])
    for row in statistics:
        md_lines.append(
            f"| {row['branch']} | `{row['metric']}` | {row['n']} | {row['mean']:.5g} | {row['sd']:.5g} | "
            f"{row['median']:.5g} | {row['p05']:.5g} | {row['p95']:.5g} |"
        )
    md_lines.extend(["", "## Warnings", ""])
    md_lines.extend([f"- {warning}" for warning in warnings] or ["- None recorded."])
    if roi_summary:
        md_lines.extend([
            "", "## Anatomical overlap", "",
            "Nominal overlap is computed from existing image-world affines; registration is not performed or validated. MRSI spatial response may extend beyond nominal voxel boundaries. ROI overlap does not restrict fitting.",
            "", f"Mean segmentation field-of-view coverage: {roi_summary['mean_segmentation_fov_coverage_pct']:.2f}%.",
            "Overlap-weighted estimates are descriptive voxel-weighted summaries, not pure-ROI concentrations or fits of an averaged ROI spectrum.",
        ])
    md_lines.extend(["", "## Reproducibility", "", "Full commands, checksums, software versions, resolved configuration, and processing status are in the linked provenance files.", ""])
    report_md.parent.mkdir(parents=True, exist_ok=True)
    report_md.write_text("\n".join(md_lines) + "\n")

    quick_links = [
        ("QC folder", "qc/"),
        ("Map folder", "maps/"),
        ("Figures", "figures/"),
        ("Native FSL-MRS outputs", "processing/fsl-mrs/"),
        ("Logs", "logs/"),
        ("Provenance and resolved configuration", "provenance/"),
    ]
    link_html = " ".join(f'<a class="quick-link" href="{href}">{label}</a>' for label, href in quick_links)
    qc_links = "".join(
        f'<li><a href="{html.escape(rel(path))}">{html.escape(branch)} voxel QC table</a></li>'
        for branch, path in branch_qc_paths.items()
    )
    map_links = "".join(
        f'<li><a href="{html.escape(rel(entry["path"]))}">{html.escape(entry["title"])}</a></li>'
        for entry in map_outputs
    )
    figure_html = []
    for entry in qc_figures:
        figure_html.append(
            f'<figure><a href="{html.escape(rel(entry["path"]))}"><img loading="lazy" src="{html.escape(rel(entry["path"]))}" alt="{html.escape(entry["title"])}"></a>'
            f'<figcaption>{html.escape(entry["title"])}</figcaption></figure>'
        )
    heatmap_html = []
    for entry in map_figures:
        heatmap_html.append(
            f'<figure><a href="{html.escape(rel(entry["path"]))}"><img loading="lazy" src="{html.escape(rel(entry["path"]))}" alt="{html.escape(entry["title"])}"></a>'
            f'<figcaption><a href="{html.escape(rel(entry["map_path"]))}">{html.escape(entry["title"])} map</a></figcaption></figure>'
        )
    voxel_pdf_html = "".join(
        f'<li><a href="{html.escape(rel(entry["path"]))}">{html.escape(entry["title"])}</a></li>'
        for entry in figures if entry["kind"] == "voxel-diagnostics"
    )
    comparison_html = ""
    for path in manifest.get("water_removal_comparison", []):
        comparison_html += f'<li><a href="{html.escape(rel(path))}">{html.escape(Path(path).name)}</a></li>'
    stats_rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row[key]))}</td>" for key in ("branch", "metric", "n")) +
        "".join(f"<td>{row[key]:.5g}</td>" for key in ("mean", "sd", "median", "p05", "p95")) + "</tr>"
        for row in statistics
    )
    warning_html = "".join(f"<li>{html.escape(warning)}</li>" for warning in warnings) or "<li>None recorded.</li>"
    roi_html = ""
    if roi_summary:
        roi_links = []
        for label, key in (("Per-voxel overlap table", "overlap_table"), ("Overlap-weighted summaries", "weighted_estimates_table")):
            if roi_summary.get(key):
                roi_links.append(f'<li><a href="{html.escape(rel(roi_summary[key]))}">{label}</a></li>')
        for roi_map in sorted(layout.roi_dir.glob("*_statmap.nii.gz")):
            roi_links.append(f'<li><a href="{html.escape(rel(roi_map))}">{html.escape(roi_map.name)}</a></li>')
        roi_html = (
            "<section><h2>Anatomical overlap</h2><p>Nominal overlap from existing image-world affines; registration is not performed or validated. "
            "MRSI spatial response may extend beyond nominal voxel boundaries. ROI overlap does not restrict fitting.</p>"
            f"<p>Mean segmentation field-of-view coverage: {roi_summary['mean_segmentation_fov_coverage_pct']:.2f}%.</p>"
            "<p>Overlap-weighted estimates are descriptive, not pure-ROI concentrations or fits of an averaged ROI spectrum.</p>"
            f"<ul>{''.join(roi_links)}</ul></section>"
        )
    html_document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Metabolix report {html.escape(layout.stem)}</title>
<style>
:root{{--ink:#192b32;--muted:#5b6b70;--line:#d8e0de;--wash:#f3f7f5;--accent:#176b68;--signal:#bd5b35}}
*{{box-sizing:border-box}}body{{margin:0;color:var(--ink);font:15px/1.55 system-ui,sans-serif;background:white}}
header{{padding:2rem clamp(1rem,4vw,4rem);background:var(--wash);border-bottom:1px solid var(--line)}}main{{max-width:1500px;margin:auto;padding:1rem clamp(1rem,4vw,4rem) 4rem}}
h1{{font-size:1.8rem;margin:.15rem 0}}h2{{font-size:1.3rem;margin:2rem 0 .6rem}}h3{{font-size:1rem;margin:.5rem 0;color:var(--muted)}}p,.muted{{color:var(--muted)}}
.notice{{border-left:4px solid var(--signal);padding:.8rem 1rem;background:#fff7f2}}.quick-link{{display:inline-block;margin:.25rem .35rem .25rem 0;padding:.45rem .7rem;border:1px solid var(--line);color:var(--accent);text-decoration:none;font-weight:650}}
.facts{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:.5rem 1.5rem;padding:1rem 0}}.facts div{{border-bottom:1px solid var(--line);padding:.5rem 0}}.facts b{{display:block;font-size:.78rem;text-transform:uppercase;color:var(--muted)}}
.gallery{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,390px),1fr));gap:1rem}}figure{{margin:0;border:1px solid var(--line);background:white}}figure img{{display:block;width:100%;height:auto}}figcaption{{padding:.65rem .8rem;color:var(--muted)}}
.maps{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,260px),1fr));gap:.75rem}}.maps figure{{font-size:.9rem}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;display:block;overflow-x:auto}}th,td{{text-align:left;border-bottom:1px solid var(--line);padding:.45rem .6rem;white-space:nowrap}}th{{background:var(--wash)}}
details{{margin:1rem 0;border:1px solid var(--line);padding:.8rem}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f6f8f7;padding:1rem;max-height:34rem;overflow:auto}}
@media(max-width:600px){{header{{padding:1.25rem 1rem}}main{{padding:0 1rem 2rem}}h1{{font-size:1.4rem}}}}
</style></head><body>
<header><p class="muted">METABOLIX / MRS DERIVATIVE QC</p><h1>{html.escape(layout.stem)}</h1>
<p class="notice">Processing completion is not evidence of validated metabolite measurements. Interpret spectral QC, fit uncertainty, anatomical overlap, and quantification units separately.</p>
{link_html}</header><main>
<section><h2>Acquisition and fitting</h2><div class="facts">
<div><b>Source</b>{html.escape(str(item.path))}</div><div><b>Discovery</b>{'Standard BIDS MRS' if item.standard_bids else 'Legacy *_mrsi compatibility layout'}</div>
<div><b>Shape</b>{html.escape(str(validated.image.shape))}</div><div><b>Valid spatial voxels</b>{valid_count} / {total_count}</div>
<div><b>Basis</b>{html.escape(basis.name)}</div><div><b>Branches</b>{html.escape(', '.join(branches))}</div>
<div><b>Water removal</b>{html.escape(config['processing']['water_removal'])}</div><div><b>Reporting focus</b>{html.escape(', '.join(config['report']['metabolites']) or 'all fitted components')}</div>
</div><p>Reporting focus changes report emphasis only; it does not remove model components. FSL-MRS internal values are relative estimates. Raw fit scaling is arbitrary and is not a molar concentration.</p></section>
<section><h2>Preprocessing and spectral QC</h2><p>Curves summarize voxel/coil magnitude spectra with the 10th-90th percentile band. The time-domain panel shows mean FID magnitude, not a coherent complex mean across spatial voxels.</p><div class="gallery">{''.join(figure_html)}</div></section>
<section><h3>All-voxel fit diagnostics</h3><p>Each PDF contains every fitted voxel's spectrum, fit, baseline, and residual, paginated for browsing.</p><ul>{voxel_pdf_html}</ul></section>
<section><h2>Basis used</h2><p>{html.escape('; '.join(basis_facts))}</p><p>The normalized component spectra are shown for inspection; their presence does not establish sequence equivalence or quantitative validity.</p><div class="gallery">{''.join(f'<figure><a href="{html.escape(rel(entry["path"]))}"><img loading="lazy" src="{html.escape(rel(entry["path"]))}" alt="{html.escape(entry["title"])}"></a><figcaption>{html.escape(entry["title"])}</figcaption></figure>' for entry in figures if "basis" in entry["kind"])}</div></section>
<section><h2>Metabolite maps</h2><p>Heatmaps and downloadable NIfTI maps use the BIDS entity prefix. Review the scale and units in each map filename/report entry.</p><div class="maps">{''.join(heatmap_html)}</div></section>
<section><h3>Downloadable NIfTI maps</h3><ul>{map_links}</ul></section>
<section><h2>Voxel QC tables</h2><ul>{qc_links}</ul></section>
{('<section><h2>Water-removal comparison</h2><p>Paired differences are water-removed minus original. No branch is selected automatically.</p><ul>' + comparison_html + '</ul></section>') if comparison_html else ''}
<section><h2>Voxel-wise distributions</h2><table><thead><tr><th>Branch</th><th>Metric</th><th>N</th><th>Mean</th><th>SD</th><th>Median</th><th>P05</th><th>P95</th></tr></thead><tbody>{stats_rows}</tbody></table></section>
{roi_html}
<section><h2>Processing warnings</h2><ul>{warning_html}</ul></section>
<section><h2>Reproducibility</h2><p>Full per-stage provenance, checksums, resolved configuration, and exact command argument arrays are available in <a href="provenance/">provenance/</a>. Native FSL outputs and process logs are retained under <a href="processing/fsl-mrs/">processing/fsl-mrs/</a> and <a href="logs/">logs/</a>.</p>
<details><summary>Configuration summary</summary><pre>{html.escape(config_text)}</pre></details>
<details><summary>Run provenance summary</summary><pre>{html.escape(json.dumps(repro, indent=2, default=str))}</pre></details></section>
</main></body></html>'''
    report_html.write_text(html_document + "\n")
    return report_md, report_html


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


def _write_roi_weighted_estimates(layout: OutputLayout, branches: list[str], config: dict[str, Any]) -> Path:
    overlap_path = layout.roi_dir / bids_filename(layout.item, "metabolix-roi-overlap", "voxelqc", ".tsv")
    with overlap_path.open(newline="") as stream:
        overlap = {
            (int(row["i"]), int(row["j"]), int(row["k"])): row
            for row in csv.DictReader(stream, delimiter="\t")
        }
    threshold = config["roi"].get("minimum_overlap_percent")
    roi_columns = [key for key in next(iter(overlap.values()), {}) if key.endswith("_overlap_pct")]
    rows = []
    for branch in branches:
        with layout.branch_qc_table(branch).open(newline="") as stream:
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
    output_path = layout.roi_dir / bids_filename(layout.item, "metabolix-roi-weighted", "summary", ".tsv")
    with output_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _process_one(
    bids_dir: Path,
    output_root: Path,
    item: MRSInput,
    config: dict[str, Any],
    dry_run: bool,
    logger: logging.Logger,
) -> None:
    validated = validate_mrs(item, logger)
    layout = OutputLayout(output_root, item)
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
        run_dir = layout.bundle_dir
        clean_input = layout.work_dir / bids_filename(item, "metabolix-validInput", "mrs", ".nii.gz")
        planned_input = clean_input
        coil_mode = config["processing"]["coil_combine"]
        if validated.coil_dimension is not None and coil_mode in {"auto", "on"}:
            output_name = f"{layout.stem}_desc-coilcombined_mrs"
            planned_input = layout.stage_dir("coilcombine") / f"{output_name}.nii.gz"
            logger.info(
                "DRY RUN argv: %s",
                json.dumps([find_executable("fsl_mrs_proc", config["basis"].get("binary_dir")) or "fsl_mrs_proc", "coilcombine", "--file", str(clean_input), *( ["--reference", str(coil_reference_path)] if coil_reference_path else [] ), "--output", str(layout.stage_dir("coilcombine")), "--filename", output_name, "--generateReports"]),
            )
        water_mode = config["processing"]["water_removal"]
        fit_inputs = {"original": planned_input}
        if water_mode in {"on", "compare"}:
            output_name = f"{layout.stem}_desc-waterremoved_mrs"
            water_input = layout.stage_dir("waterremove") / f"{output_name}.nii.gz"
            logger.info(
                "DRY RUN argv: %s",
                json.dumps([find_executable("fsl_mrs_proc", config["basis"].get("binary_dir")) or "fsl_mrs_proc", "remove", "--file", str(planned_input), "--ppm", *map(str, config["processing"]["water_ppm"]), "--output", str(layout.stage_dir("waterremove")), "--filename", output_name, "--generateReports"]),
            )
            fit_inputs = {"waterremoved": water_input} if water_mode == "on" else {"original": planned_input, "waterremoved": water_input}
        if basis is not None and fsl_mrsi:
            planned_mask = layout.work_dir / bids_filename(item, "metabolix-fitMask", "mask", ".nii.gz")
            for branch, data_path in fit_inputs.items():
                fit_dir = layout.branch_dir(branch)
                logger.info("DRY RUN argv: %s", json.dumps(_fit_command(fsl_mrsi, data_path, basis, planned_mask, fit_dir, config)))
        if basis_metadata:
            logger.info("DRY RUN basis metadata: %s", json.dumps(basis_metadata, sort_keys=True))
        return
    assert basis is not None
    fsl_version = executable_version("fsl_mrsi", config["basis"].get("binary_dir"))
    fingerprint = _fingerprint(item, basis, config, fsl_version)
    run_dir = layout.bundle_dir
    manifest_path = layout.manifest_path()
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
    layout.provenance_dir.mkdir(parents=True, exist_ok=True)
    layout.log_dir.mkdir(parents=True, exist_ok=True)
    log_path = layout.log_path()
    handler = logging.FileHandler(log_path, mode="a" if resume_manifest else "w")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    dump_yaml(config, layout.config_path())
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
        (layout.provenance_dir / bids_filename(item, "metabolix-basis", "metadata", ".json")).write_text(json.dumps(basis_metadata, indent=2) + "\n")
    try:
        work = layout.work_dir
        work.mkdir(parents=True, exist_ok=True)
        clean_input = work / bids_filename(item, "metabolix-validInput", "mrs", ".nii.gz")
        prepare_finite_copy(validated, clean_input)
        mask = work / bids_filename(item, "metabolix-fitMask", "mask", ".nii.gz")
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
            stage = layout.stage_dir("coilcombine")
            stem = f"{layout.stem}_desc-coilcombined_mrs"
            if "coilcombine" not in completed_stages:
                if stage.exists():
                    shutil.rmtree(stage)
                stage.parent.mkdir(parents=True, exist_ok=True)
                argv = [_resolve_fsl("fsl_mrs_proc", config), "coilcombine", "--file", str(clean_input)]
                if coil_reference_path:
                    argv.extend(["--reference", str(coil_reference_path)])
                argv.extend(["--output", str(stage), "--filename", stem, "--generateReports"])
                command_record = run_command(argv, stage, "coilcombine", logger, log_dir=layout.log_dir, log_stem=layout.stem)
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
            stage = layout.stage_dir("waterremove")
            stem = f"{layout.stem}_desc-waterremoved_mrs"
            if "waterremove" not in completed_stages:
                if stage.exists():
                    shutil.rmtree(stage)
                stage.parent.mkdir(parents=True, exist_ok=True)
                argv = [_resolve_fsl("fsl_mrs_proc", config), "remove", "--file", str(processed_input), "--ppm", *map(str, config["processing"]["water_ppm"]), "--output", str(stage), "--filename", stem, "--generateReports"]
                command_record = run_command(argv, stage, "waterremove", logger, log_dir=layout.log_dir, log_stem=layout.stem)
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
        branch_qc_paths: dict[str, Path] = {}
        map_outputs: list[dict[str, str]] = []
        diagnostic_paths: dict[str, Path] = {}
        for branch, data_path in branches.items():
            stage = layout.branch_dir(branch)
            stage_name = f"fit_{branch}"
            if stage_name not in completed_stages:
                if stage.exists():
                    shutil.rmtree(stage)
                stage.parent.mkdir(parents=True, exist_ok=True)
                argv = _fit_command(fit_executable, data_path, basis, mask, stage, config)
                command_record = run_command(argv, stage, stage_name, logger, log_dir=layout.log_dir, log_stem=layout.stem)
                record_stage(stage_name, command_record)
                require_success(stage_name, command_record)
                fit_reasons = validated.invalid_reasons.copy()
                fit_reasons[validated.valid_voxels & ~fit_mask] = "outside_fit_mask"
                metric_maps = _load_metric_maps(stage, tuple(validated.image.shape[:3]), data_path, config, logger)
                qc_path = layout.branch_qc_table(branch)
                _write_qc_table(qc_path, item, fit_mask, fit_reasons, branch, "fit_completed", metric_maps, config)
            else:
                qc_path = layout.branch_qc_table(branch)
                if not qc_path.is_file():
                    raise ValueError(f"Manifest marks {stage_name} complete but its BIDS-named voxel QC table is missing.")
            branch_qc_paths[branch] = qc_path
            map_outputs.extend(_publish_fsl_maps(layout, branch, stage))
            fsl_mrsi_runtime = find_executable("fsl_mrsi", config["basis"].get("binary_dir"))
            diagnostic_pdf = layout.figure_dir / bids_filename(item, f"metabolix-{branch}-voxelDiagnostics", "report", ".pdf")
            diagnostic_paths[branch] = diagnostic_pdf
            if fsl_mrsi_runtime and not diagnostic_pdf.is_file():
                make_voxel_diagnostic_pdf(
                    data_path,
                    stage,
                    diagnostic_pdf,
                    mask,
                    validated.metadata,
                    fsl_mrsi_runtime,
                )

        roi_report = None
        if config["roi"].get("freesurfer_dir") or config["roi"].get("map"):
            roi_report = calculate_roi_overlap(item, validated, run_dir, config, logger)
            for branch in branches:
                _add_roi_columns(branch_qc_paths[branch], Path(roi_report["overlap_table"]))
            roi_report["weighted_estimates_table"] = str(_write_roi_weighted_estimates(layout, list(branches), config))
            if config["roi"].get("minimum_overlap_percent") is not None:
                threshold = float(config["roi"]["minimum_overlap_percent"])
                manifest["warnings"].append(f"ROI summary threshold {threshold:g}% is reporting-only and did not restrict fitting.")

        if "original" in branches and "waterremoved" in branches:
            comparison_paths = _write_branch_comparison(
                branch_qc_paths["original"],
                branch_qc_paths["waterremoved"],
                run_dir,
                item,
            )
            manifest["water_removal_comparison"] = [str(path) for path in comparison_paths]

        spectral_inputs = {}
        if validated.coil_dimension is not None and do_coil:
            spectral_inputs["Before coil combination"] = str(item.path)
            spectral_inputs["After coil combination"] = str(processed_input)
        else:
            spectral_inputs["Input / fitting data"] = str(processed_input)
        if water_mode in {"on", "compare"}:
            spectral_inputs["Before water removal"] = str(processed_input)
            spectral_inputs["After water removal"] = str(next(path for name, path in branches.items() if name == "waterremoved"))
        figures = []
        fsl_mrsi_runtime = find_executable("fsl_mrsi", config["basis"].get("binary_dir"))
        if fsl_mrsi_runtime:
            try:
                heatmap_candidates = [entry for entry in map_outputs if entry["scale"] == "concs/internal"]
                if not heatmap_candidates:
                    heatmap_candidates = [entry for entry in map_outputs if entry["scale"] == "concs/raw"]
                figures = generate_report_figures(
                    item,
                    layout.figure_dir,
                    layout.work_dir,
                    fsl_mrsi_runtime,
                    validated.metadata,
                    spectral_inputs,
                    str(basis),
                    heatmap_candidates,
                )
            except Exception as exc:
                warning = f"Application QC figure generation failed: {exc}"
                manifest["warnings"].append(warning)
                logger.warning(warning)
        for branch, path in diagnostic_paths.items():
            if path.is_file():
                figures.append({"path": str(path), "title": f"{branch} all-voxel fit diagnostics", "kind": "voxel-diagnostics"})
        report_md, report_html = _write_report(
            layout,
            item,
            config,
            validated,
            basis,
            list(branches),
            branch_qc_paths,
            map_outputs,
            figures,
            roi_report,
            manifest["warnings"],
            manifest,
        )
        manifest["reports"] = {"html": str(report_html), "markdown": str(report_md)}
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
    tables = sorted(
        path for path in derivative_root.rglob("*_voxelqc.tsv")
        if "waterRemovalComparison" not in path.name
    )
    if not tables:
        raise ValueError(f"No BIDS-named *_voxelqc.tsv participant tables found under {derivative_root}.")
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