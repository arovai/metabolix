"""NIfTI-MRS and basis validation helpers."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from metabolix.discovery import MRSInput


@dataclass
class ValidatedInput:
    image: nib.spatialimages.SpatialImage
    metadata: dict[str, Any]
    data: np.ndarray
    valid_voxels: np.ndarray
    invalid_reasons: np.ndarray
    coil_dimension: int | None
    warnings: list[str]


def _embedded_metadata(image: nib.spatialimages.SpatialImage) -> dict[str, Any]:
    for extension in image.header.extensions:
        try:
            content = extension.get_content()
            if isinstance(content, bytes):
                content = content.decode("utf-8")
            value = json.loads(content)
            if isinstance(value, dict):
                return value
        except (UnicodeDecodeError, json.JSONDecodeError, AttributeError):
            continue
    return {}


def _dimension_tag(metadata: dict[str, Any], axis: int) -> str | None:
    tag = metadata.get(f"dim_{axis + 1}")
    if isinstance(tag, dict):
        tag = tag.get("tag")
    return str(tag).upper() if tag is not None else None


def validate_mrs(item: MRSInput, logger: logging.Logger) -> ValidatedInput:
    try:
        image = nib.load(str(item.path))
    except Exception as exc:
        raise ValueError(f"Could not read NIfTI-MRS file {item.path}: {exc}") from exc
    metadata = _embedded_metadata(image)
    if item.sidecar:
        try:
            metadata.update(json.loads(item.sidecar.read_text()))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Invalid JSON sidecar {item.sidecar}: {exc}") from exc
    if len(image.shape) < 4 or len(image.shape) > 7:
        raise ValueError(f"Expected 4-7 NIfTI-MRS dimensions; got shape {image.shape} for {item.path}.")
    if not np.issubdtype(image.get_data_dtype(), np.complexfloating):
        raise ValueError(f"NIfTI-MRS input must contain complex data: {item.path}")
    affine = np.asarray(image.affine, dtype=float)
    if affine.shape != (4, 4) or not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3, :3])) < 1e-8:
        raise ValueError(f"Invalid or singular image affine: {item.path}")
    if image.header.get_zooms()[3] <= 0:
        raise ValueError(f"Spectral dwell time must be positive: {item.path}")
    nucleus = metadata.get("Nucleus", metadata.get("ResonantNucleus"))
    frequency = metadata.get("SpectrometerFrequency")
    if not nucleus or not frequency:
        raise ValueError(f"NIfTI-MRS metadata must include ResonantNucleus and SpectrometerFrequency: {item.path}")
    if isinstance(frequency, (list, tuple)):
        frequency = frequency[0] if frequency else None
    if frequency is None or not np.isfinite(float(frequency)) or float(frequency) <= 0:
        raise ValueError(f"SpectrometerFrequency must be finite and positive: {item.path}")
    if isinstance(nucleus, (list, tuple)):
        nucleus = nucleus[0] if nucleus else ""
    if str(nucleus).replace("^", "").strip() not in {"1H", "H1"}:
        raise ValueError(f"Only proton MRS (1H) is currently supported; found {nucleus!r}.")
    kspace = metadata.get("kSpace", metadata.get("KSpace", False))
    kspace_encoded = any(kspace) if isinstance(kspace, (list, tuple)) else bool(kspace)
    repr_value = str(metadata.get("SpatialRepresentation", "")).lower()
    if "k-space" in repr_value or "kspace" in repr_value or kspace_encoded:
        raise ValueError(f"K-space input is unsupported; provide image-space NIfTI-MRS: {item.path}")

    coil_axes = []
    for axis in range(4, len(image.shape)):
        tag = _dimension_tag(metadata, axis)
        if tag == "DIM_COIL":
            coil_axes.append(axis)
        elif tag is None:
            raise ValueError(f"Higher dimension {axis + 1} has no NIfTI-MRS tag; cannot safely infer its meaning.")
        else:
            raise ValueError(f"Unsupported NIfTI-MRS dimension tag {tag!r} on dimension {axis + 1}; dynamics/editing are not averaged.")
    if len(coil_axes) > 1:
        raise ValueError("Multiple coil dimensions are unsupported.")
    coil_dimension = coil_axes[0] if coil_axes else None

    qform, qcode = image.get_qform(coded=True)
    sform, scode = image.get_sform(coded=True)
    warnings: list[str] = []
    if qcode and scode and not np.allclose(qform, sform, atol=0.01):
        warnings.append("qform and sform disagree; selected image affine is used, anatomical geometry needs review.")
    if not metadata.get("EchoTime"):
        warnings.append("EchoTime is absent from NIfTI-MRS metadata.")
    if not metadata.get("RepetitionTime"):
        warnings.append("RepetitionTime is absent from NIfTI-MRS metadata.")
    for key in ("EchoTime", "RepetitionTime"):
        value = metadata.get(key)
        if value is not None and (not np.isfinite(float(value)) or float(value) <= 0):
            raise ValueError(f"{key} must be finite and positive when present: {item.path}")
    spectral_width = metadata.get("SpectralWidth")
    if spectral_width is not None:
        width = float(spectral_width)
        measured = 1.0 / float(image.header.get_zooms()[3])
        if not np.isfinite(width) or width <= 0:
            raise ValueError(f"SpectralWidth must be finite and positive: {item.path}")
        if not np.isclose(width, measured, rtol=1e-3):
            warnings.append(f"SpectralWidth metadata ({width:g} Hz) differs from NIfTI dwell-derived width ({measured:g} Hz).")

    data = np.asanyarray(image.dataobj)
    finite = np.isfinite(data)
    voxel_finite = finite.all(axis=tuple(range(3, data.ndim)))
    voxel_nonzero = np.any(np.where(finite, data, 0) != 0, axis=tuple(range(3, data.ndim)))
    valid = voxel_finite & voxel_nonzero
    reasons = np.full(data.shape[:3], "", dtype=object)
    reasons[~voxel_finite] = "nonfinite_signal"
    reasons[voxel_finite & ~voxel_nonzero] = "zero_signal"
    logger.info("Validated %s: %s; valid spatial voxels %d/%d", item.path, data.shape, int(valid.sum()), valid.size)
    for warning in warnings:
        logger.warning("%s: %s", item.path.name, warning)
    return ValidatedInput(image, metadata, data, valid, reasons, coil_dimension, warnings)


def prepare_finite_copy(validated: ValidatedInput, destination: Path) -> None:
    """Write a derived copy with invalid spatial FIDs zeroed for stable external processing."""
    data = np.array(validated.data, copy=True)
    data[~validated.valid_voxels, ...] = 0
    header = validated.image.header.copy()
    header.extensions.clear()
    metadata_json = json.dumps(validated.metadata, separators=(",", ":"), allow_nan=False).encode("utf-8")
    header.extensions.append(nib.nifti1.Nifti1Extension(44, metadata_json))
    clean = nib.Nifti1Image(data, validated.image.affine, header=header)
    clean.header.set_qform(*validated.image.get_qform(coded=True))
    clean.header.set_sform(*validated.image.get_sform(coded=True))
    destination.parent.mkdir(parents=True, exist_ok=True)
    nib.save(clean, str(destination))


def write_validity_mask(valid: np.ndarray, affine: np.ndarray, path: Path) -> None:
    image = nib.Nifti1Image(valid.astype(np.uint8), affine)
    image.header.set_xyzt_units("mm")
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(image, str(path))


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    if path.is_dir():
        for member in sorted(candidate for candidate in path.rglob("*") if candidate.is_file()):
            digest.update(str(member.relative_to(path)).encode())
            with member.open("rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
        return digest.hexdigest()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_basis(path: Path, bandwidth_hz: float, points: int, fsl_executable: str) -> dict[str, Any]:
    """Read and format a basis using the Python runtime paired with FSL-MRS."""
    executable_dir = Path(fsl_executable).resolve().parent
    python_candidates = [executable_dir / "python3", executable_dir / "python", executable_dir / "python3.12"]
    interpreter = next((candidate for candidate in python_candidates if candidate.is_file() and os.access(candidate, os.X_OK)), None)
    if interpreter is None:
        raise ValueError(f"Could not locate a Python interpreter alongside FSL-MRS executable {fsl_executable}.")
    program = r"""
import json, sys
import numpy as np
from fsl_mrs.utils import mrs_io
basis = mrs_io.read_basis(sys.argv[1])
original = np.asarray(basis.original_basis_array)
formatted = np.asarray(basis.get_formatted_basis(float(sys.argv[2]), int(sys.argv[3])))
if not np.isfinite(original).all() or not np.isfinite(formatted).all():
    raise ValueError('Basis contains nonfinite values after target sampling conversion')
print(json.dumps({
    'component_names': list(basis.names),
    'central_frequency_mhz': float(basis.cf),
    'original_points': int(original.shape[0]),
    'component_count': int(original.shape[1]),
    'original_array_shape': list(original.shape),
    'original_bandwidth_hz': float(basis.original_bw),
    'original_dwell_time_s': float(basis.original_dwell),
    'target_bandwidth_hz': float(sys.argv[2]),
    'target_points': int(sys.argv[3]),
    'formatted_array_shape': list(formatted.shape),
    'formatted_finite': bool(np.isfinite(formatted).all()),
}))
"""
    completed = subprocess.run(
        [str(interpreter), "-c", program, str(path), str(float(bandwidth_hz)), str(int(points))],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise ValueError(f"FSL-MRS could not read/format basis {path}: {completed.stderr.strip() or completed.stdout.strip()}")
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValueError(f"FSL-MRS basis inspection returned invalid metadata for {path}.") from exc


def make_voxel_diagnostic_pdf(
    input_path: Path,
    fit_dir: Path,
    output_path: Path,
    mask_path: Path,
    metadata: dict[str, Any],
    fsl_executable: str,
) -> None:
    """Create paginated real-spectrum, fit, baseline, and residual plots for every fitted voxel."""
    executable_dir = Path(fsl_executable).resolve().parent
    python_candidates = [executable_dir / "python3", executable_dir / "python", executable_dir / "python3.12"]
    interpreter = next((candidate for candidate in python_candidates if candidate.is_file() and os.access(candidate, os.X_OK)), None)
    if interpreter is None:
        raise ValueError(f"Could not locate a Python interpreter alongside FSL-MRS executable {fsl_executable}.")
    program = r"""
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import nibabel as nib
import numpy as np
from fsl_mrs.core import MRS
from fsl_mrs.utils.misc import FIDToSpec

input_path, fit_dir, output_path, mask_path, metadata_path = map(Path, sys.argv[1:])
metadata = json.loads(metadata_path.read_text())
source = nib.load(str(input_path))
fit = nib.load(str(fit_dir / 'fit' / 'fit.nii.gz'))
baseline = nib.load(str(fit_dir / 'fit' / 'baseline.nii.gz'))
residual = nib.load(str(fit_dir / 'fit' / 'residual.nii.gz'))
mask = np.asanyarray(nib.load(str(mask_path)).dataobj) > 0
data_values = np.asanyarray(source.dataobj)
fit_values = np.asanyarray(fit.dataobj)
baseline_values = np.asanyarray(baseline.dataobj)
residual_values = np.asanyarray(residual.dataobj)
frequency = metadata.get('SpectrometerFrequency', metadata.get('SpectrometerFrequencyMHz'))
if isinstance(frequency, (list, tuple)):
    frequency = frequency[0]
nucleus = metadata.get('ResonantNucleus', metadata.get('Nucleus', '1H'))
if isinstance(nucleus, (list, tuple)):
    nucleus = nucleus[0]
bandwidth = 1.0 / float(source.header.get_zooms()[3])
voxels = list(zip(*np.where(mask)))
output_path.parent.mkdir(parents=True, exist_ok=True)
with PdfPages(output_path) as pdf:
    for first in range(0, len(voxels), 12):
        page_voxels = voxels[first:first + 12]
        figure, axes = plt.subplots(4, 3, figsize=(14, 12), squeeze=False)
        for axis, index in zip(axes.flat, page_voxels):
            mrs = MRS(FID=np.asanyarray(data_values[index]).conj(), cf=float(frequency), bw=bandwidth, nucleus=str(nucleus))
            ppm = mrs.getAxes()
            region = (ppm >= 0.0) & (ppm <= 6.0)
            signals = [
                (np.asanyarray(data_values[index]).conj(), 'data', 'black'),
                (np.asanyarray(fit_values[index]).conj(), 'fit', 'tab:red'),
                (np.asanyarray(baseline_values[index]).conj(), 'baseline', 'tab:blue'),
                (np.asanyarray(residual_values[index]).conj(), 'residual', 'tab:green'),
            ]
            for signal, label, color in signals:
                spectrum = FIDToSpec(signal)
                axis.plot(ppm[region], np.real(spectrum)[region], color=color, linewidth=0.8, label=label)
            axis.set_xlim(4.0, 1.8)
            axis.set_title(f'Voxel {index}', fontsize=8)
            axis.set_xlabel('ppm')
            axis.tick_params(labelsize=7)
        for axis in axes.flat[len(page_voxels):]:
            axis.set_visible(False)
        axes.flat[0].legend(fontsize=7, loc='best')
        figure.suptitle(f'Voxel-wise FSL-MRS diagnostics: voxels {first + 1}-{first + len(page_voxels)} / {len(voxels)}')
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)
"""
    metadata_path = fit_dir / "voxel_diagnostics_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    completed = subprocess.run(
        [str(interpreter), "-c", program, str(input_path), str(fit_dir), str(output_path), str(mask_path), str(metadata_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise ValueError(f"Could not generate voxel diagnostic PDF: {completed.stderr.strip() or completed.stdout.strip()}")