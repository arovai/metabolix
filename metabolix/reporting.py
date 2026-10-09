"""Figure generation for application-owned MRSI reports."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

from metabolix.discovery import MRSInput
from metabolix.output_layout import bids_filename


def _fsl_python(fsl_executable: str) -> Path:
    directory = Path(fsl_executable).resolve().parent
    candidates = (directory / "python3", directory / "python", directory / "python3.12")
    interpreter = next((path for path in candidates if path.is_file() and os.access(path, os.X_OK)), None)
    if interpreter is None:
        raise RuntimeError(f"Could not find the FSL-MRS Python interpreter beside {fsl_executable}.")
    return interpreter


def generate_report_figures(
    item: MRSInput,
    figure_dir: Path,
    work_dir: Path,
    fsl_executable: str,
    metadata: dict[str, Any],
    spectral_inputs: dict[str, str],
    basis_path: str,
    map_inputs: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Create QC, basis, and metabolite map figures in an acquisition figure folder."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    payload_path = work_dir / bids_filename(item, "metabolix-reportInputs", "json", ".json")
    payload_path.write_text(
        json.dumps(
            {
                "spectral_inputs": spectral_inputs,
                "metadata": metadata,
                "basis_path": basis_path,
                "map_inputs": map_inputs,
                "figure_dir": str(figure_dir),
                "entity_stem": item_stem(item),
            },
            indent=2,
        )
        + "\n"
    )
    program = r'''
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from fsl_mrs.core import MRS
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils.misc import FIDToSpec

payload = json.loads(Path(sys.argv[1]).read_text())
metadata = payload['metadata']
frequency = metadata.get('SpectrometerFrequency', metadata.get('SpectrometerFrequencyMHz'))
if isinstance(frequency, (list, tuple)):
    frequency = frequency[0]
nucleus = metadata.get('ResonantNucleus', metadata.get('Nucleus', '1H'))
if isinstance(nucleus, (list, tuple)):
    nucleus = nucleus[0]
figure_dir = Path(payload['figure_dir'])
stem = payload['entity_stem']
spectral_width = metadata.get('SpectralWidth')

def load_fids(path):
    image = nib.load(path)
    values = np.asanyarray(image.dataobj)
    points = values.shape[3]
    bandwidth = float(spectral_width or 1.0 / image.header.get_zooms()[3])
    if values.ndim == 5:
        fids = values.reshape(-1, points, values.shape[4]).transpose(0, 2, 1)
    elif values.ndim == 4:
        fids = values.reshape(-1, points)[:, None, :]
    else:
        raise ValueError(f'Unexpected MRS shape for report figure: {values.shape}')
    valid = np.isfinite(fids).all(axis=(1, 2)) & np.any(fids != 0, axis=(1, 2))
    return fids[valid].conj(), bandwidth

def summarize_spectra(path):
    fids, bandwidth = load_fids(path)
    if not len(fids):
        raise ValueError(f'No finite nonzero FIDs available for {path}')
    # Vectorized FSL-MRS Fourier convention: FIDToSpec applies the special first-point rule.
    spectra = FIDToSpec(fids.copy(), axis=-1)
    magnitude = np.abs(spectra).reshape(-1, spectra.shape[-1])
    mean = np.mean(magnitude, axis=0)
    p10, p90 = np.percentile(magnitude, [10, 90], axis=0)
    times_ms = np.arange(fids.shape[-1]) * (1000.0 / bandwidth)
    mean_fid_magnitude = np.mean(np.abs(fids), axis=(0, 1))
    reference = MRS(FID=fids[0, 0], cf=float(frequency), bw=bandwidth, nucleus=str(nucleus))
    ppm = reference.getAxes()
    return ppm, mean, p10, p90, times_ms, mean_fid_magnitude, len(magnitude)

series = []
for label, path in payload['spectral_inputs'].items():
    if path and Path(path).is_file():
        values = summarize_spectra(path)
        series.append((label, values))

if series:
    figure, axes = plt.subplots(1, 2, figsize=(13, 5))
    for label, (ppm, mean, p10, p90, times_ms, fid_magnitude, count) in series:
        region = (ppm >= 0.0) & (ppm <= 6.0)
        axes[0].plot(ppm[region], mean[region], linewidth=1.2, label=f'{label} (n={count})')
        axes[0].fill_between(ppm[region], p10[region], p90[region], alpha=0.10)
        axes[1].plot(times_ms, fid_magnitude, linewidth=1.2, label=label)
    axes[0].set_xlim(6.0, 0.0)
    axes[0].set_xlabel('Chemical shift (ppm)')
    axes[0].set_ylabel('Mean voxel-wise spectral magnitude (a.u.)')
    axes[0].set_title('Spatial distribution summary: mean and 10th-90th percentile')
    axes[1].set_xlabel('Time (ms)')
    axes[1].set_ylabel('Mean voxel/coil FID magnitude (a.u.)')
    axes[1].set_title('FID magnitude QC (not a coherent spatial mean)')
    for axis in axes:
        axis.legend(fontsize=8)
        axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(figure_dir / f'{stem}_desc-metabolix_qc-spatialSpectraAndFID.png', dpi=160)
    plt.close(figure)

basis = mrs_io.read_basis(payload['basis_path'])
basis_array = np.asarray(basis.original_basis_array).conj()
basis_spectra = np.abs(FIDToSpec(basis_array.copy(), axis=0))
ppm = basis.original_ppm_shift_axis
names = list(basis.names)
columns = 4
rows = int(np.ceil(len(names) / columns))
figure, axes = plt.subplots(rows, columns, figsize=(14, max(5, rows * 2.7)), squeeze=False)
for index, axis in enumerate(axes.flat):
    if index >= len(names):
        axis.set_visible(False)
        continue
    region = (ppm >= 0.0) & (ppm <= 6.0)
    component = basis_spectra[:, index]
    scale = np.max(component[region]) if np.any(region) else np.max(component)
    axis.plot(ppm[region], component[region] / scale if scale else component[region], color='#28666e', linewidth=0.9)
    axis.set_xlim(6.0, 0.0)
    axis.set_title(names[index], fontsize=9)
    axis.set_xlabel('ppm', fontsize=8)
    axis.tick_params(labelsize=7)
figure.suptitle(f'Basis components used by FSL-MRS: {Path(payload["basis_path"]).name}\nNormalized magnitude spectra; components are displayed, not validation evidence')
figure.tight_layout()
figure.savefig(figure_dir / f'{stem}_desc-metabolix_basisComponents.png', dpi=160)
plt.close(figure)

for entry in payload['map_inputs']:
    image = nib.load(entry['path'])
    values = np.squeeze(np.asanyarray(image.dataobj))
    if values.ndim == 2:
        values = values[:, :, None]
    if values.ndim != 3:
        continue
    # MRSI data are often single-slice; render every slice when there are multiple.
    slice_count = values.shape[2]
    columns = min(4, slice_count)
    rows = int(np.ceil(slice_count / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(3.4 * columns, 3.2 * rows), squeeze=False)
    finite = values[np.isfinite(values)]
    if not finite.size:
        plt.close(figure)
        continue
    vmin, vmax = np.percentile(finite, [2, 98])
    if vmin == vmax:
        vmin, vmax = float(finite.min()), float(finite.max() + 1e-6)
    image_artist = None
    for z, axis in enumerate(axes.flat):
        if z >= slice_count:
            axis.set_visible(False)
            continue
        image_artist = axis.imshow(values[:, :, z].T, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax, interpolation='nearest')
        axis.set_title(f'z={z}', fontsize=8)
        axis.set_xlabel('i')
        axis.set_ylabel('j')
    figure.colorbar(image_artist, ax=axes.ravel().tolist(), shrink=0.75, label=entry.get('units', 'map value'))
    figure.suptitle(entry['title'])
    figure.tight_layout()
    figure.savefig(entry['figure_path'], dpi=160)
    plt.close(figure)
'''
    completed = subprocess.run(
        [str(_fsl_python(fsl_executable)), "-c", program, str(payload_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"FSL-MRS report figure generation failed: {completed.stderr.strip() or completed.stdout.strip()}")

    stem = item_stem(item)
    figures = []
    preprocessing_figure = figure_dir / f"{stem}_desc-metabolix_qc-spatialSpectraAndFID.png"
    if preprocessing_figure.is_file():
        figures.append({
            "path": str(preprocessing_figure),
            "title": "Pre/post-processing mean voxel-wise spectra and FID magnitude",
            "kind": "preprocessing-qc",
        })
    basis_figure = figure_dir / f"{stem}_desc-metabolix_basisComponents.png"
    if basis_figure.is_file():
        figures.append({"path": str(basis_figure), "title": "Basis component spectra", "kind": "basis"})
    for entry in map_inputs:
        figure_path = Path(entry["figure_path"])
        if figure_path.is_file():
            figures.append({
                "path": str(figure_path),
                "map_path": entry["path"],
                "title": entry["title"],
                "kind": "metabolite-heatmap",
            })
    return figures


def item_stem(item: MRSInput) -> str:
    from metabolix.output_layout import entity_stem

    return entity_stem(item)