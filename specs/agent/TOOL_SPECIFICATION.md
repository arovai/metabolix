# Metabolix Tool Specification

This specification was synthesized from the repository's unstructured notes using the structure in the ln2t BIDS application template. The cloned template's actual instruction paths are under `templates/agentic_coding_instructions/` rather than the shortened path in the prompt.

## 1. Purpose and Overview

Metabolix is a CLI-first BIDS application for processing already-converted, image-space multi-voxel proton MRS (MRSI) data using FSL-MRS. It discovers BIDS `mrs/` NIfTI-MRS acquisitions, validates spatial/spectral metadata and unsupported dimensions, combines coils when needed, optionally creates residual-water removal branches, fits all eligible spatial voxels, exports application-owned QC/provenance, and can estimate nominal anatomical ROI overlap using image-world affines.

It supports `participant` processing and `group` aggregation of participant outputs. It does not reconstruct vendor raw data, process k-space, register anatomy, provide absolute quantification, or establish validated scientific measurements. Input files are read-only. Normal fitting requires a compatible explicit basis or unambiguous configured mapping.

Supported layout is BIDS MRS `sub-*/[ses-*/]mrs/*_mrs.nii[.gz]` with sidecar JSON, plus a documented compatibility discovery for the supplied legacy-style `sub-*/mrs/sub-*_mrsi.nii.gz`. `OUTPUT_DIR` is the derivative dataset root (for example `/data/derivatives/metabolix`) and contains subject/session/acquisition/run-specific directories and BIDS-like entities.

## 2. Core Libraries and Dependencies

- Python 3.10+.
- Runtime package: NumPy for validation and numerical geometry, nibabel for NIfTI/NIfTI-MRS and affine operations, PyYAML for safe YAML configuration.
- FSL-MRS supplies `fsl_mrs_proc` and `fsl_mrsi` executables and its scientific fitting implementation. Metabolix invokes these through argument arrays and does not substitute `fslpy` for FSL binaries.
- Optional HTML/report plotting may use Matplotlib when installed; application TSV, JSON, Markdown and NIfTI outputs must not depend on FSL-generated HTML reports.
- The wrapper itself is pip-installable. Full FSL-MRS runtime setup is documented through the supported FSL installer or FSL's Conda channel; installation must never download FSL implicitly.

## 3. Input Data Specifications

### BIDS entities and discovery

`sub` is required. `ses`, `acq`, `run`, `task`, and `rec` are optional and preserved in output identity. Participant/session/acquisition/run filters select matching files. Do not silently collapse two acquisitions to the same output path.

Primary discovery: image-space NIfTI-MRS in a BIDS `mrs` directory with matching JSON sidecar. Compatibility discovery: the example dataset's `sub-*_mrsi.nii[.gz]` under `sub-*/mrs/`, even when its suffix is not standard BIDS MRS. The compatibility path is reported in provenance and does not imply full BIDS validation.

Anatomy/segmentation is optional. ROI sources may be FreeSurfer `aseg.mgz`, a 3D integer label image, or a binary mask, with explicit labels for integer maps. FreeSurfer's scanner/world affine is used, not tkregister/surface coordinates. Matching T1/anatomy can be checked where available; alignment still requires visual/external validation.

### Validation

- Load a 3D spatial grid plus a spectral dimension and optional tagged higher dimensions using NIfTI-MRS metadata; never assume dimension 5 is a coil dimension.
- Require complex samples, finite nonzero signal checks, a finite nonsingular affine, positive spectral dwell/sampling, and readable metadata (including nucleus and spectrometer frequency where required by FSL-MRS). Record optional TE/TR and warn when absent.
- Reject k-space representations and unsupported dynamic/editing dimensions with actionable messages. Accept known coil dimension tags for coil combination and accept already coil-combined data.
- Verify qform/sform agreement and geometry compatibility for masks/segmentations; flag disagreement rather than silently choosing anatomical interpretation.
- Validate basis readability and components, frequency, dwell/bandwidth, point count and available sequence metadata. Formatting/resampling does not establish sequence compatibility.
- Build the fitting mask from every spatial voxel with finite, nonzero FID unless an explicit compatible mask is provided. Invalid voxels stay represented in QC as invalid, never as measured zeros.

## 4. Command-Line Options and Configuration

Required BIDS-App form:

```text
metabolix BIDS_DIR OUTPUT_DIR {participant|group} [options]
```

Common options: `--participant-label LABEL [LABEL ...]`, `--session LABEL`, `--acquisition LABEL`, `--run LABEL`, `--config YAML`, `--basis PATH`, optional validated `--coil-reference PATH`, `--freesurfer-dir PATH`, `--roi-labels INT [INT ...]`, `--water-removal {off,on,compare}`, `--n-jobs N`, `--dry-run`, `--resume`, `--overwrite`, `--verbose`, `--check-deps`, `--write-default-config PATH`, `--write-example-config PATH`, and `--export-config PATH`.

Fitting options include `--ppmlim LOW HIGH`, `--ignore NAME ...`, repeatable `--combine NAME NAME`, `--internal-ref NAME ...`, `--baseline VALUE`, `--algorithm NAME`, optional `--fit-mask PATH`, and `--roi-step-mm VALUE`. CLI values override YAML values only when explicitly supplied. Precedence is built-in neutral defaults < YAML < explicit CLI. CLI list values replace the YAML list. Boolean switches support explicit positive/negative forms where applicable. Unknown YAML keys and invalid values are errors. Relative config paths are resolved relative to the YAML file.

Configuration is grouped into `processing`, `basis`, `fit`, `roi`, `execution`, `qc`, and `report`. See `config/default.yaml` and `config/ge_press_te144.yaml`. The latter is a clearly acquisition-specific example, not a universal model. Optional QC limits only add flags; they never reject estimates.

Neutral defaults: all valid spatial voxels, automatic coil handling from dimension tags, water removal off, no anatomical restriction, one worker, QC/report generation on, no absolute quantification, no arbitrary rejection threshold. Missing/ambiguous basis mapping is a clear preflight error before fitting.

## 5. Processing Workflow

1. Parse command/config, validate keys and values, resolve relative paths, and persist full resolved configuration.
2. Discover and select MRS inputs and record normalized BIDS entities. Validate image-space representation, geometry, signal, dimensions and metadata without modifying source data.
3. For tagged coil data, call `fsl_mrs_proc coilcombine` via `subprocess.run(argv, ...)`; use a reference only when explicitly provided and validated. Verify output spectral/spatial shape and metadata. Record the covariance sample warning if present; it is cautionary, not an automatic failure. Do not assume an HTML report exists.
4. `off`: fit coil-combined/original input. `on`: preserve input and fit water-removed branch. `compare`: generate and fit both branches with identical basis, mask and model. Water ppm limits are configured; no automatic water-peak threshold selects a branch.
5. Validate and record the selected basis: checksum, components, metadata, path, and compatibility caveats. Run `fsl_mrsi` over all valid spatial voxels, passing configured model settings and bounded worker count. Combined components report sums, not fixed component ratios. Reporting selection does not prune fitted components.
6. Preserve native FSL outputs. Export application-owned voxel status/coordinates, summaries where extraction is reliable, warnings, exact commands, captured output, timestamps, versions and basis details. Treat FSL `*_sd` maps as CRLB-derived percentage uncertainty estimates (possibly capped), not absolute SD or fully propagated ratio uncertainty. Undefined ratios from zero references are missing with flags; native zero maps remain untouched. Outside-mask zeros are not observations.
7. If requested, calculate nominal ROI overlap by midpoint sampling within each MRS voxel, mapping MRS world coordinates through inverse segmentation affine. Use nearest-neighbour label lookup. Default maximum sampling step is 0.5 mm; repeat at 1 mm and report sensitivity. The entire MRS voxel is the denominator; segmentation field-of-view coverage is a separate metric. Save per-voxel overlap, maps, ROI volume/coverage, and optional thresholded summaries. Do not use ROI to restrict fitting.
8. Produce accessible TSV/JSON/Markdown and an HTML report when report dependencies are available; report all voxels using lazy/paginated diagnostics. Compare water branches with paired differences and invalid-fit counts, never auto-select a winner from fit error alone.
9. Group level aggregates existing participant QC tables only; it does not refit participant data. Preserve source references and subject/run/acquisition identity.

FSL command arguments must be safe arrays, never shell interpolation. Support paths containing spaces. Resume is permitted only when source checksums, resolved relevant configuration and tool version match a recorded manifest. Existing outputs are not overwritten by default; partial results remain discoverable after failures.

## 6. Output Specifications

Suggested layout:

```text
OUTPUT_DIR/
  dataset_description.json
  sub-XX/[ses-YY/]<acquisition-specific-run>/
    config_resolved.yaml
    manifest.json
    processing.log
    01_coilcombine/...
    02_waterremove/...
    03_fit_original/...
    03_fit_waterremoved/...
    03_fit_*/voxel_qc.tsv, voxel_qc_definitions.json
    qc/water_removal_paired_differences.tsv (compare mode)
    roi/voxel_overlap.tsv, overlap maps, roi_summary.tsv
    report.md, report.html, voxel_diagnostics.pdf
  group/group_summary.tsv
  logs/
```

Include derivative dataset metadata (`GeneratedBy`, source dataset references), source checksums, all selected entities, package/Python/FSL-MRS versions, basis checksum/metadata, resolved configuration, and per-command argv, timestamp, exit status, stdout/stderr locations. Never alter native FSL products to sanitize invalid outputs.

QC rows include entities, zero-based `i,j,k`, branch, validity/status/reason, available components and combined results, raw amplitudes versus internally referenced values, CRLB percentage estimates, component SNR/FWHM, residual summaries, fit parameters and ROI overlap/QC flags. Missing values stay missing. If an installed FSL output format does not provide a field safely, report it unavailable rather than infer it.

## 7. Errors and Edge Cases

- No matching input, ambiguous basis, invalid/unsupported dimensions, k-space input, malformed metadata, incompatible mask/affine, or missing required executable: fail early with a corrective message.
- Missing optional TE/TR, qform/sform disagreement, less than recommended noise samples, partial segmentation coverage, or absent FSL report: retain outputs and report a warning where safe.
- Invalid/zero-signal voxel: preserve row with status/reason and exclude it from fitting mask.
- A zero internal reference: export ratio as missing plus an explicit flag; do not rewrite native FSL map.
- Per-acquisition output collision: refuse to merge/overwrite.
- Partial subprocess failure: preserve completed stage directories and exact captured logs; return nonzero and identify affected input.
- Group mode with no participant tables: clear error. No subject-level refitting at group level.

## 8. Validation Notes

No test suite is to be created. Validate using package compilation/syntax checks, installed CLI `--help`/`--version`, config export and invalid-config behavior, dependency diagnostics, dry-run discovery/config/provenance, and geometry/ROI sanity only on available small data if non-expensive. Do not run full-data fitting without explicit authorization. The provided reference scripts/log describe prior user/Codex exploration; their reported processing results are not newly reproduced by this implementation.

## 9. References

- FSL-MRS installation: https://pages.fmrib.ox.ac.uk/fsl/fsl_mrs/install.html
- FSL-MRS processing: https://pages.fmrib.ox.ac.uk/fsl/fsl_mrs/processing.html
- FSL-MRS fitting: https://pages.fmrib.ox.ac.uk/fsl/fsl_mrs/fitting.html
- FSL-MRS basis simulation: https://pages.fmrib.ox.ac.uk/fsl/fsl_mrs/simulation.html
- BIDS MRS specification: https://bids-specification.readthedocs.io/
- MRSCloud: https://github.com/shui5/MRSCloud
- Consensus basis-set considerations: https://pmc.ncbi.nlm.nih.gov/articles/PMC7442593/
- Local exploratory references: dataset `code/PROCESSING_LOG.md`, `code/report_thalamus_overlap.py`, and `code/check_mrsi_pilot.py`.

## 10. Implementation Notes

- Package name and executable: `metabolix`; use `pyproject.toml` entry point.
- Keep CLI/config/discovery/validation/subprocess/QC/ROI/report code modular and dependency-light.
- Store no pilot-four-voxel behavior; all-valid-spatial-voxel mask is default.
- Distinguish raw data, coil-combined and water branches; retain all branch inputs and native outputs.
- Do not infer concentration units, SNR from peak-height ratios, spectral quality from successful execution, or registration validity from matching headers.
- Document anatomy overlap as nominal geometric overlap and note finite MRSI point-spread response.
- Include practical tutorial using the provided `/home/arovai/datasets_shortcuts/2021-Thankful_Alligator-640c8fba6b3e` dataset paths. A sample real-data run is prepared only; no costly fitting during development.