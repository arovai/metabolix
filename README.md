# Metabolix

Metabolix is a pip-installable, CLI-first BIDS application for all-spatial-voxel proton MRSI processing with FSL-MRS. It validates already-converted image-space NIfTI-MRS, combines tagged receiver coils when needed, optionally compares residual-water removal, fits every valid spatial voxel, and can report nominal anatomical ROI overlap. Input data are never modified.

Successful execution does not establish validated metabolite measurements. Basis compatibility, spectral quality, fit uncertainty, nominal anatomical overlap, and quantification units must be evaluated separately.

## Installation

The Metabolix wrapper is pip-installable:

```bash
python -m pip install .
metabolix --version
```

FSL-MRS is a separate scientific runtime and is not installed by pip or by Metabolix. The official installation guidance recommends the FSL installer or its Conda package. For a minimal Conda environment:

```bash
conda create -n metabolix-fsl python=3.11
conda activate metabolix-fsl
conda install -c conda-forge \
	-c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/ \
	fsl_mrs
python -m pip install metabolix-bids
fsl_mrs --version
fsl_mrs_verify
metabolix --check-deps
```

Micromamba can replace Conda. Do not install FSL-MRS into the base environment. FSL-MRS may also be installed through the main FSL installer. Metabolix searches executables on `PATH`; `basis.binary_dir` / `--binary-dir` can point to another executable directory. `fslpy` alone is not a substitute for FSL-MRS binaries. See [official FSL-MRS installation instructions](https://pages.fmrib.ox.ac.uk/fsl/fsl_mrs/install.html).

Runtime log levels are ANSI-colored when writing to an interactive terminal. Set `NO_COLOR=1` to disable colors. Missing or ambiguous basis errors include suggested YAML/CLI remedies and point to the relevant sections below.

## Quick Start

```bash
metabolix --write-default-config mrsi.yaml
metabolix --check-deps
metabolix /data/bids /data/derivatives/metabolix participant \
	--participant-label 11 \
	--config mrsi.yaml \
	--basis /data/basis/compatible_basis.BASIS \
	--dry-run
```

Actual processing:

```bash
metabolix /data/bids /data/derivatives/metabolix participant \
	--participant-label 11 \
	--config mrsi.yaml \
	--basis /data/basis/compatible_basis.BASIS \
	--n-jobs 4
```

Water-removal sensitivity comparison and optional FreeSurfer overlap:

```bash
metabolix /data/bids /data/derivatives/metabolix participant \
	--config mrsi.yaml \
	--water-removal compare \
	--freesurfer-dir /data/derivatives/freesurfer_7.3.2 \
	--roi-labels 10 49 \
	--n-jobs 4
```

Group mode aggregates existing participant QC tables and does not refit:

```bash
metabolix /data/bids /data/derivatives/metabolix group
```

Generate neutral and illustrative GE PRESS TE144 configurations with `--write-default-config PATH` and `--write-example-config PATH`. The acquisition profile in `config/ge_press_te144.yaml` contains example-specific model choices and placeholder paths, not universal defaults.

## Configuration

Precedence is built-in defaults < YAML configuration < explicitly provided CLI options. Omitted CLI options do not overwrite YAML. A CLI list replaces the configured list. Relative paths in YAML are resolved against the YAML file. Unknown keys and invalid values fail before processing. `--export-config PATH` writes the resolved configuration.

Configuration sections:

- `processing`: `coil_combine` (`auto|on|off`), optional validated `coil_reference`, `water_removal` (`off|on|compare`), and `water_ppm`.
- `basis`: explicit `path`, optional `mapping` entries with `match` metadata/entity values and `path`, and optional `binary_dir`.
- `fit`: `ppmlim`, `ignore`, repeated component `combine` groups, `internal_reference`, FSL-MRS `baseline`, `algorithm` (`Newton|MH`), and explicit `mask`.
- `roi`: `freesurfer_dir` or label `map`, integer `labels`, `sampling_step_mm` (maximum 0.5), and optional reporting-only `minimum_overlap_percent`.
- `execution`: bounded `n_jobs` (1-256), `resume`, and `overwrite`.
- `qc`: optional `max_crlb_percent`, `min_snr`, and `max_linewidth_hz` operational flags; unset by default and never used to reject estimates.
- `report`: `enabled` and reporting-focused metabolite names. This does not prune model components.
- `selection`: participant, session, acquisition and run lists.

Neutral defaults process all valid spatial voxels; detect coil dimensions from tags; leave water removal off; impose no anatomical restriction; use one worker; enable QC/reporting; and do not perform absolute quantification or uncertainty-threshold rejection. Fitting requires a compatible basis. The example GE PRESS TE144 profile is not a universal scientific protocol.

### Obtaining a basis

If you do not already have a suitable basis, consider submitting a basis-generation job through [MRSCloud](https://mricloud.org/) (project information and source: [MRSCloud](https://github.com/shui5/MRSCloud)). Enter the acquisition's actual scanner/vendor, sequence and localization, field strength, echo time, and metabolite selection. Where offered, request an LCModel `.BASIS` export, which FSL-MRS can read. Then point `basis.path` or `--basis` to the downloaded file. MRSCloud availability, supported sequences, and job options depend on the service; its generated basis is a candidate, not proof of equivalence to the scanner sequence. Check component names, sequence assumptions, field strength, TE, frequency, and sampling before fitting. Readable/resampled output alone does not validate the model.

The report metabolite list controls which metabolite maps are summarized in report tables only. Full per-voxel maps and QC values remain available.

## Inputs and Outputs

Primary input is already-converted image-space complex NIfTI-MRS under BIDS `sub-*/[ses-*/]mrs/*_mrs.nii[.gz]` with metadata. The documented compatibility pattern `sub-*/mrs/sub-*_mrsi.nii[.gz]` supports the supplied dataset, but does not claim full BIDS conformance. K-space and unsupported higher dimensions are rejected. Dynamic/editing dimensions are not silently averaged.

The positional `OUTPUT_DIR` is the derivative dataset root and contains `dataset_description.json`, subject/session/source-specific folders, processing stages, resolved configuration, manifest, logs, branch-specific voxel QC tables, optional ROI overlap maps/tables, `report.md`/`report.html`, and paginated `voxel_diagnostics.pdf` files for all fitted voxels. Native FSL-MRS outputs are retained. Group summary is `OUTPUT_DIR/group/group_summary.tsv`.

ROI overlap uses scanner/world affines, fine sampling at <=0.5 mm with a 1 mm sensitivity comparison, whole-MRS-voxel denominator, and separate segmentation field-of-view coverage. It never restricts the default fit mask. Existing header alignment is assumed, not validated registration.

## Scientific Notes and Tutorial

See [docs/TUTORIAL.md](docs/TUTORIAL.md) for the actual dataset paths and walkthrough, and [docs/LIMITATIONS.md](docs/LIMITATIONS.md) for interpretation limits and troubleshooting. The notes' exploratory pilot and processing log live outside this package and were not rerun during implementation. No expensive real-data fitting was run here.

Structured audit artifacts are in `specs/agent/`: the tool specification, finalized implementation instructions, and source-notes summary.

## Development Validation

The project intentionally includes no test suite, pytest/coverage dependencies, fixtures, or testing CI, as requested. Syntax, CLI help, configuration handling, dependency diagnostics and dry-run checks are the lightweight validation surface.
