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

The positional `OUTPUT_DIR` is the derivative dataset root. Participant products follow `sub-*/[ses-*/]mrs/`; each acquisition has an entity-named `*_desc-metabolix/` bundle. The report is named `sub-..._desc-metabolix_report.html` (with a matching Markdown report). Downloadable metabolite/QC maps are BIDS-named `*_statmap.nii.gz` files with JSON sidecars; branch voxel tables use `*_voxelqc.tsv`. Detailed material is grouped into `figures/`, `maps/`, `qc/`, `roi/`, `logs/`, `provenance/`, `work/`, and `processing/fsl-mrs/`. Native FSL-MRS outputs are retained. Group summary is `OUTPUT_DIR/group/group_summary.tsv`.

Typical layout:

```text
OUTPUT_DIR/
  dataset_description.json
  sub-11/[ses-01/]mrs/
    sub-11[_ses-01][_task-...][_acq-...][_run-...]_desc-metabolix/
      sub-11[_ses-01][_task-...][_acq-...][_run-...]_desc-metabolix_report.html
      figures/       preprocessing QC, basis and per-metabolite heatmaps
	maps/          BIDS-named raw/internal/QC NIfTI statmaps and JSON sidecars
      qc/            branch voxel tables and paired comparison tables
	processing/fsl-mrs/  native FSL-MRS products by branch/stage
      roi/ logs/ provenance/ work/
```

The optional bracketed entities appear only when present in the input; MRS-specific `voi` and `nuc` entities are retained too. The report's first sections provide links to QC tables, downloadable raw/internal/QC maps, figures, native outputs, logs and provenance. Heatmaps emphasize the internally referenced maps when available; the detailed map directory retains all published scales and component QC maps.

See [ROI estimates and metabolite units](docs/ROI_ANALYSIS.md) for the averaging formula, overlap thresholds, small-ROI limitations, output table columns, and the meaning of `internal` maps.

ROI overlap uses scanner/world affines, fine sampling at <=0.5 mm with a 1 mm sensitivity comparison, whole-MRS-voxel denominator, and separate segmentation field-of-view coverage. It never restricts the default fit mask. Existing header alignment is assumed, not validated registration.

## Scientific Notes and Tutorial

See [docs/TUTORIAL.md](docs/TUTORIAL.md) for the actual dataset paths and walkthrough, and [docs/LIMITATIONS.md](docs/LIMITATIONS.md) for interpretation limits and troubleshooting. The notes' exploratory pilot and processing log live outside this package and were not rerun during implementation. No expensive real-data fitting was run here.

Structured audit artifacts are in `specs/agent/`: the tool specification, finalized implementation instructions, and source-notes summary.

## Development Validation

The project intentionally includes no test suite, pytest/coverage dependencies, fixtures, or testing CI, as requested. Syntax, CLI help, configuration handling, dependency diagnostics and dry-run checks are the lightweight validation surface.


### Resuming a subset of participants

You may change `--participant-label` when resuming with `--resume`. Participant
selection is excluded from each input's fingerprint; input contents, basis,
processing settings, and software versions remain checked. Existing manifests
created with participant selection in their fingerprint are supported: the
application verifies them against current files and settings while allowing only
the participant list to differ. Completed processing stages are reused. Do not
use `--overwrite` to resume an interrupted run.


Group output is a participant-wide table: one row per analysed participant, separate
columns for acquisition/branch/ROI/measure combinations, and an automatic join
with input `participants.tsv` metadata. See [group table details](docs/TUTORIAL.md#participant-wide-group-table).


## Output labels and units

These definitions describe estimates and QC, not acceptance criteria. Blank/n/a estimates are missing; map zeros may also encode unfitted or invalid results.

| Label | Units | Meaning |
|---|---|---|
| `crlb_percent` — CRLB / relative fit uncertainty | percent | CRLB means Cramér–Rao lower bound: a model-based lower bound on estimator variance under statistical assumptions. For Newton fitting, FSL-MRS reports 100 × sqrt(CRLB variance) / fitted estimate, i.e. relative SD, not variance. 20 means estimated SD about 20% of the fitted value; it is not a concentration, percent of reference, or 95% confidence interval. For MH (Metropolis–Hastings) fitting, the same legacy crlb_percent label holds posterior-sample SD divided by the posterior mean × 100, not a CRLB. The inspected FSL-MRS implementation caps percentages above 999 and maps undefined NaN percentages to 999. 999 is a capped/undefined result, not a precise uncertainty. Near-zero estimates can have very large relative uncertainty. Low uncertainty does not exclude artifacts, incorrect basis/sequence assumptions, bias or partial-volume mixing. There is no universal acceptance threshold; uncertainty of an internal ratio or ROI mean is not automatically propagated. |
| `raw_fit` — Raw fitted amplitude | arbitrary units | FSL-MRS fitted amplitude before internal-reference scaling. Not a molar concentration; scaling can differ between spectra. |
| `internal_reference_ratio` — Internal-reference ratio | dimensionless | Metabolite estimate divided by the sum of configured reference metabolites within the same voxel (default Cr+PCr). 1.5 means 150% of reference, not 150% of the total spectrum. A zero/nonfinite reference is invalid. Native/published FSL maps may encode invalid or unfitted values as zero; inspect voxel-table validity. |
| `snr` — Signal-to-noise ratio (SNR) | dimensionless | FSL-MRS component-specific spectral QC SNR, using its noise estimation and matched-filter conventions. Not the residual-water/metabolite peak ratio; values depend on implementation and component and are not a universal quality score. |
| `linewidth_fwhm_hz` — Full width at half maximum (FWHM) | Hz | FSL-MRS component spectral width at half maximum. Broader peaks can impair separation, but linewidth alone does not establish fit quality. Distinct from added nuisance broadening. |
| `residual_l2_over_input_l2` — Relative residual norm | dimensionless | Euclidean norm of the complex time-domain residual divided by the complex input FID norm. Descriptive fit discrepancy; not a percentage of metabolite loss, a concentration uncertainty or a standalone quality criterion. |
| `overlap_pct` — Nominal ROI overlap | percent | 100 × ROI volume inside an MRS voxel / entire nominal MRS voxel volume. Fine anatomical sampling does not upsample MRS signal or model its spatial response. Existing image alignment is assumed. |
| `overlap_mm3` — Nominal ROI volume within voxel | mm^3 | Estimated intersection volume of the ROI with the nominal MRS voxel. |
| `sampling_difference_pp` — Overlap sampling sensitivity | percentage points | Absolute difference between fine and 1 mm overlap estimates. A numerical sensitivity check, not a rigorous error bound or a relative percent change. |
| `segmentation_fov_pct` — Segmentation field-of-view (FOV) coverage | percent | Percentage of the entire MRS voxel lying inside the segmentation image field of view, not necessarily inside an ROI. |
| `weighted_mean` — Overlap-weighted ROI mean | same as underlying metric | Sum(overlap fraction × finite voxel estimate) / sum(overlap fraction) over eligible voxels. Fits remain on the original grid; spectra are not averaged/refitted. All positive overlaps contribute by default. QC flags alone do not exclude estimates. A mean of voxel ratios is not generally a ratio of mean amplitudes, and is not a pure-ROI concentration. |
| `contributing_voxels` — Contributing voxel count | count | Integer count of contributing voxels (also called weighted_voxel_count); not the sum of weights or an effective sample size. Counts can differ by metabolite and scale. |
| `minimum_overlap_percent` — Minimum ROI overlap threshold | percent | A reporting-only inclusion threshold relative to whole MRS voxel volume. Equality qualifies; zero overlap never contributes. Blank/n/a means no configured minimum. It is not a universal recommended quality threshold. |
| `mean` — Arithmetic mean | same as underlying metric | Mean across finite estimates in completed voxels. A mean of voxel uncertainty percentages is descriptive, not uncertainty of the mean. |
| `sd` — Across-voxel standard deviation (SD) | same as underlying metric | Population SD (ddof=0) across included voxel values, not standard error and not the per-metabolite fit SD. Native FSL uncertainties/*_sd maps instead contain percentage fit uncertainty. |
| `median` — Median | same as underlying metric | 50th percentile across included finite voxel values. |
| `p05` — 5th percentile | same as underlying metric | 5th percentile across included finite voxel values; not a confidence limit. |
| `p95` — 95th percentile | same as underlying metric | 95th percentile across included finite voxel values; not a confidence limit. |
| `n` — Number of finite values | count | Number of finite metric values included in the summary, not necessarily all fitted voxels. |
| `status` — Processing status | categorical | fit_completed means processing returned a fit, not scientifically validated quantification. zero_signal, nonfinite_signal and outside_fit_mask describe exclusion. |
| `reason` — Exclusion reason | categorical | Reason a voxel was not included in fitting; blank for eligible voxels. |
| `qc_flags` — Quality-control (QC) flags | categorical | Semicolon-separated operational flags: crlb_above_operational_threshold, snr_below_operational_threshold, linewidth_above_operational_threshold. Flags retain estimates and do not automatically exclude them from ROI/group summaries. Thresholds are user choices, not universal acceptance rules. |
| `branch` — Processing branch | categorical | original is the fitted branch without additional residual-water removal (normally after coil combination); waterremoved includes that optional removal. No preferred branch is automatically selected. |
| `difference` — Paired branch difference | same as underlying metric | Water-removed value minus original value in matched eligible voxels. Differences of percentage metrics are percentage points, not relative percent changes. No branch is automatically selected. |
| `fid` — Free induction decay (FID) | input signal units | Complex time-domain MRS signal; Fourier transformation produces a spectrum. Real, imaginary and magnitude displays are different representations. |
| `ppm` — Chemical shift | ppm | Parts-per-million frequency coordinate; spectral window boundaries refer to chemical shift, not concentration. |
| `te_tr` — Echo time (TE) and repetition time (TR) | seconds in BIDS metadata | Acquisition timing relevant to basis compatibility and relaxation weighting; often displayed in milliseconds elsewhere. |
| `basis` — Spectral basis | component-dependent | Simulated or measured component spectra used by the fitting model. Successful loading/resampling does not establish acquisition compatibility. |
| `roi` — Region of interest (ROI) | label | Anatomical region defined by segmentation or mask. Coverage of the ROI by the grid and purity of an MRS voxel have different denominators. |

Model nuisance parameters: `p0` is phase in degrees; `p1` is first-order phase in seconds; `shift_group*` is in ppm; `gamma_group*`, `sigma_group*` and `combined_lw_group*` are added broadening in Hz, excluding intrinsic basis linewidth.

Metabolite abbreviations: Ala: Alanine; Cr: Creatine; PCr: Phosphocreatine; NAA: N-acetylaspartate; NAAG: N-acetylaspartylglutamate; GPC: Glycerophosphocholine; PCh: Phosphocholine; Glu: Glutamate; Gln: Glutamine; GSH: Glutathione; Lac: Lactate; mI: Myo-inositol; sI: Scyllo-inositol; Tau: Taurine; H2O: Water; CrCH2: Auxiliary creatine methylene basis component. Combined labels joined by + are sums of fitted components, not constraints on their ratios.

Sources: https://open.oxcin.ox.ac.uk/pages/fslcourse/lectures/MRS4_data_quality.pdf, https://open.oxcin.ox.ac.uk/pages/fsl/fsl_mrs/fitting.html

Definitions are generated from `metabolix/output_definitions.py` for participant reports and output sidecars. Group mode writes a readable `group_summary.md` and per-column definitions in `group_summary.json`; the TSV remains one row per participant. Existing `crlb_percent_*` names are retained for compatibility, including MH-derived uncertainties whose method is identified in metadata.
