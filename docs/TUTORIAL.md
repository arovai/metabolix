# Metabolix Tutorial

This walkthrough uses the example dataset and basis paths recorded in the project notes. It prepares commands but does not imply that a full-dataset run has been executed or scientifically approved.

## Install and inspect dependencies

FSL-MRS is expected in an activated FSL/Conda environment. The wrapper can be installed into that environment with pip, or installed separately if the FSL-MRS executables are on `PATH`:

```bash
python -m pip install metabolix-bids
metabolix --version
metabolix --check-deps
```

The application is a pip package; FSL-MRS and FSL binaries are not. See the README installation section for official setup.

## Example paths

```bash
DATASET=/home/arovai/datasets_shortcuts/2021-Thankful_Alligator-640c8fba6b3e
BIDS_DIR="$DATASET/rawdata"
OUTPUT_DIR="$DATASET/derivatives/metabolix"
BASIS="$DATASET/derivatives/fsl-mrsi/mrscloud/20261005/LCModel_GE_UnEdited_PRESS_144_.BASIS"
FS_DIR="$DATASET/derivatives/freesurfer_7.3.2"
```

The example MRS file is `rawdata/sub-11/mrs/sub-11_mrsi.nii.gz`. It uses the documented compatibility naming pattern rather than the standard `_mrs` suffix. Metabolix records that distinction and does not claim the compatibility path has been fully BIDS-validated.

If generating a new basis for a different acquisition, MRSCloud is one possible route: submit the scanner/sequence, field strength, TE, and metabolite selection, and request an LCModel `.BASIS` export if available. The existing example basis was generated for GE unedited PRESS TE144, but its sequence equivalence still requires scientific review. See [Obtaining a basis](../README.md#obtaining-a-basis).

## Dry-run all spatial voxels

Generate an editable acquisition-profile copy, then set paths for this dataset:

```bash
metabolix --write-example-config mrsi.yaml
```

Edit `basis.path` to `$BASIS` and `roi.freesurfer_dir` to `$FS_DIR`, or provide those values as CLI overrides. Inspect the selected input, metadata validation, and resolved configuration without running FSL-MRS:

```bash
metabolix "$BIDS_DIR" "$OUTPUT_DIR" participant \
  --participant-label 11 \
  --config mrsi.yaml \
  --basis "$BASIS" \
  --freesurfer-dir "$FS_DIR" \
  --dry-run \
  --export-config resolved.yaml
```

The default fit mask is every finite, nonzero spatial voxel. No four-voxel pilot mask or thalamic restriction is used. Dry-run validates discovery and input metadata only; it does not establish basis compatibility or fitting success.

## Run the participant workflow

This command performs processing and fitting across the selected participant's MRS acquisitions. Run it only when an actual-data execution has been authorized and the chosen basis/model has been reviewed:

```bash
metabolix "$BIDS_DIR" "$OUTPUT_DIR" participant \
  --participant-label 11 \
  --config mrsi.yaml \
  --basis "$BASIS" \
  --freesurfer-dir "$FS_DIR" \
  --roi-labels 10 49 \
  --n-jobs 4
```

The example profile uses GE unedited PRESS TE144 component selection and an internal Cr+PCr reference. Combined metabolites are display/quantification sums; `--combine` does not constrain component ratios. These choices are not universal or independently validated for every acquisition.

To prepare a matched water-removal sensitivity comparison, both branches use the same basis, mask, fitting limits and component model:

```bash
metabolix "$BIDS_DIR" "$OUTPUT_DIR" participant \
  --participant-label 11 \
  --config mrsi.yaml \
  --basis "$BASIS" \
  --freesurfer-dir "$FS_DIR" \
  --roi-labels 10 49 \
  --water-removal compare \
  --water-ppm 4.5 4.9 \
  --n-jobs 4
```

Water-removal comparison is optional. A prominent water peak alone does not establish that removal improves metabolite estimates. The tool does not choose a winner by residual size or uncertainty.

## Outputs and interpretation

Results are written under `OUTPUT_DIR/sub-11/[ses-...]/mrs/<BIDS-entity-stem>_desc-metabolix/`. Open the matching `sub-..._desc-metabolix_report.html` first; it links the BIDS-named voxel QC tables, downloadable statistic maps, figures, native FSL-MRS products, logs, ROI outputs, and provenance. The acquisition bundle has tidy `figures/`, `maps/`, `qc/`, `roi/`, `logs/`, `provenance/`, `work/`, and `processing/fsl-mrs/` folders. Compare-mode paired differences are in `qc/`; each branch's voxel table includes all spatial coordinates and invalid-input statuses.

The report includes pre/post coil-combination spectra and mean FID-magnitude QC, before/after water-removal QC when requested, a plot of the selected basis components, and heatmaps for each available internally referenced metabolite map. Each fit branch also includes a paginated BIDS-named PDF showing every fitted voxel's spectrum, fit, baseline, and residual. The plots use FSL-MRS's reader/spectral conventions and are generated without requiring a manually selected diagnostic voxel.

FreeSurfer `aseg.mgz` labels 10 and 49 denote left/right thalamus in this preset. Overlap is estimated by fine world-affine sampling, uses the entire nominal MRS voxel as denominator, and reports segmentation field-of-view coverage separately. It assumes existing scanner/world affines are aligned; no registration is performed. MRSI spatial response may spread signal beyond nominal boundaries.

The ROI weighted summary is `roi/sub-11_desc-metabolix-roi-weighted_summary.tsv`
(with additional source entities when applicable). It contains separate raw and
internal-reference weighted means per ROI, metabolite and processing branch,
along with overlap-weighted means of available component CRLB-derived uncertainty,
SNR and FWHM linewidth.
Default ROI averaging includes all positive overlaps; `--minimum-overlap-percent`
can restrict reporting, without restricting fitting. This averages fitted values,
not spectra, and does not upsample MRS data. Internal values are ratios to the
configured reference (Cr+PCr by default); multiply by 100 for percent of reference.
Read [ROI estimates and metabolite units](ROI_ANALYSIS.md) before interpreting
small ROIs, contributing counts, QC flags or uncertainty.

Newton `*_sd` products are CRLB-derived percentage estimates, not absolute standard deviations and not fully propagated uncertainty on ratios. Ratios with zero reference are undefined. Never describe arbitrary fit amplitudes or creatine ratios as molar concentrations. Fit completion, anatomical overlap, spectral quality, fitting uncertainty and quantification units are distinct claims.

To aggregate existing participant QC tables without refitting:

```bash
metabolix "$BIDS_DIR" "$OUTPUT_DIR" group
```

Group output is `OUTPUT_DIR/group/group_summary.tsv`, with one row per participant and a JSON column dictionary.

## Participant-wide group table

Group mode reads existing participant results without running any fitting:

```bash
metabolix "$BIDS_DIR" "$OUTPUT_DIR" group
```

`group/group_summary.tsv` contains **one row per analysed participant**, keyed by
`participant_id` (for example `sub-12`). It includes whole-grid statistics and,
when available, participant ROI metabolite and QC summaries. There is no pooled
`GROUP` row. Session, acquisition, run and other non-subject entities, processing
branches and ROIs are kept in separate columns, never silently averaged together.

Example columns for an acquisition without additional entities:

```text
participant_id
age
sex
group
metabolix__branch-original__scope-grid__internal_reference_ratio_NAA__mean
metabolix__branch-original__scope-roi__roi-left_thalamus__internal_reference_ratio_NAA__weighted_mean
metabolix__branch-waterremoved__scope-roi__roi-right_thalamus__raw_fit_NAA__weighted_mean
```

Grid metrics retain finite-value count (`n`), mean, population SD and median,
plus fit-status counts. ROI metrics retain the overlap-weighted mean,
contributing voxel count and configured overlap threshold. ROI estimates are
copied from participant summaries, using their original reporting settings.
Internal values remain ratios, not percentages; multiply by 100 for percent of
reference. QC flags alone do not exclude otherwise eligible estimates.

If `BIDS_DIR/participants.tsv` exists, all its columns are automatically joined
by `participant_id`. An optional `sub-` prefix is normalized on both sides;
leading zeros are preserved (`sub-01` is distinct from `sub-1`). Metadata values
are read as text, preserving categories and missing-value markers. Analysed
participants without metadata remain in the table with `n/a` values; participants
with no selected analysis are not added. Duplicate normalized IDs, malformed
metadata, conflicting column names or duplicate measurements cause a clear error
instead of silently dropping or averaging records. Missing `participant_id` in
an existing participants file is an error. Without a participants file, only
analysis columns and participant IDs are written.

`group/group_summary.json` documents each measure column's branch, ROI, metric,
statistic, acquisition entities and source tables. Double underscores separate
column components; special characters and embedded separators are URL-encoded.
Missing estimates are `n/a`. For the overlap-threshold field, `n/a` denotes no
configured minimum, as in the participant ROI summary.

Use `--participant-label` to select a subset; session/acquisition/run filters are
also respected. After preserving an older summary if needed, use `--overwrite`
in **group mode** to regenerate the group TSV and JSON. This does not overwrite
participant fits. Changing the group format does not require rerunning participant
processing. `--dry-run` validates aggregation and the metadata join without
writing outputs.


Participant HTML/Markdown reports include an output-label glossary. NIfTI maps and TSV tables carry JSON descriptions and units; voxel tables also retain their legacy `_definitions.json` sidecar. `group_summary.md` explains group labels, while `group_summary.json` describes every measure column. See the [README glossary](../README.md#output-labels-and-units) for CRLB, SNR, linewidth, internal ratios, ROI overlap and distribution statistics.
