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

Results are written under `OUTPUT_DIR/sub-11/.../` with `dataset_description.json` at the output root. Stage folders retain original/coil-combined/water-removed inputs and native FSL-MRS outputs. `manifest.json` records source/basis checksums, resolved configuration, command argument arrays, timestamps and exit states; stdout/stderr are in `logs/`. Each `03_fit_*/voxel_qc.tsv` retains voxel coordinates and invalid-input statuses; compare-mode paired differences are in the run-level `qc/` folder. ROI overlap tables and percentage maps are separate from the fitting mask.

Each fit branch includes an application Markdown/HTML report and a paginated PDF showing every fitted voxel's spectrum, fit, baseline, and residual. The plots use FSL-MRS's reader/spectral conventions and are generated without requiring a manually selected diagnostic voxel.

FreeSurfer `aseg.mgz` labels 10 and 49 denote left/right thalamus in this preset. Overlap is estimated by fine world-affine sampling, uses the entire nominal MRS voxel as denominator, and reports segmentation field-of-view coverage separately. It assumes existing scanner/world affines are aligned; no registration is performed. MRSI spatial response may spread signal beyond nominal boundaries.

Newton `*_sd` products are CRLB-derived percentage estimates, not absolute standard deviations and not fully propagated uncertainty on ratios. Ratios with zero reference are undefined. Never describe arbitrary fit amplitudes or creatine ratios as molar concentrations. Fit completion, anatomical overlap, spectral quality, fitting uncertainty and quantification units are distinct claims.

To aggregate existing participant QC tables without refitting:

```bash
metabolix "$BIDS_DIR" "$OUTPUT_DIR" group
```

Group output is `OUTPUT_DIR/group/group_summary.tsv`.