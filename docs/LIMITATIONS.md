# Scientific Limitations and Troubleshooting

## Scientific limits

- Inputs must already be reconstructed image-space NIfTI-MRS. Vendor raw data and k-space are not processed or implicitly reinterpreted.
- Coil handling is based on NIfTI-MRS dimension tags, not a fixed dimension index. Untagged or unsupported higher dimensions fail rather than being averaged.
- An optional coil reference must match the MRS grid, tagged coil count, nucleus, and spectrometer frequency. The example dataset has no identified unsuppressed-water reference; none is inferred or used by default.
- A readable/resampled basis is not necessarily sequence-compatible. PRESS subecho timings, RF waveforms, localization, field strength, echo time and scanner implementation may matter. The GE MRSCloud TE144 basis in the example notes is a candidate, not a validated universal basis; its `.BASIS` export does not include the MM/lipid components from the separate `_wMM.mat` file.
- No separate unsuppressed-water reference is known for the example data. Metabolix does not perform water-reference absolute quantification or reference-based corrections without explicitly configured, validated reference data.
- `--combine` forms sums for reporting; it does not fix component ratios. Restricting the reporting list does not remove model components.
- Residual-water removal is an optional sensitivity branch. Peak-height ratios are not SNR, and a spectral-change norm is not metabolite loss. No branch is automatically chosen.
- FSL-MRS's averaged-FID HTML report is not a voxel-wise fit summary. Native per-voxel products should be inspected in the appropriate viewer.
- Newton `*_sd` maps are CRLB-derived percentage uncertainty estimates, can be capped, and are not absolute SD or necessarily propagated uncertainty for ratios. Zero reference values make ratios undefined. Outside-mask zeros are not measured absence.
- Anatomical overlap is nominal geometry calculated from existing image-world affines. Matching headers do not prove registration. No FSL FLIRT transform is treated as a generic world transform. Finite MRSI point-spread response is not modeled.
- Successful processing does not establish validated metabolite measurements. Report spectral quality, fit uncertainty, nominal ROI overlap and quantification units separately.

## Troubleshooting

**`fsl_mrs_proc` or `fsl_mrsi` not found**: install FSL-MRS using the official FSL installer or its Conda channel, activate that environment, verify with `fsl_mrs --version` and `fsl_mrs_verify`, then run `metabolix --check-deps`. Alternatively set `basis.binary_dir` or `--binary-dir`. Installing `fslpy` alone does not install FSL executables.

**No MRS files discovered**: verify that input is the BIDS root (commonly `rawdata/`) and that files match `sub-*/[ses-*/]mrs/*_mrs.nii[.gz]`. Legacy `*_mrsi.nii[.gz]` names under subject `mrs/` are also supported. Other non-standard layouts need conversion or explicit normalization before use.

**Unsupported dimension tag**: inspect NIfTI-MRS metadata. `DIM_COIL` is handled; dynamic/editing/unknown dimensions are rejected. Do not relabel data without understanding the acquisition dimension.

**No basis or ambiguous basis mapping**: provide `--basis` or define one unambiguous `basis.mapping` match. Check the basis sequence, TE and frequency against the data. Successful resampling alone is not evidence of scientific compatibility.

If a suitable basis is unavailable, one practical option is to submit a sequence-specific basis-generation job through [MRSCloud](https://mricloud.org/) (see [MRSCloud project information](https://github.com/shui5/MRSCloud)). Supply the actual vendor/scanner, sequence/localization, field strength, TE, and metabolite set; request an LCModel `.BASIS` export when available for FSL-MRS. MRSCloud's supported choices and access may vary. Treat the result as a candidate: verify the exported components and metadata and assess sequence agreement before using it. A successful download, parse, or resampling does not validate the basis scientifically.

**Existing output directory**: default behavior refuses overwrite. Use `--resume` only when the run manifest matches the current inputs/configuration, or explicitly select `--overwrite` after preserving any results that should be retained.

**ROI overlap appears implausible**: inspect qform/sform warnings, segmentation coverage, coordinate conventions, and anatomical overlays. Metabolix assumes existing world affines are aligned; it does not register the images.

**Water-removal changes look large**: compare spectra at matched voxels and inspect the relevant spectral intervals. Retain both branches; do not infer metabolite loss or select a branch from a peak-height ratio.