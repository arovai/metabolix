# Metabolix Implementation Instructions

This is the finalized coding contract, based on the ln2t `INSTRUCTIONS_TEMPLATE.md` and `specs/agent/TOOL_SPECIFICATION.md`.

## Application contract

- Deliver a pip-installable, CLI-only BIDS app named `metabolix` with positional syntax `BIDS_DIR OUTPUT_DIR {participant|group}`.
- Read only input data. Discover standard BIDS MRS `sub-*/[ses-*/]mrs/*_mrs.nii[.gz]` and explicitly support the supplied legacy `*_mrsi.nii.gz` convention. State the accommodation in output provenance and docs.
- Support participant/session/acquisition/run selection, participant processing, and group aggregation of existing participant results.
- Defaults < YAML < explicitly supplied CLI arguments. Omitted CLI options must not override YAML. Boolean flags must support explicit on/off. Lists replace, not append to, YAML lists. Reject unknown config keys. Resolve relative config paths against the YAML's directory.
- Include default/example config generators, dependency diagnostics, resolved-config export, dry-run, and checksum/config-sensitive resume. Refuse accidental overwrite unless explicitly requested.
- Validate image-space complex NIfTI-MRS dimensions/tags, sampling, nucleus/frequency and available TE/TR, affine/qform/sform, finite/nonzero data, basis and optional masks. Never mistake spatial voxels for k-space or assume a fixed dimension is coils. Reject unsupported dynamic/editing dimensions.
- All valid spatial voxels are fitted by default. Preserve invalid voxels as QC records with reasons. No pilot-voxel restriction or anatomy-driven fitting mask by default.
- Invoke `fsl_mrs_proc` and `fsl_mrsi` with `subprocess` argument arrays, capturing exact argv, timestamps, output, exit status and logs. Find commands on PATH or through explicit configured binary directory. Do not shell-interpolate or hard-code `/opt/fsl/bin`.
- Retain native FSL outputs and inputs per stage. Coil-combine only when metadata indicates a coil dimension; preserve already-combined input. Optional water modes are off/on/compare and share fitting parameters/mask/basis in compare mode. Never infer an optimal branch.
- Use a coil reference only when explicitly configured and validated against the input's spatial grid, tagged coil dimension/count, nucleus, and frequency; never infer or download a reference.
- Basis is explicit or matched unambiguously from configured mapping. Record checksum/components/frequency/bandwidth/points/sequence fields and compatibility cautions. No silent basis substitution. No H2O reference correction without an explicitly validated measured reference.
- Export accessible application-owned QC/provenance. Mark CRLB-derived `*_sd` as percentage uncertainty, not absolute SD; handle capped values cautiously and do not claim ratio uncertainty is propagated. Undefined reference ratios are missing with flags. Outside-mask zeros are not observations.
- ROI overlap is optional and must not change the fit mask. Support FreeSurfer/int-label/binary images. Map world coordinates with affines, sample nominal voxel volumes at <=0.5 mm, compare 1 mm, use full voxel denominator, and report segmentation FOV coverage separately. Preserve qform/sform warnings and registration caveat. No FLIRT transform may be treated as generic world affine.
- Group mode aggregates existing participant QC summaries only.
- Organize outputs as `OUTPUT_DIR/sub-*/[ses-*/]mrs/<entity-stem>_desc-metabolix/`; use entity-rich BIDS names for reports, QC tables, statistic maps, ROI products and logs. Keep native FSL products in `processing/fsl-mrs/` and optional detail in tidy `figures/`, `maps/`, `qc/`, `roi/`, `logs/`, `provenance/`, and `work/` subfolders.
- The HTML report is the user entry point: add pre/post coil combination and optional water-removal spectral/FID QC, basis-component spectra, per-metabolite spatial heatmaps, direct links to QC/maps/ROI/native outputs, and configuration/reproducibility summaries. Use FSL-MRS spectral conventions.
- Provide README, neutral config, acquisition-specific GE PRESS TE144 example config, complete tutorial, scientific limitations/troubleshooting, and real command examples.
- Do not add unit tests, test directories, pytest/coverage dependencies/config or test CI. Do not run expensive whole-dataset processing without explicit authorization.

## Implementation quality and failure behavior

- Keep the runtime small (NumPy, nibabel, PyYAML); do not claim pip installs FSL. Document official FSL installer and Conda channels and actionable diagnostics.
- Paths with spaces must work. Output identities must retain entities and avoid collisions.
- On stage failure, preserve successful earlier outputs and write failure details; exit nonzero.
- A successful fit or visually plausible report is not evidence of validated quantification. Do not assign molar units to arbitrary amplitudes/ratios.
- Validate with syntax/compile, CLI help/version/config checks, dependency diagnostics and dry-run only; do not execute real full-data fitting in this implementation session.