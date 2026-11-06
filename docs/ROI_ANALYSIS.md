# ROI estimates and metabolite units

## What is averaged?

Metabolix fits each eligible original MRS voxel independently, then averages its
fitted metabolite estimates using anatomical overlap weights. It does not average
ROI spectra and refit them, and does not upsample the spectroscopic signal.

For a given ROI, metabolite, scale and processing branch:

```text
ROI estimate = sum(overlap_fraction[v] * estimate[v]) / sum(overlap_fraction[v])
overlap_fraction[v] = ROI volume inside voxel v / entire nominal MRS voxel volume
```

The sums include eligible voxels only. A voxel must have `status=fit_completed`, a
finite estimate for that metric, positive overlap, and satisfy any configured
overlap threshold. Missing estimates are excluded separately for each metric, so
contributing voxel counts can differ between metabolites and scales.

**QC flags do not currently exclude a voxel from this average.** Configured SNR,
linewidth or uncertainty flags remain available for review in the voxel-QC table.
A completed fit is not necessarily an acceptable measurement. ROI summaries do
not currently include propagated uncertainty.

## Overlap sampling and thresholds

Anatomical overlap is estimated with midpoint samples spaced at most 0.5 mm
within each original MRS voxel by default (`--roi-step-mm`). These samples locate
segmentation labels using image-world affines; they do not create additional MRS
measurements or improve MRS resolution. A 1 mm sampling comparison measures
numerical sensitivity, not a rigorous error bound. Segmentation-field-of-view
coverage is reported separately; the overlap denominator remains the whole voxel.

The default `roi.minimum_overlap_percent: null` includes every positive overlap.
To require at least 20% of each contributing MRS voxel to belong to the ROI, use:

```bash
--minimum-overlap-percent 20
```

Or configure:

```yaml
roi:
  minimum_overlap_percent: 20
```

Twenty percent is an illustrative setting, not a validated universal threshold.
The threshold affects ROI reporting only, not the fitting mask. Equality to the
threshold is included; zero-overlap voxels are always excluded. If no eligible
estimate remains, the summary value is blank and the contributing count is zero.
Choose reporting settings before a run: changing this setting currently changes
the resume fingerprint; it is not a reporting-only resume override. Preserve
existing results rather than using `--overwrite` just to change a threshold.

## Small ROIs and partial-volume mixing

If three voxels each overlap an ROI by 3%, their normalized weights are equal and
the ROI estimate is their ordinary mean. Each estimate still reflects the entire
MRS voxel, including the approximately 97% outside the ROI. Weighting cannot
isolate a small structure's metabolite concentration. Increasing a threshold can
leave no usable voxels and does not increase acquisition resolution.

Inspect both voxel overlap percentages and the fraction of the whole ROI covered
by the grid. These have different denominators: high ROI coverage does not imply
high ROI purity within an MRS voxel. Geometry assumes correct existing alignment;
no registration is performed. MRSI spatial response can also spread signal beyond
nominal voxel boundaries.

## Raw values, internal ratios and percentages

| Scale or prefix | Meaning |
|---|---|
| `raw_fit_<metabolite>` / native `concs/raw` | Fitted amplitude in arbitrary units; not a molar concentration. |
| `internal_reference_ratio_<metabolite>` / native `concs/internal` | Voxel-wise metabolite estimate divided by the internal metabolite reference. |
| Newton `*_sd` / `crlb_percent_*` | CRLB-derived percentage uncertainty, not the metabolite's percentage of its reference. |

The reference defaults to total creatine, `Cr + PCr`, unless changed with
`--internal-ref` or `fit.internal_reference`. For example:

```text
internal NAA = fitted NAA / (fitted Cr + fitted PCr)
1.5 internal units = 150% of the reference
```

This is neither a percentage of the total spectrum nor an absolute concentration.
Multiply an internal ratio by 100 to express it as a percentage of the reference;
Metabolix does not currently export a separate percentage column. The reference
components are summed; a combined metabolite numerator is available only if
produced by the fit's reporting configuration.

The ROI internal estimate is the overlap-weighted mean of voxel-wise ratios. It
is generally **not** the ratio of the weighted mean metabolite and reference
amplitudes. The ROI raw estimate is a separate weighted average of raw amplitudes.

A zero or nonfinite reference makes internal ratios invalid. Metabolix marks
these as missing in its voxel table when the reference component maps are
available, and excludes them from the corresponding ROI mean. Native/published
FSL maps may still contain zeros used for invalid values; consult the table and
QC rather than interpreting these as biological absence. Outside-mask zeros
also do not represent measurements. Percentage uncertainty may be capped at
999 and is not automatically propagated ratio or ROI uncertainty.

## Finding the outputs

Within each acquisition bundle, for example
`OUTPUT_DIR/sub-12/mrs/sub-12_desc-metabolix/`, inspect:

| File in `roi/` | Contents |
|---|---|
| `sub-12_desc-metabolix-roi-weighted_summary.tsv` | One row per branch, ROI and available raw/internal metabolite metric. |
| `sub-12_desc-metabolix-roi-overlap_voxelqc.tsv` | Per-voxel ROI percentages, overlap volumes, sampling differences and segmentation FOV coverage. |
| `sub-12_desc-metabolix-roi_summary.json` | Geometry, labels, sampling settings and ROI grid-coverage summaries. |
| `*_statmap.nii.gz` and sidecars | ROI overlap and segmentation FOV coverage maps, in percent. |

Session, acquisition and other applicable source entities are retained in names.
The weighted summary has these columns:

| Column | Meaning |
|---|---|
| `branch` | `original` or `waterremoved`. |
| `roi` | Label name, such as `left_thalamus`. |
| `metric` | Scale plus metabolite, such as `internal_reference_ratio_NAA`. |
| `overlap_weighted_voxel_estimate` | Weighted mean defined above; blank if no eligible estimates. |
| `weighted_voxel_count` | Integer number of contributing voxels, despite its name; not a sum of weights or effective sample size. |
| `minimum_overlap_threshold_percent` | Applied threshold; blank means no minimum. |
| `interpretation` | Reminder that this is not a pure-ROI concentration or an averaged-spectrum fit. |

Use `qc/` voxel tables alongside these summaries to examine individual fits and
flags. The ROI table currently includes raw and internal estimates, not all QC
map types or propagated ROI uncertainty. See the [tutorial](TUTORIAL.md) for
processing commands and [limitations](LIMITATIONS.md) for broader QC guidance.
