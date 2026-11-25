"""Shared, human-readable definitions for Metabolix outputs."""
import html

SOURCES = [
    'https://open.oxcin.ox.ac.uk/pages/fslcourse/lectures/MRS4_data_quality.pdf',
    'https://open.oxcin.ox.ac.uk/pages/fsl/fsl_mrs/fitting.html',
]
METABOLITES = dict(Ala='Alanine', Cr='Creatine', PCr='Phosphocreatine', NAA='N-acetylaspartate',
                  NAAG='N-acetylaspartylglutamate', GPC='Glycerophosphocholine', PCh='Phosphocholine',
                  Glu='Glutamate', Gln='Glutamine', GSH='Glutathione', Lac='Lactate',
                  mI='Myo-inositol', sI='Scyllo-inositol', Tau='Taurine', H2O='Water',
                  CrCH2='Auxiliary creatine methylene basis component')

# Label, units, description. These descriptions are shared by reports and sidecars.
TERMS = {
 'crlb_percent': ('CRLB / relative fit uncertainty', 'percent',
  'CRLB means Cramér–Rao lower bound: a model-based lower bound on estimator variance under statistical assumptions. '
  'For Newton fitting, FSL-MRS reports 100 × sqrt(CRLB variance) / fitted estimate, i.e. relative SD, not variance. '
  '20 means estimated SD about 20% of the fitted value; it is not a concentration, percent of reference, or 95% confidence interval. '
  'For MH (Metropolis–Hastings) fitting, the same legacy crlb_percent label holds posterior-sample SD divided by the posterior mean × 100, not a CRLB. '
  'The inspected FSL-MRS implementation caps percentages above 999 and maps undefined NaN percentages to 999. '
  '999 is a capped/undefined result, not a precise uncertainty. Near-zero estimates can have very large relative uncertainty. '
  'Low uncertainty does not exclude artifacts, incorrect basis/sequence assumptions, bias or partial-volume mixing. '
  'There is no universal acceptance threshold; uncertainty of an internal ratio or ROI mean is not automatically propagated.'),
 'raw_fit': ('Raw fitted amplitude', 'arbitrary units', 'FSL-MRS fitted amplitude before internal-reference scaling. Not a molar concentration; scaling can differ between spectra.'),
 'internal_reference_ratio': ('Internal-reference ratio', 'dimensionless', 'Metabolite estimate divided by the sum of configured reference metabolites within the same voxel (default Cr+PCr). 1.5 means 150% of reference, not 150% of the total spectrum. A zero/nonfinite reference is invalid. Native/published FSL maps may encode invalid or unfitted values as zero; inspect voxel-table validity.'),
 'snr': ('Signal-to-noise ratio (SNR)', 'dimensionless', 'FSL-MRS component-specific spectral QC SNR, using its noise estimation and matched-filter conventions. Not the residual-water/metabolite peak ratio; values depend on implementation and component and are not a universal quality score.'),
 'linewidth_fwhm_hz': ('Full width at half maximum (FWHM)', 'Hz', 'FSL-MRS component spectral width at half maximum. Broader peaks can impair separation, but linewidth alone does not establish fit quality. Distinct from added nuisance broadening.'),
 'residual_l2_over_input_l2': ('Relative residual norm', 'dimensionless', 'Euclidean norm of the complex time-domain residual divided by the complex input FID norm. Descriptive fit discrepancy; not a percentage of metabolite loss, a concentration uncertainty or a standalone quality criterion.'),
 'overlap_pct': ('Nominal ROI overlap', 'percent', '100 × ROI volume inside an MRS voxel / entire nominal MRS voxel volume. Fine anatomical sampling does not upsample MRS signal or model its spatial response. Existing image alignment is assumed.'),
 'overlap_mm3': ('Nominal ROI volume within voxel', 'mm^3', 'Estimated intersection volume of the ROI with the nominal MRS voxel.'),
 'sampling_difference_pp': ('Overlap sampling sensitivity', 'percentage points', 'Absolute difference between fine and 1 mm overlap estimates. A numerical sensitivity check, not a rigorous error bound or a relative percent change.'),
 'segmentation_fov_pct': ('Segmentation field-of-view (FOV) coverage', 'percent', 'Percentage of the entire MRS voxel lying inside the segmentation image field of view, not necessarily inside an ROI.'),
 'weighted_mean': ('Overlap-weighted ROI mean', 'same as underlying metric', 'Sum(overlap fraction × finite voxel value) / sum(overlap fraction) over eligible voxels. Fits remain on the original grid; spectra are not averaged/refitted. All positive overlaps contribute by default. QC flags alone do not exclude values. A mean of voxel ratios is not generally a ratio of mean amplitudes, and is not a pure-ROI concentration. QC metric means are descriptive averages, not QC recomputed from an ROI spectrum; a mean of CRLB-derived values is not uncertainty of the ROI mean.'),
 'contributing_voxels': ('Contributing voxel count', 'count', 'Integer count of contributing voxels (also called weighted_voxel_count); not the sum of weights or an effective sample size. Counts can differ by metabolite and scale.'),
 'minimum_overlap_percent': ('Minimum ROI overlap threshold', 'percent', 'A reporting-only inclusion threshold relative to whole MRS voxel volume. Equality qualifies; zero overlap never contributes. Blank/n/a means no configured minimum. It is not a universal recommended quality threshold.'),
 'mean': ('Arithmetic mean', 'same as underlying metric', 'Mean across finite estimates in completed voxels. A mean of voxel uncertainty percentages is descriptive, not uncertainty of the mean.'),
 'sd': ('Across-voxel standard deviation (SD)', 'same as underlying metric', 'Population SD (ddof=0) across included voxel values, not standard error and not the per-metabolite fit SD. Native FSL uncertainties/*_sd maps instead contain percentage fit uncertainty.'),
 'median': ('Median', 'same as underlying metric', '50th percentile across included finite voxel values.'),
 'p05': ('5th percentile', 'same as underlying metric', '5th percentile across included finite voxel values; not a confidence limit.'),
 'p95': ('95th percentile', 'same as underlying metric', '95th percentile across included finite voxel values; not a confidence limit.'),
 'n': ('Number of finite values', 'count', 'Number of finite metric values included in the summary, not necessarily all fitted voxels.'),
 'status': ('Processing status', 'categorical', 'fit_completed means processing returned a fit, not scientifically validated quantification. zero_signal, nonfinite_signal and outside_fit_mask describe exclusion.'),
 'reason': ('Exclusion reason', 'categorical', 'Reason a voxel was not included in fitting; blank for eligible voxels.'),
 'qc_flags': ('Quality-control (QC) flags', 'categorical', 'Semicolon-separated operational flags: crlb_above_operational_threshold, snr_below_operational_threshold, linewidth_above_operational_threshold. Flags retain estimates and do not automatically exclude them from ROI/group summaries. Thresholds are user choices, not universal acceptance rules.'),
 'branch': ('Processing branch', 'categorical', 'original is the fitted branch without additional residual-water removal (normally after coil combination); waterremoved includes that optional removal. No preferred branch is automatically selected.'),
 'difference': ('Paired branch difference', 'same as underlying metric', 'Water-removed value minus original value in matched eligible voxels. Differences of percentage metrics are percentage points, not relative percent changes. No branch is automatically selected.'),
 'fid': ('Free induction decay (FID)', 'input signal units', 'Complex time-domain MRS signal; Fourier transformation produces a spectrum. Real, imaginary and magnitude displays are different representations.'),
 'ppm': ('Chemical shift', 'ppm', 'Parts-per-million frequency coordinate; spectral window boundaries refer to chemical shift, not concentration.'),
 'te_tr': ('Echo time (TE) and repetition time (TR)', 'seconds in BIDS metadata', 'Acquisition timing relevant to basis compatibility and relaxation weighting; often displayed in milliseconds elsewhere.'),
 'basis': ('Spectral basis', 'component-dependent', 'Simulated or measured component spectra used by the fitting model. Successful loading/resampling does not establish acquisition compatibility.'),
 'roi': ('Region of interest (ROI)', 'label', 'Anatomical region defined by segmentation or mask. Coverage of the ROI by the grid and purity of an MRS voxel have different denominators.'),
}


def describe(metric, config=None):
    """Return a definition for an exact output column, with acquisition context."""
    if metric.startswith('waterremoved_minus_original_'):
        info = describe(metric.removeprefix('waterremoved_minus_original_'), config)
        return {**info, 'Description': TERMS['difference'][2] + ' Underlying metric: ' + info['Description'],
                'Units': 'percentage points' if info['Units'] == 'percent' else info['Units']}
    key = metric
    metabolite = None
    for prefix in ('raw_fit', 'internal_reference_ratio', 'crlb_percent', 'snr', 'linewidth_fwhm_hz'):
        if metric.startswith(prefix + '_'):
            key, metabolite = prefix, metric[len(prefix)+1:]
            break
    for suffix in ('sampling_difference_pp', 'overlap_pct', 'overlap_mm3'):
        if metric.endswith('_' + suffix): key = suffix
    aliases = {'overlap_weighted_voxel_estimate':'weighted_mean', 'weighted_voxel_count':'contributing_voxels',
               'minimum_overlap_threshold_percent':'minimum_overlap_percent'}
    key = aliases.get(key, key)
    if key in TERMS:
        title, units, description = TERMS[key]
    elif metric.startswith('nuisance_'):
        parameter = metric.removeprefix('nuisance_')
        if parameter == 'p0':
            units, description = 'degrees', 'Fitted zero-order spectral phase.'
        elif parameter == 'p1':
            units, description = 'seconds', 'FSL-MRS first-order phase expressed as an equivalent time shift.'
        elif parameter.startswith('shift_group'):
            units, description = 'ppm', 'Fitted frequency shift for a model metabolite group.'
        elif parameter.startswith(('gamma_group','sigma_group','combined_lw_group')):
            units, description = 'Hz', 'Additional model broadening: gamma is Lorentzian, sigma Gaussian, combined_lw their combined width. Excludes linewidth already implicit in the basis; not a concentration.'
        else:
            units, description = 'see native FSL metadata', 'Fitted nuisance parameter; not a metabolite concentration.'
        title = 'Model nuisance parameter'
    elif metric in ('i','j','k'):
        title, units, description = 'Voxel index', 'index', 'Zero-based spatial voxel index on the original MRS grid.'
    elif metric.startswith('entity_') or metric == 'participant_id':
        title, units, description = 'BIDS identifier', 'identifier', 'Source participant/acquisition identifier; preserve leading zeros.'
    elif metric in ('voxel_rows','fit_completed_rows','invalid_or_excluded_rows'):
        title, units, description = 'Voxel status count', 'count', {'voxel_rows':'All voxel rows in the branch.', 'fit_completed_rows':'Rows with completed fitting; not a count of scientifically acceptable spectra.', 'invalid_or_excluded_rows':'Rows not marked fit_completed.'}[metric]
    elif metric in ('metric','interpretation'):
        title, units, description = 'Metric description', 'text', 'Metric identifier or interpretation note; consult the associated metric definition.'
    else:
        title, units, description = metric, 'not specified', 'Source-defined label; consult its source metadata.'
    result = dict(LongName=title, Units=units, Description=description)
    if metabolite:
        result['Metabolite'] = metabolite
        result['MetaboliteName'] = ' + '.join(METABOLITES.get(m,m) for m in metabolite.split('+'))
    if key in ('crlb_percent','internal_reference_ratio'):
        fit = (config or {}).get('fit', {})
        result['FitAlgorithm'] = (fit.get('algorithm') or 'Newton') if config is not None else 'See source run provenance'
        result['InternalReference'] = (fit.get('internal_reference') or ['Cr','PCr']) if config is not None else 'See source run provenance (FSL default Cr+PCr)'
    if key == 'crlb_percent':
        result['UncertaintyMethod'] = ('CRLB-derived relative SD' if result['FitAlgorithm']=='Newton' else
                                       'Posterior-sample relative SD, not CRLB' if result['FitAlgorithm']=='MH' else
                                       'Newton: CRLB-derived; MH: posterior-sample-derived')
        result['DisplayCapPercent'] = 999
    return result


def glossary():
    return {key: describe(key) for key in TERMS}


def glossary_markdown(config=None):
    lines = ['## Output labels and units', '',
             'These definitions describe estimates and QC, not acceptance criteria. Blank/n/a estimates are missing; map zeros may also encode unfitted or invalid results.', '']
    if config is not None:
        fit = config.get('fit', {})
        lines += [f"Fit algorithm: {fit.get('algorithm') or 'Newton'}. Internal reference: {' + '.join(fit.get('internal_reference') or ['Cr','PCr'])}.", '']
    lines += ['| Label | Units | Meaning |', '|---|---|---|']
    for key in TERMS:
        info = describe(key, config)
        lines.append(f"| `{key}` — {info['LongName']} | {info['Units']} | {info['Description']} |")
    lines += ['', 'Model nuisance parameters: `p0` is phase in degrees; `p1` is first-order phase in seconds; `shift_group*` is in ppm; `gamma_group*`, `sigma_group*` and `combined_lw_group*` are added broadening in Hz, excluding intrinsic basis linewidth.', '',
              'Metabolite abbreviations: ' + '; '.join(f'{key}: {value}' for key,value in METABOLITES.items()) + '. Combined labels joined by + are sums of fitted components, not constraints on their ratios.', '',
              'Sources: ' + ', '.join(SOURCES), '']
    return '\n'.join(lines)


def glossary_html(config=None):
    rows = ''.join('<tr><td>'+html.escape(key)+' — '+html.escape(describe(key,config)['LongName'])+'</td><td>'+html.escape(describe(key,config)['Units'])+'</td><td>'+html.escape(describe(key,config)['Description'])+'</td></tr>' for key in TERMS)
    context = ''
    if config is not None:
        fit=config.get('fit',{})
        context='<p>Fit algorithm: '+html.escape(fit.get('algorithm') or 'Newton')+'. Internal reference: '+html.escape(' + '.join(fit.get('internal_reference') or ['Cr','PCr']))+'.</p>'
    return '<section id="output-definitions"><h2>Output labels and units</h2>'+context+'<p>Missing estimates are blank/n/a; map zeros can encode unfitted or invalid results. QC flags are not automatic exclusions.</p><table><thead><tr><th>Label</th><th>Units</th><th>Meaning</th></tr></thead><tbody>'+rows+'</tbody></table><p>Nuisance parameters: p0 phase (degrees), p1 first-order phase (seconds), shift_group* (ppm), gamma/sigma/combined_lw_group* added broadening (Hz). These are not metabolite concentrations.</p><p>'+html.escape('; '.join(k+': '+v for k,v in METABOLITES.items()))+'. Labels joined by + denote sums.</p><p>Sources: '+', '.join('<a href="'+u+'">FSL-MRS documentation</a>' for u in SOURCES)+'</p></section>'
