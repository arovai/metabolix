"""Participant-wide group tables from existing voxel and ROI summaries."""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any
from urllib.parse import quote

import numpy as np

from metabolix.output_definitions import describe, glossary, glossary_markdown


def _read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f'Duplicate column names in {path}.')
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f'Malformed TSV row in {path}.')
        return fields, rows


def _participant_id(value: str) -> str:
    value = value.strip()
    label = value.removeprefix('sub-')
    if not label or label.lower() in {'n/a', 'nan'}:
        raise ValueError('Missing participant_id in group inputs.')
    return 'sub-' + label


def _part(value: str) -> str:
    # Encode separators inside values so different tuples cannot share a column.
    return quote(value, safe='-+').replace('__', '%5F%5F')


def run_group(bids_dir: Path, output_dir: Path, config: dict[str, Any],
              dry_run: bool, logger: logging.Logger) -> int:
    tables = sorted(path for path in output_dir.rglob('*_voxelqc.tsv')
                    if path.parent.name == 'qc' and 'waterRemovalComparison' not in path.name)
    if not tables:
        raise ValueError(f'No participant voxel-QC tables found under {output_dir}.')
    summary_path = output_dir / 'group/group_summary.tsv'
    dictionary_path = summary_path.with_suffix('.json')
    if not dry_run and (summary_path.exists() or dictionary_path.exists()) and not config['execution']['overwrite']:
        raise ValueError(f'Group summary already exists: {summary_path}. Use --overwrite to regenerate group outputs only.')

    participant_path = bids_dir / 'participants.tsv'
    metadata_fields: list[str] = []
    metadata: dict[str, dict[str, str]] = {}
    if participant_path.exists():
        fields, rows = _read_table(participant_path)
        if 'participant_id' not in fields:
            raise ValueError(f'{participant_path} requires a participant_id column.')
        metadata_fields = [key for key in fields if key != 'participant_id']
        for row in rows:
            participant = _participant_id(row['participant_id'])
            if participant in metadata:
                raise ValueError(f'Duplicate participant_id {participant} in {participant_path}.')
            metadata[participant] = {key: row[key] for key in metadata_fields}

    participants: dict[str, dict[str, Any]] = {}
    columns: dict[str, dict[str, Any]] = {}
    source_contexts: dict[str, Any] = {}
    selected = {_participant_id(value) for value in config['selection']['participant_label']}
    selection_keys = {'ses': 'session', 'acq': 'acquisition', 'run': 'run'}

    def add(participant: str, entities: dict[str, str], branch: str, scope: str,
            metric: str, statistic: str, value: Any, source: Path, roi: str | None = None) -> None:
        parts = ['metabolix'] + [f'{key}-{_part(val)}' for key, val in sorted(entities.items()) if key != 'sub']
        parts += [f'branch-{_part(branch)}', f'scope-{scope}']
        if roi is not None:
            parts.append('roi-' + _part(roi))
        # Escape embedded delimiters while keeping ordinary metric names readable.
        parts += [_part(metric), statistic]
        name = '__'.join(parts)
        if name in metadata_fields:
            raise ValueError(f'Participant metadata column conflicts with generated measure: {name}')
        row = participants.setdefault(participant, {'participant_id': participant})
        if name in row:
            raise ValueError(f'Duplicate measurement for {participant}, {name}; check duplicate runs/tables. No averaging performed.')
        row[name] = value
        entry = columns.setdefault(name, dict(entities={k: v for k, v in entities.items() if k != 'sub'},
                                              branch=branch, scope=scope, roi=roi, metric=metric,
                                              statistic=statistic, source_tables=[]))
        entry['MetricDefinition'] = describe(metric)
        entry['StatisticDefinition'] = describe(statistic)
        entry['Units'] = (entry['StatisticDefinition']['Units']
                          if statistic in {'n', 'count', 'contributing_voxels', 'minimum_overlap_percent'}
                          else entry['MetricDefinition']['Units'])
        if statistic == 'count':
            entry['StatisticDefinition'] = dict(LongName='Count', Units='count', Description='Number of voxel rows with the indicated processing status.')
            entry['Units'] = 'count'
        relative = str(source.relative_to(output_dir))
        if relative not in entry['source_tables']:
            entry['source_tables'].append(relative)

    metric_prefixes = ('raw_fit_', 'internal_reference_ratio_', 'crlb_percent_', 'snr_', 'linewidth_fwhm_hz_')
    table_count = 0
    for table in tables:
        fields, rows = _read_table(table)
        if not rows:
            continue
        groups: dict[tuple[tuple[tuple[str, str], ...], str], list[dict[str, str]]] = {}
        for row in rows:
            entities = {key.removeprefix('entity_'): val for key, val in row.items()
                        if key.startswith('entity_') and val}
            participant = _participant_id(entities.get('sub', ''))
            if selected and participant not in selected:
                continue
            if any(config['selection'].get(option) and entities.get(key) not in config['selection'][option]
                   for key, option in selection_keys.items()):
                continue
            branch = row.get('branch', '')
            if not branch:
                raise ValueError(f'Missing branch in {table}.')
            groups.setdefault((tuple(sorted(entities.items())), branch), []).append(row)
        if not groups:
            continue
        table_count += 1
        manifests = list((table.parent.parent / 'provenance').glob('*_provenance.json'))
        context = None
        if len(manifests) == 1:
            manifest = json.loads(manifests[0].read_text())
            saved_config = manifest.get('resolved_config')
            if isinstance(saved_config, dict):
                fit = saved_config.get('fit', {})
                context = {'FitAlgorithm': fit.get('algorithm') or 'Newton',
                           'InternalReference': fit.get('internal_reference') or ['Cr','PCr'],
                           'FSLMRSVersion': manifest.get('fsl_mrs_version'),
                           'Provenance': str(manifests[0].relative_to(output_dir))}
        source_contexts[str(table.relative_to(output_dir))] = context
        metrics = [key for key in fields if key.startswith(metric_prefixes) or key == 'residual_l2_over_input_l2']
        roi_tables = sorted((table.parent.parent / 'roi').glob('*_desc-metabolix-roi-weighted_summary.tsv'))
        if len(roi_tables) > 1:
            raise ValueError(f'Ambiguous ROI summary tables for {table}.')
        roi_rows = _read_table(roi_tables[0])[1] if roi_tables else []
        for roi_table in roi_tables:
            source_contexts[str(roi_table.relative_to(output_dir))] = context
        for (entity_tuple, branch), voxel_rows in groups.items():
            entities = dict(entity_tuple)
            participant = _participant_id(entities['sub'])
            completed = [row for row in voxel_rows if row.get('status') == 'fit_completed']
            for metric, value in [('voxel_rows', len(voxel_rows)), ('fit_completed_rows', len(completed)),
                                  ('invalid_or_excluded_rows', len(voxel_rows) - len(completed))]:
                add(participant, entities, branch, 'grid', metric, 'count', value, table)
            for metric in metrics:
                values = []
                for row in completed:
                    try:
                        value = float(row[metric])
                    except (KeyError, ValueError):
                        continue
                    if np.isfinite(value):
                        values.append(value)
                array = np.asarray(values, dtype=float)
                statistics = {'n': len(values), 'mean': float(array.mean()) if values else 'n/a',
                              'sd': float(array.std(ddof=0)) if values else 'n/a',
                              'median': float(np.median(array)) if values else 'n/a'}
                for statistic, value in statistics.items():
                    add(participant, entities, branch, 'grid', metric, statistic, value, table)
            for row in roi_rows:
                if row['branch'] != branch:
                    continue
                for field, statistic in [('overlap_weighted_voxel_estimate', 'weighted_mean'),
                                         ('weighted_voxel_count', 'contributing_voxels'),
                                         ('minimum_overlap_threshold_percent', 'minimum_overlap_percent')]:
                    value = row[field] or 'n/a'
                    if value != 'n/a':
                        try:
                            if not np.isfinite(float(value)):
                                value = 'n/a'
                        except ValueError as exc:
                            raise ValueError(f'Invalid numeric ROI value in {roi_tables[0]}: {field}') from exc
                    add(participant, entities, branch, 'roi', row['metric'], statistic,
                        value, roi_tables[0], roi=row['roi'])
    if not participants:
        raise ValueError('No participant rows match the selection.')
    missing = sorted(set(participants) - set(metadata)) if participant_path.exists() else []
    if missing:
        logger.warning('%d analysed participants have no participants.tsv match; metadata will be n/a.', len(missing))
    if dry_run:
        logger.info('DRY RUN group: %d participant rows, %d measure columns, %d metadata columns.',
                    len(participants), len(columns), len(metadata_fields))
        return 0
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ['participant_id', *metadata_fields, *sorted(columns)]
    with summary_path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', restval='n/a')
        writer.writeheader()
        for participant in sorted(participants):
            writer.writerow({**metadata.get(participant, {}), **participants[participant]})
    dictionary_path.write_text(json.dumps({
        'Description': 'One row per analysed participant; sessions/acquisitions/runs, branches and ROIs remain separate columns.',
        'ParticipantMetadataSource': str(participant_path) if participant_path.exists() else None,
        'ParticipantJoin': 'Left join of analysed participants on participant_id; optional sub- prefix normalized, leading zeros preserved.',
        'ParticipantsWithoutMetadata': missing,
        'MissingValue': 'n/a; for minimum_overlap_percent this means no configured minimum.',
        'GridStatistics': 'Finite estimates in fit_completed voxels; population SD (ddof=0). QC flags alone do not exclude voxels.',
        'ROIEstimates': 'Copied from participant overlap-weighted summaries, not recomputed or refitted. No upsampling.',
        'InternalValues': 'Ratios to configured internal reference, not percentages; multiply by 100 for percent of reference.',
        'ColumnNaming': 'Double underscores separate fields. Values are URL-encoded, including embedded double underscores; decode using column metadata.',
        'Glossary': glossary(),
        'SourceFitContexts': source_contexts,
        'Columns': columns,
    }, indent=2) + '\n')
    summary_path.with_suffix('.md').write_text('# Metabolix group summary\n\nOne row per analysed participant. See group_summary.json for each column’s source, metric, statistic and units. Grid uncertainty means/SDs describe voxel distributions, not uncertainty of an ROI/group mean.\n\n' + glossary_markdown())
    logger.info('Wrote %d participant rows from %d voxel-QC tables to %s; no refitting performed.',
                len(participants), table_count, summary_path)
    return 0
