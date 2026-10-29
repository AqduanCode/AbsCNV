#!/usr/bin/env python3
import csv
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
RAW = Path(os.environ.get('CNV_RAW_DIR', str(SCRIPT_DIR / 'raw')))
TRUTH_PATH = RAW / 'truth.tsv'
OUT = Path(os.environ.get('CNV_EVAL_OUT', str(SCRIPT_DIR / 'output')))
EVENT_STANDARD = os.environ.get('CNV_EVENT_STANDARD', 'my_method_greedy')
GAIN_LOG2 = math.log2(1.1)
LOSS_LOG2 = math.log2(0.9)

INPUTS = {
    'CNVkit': RAW / 'CNVkit.call.cns',
    'ControlFREEC': RAW / 'ControlFREEC.CNVs',
    'ichorCNA': RAW / 'ichorCNA.seg',
    'BICseq2': RAW / 'BICseq2.tsv',
}

REMOTE_BASE = '/home/data2/DAQ/HCC1395_like_CN_50events_benchmark_20260921'
DISPLAY_PATHS = {
    'truth': f'{REMOTE_BASE}/benchmark/Sample_B_vs_Sample_A.benchmark.tsv',
    'CNVkit': f'{REMOTE_BASE}/Tools/CNVkit_10kb/Sample_B.sorted.call.cns',
    'ControlFREEC': f'{REMOTE_BASE}/Tools/ControlFREEC_10kb/Sample_B.sorted.bam_CNVs',
    'ichorCNA': f'{REMOTE_BASE}/Tools/ichorCNA_10kb/Sample_B_vs_Sample_A_10kb.seg',
    'BICseq2': f'{REMOTE_BASE}/Tools/BICseq2_10kb/results/Sample_B_vs_Sample_A.10kb.bicseq2_segments.tsv',
}


def chrom_norm(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith('chr') else 'chr' + chrom


def chrom_key(chrom):
    value = chrom_norm(chrom)[3:]
    return (0, int(value)) if value.isdigit() else (1, value)


def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def parse_number(value):
    """Return a finite float, or None for a blank/invalid numeric field."""
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def direction_from_cn(value, neutral_cn=2.0):
    number = parse_number(value)
    if number is None:
        return None
    if number > neutral_cn:
        return 'gain'
    if number < neutral_cn:
        return 'loss'
    return None


def direction_from_label(value):
    label = str(value or '').strip().upper()
    if not label or label in {'NEUT', 'NEUTRAL', 'DIPLOID', 'NORMAL'}:
        return None
    loss_tokens = ('LOSS', 'DEL', 'DELETION', 'HETD', 'HOMD')
    gain_tokens = ('GAIN', 'DUP', 'DUPLICATION', 'AMP', 'HLAMP')
    if any(token in label for token in loss_tokens):
        return 'loss'
    if any(token in label for token in gain_tokens):
        return 'gain'
    return None


def load_truth():
    result = []
    for row in read_tsv(TRUTH_PATH):
        delta = int(row['delta_B_minus_A'])
        if delta == 0:
            continue
        result.append({
            'id': row['event_id'], 'chrom': chrom_norm(row['chrom']),
            'start': int(row['start']), 'end': int(row['end']),
            'direction': 'gain' if delta > 0 else 'loss',
            'source': row['pair_id'], 'value': str(delta),
        })
    return result


def load_cnvkit(path):
    calls = []
    for i, row in enumerate(read_tsv(path), 1):
        # Native absolute CN has priority. Continuous log2 is only a fallback
        # when CN is absent or invalid.
        cn = parse_number(row.get('cn'))
        if cn is not None:
            direction = direction_from_cn(cn)
            source, value = 'absolute_cn', str(row.get('cn')).strip()
        else:
            log2_value = parse_number(row.get('log2'))
            direction = ('gain' if log2_value is not None and log2_value > GAIN_LOG2
                         else ('loss' if log2_value is not None and log2_value < LOSS_LOG2 else None))
            source, value = 'log2_threshold_fallback', str(row.get('log2', '')).strip()
        if direction:
            calls.append({'id': f'CNVkit_{i}', 'chrom': chrom_norm(row['chromosome']),
                          'start': int(row['start']), 'end': int(row['end']),
                          'direction': direction, 'source': source, 'value': value})
    return calls


def load_freec(path):
    calls = []
    with open(path) as handle:
        for i, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                continue
            direction = direction_from_label(fields[4])
            # The native file supplies gain/loss. If that label is ever absent,
            # use the reported absolute copy number before any signal fallback.
            source = 'native_event_label'
            if direction is None and str(fields[4]).strip().upper() not in {'NEUT', 'NEUTRAL'}:
                direction = direction_from_cn(fields[3])
                source = 'absolute_cn_fallback'
            if direction is None:
                continue
            calls.append({'id': f'ControlFREEC_{i}', 'chrom': chrom_norm(fields[0]),
                          'start': int(fields[1]), 'end': int(fields[2]),
                          'direction': direction, 'source': source, 'value': fields[3]})
    return calls


def load_ichor(path):
    calls = []
    for i, row in enumerate(read_tsv(path), 1):
        # Absolute copy.number has priority; event is the second native source;
        # median logR is used only if neither native field gives a direction.
        cn = parse_number(row.get('copy.number'))
        event = str(row.get('event', '')).strip()
        if cn is not None:
            direction = direction_from_cn(cn)
            source, value = f'absolute_cn;event={event}', str(row.get('copy.number')).strip()
        else:
            direction = direction_from_label(event)
            source, value = f'native_event={event}', event
            if direction is None and event.upper() not in {'NEUT', 'NEUTRAL'}:
                median = parse_number(row.get('median'))
                direction = ('gain' if median is not None and median > GAIN_LOG2
                             else ('loss' if median is not None and median < LOSS_LOG2 else None))
                source, value = 'median_logR_threshold_fallback', str(row.get('median', '')).strip()
        if direction:
            # ichorCNA .seg coordinates are 1-based inclusive; convert to 0-based half-open.
            calls.append({'id': f'ichorCNA_{i}', 'chrom': chrom_norm(row['chr']),
                          'start': int(row['start']) - 1, 'end': int(row['end']),
                          'direction': direction, 'source': source, 'value': value})
    return calls


def load_bicseq(path):
    calls = []
    rows = read_tsv(path)
    if not rows:
        return calls
    log_col = 'log2.copyRatio'
    if log_col not in rows[0]:
        raise ValueError(f'BICseq2 column {log_col!r} not found; columns={list(rows[0])}')
    chrom_col = 'chrom' if 'chrom' in rows[0] else 'chromosome'
    start_col = 'start'
    end_col = 'end'
    for i, row in enumerate(rows, 1):
        value = float(row[log_col])
        direction = 'gain' if value > GAIN_LOG2 else ('loss' if value < LOSS_LOG2 else None)
        if direction:
            # BICseq2 report.R emits 1-based inclusive coordinates.
            calls.append({'id': f'BICseq2_{i}', 'chrom': chrom_norm(row[chrom_col]),
                          'start': int(float(row[start_col])) - 1, 'end': int(float(row[end_col])),
                          'direction': direction, 'source': 'log2.copyRatio', 'value': str(value)})
    return calls


LOADERS = {'CNVkit': load_cnvkit, 'ControlFREEC': load_freec,
           'ichorCNA': load_ichor, 'BICseq2': load_bicseq}


def merge_calls(calls, tool):
    grouped = defaultdict(list)
    for call in calls:
        if call['end'] > call['start']:
            grouped[(call['chrom'], call['direction'])].append(call)
    merged = []
    counter = 0
    for (chrom, direction), items in grouped.items():
        items.sort(key=lambda x: (x['start'], x['end']))
        current_start, current_end = items[0]['start'], items[0]['end']
        members = [items[0]['id']]
        for item in items[1:]:
            if item['start'] <= current_end:
                current_end = max(current_end, item['end'])
                members.append(item['id'])
            else:
                counter += 1
                merged.append({'id': f'{tool}_merged_{counter}', 'chrom': chrom,
                               'start': current_start, 'end': current_end,
                               'direction': direction, 'source': ','.join(members), 'value': ''})
                current_start, current_end, members = item['start'], item['end'], [item['id']]
        counter += 1
        merged.append({'id': f'{tool}_merged_{counter}', 'chrom': chrom,
                       'start': current_start, 'end': current_end,
                       'direction': direction, 'source': ','.join(members), 'value': ''})
    return sorted(merged, key=lambda x: (chrom_key(x['chrom']), x['start'], x['end'], x['direction']))


def overlap(a, b):
    if a['chrom'] != b['chrom'] or a['direction'] != b['direction']:
        return 0
    return max(0, min(a['end'], b['end']) - max(a['start'], b['start']))


def maximum_event_matching(preds, truths):
    candidates = {}
    for pi, pred in enumerate(preds):
        edges = []
        plen = pred['end'] - pred['start']
        for ti, truth in enumerate(truths):
            ov = overlap(pred, truth)
            if not ov:
                continue
            tlen = truth['end'] - truth['start']
            pf, tf = ov / plen, ov / tlen
            if pf >= 0.5 and tf >= 0.5:
                edges.append((ti, min(pf, tf), ov, pf, tf))
        candidates[pi] = sorted(edges, key=lambda x: (-x[1], -x[2], x[0]))
    truth_to_pred = {}

    def augment(pi, seen):
        for ti, *_ in candidates[pi]:
            if ti in seen:
                continue
            seen.add(ti)
            if ti not in truth_to_pred or augment(truth_to_pred[ti], seen):
                truth_to_pred[ti] = pi
                return True
        return False

    for pi in sorted(candidates, key=lambda x: (len(candidates[x]), x)):
        augment(pi, set())
    matches = []
    for ti, pi in sorted(truth_to_pred.items()):
        pred, truth = preds[pi], truths[ti]
        ov = overlap(pred, truth)
        matches.append((pi, ti, ov, ov / (pred['end'] - pred['start']),
                        ov / (truth['end'] - truth['start'])))
    return matches


def greedy_event_matching(preds, truths):
    """Match exactly as compare_cnv_accuracy.py used for the user's method."""
    candidates = []
    for pi, pred in enumerate(preds):
        plen = pred['end'] - pred['start']
        for ti, truth in enumerate(truths):
            ov = overlap(pred, truth)
            if not ov:
                continue
            tlen = truth['end'] - truth['start']
            pf, tf = ov / plen, ov / tlen
            if pf >= 0.5 and tf >= 0.5:
                union = plen + tlen - ov
                candidates.append((min(pf, tf), ov / union, ov, pf, tf, pi, ti))
    candidates.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
    used_p, used_t, matches = set(), set(), []
    for score, iou, ov, pf, tf, pi, ti in candidates:
        if pi in used_p or ti in used_t:
            continue
        used_p.add(pi); used_t.add(ti)
        matches.append((pi, ti, ov, pf, tf))
    return matches


def union_intervals(records):
    grouped = defaultdict(list)
    for row in records:
        grouped[(row['chrom'], row['direction'])].append((row['start'], row['end']))
    result = {}
    for key, intervals in grouped.items():
        intervals.sort()
        merged = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        result[key] = merged
    return result


def union_size(unions, direction='all'):
    return sum(end - start for (chrom, direc), intervals in unions.items()
               if direction == 'all' or direc == direction for start, end in intervals)


def intersection_size(pred_union, truth_union, direction='all'):
    total = 0
    for key, pred_intervals in pred_union.items():
        if direction != 'all' and key[1] != direction:
            continue
        truth_intervals = truth_union.get(key, [])
        i = j = 0
        while i < len(pred_intervals) and j < len(truth_intervals):
            ps, pe = pred_intervals[i]; ts, te = truth_intervals[j]
            total += max(0, min(pe, te) - max(ps, ts))
            if pe <= te: i += 1
            else: j += 1
    return total


def scores(tp, predicted, actual):
    precision = tp / predicted if predicted else (1.0 if actual == 0 else 0.0)
    recall = tp / actual if actual else 1.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def main():
    missing = [str(path) for path in [TRUTH_PATH, *INPUTS.values()] if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError('Missing/empty required inputs:\n' + '\n'.join(missing))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'normalized_calls').mkdir(exist_ok=True)
    truth = load_truth()
    event_metrics, base_metrics, match_rows, unmatched_truth, unmatched_pred = [], [], [], [], []
    manifest = [{'name': 'truth', 'path': DISPLAY_PATHS['truth'], 'bytes': TRUTH_PATH.stat().st_size}]
    truth_union = union_intervals(truth)

    for tool, path in INPUTS.items():
        raw = LOADERS[tool](path)
        if EVENT_STANDARD == 'my_method_greedy':
            # The user's original event-level evaluator did not merge prediction
            # fragments before matching. Keep raw calls for an apples-to-apples comparison.
            calls = sorted(raw, key=lambda x: (chrom_key(x['chrom']), x['start'], x['end'], x['direction']))
        elif EVENT_STANDARD == 'merged_maximum':
            calls = merge_calls(raw, tool)
        else:
            raise ValueError('CNV_EVENT_STANDARD must be my_method_greedy or merged_maximum')
        manifest.append({'name': tool, 'path': DISPLAY_PATHS[tool], 'bytes': path.stat().st_size})
        write_tsv(OUT / 'normalized_calls' / f'{tool}.tsv', calls,
                  ['id', 'chrom', 'start', 'end', 'direction', 'source', 'value'])
        for direction in ('all', 'gain', 'loss'):
            pidx = [i for i, x in enumerate(calls) if direction == 'all' or x['direction'] == direction]
            tidx = [i for i, x in enumerate(truth) if direction == 'all' or x['direction'] == direction]
            subpred, subtruth = [calls[i] for i in pidx], [truth[i] for i in tidx]
            matches = (greedy_event_matching(subpred, subtruth)
                       if EVENT_STANDARD == 'my_method_greedy'
                       else maximum_event_matching(subpred, subtruth))
            tp = len(matches); precision, recall, f1 = scores(tp, len(subpred), len(subtruth))
            event_metrics.append({'tool': tool, 'direction': direction, 'tp': tp,
                                  'fp': len(subpred)-tp, 'fn': len(subtruth)-tp,
                                  'predicted_events': len(subpred), 'truth_events': len(subtruth),
                                  'precision': f'{precision:.9f}', 'recall': f'{recall:.9f}', 'f1': f'{f1:.9f}'})
            if direction == 'all':
                matched_p, matched_t = set(), set()
                for pi, ti, ov, pf, tf in matches:
                    pred, tru = subpred[pi], subtruth[ti]
                    matched_p.add(pi); matched_t.add(ti)
                    match_rows.append({'tool': tool, 'prediction_id': pred['id'], 'truth_id': tru['id'],
                                       'chrom': pred['chrom'], 'direction': pred['direction'],
                                       'pred_start': pred['start'], 'pred_end': pred['end'],
                                       'truth_start': tru['start'], 'truth_end': tru['end'],
                                       'overlap_bp': ov, 'prediction_overlap_fraction': f'{pf:.9f}',
                                       'truth_overlap_fraction': f'{tf:.9f}'})
                for i, tru in enumerate(subtruth):
                    if i not in matched_t:
                        unmatched_truth.append({'tool': tool, **tru})
                for i, pred in enumerate(subpred):
                    if i not in matched_p:
                        unmatched_pred.append({'tool': tool, **pred})
        pred_union = union_intervals(calls)
        for direction in ('all', 'gain', 'loss'):
            pb = union_size(pred_union, direction); tb = union_size(truth_union, direction)
            matched = intersection_size(pred_union, truth_union, direction)
            precision, recall, f1 = scores(matched, pb, tb)
            base_metrics.append({'tool': tool, 'direction': direction, 'matched_bp': matched,
                                 'predicted_cnv_bp': pb, 'truth_cnv_bp': tb,
                                 'false_positive_bp': pb-matched, 'false_negative_bp': tb-matched,
                                 'precision': f'{precision:.9f}', 'recall': f'{recall:.9f}', 'f1': f'{f1:.9f}'})

    write_tsv(OUT / 'event_level_metrics.tsv', event_metrics,
              ['tool','direction','tp','fp','fn','predicted_events','truth_events','precision','recall','f1'])
    write_tsv(OUT / 'base_level_metrics.tsv', base_metrics,
              ['tool','direction','matched_bp','predicted_cnv_bp','truth_cnv_bp','false_positive_bp','false_negative_bp','precision','recall','f1'])
    write_tsv(OUT / 'event_matches.tsv', match_rows,
              ['tool','prediction_id','truth_id','chrom','direction','pred_start','pred_end','truth_start','truth_end','overlap_bp','prediction_overlap_fraction','truth_overlap_fraction'])
    write_tsv(OUT / 'unmatched_truth_events.tsv', unmatched_truth,
              ['tool','id','chrom','start','end','direction','source','value'])
    write_tsv(OUT / 'unmatched_predictions.tsv', unmatched_pred,
              ['tool','id','chrom','start','end','direction','source','value'])
    write_tsv(OUT / 'input_manifest.tsv', manifest, ['name','path','bytes'])
    if EVENT_STANDARD == 'my_method_greedy':
        event_note = ('- Event calls are not merged before evaluation.\n'
                      '- Event-level matching exactly follows the user-method evaluator: '
                      'greedy one-to-one assignment ranked by minimum reciprocal overlap, then IoU, then overlap bp.')
    else:
        event_note = ('- Same-tool overlapping or directly adjacent calls of the same direction are merged.\n'
                      '- Event-level matching uses maximum-cardinality one-to-one assignment.')
    definition = f'''# Direction-aware CNV evaluation\n\nTruth: `{DISPLAY_PATHS['truth']}`\nEvent standard: `{EVENT_STANDARD}`\n\n## CNV classification priority\n\n- Coordinate system: all records normalized to 0-based half-open.\n- Direction: Sample_B relative to Sample_A; gain and loss must agree.\n- Absolute copy-number rule: CN > 2 is gain, CN < 2 is loss, and CN = 2 is neutral.\n- CNVkit: use native `cn`; use `log2` threshold only when `cn` is missing or invalid.\n- Control-FREEC: use native gain/loss label; use reported absolute CN only if that label is missing or unrecognized.\n- ichorCNA: use native `copy.number`; then native `event`; use `median` threshold only when both native fields are unavailable or unrecognized.\n- BIC-seq2: no absolute CN or categorical gain/loss field is present, so use `log2.copyRatio`.\n- Continuous-signal fallback thresholds: gain if ratio > 1.1 (log2 > {GAIN_LOG2:.9f}); loss if ratio < 0.9 (log2 < {LOSS_LOG2:.9f}).\n- Native CN/label and continuous signal can disagree; native CN/label takes priority by design in this reorganization.\n\n## Evaluation\n\n{event_note}\n- Event-level TP requires both prediction overlap fraction >= 0.5 and truth overlap fraction >= 0.5.\n- Event precision = matched pairs / predicted events; recall = matched pairs / truth events.\n- Base-level: direction-specific interval unions; matched bases are same-direction intersections.\n- Base precision = matched bp / predicted CNV bp; recall = matched bp / truth CNV bp.\n- F1 is the harmonic mean of precision and recall.\n'''
    (OUT / 'README.md').write_text(definition)
    print(OUT)
    for row in event_metrics:
        if row['direction'] == 'all': print('EVENT', row)
    for row in base_metrics:
        if row['direction'] == 'all': print('BASE', row)


if __name__ == '__main__':
    main()
