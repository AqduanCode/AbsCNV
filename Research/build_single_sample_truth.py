"""Build the single-sample event truth from its absolute-CN simulation profile.

Uses the original truth loader and same-direction merge definitions. Coordinates
are preserved; adjacent/overlapping non-CN2 intervals of the same direction merge.
"""
import argparse, csv
from collections import defaultdict
from pathlib import Path
TRUTH_PATH = None

def chrom_norm(value):
    value = str(value)
    return value if value.startswith('chr') else 'chr' + value

def chrom_key(value):
    value = chrom_norm(value)[3:]
    return (0, int(value)) if value.isdigit() else (1, value)

def direction(cn):
    if cn > 2:
        return 'gain'
    if cn < 2:
        return 'loss'
    return 'neutral'

def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def load_truth_segments():
    result = []
    for i, row in enumerate(read_tsv(TRUTH_PATH), 1):
        cn = float(row['absolute_CN'])
        result.append({
            'id': f'truth_{i}', 'chrom': chrom_norm(row['chrom']),
            'start': int(row['start']), 'end': int(row['end']),
            'direction': direction(cn), 'cn': cn, 'source': row['state'],
            'value': str(cn),
        })
    return result

def merge_direction_calls(calls, prefix):
    grouped = defaultdict(list)
    for call in calls:
        if call['direction'] != 'neutral' and call['end'] > call['start']:
            grouped[(call['chrom'], call['direction'])].append(call)
    merged = []
    counter = 0
    for (chrom, direc), records in grouped.items():
        records.sort(key=lambda x: (x['start'], x['end']))
        start, end = records[0]['start'], records[0]['end']
        members = [records[0]['id']]
        for record in records[1:]:
            if record['start'] <= end:
                end = max(end, record['end'])
                members.append(record['id'])
            else:
                counter += 1
                merged.append({'id': f'{prefix}_{counter}', 'chrom': chrom,
                               'start': start, 'end': end, 'direction': direc,
                               'source': ','.join(members), 'value': '', 'cn': ''})
                start, end, members = record['start'], record['end'], [record['id']]
        counter += 1
        merged.append({'id': f'{prefix}_{counter}', 'chrom': chrom,
                       'start': start, 'end': end, 'direction': direc,
                       'source': ','.join(members), 'value': '', 'cn': ''})
    return sorted(merged, key=lambda x: (chrom_key(x['chrom']), x['start'], x['end']))

def main():
    global TRUTH_PATH
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args=p.parse_args(); TRUTH_PATH=args.profile
    rows=merge_direction_calls([x for x in load_truth_segments() if x['direction']!='neutral'], 'truth_event')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['id','chrom','start','end','direction','source','value','cn'],delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)
    print(f'Wrote {len(rows)} truth events to {args.output}')

if __name__=='__main__': main()
