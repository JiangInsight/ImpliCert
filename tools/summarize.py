import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import re
import statistics
from common import save

def percentile(values, probability):
    return sorted(values)[max(0, math.ceil(probability * len(values)) - 1)]

def parse(folder):
    record = json.loads((folder / 'record.json').read_text())
    rows = (folder / 'run.log').read_text(errors='replace').splitlines()
    batch = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in rows if line.startswith('batch=')]
    accepted = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in rows if line.startswith('accepted ')]
    work = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in rows if line.startswith('work ')]
    latency = [float(row['total_s']) for row in accepted if int(row['batch']) > 0] if accepted else [float(row['upd'].rstrip('s')) + float(row['repair'].rstrip('s')) for row in batch]
    baseline = record['identity'].get('baseline')
    if baseline in ('cap', 'cq') and (folder / 'result.query_times.tsv').exists():
        with (folder / 'result.query_times.tsv').open() as handle:
            queries = list(csv.DictReader(handle, delimiter='\t'))
        scale, update, query = (1e-9, 'update_ns', 'query_ns') if baseline == 'cap' else (1, 'update_s', 'query_s')
        cumulative = [(float(row[update]) + float(row[query])) * scale for row in queries]
        latency = [right - left for left, right in zip(cumulative, cumulative[1:])]
    elif baseline == 'udshp' and (folder / 'result.csv.gz').exists():
        with gzip.open(folder / 'result.csv.gz', 'rt') as handle:
            queries = list(csv.DictReader(handle))
        with gzip.open(folder / 'result.updates.csv.gz', 'rt') as handle:
            updates = list(csv.DictReader(handle))
        cumulative = []
        total = 0
        for row in updates:
            total += int(row['service_update_ns'])
            cumulative.append(total)
        latency = [(cumulative[int(right['step'])] - cumulative[int(left['step'])] + int(right['query_ns'])) * 1e-9 for left, right in zip(queries, queries[1:])]
    result = {'case': folder.name, 'status': record['status'], 'contract': record.get('validation', {}).get('contract', 'unvalidated'), 'batches': len(latency), 'wall_s': record['wall_s'], 'version': record['identity'].get('job', {}).get('version', 'baseline')}
    if latency:
        result.update(mean_s=statistics.mean(latency), median_s=statistics.median(latency), p95_s=percentile(latency, .95), max_s=max(latency), stdev_s=statistics.stdev(latency) if len(latency) > 1 else 0)
    if work:
        for key in ['full_passes', 'full_visits', 'prefix_visits', 'local_visits']:
            if all(key in row for row in work):
                result[key] = sum(int(row[key]) for row in work)
    result['affected_cliques'] = sum(int(row.get('cliques_upd', 0)) for row in batch)
    result['reported_local_cliques'] = sum(int(row.get('cliques_rep', 0)) for row in batch)
    result['complete_work_including_fallback'] = bool(work and all(all(key in row for key in ['full_visits', 'prefix_visits', 'local_visits']) for row in work))
    checks = [line for line in rows if line.startswith('check ')]
    warm, cold = [], []
    for line in checks:
        w = re.search(r'warm_global iters=\d+ prep=([\d.]+)s iter=([\d.]+)s', line)
        c = re.search(r'cold_global iters=\d+ iter=([\d.]+)s', line)
        p = re.search(r'cold_prep=([\d.]+)s', line)
        if w:
            warm.append(float(w[1]) + float(w[2]))
        if c and p:
            cold.append(float(c[1]) + float(p[1]))
    if warm:
        result['warm_checkpoint_mean_s'] = statistics.mean(warm)
    if cold:
        result['cold_checkpoint_mean_s'] = statistics.mean(cold)
        if latency and statistics.mean(latency) > 0:
            result['cold_checkpoint_over_batch_mean'] = statistics.mean(cold) / statistics.mean(latency)
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    rows = [parse(path.parent) for path in sorted(args.directory.rglob('record.json'))]
    if not rows:
        raise ValueError('No records found')
    save(args.directory / 'summary.json', rows)
    grouped = {}
    for path, row in zip(sorted(args.directory.rglob('record.json')), rows):
        record = json.loads(path.read_text())
        job = record['identity'].get('job')
        if not job or row['status'] != 'ok' or 'mean_s' not in row:
            continue
        key = tuple((field, str(job.get(field, ''))) for field in ['dataset', 'version', 'arm', 'strategy', 'mode', 'k', 'epsilon', 'batch', 'batches', 'threads']) + (('graph_sha256', record['identity']['graph_sha256']),)
        grouped.setdefault(key, []).append(row['mean_s'])
    repeats = [dict(key, runs=len(values), mean_of_run_means_s=statistics.mean(values), median_run_mean_s=statistics.median(values), stdev_run_means_s=statistics.stdev(values) if len(values) > 1 else None, min_run_mean_s=min(values), max_run_mean_s=max(values)) for key, values in grouped.items()]
    save(args.directory / 'repeat-summary.json', repeats)
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (args.directory / 'summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, keys)
        writer.writeheader()
        writer.writerows(rows)
    print(len(rows), 'records summarized')

if __name__ == '__main__':
    main()
