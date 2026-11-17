import argparse
import json
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / 'tools'))
from common import ROOT, digest, execute, save
from plans import DATASETS, SUITES, plan
from data import sample

def validate(folder, row):
    folder = Path(folder)
    log = (folder / 'run.log').read_text(errors='replace')
    if re.search(r'^hdump done ', log, re.M):
        if not (folder / 'hypergraph.updates').is_file() or not (folder / 'hypergraph.queries.tsv').is_file():
            raise ValueError('Missing exported hypergraph trace')
        return {'contract': 'hypergraph-export-only', 'numerical_maintenance_timed': False}
    batches = [int(v) for v in re.findall(r'^batch=(\d+)', log, re.M)]
    if row['version'] == 'trace':
        accepted = re.findall(r'^accepted batch=(\d+)', log, re.M)
        certificates = [json.loads(line) for line in log.splitlines() if line.startswith('{') and '"rigorous"' in line]
        if list(map(int, accepted)) != list(range(row['batches'] + 1)) or len(certificates) != row['batches'] + 1 or not all(x.get('rigorous') for x in certificates):
            raise ValueError('Missing rigorous per-batch acceptance')
        return {'accepted_checkpoints': len(certificates), 'contract': 'independent-outward-k4-epsilon-1/100'}
    if batches != list(range(1, row['batches'] + 1)) or not re.search(r'^memory ', log, re.M):
        raise ValueError('Incomplete numerical trajectory')
    checks = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in log.splitlines() if line.startswith('check ')]
    for check in checks:
        if 'exact' in check and check['exact'] != check['cS'] and float(check['exact']) != float(check['cS']):
            raise ValueError('Witness recount mismatch')
        if float(check.get('load_relerr', 0)) > 1e-6:
            raise ValueError('Load recount mismatch')
    return {'batches': len(batches), 'numerical_checks': len(checks), 'contract': 'binary64-numerical-closure'}

def run_job(row, args):
    row = dict(row, env=dict(row['env']))
    if args.version:
        if row['version'] != 'guarded':
            raise ValueError('Version override applies only to the ordinary numerical suites')
        row['version'] = args.version
    for name in ['batch', 'batches', 'threads', 'hold', 'check']:
        value = getattr(args, name)
        if value is not None:
            row[name] = value
    if args.limit:
        row['limit'] = args.limit
    if args.graph:
        graph = args.graph.resolve()
        if args.mode:
            row['mode'] = args.mode
    else:
        graph = args.data.resolve() / DATASETS[row['dataset']]['file']
    if row.get('sample') and row['sample'] != 100:
        sampled = args.data.resolve() / f'{row["dataset"]}_v{row["sample"]}.bin'
        if not sampled.exists():
            sample(graph, sampled, row['sample'] / 100, 47)
        graph = sampled
    if not graph.is_file():
        raise FileNotFoundError(f'{graph}; run tools/data.py fetch {row["dataset"]}')
    output = (args.output / row['name']).resolve()
    output.mkdir(parents=True, exist_ok=True)
    version = row['version']
    executable = ROOT / 'versions/v7/implicert-v7.exe' if version == 'v7' else ROOT / 'build' / (version + '.exe')
    if version == 'v7' and digest(executable) != json.loads((ROOT / 'versions/v7/binary.json').read_text())['sha256']:
        raise ValueError('Historical v7 executable identity mismatch')
    if not executable.exists():
        raise FileNotFoundError(f'{executable}; run tools/build.py {version}')
    env = dict(row['env'], OMP_NUM_THREADS=str(row['threads']))
    if row.get('fixed_edges'):
        trace = args.output.resolve() / 'fixed-edges.txt'
        key = 'KCDS_WRITE_EDGES' if row['arm'] == 'hybrid' else 'KCDS_READ_EDGES'
        if key == 'KCDS_READ_EDGES' and not trace.exists():
            raise FileNotFoundError('Run the hybrid policy first to record the shared adversarial trace')
        env[key] = str(trace)
    if args.affinity:
        env['MEMCAP_AFFINITY'] = args.affinity
    if version == 'trace':
        if row['k'] != 4 or float(row['epsilon']) != .01:
            raise ValueError('Published strict barrier supports k=4 and epsilon=0.01')
        checkpoint = output / 'checkpoints'
        checkpoint.mkdir(exist_ok=True)
        env.update(TRACE_ARM=row.get('arm', 'maintained'), TRACE_DIR=str(checkpoint), TRACE_PYTHON=sys.executable, TRACE_CHECKER=str(ROOT / 'verify_checkpoint.py'), TRACE_NATIVE=str(ROOT / 'build/verify_k4.exe'), KCDS_NOCHECK='1')
    elif args.snapshots:
        env['KCDS_DUMP'] = str(output / 'snapshot')
    if args.hypergraph:
        if not 2 <= row['k'] <= 6:
            raise ValueError('Hyperedge export requires 2 <= k <= 6')
        env['KCDS_HDUMP'] = str(output / 'hypergraph')
    command = [ROOT / 'build/memcap.exe', args.memory_mb, executable, graph, row['k'], row['epsilon'], row['batch'], row['batches'], row['mode'], row['seed'], row['check'], row['hold'], row['strategy'], row['q']]
    record = execute(command, output, env, row['limit'], {'suite': args.suite, 'job': row, 'graph_sha256': digest(graph), 'executable_sha256': digest(executable)}, cwd=output)
    if record['status'] == 'ok' and 'validation' not in record:
        try:
            record['validation'] = validate(output, row)
        except (ValueError, KeyError) as error:
            record.update(status='invalid', validation_error=str(error))
        save(output / 'record.json', record)
    return record

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['list', 'run'])
    parser.add_argument('suite', choices=SUITES)
    parser.add_argument('--datasets', nargs='+')
    parser.add_argument('--case')
    parser.add_argument('--data', type=Path, default=ROOT / 'datasets')
    parser.add_argument('--output', type=Path, default=ROOT / 'runs')
    parser.add_argument('--version', choices=['guarded', 'v7', 'v7b', 'v7c'])
    parser.add_argument('--graph', type=Path)
    parser.add_argument('--mode', choices=['replay', 'stream', 'events', 'core', 'plant'])
    parser.add_argument('--batch', type=int)
    parser.add_argument('--batches', type=int)
    parser.add_argument('--threads', type=int)
    parser.add_argument('--hold')
    parser.add_argument('--check', type=int)
    parser.add_argument('--limit', type=float)
    parser.add_argument('--memory-mb', type=int, default=32768)
    parser.add_argument('--affinity')
    parser.add_argument('--snapshots', action='store_true')
    parser.add_argument('--hypergraph', action='store_true')
    args = parser.parse_args()
    if args.snapshots and args.hypergraph:
        parser.error('Hypergraph export skips the solver; export snapshots in a separate numerical run')
    rows = [row for row in plan(args.suite) if (not args.datasets or row['dataset'] in args.datasets) and (not args.case or args.case in row['name'])]
    if not rows:
        parser.error('No matching cases')
    if args.action == 'list':
        print(json.dumps(rows, indent=2))
        return
    results = [run_job(row, args) for row in rows]
    if any(row['status'] != 'ok' for row in results):
        raise SystemExit(1)

if __name__ == '__main__':
    main()
