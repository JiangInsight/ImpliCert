import argparse
import csv
import itertools
import json
from pathlib import Path
import struct
import subprocess
import sys
import tarfile
import io
from common import ROOT, digest, execute, save

def csr(path, n, edges):
    adjacency = [[] for _ in range(n)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    offsets = [0]
    for row in adjacency:
        row.sort()
        offsets.append(offsets[-1] + len(row))
    path.write_bytes(struct.pack('<QQ', n, len(edges)) + struct.pack(f'<{n + 1}Q', *offsets) + b''.join(struct.pack(f'<{len(row)}I', *row) for row in adjacency))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'runs/release-validation')
    parser.add_argument('--real-data', action='store_true')
    parser.add_argument('--baselines', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs = output / 'input'
    inputs.mkdir(exist_ok=True)
    graph = inputs / 'tiny.bin'
    with (ROOT / 'data/replay-edges.csv').open() as handle:
        reader = csv.reader(handle)
        next(reader)
        edges = [tuple(map(int, row)) for row in reader]
    csr(graph, 1 + max(max(edge) for edge in edges), edges)
    stress = inputs / 'stress.bin'
    csr(stress, 630, list(itertools.combinations(range(30), 2)) + [(v, v + 1) for v in range(30, 630, 2)])
    invocations = []
    def call(command):
        command = list(map(str, command))
        invocations.append(command)
        subprocess.run(command, cwd=ROOT, check=True)
    for name in ['test_contract', 'test_absorber', 'audit_micro']:
        record = execute([ROOT / 'build' / (name + '.exe')], output / name, {'OMP_NUM_THREADS': '1'}, 120, {'binary_sha256': digest(ROOT / 'build' / (name + '.exe'))})
        if record['status'] != 'ok':
            raise RuntimeError(name + ' failed')
    common = ['--datasets', 'dblp', '--graph', graph, '--batch', '2', '--batches', '4', '--hold', '.15', '--check', '1', '--threads', '1']
    for suite, extra in [('main', ['--snapshots']), ('strict', ['--case', '_r1']), ('guarded', []), ('policies', []), ('strategies', [])]:
        call([sys.executable, ROOT / 'reproduce.py', 'run', suite] + common + extra + ['--output', output / suite])
    for version in ['v7', 'v7b', 'v7c']:
        call([sys.executable, ROOT / 'reproduce.py', 'run', 'main'] + common + ['--version', version, '--output', output / version])
    for suite in ['triggered', 'work', 'ablations', 'orientation']:
        selected = ['--datasets', 'skitter', '--graph', graph, '--batch', '2', '--batches', '4', '--hold', '.15', '--check', '1', '--threads', '1', '--case', '_r1']
        call([sys.executable, ROOT / 'reproduce.py', 'run', suite] + selected + ['--output', output / suite])
    call([sys.executable, ROOT / 'reproduce.py', 'run', 'friend', '--graph', stress, '--threads', '1', '--limit', '120', '--output', output / 'friend'])
    temporal = inputs / 'temporal.txt'
    temporal.write_text(''.join(f'{u} {v} {i}\n' for i, (u, v) in enumerate(edges * 3)))
    event_text = ''.join(f'{u + 1} {v + 1} 1 {i}\n' for i, (u, v) in enumerate(edges))
    for i in range(40):
        u, v = edges[i // 2]
        event_text += f'{u + 1} {v + 1} {1 if i % 2 else -1} {len(edges) + i}\n'
    events = inputs / 'events.tar.bz2'
    with tarfile.open(events, 'w:bz2') as packed:
        content = event_text.encode()
        item = tarfile.TarInfo('fixture/out.fixture')
        item.size = len(content)
        packed.addfile(item, io.BytesIO(content))
    for kind, source, name, mode, hold, suffix in [('temporal', temporal, 'wiki', 'stream', str(len(edges)), 'tbin'), ('events', events, 'frwiki', 'events', '0', 'ebin')]:
        target = inputs / ('fixture.' + suffix)
        call([sys.executable, ROOT / 'tools/data.py', 'convert', source, target, '--kind', kind])
        call([sys.executable, ROOT / 'reproduce.py', 'run', 'main', '--datasets', name, '--graph', target, '--mode', mode, '--hold', hold, '--batch', '2', '--batches', '4', '--check', '1', '--threads', '1', '--output', output / kind])
    call([sys.executable, ROOT / 'reproduce.py', 'run', 'main'] + common + ['--hypergraph', '--output', output / 'hypergraph'])
    if args.baselines:
        snapshot = next((output / 'main').glob('*/snapshot_b4.txt'))
        call([sys.executable, ROOT / 'tools/compare.py', 'static', snapshot, '--output', output / 'static', '--threads', '1', '--repeats', '1', '--limit', '60'])
        prefix = next((output / 'hypergraph').glob('*/hypergraph.updates')).with_suffix('')
        call([sys.executable, ROOT / 'tools/compare.py', 'dynamic', prefix, '--output', output / 'dynamic', '--methods', 'cap', 'udshp', '--epsilon', '.1', '--limit', '60'])
        cq = inputs / 'cq-tiny'
        Path(str(cq) + '.updates').write_text('4\n+ 0 4 0 1 2 3\n- 0 4 0 1 2 3\n+ 1 4 0 1 2 3\n')
        Path(str(cq) + '.queries.tsv').write_text('0 init\n1 b1\n2 b2\n')
        call([sys.executable, ROOT / 'tools/compare.py', 'dynamic', cq, '--output', output / 'cq-functional', '--methods', 'cq', '--epsilon', '.5', '--limit', '60'])
        for epsilon, label in [('.1', 'e01'), ('.01', 'e001')]:
            call([sys.executable, ROOT / 'tools/compare.py', 'dynamic', cq, '--output', output / ('cq-functional-' + label), '--methods', 'cq', '--epsilon', epsilon, '--limit', '30'])
    if args.real_data:
        call([sys.executable, ROOT / 'tools/data.py', 'fetch', 'enron'])
        if digest(ROOT / 'datasets/email-Enron.bin') != '9202d8082de4fdefba3954eee4ae0d9233c3b3636575c5bd977ad570ff64cf81':
            raise ValueError('Public Enron conversion differs from the recorded experiment input')
        for suite in ['main', 'strict']:
            call([sys.executable, ROOT / 'reproduce.py', 'run', suite, '--datasets', 'enron', '--case', '_r1', '--batch', '100', '--batches', '3', '--check', '1', '--threads', '1', '--output', output / ('enron-' + suite)])
    call([sys.executable, ROOT / 'tools/summarize.py', output])
    save(output / 'validation.json', {'scope': 'release correctness and pipeline validation; not a repetition of the published timing tables', 'commands': invocations, 'real_full_enron': args.real_data, 'baseline_functions': args.baselines, 'cq_functional_epsilons': ['.5', '.1', '.01'], 'input_sha256': {'tiny': digest(graph)}})
    print('Release validation passed')

if __name__ == '__main__':
    main()
