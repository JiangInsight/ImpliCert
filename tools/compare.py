import argparse
import json
from fractions import Fraction
import math
from pathlib import Path
import re
from common import ROOT, digest, execute, save

def run(name, command, folder, args, hashes, threads=1):
    executable = ROOT / 'build' / (name + '.exe')
    env = {'OMP_NUM_THREADS': str(threads)}
    if args.affinity and threads == 1:
        env['MEMCAP_AFFINITY'] = args.affinity
    folder.mkdir(parents=True, exist_ok=True)
    return execute([ROOT / 'build/memcap.exe', args.memory_mb, executable] + command, folder, env, args.limit, {'baseline': name, 'input_sha256': hashes, 'binary_sha256': digest(executable), 'epsilon': args.epsilon, 'k': args.k}, cwd=folder)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['static', 'dynamic'])
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--methods', nargs='+')
    parser.add_argument('--epsilon', default='0.01')
    parser.add_argument('--k', type=int, default=4)
    parser.add_argument('--threads', type=int, default=24)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--limit', type=float, default=4200)
    parser.add_argument('--solver-limit', type=float, default=3600)
    parser.add_argument('--memory-mb', type=int, default=32768)
    parser.add_argument('--affinity')
    args = parser.parse_args()
    source, output = args.input.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    methods = args.methods or (['ccas', 'kcca', 'kclist', 'psctl'] if args.mode == 'static' else ['cap', 'udshp', 'cq'])
    epsilon = Fraction(args.epsilon)
    if epsilon <= 0:
        raise ValueError('epsilon must be positive')
    results = []
    if args.mode == 'static':
        hashes = {source.name: digest(source)}
        omega = None
        if 'ccas' in methods:
            maximum = run('maxclique', [source, 60], output / 'maxclique', args, hashes)
            if maximum['status'] != 'ok':
                raise RuntimeError('Maximum-clique preprocessing failed')
            text = (output / 'maxclique/run.log').read_text(errors='replace')
            omega = int(re.search(r'omega=(\d+)', text)[1])
        for name in methods:
            for repeat in range(1, (args.repeats if name == 'kclist' else 1) + 1):
                folder = output / f'{name}_r{repeat}'
                folder.mkdir(exist_ok=True)
                if name == 'ccas':
                    command = [source, args.k, 0, float(epsilon / (1 + epsilon)), omega]
                elif name == 'kcca':
                    command = [source, args.k, 0, folder / 'solution', 1, 1, float(1 + epsilon)]
                elif name == 'kclist':
                    edges = folder / 'edges.el'
                    with source.open() as handle, edges.open('w') as target:
                        next(handle)
                        for line in handle:
                            target.write(line)
                    command = [args.threads, 1000000, args.k, edges.name, 'solution', float(1 + epsilon), 3, args.solver_limit]
                elif name == 'psctl':
                    if args.k != 4 or epsilon != Fraction(1, 100):
                        raise ValueError('The matched PSCTL reconstruction uses k=4, epsilon=0.01, at most 300 passes')
                    command = [source]
                else:
                    raise ValueError('Unknown static method ' + name)
                result = run(name, command, folder, args, hashes, args.threads if name == 'kclist' else 1)
                text = (folder / 'run.log').read_text(errors='replace')
                if name == 'kclist':
                    result['quality_target_reached'] = ' certified' in text
                    result['timing_definition'] = 'cumulative iteration time plus the final report; input construction excluded'
                elif name == 'psctl':
                    passes = re.findall(r'^PASS (\d+) (\S+) (\S+) (\S+)', text, re.M)
                    result['quality_target_reached'] = bool(passes and float(passes[-1][2]) <= 1.01 * float(passes[-1][1]))
                save(folder / 'record.json', result)
                results.append(result)
    else:
        updates, queries = Path(str(source) + '.updates'), Path(str(source) + '.queries.tsv')
        hashes = {p.name: digest(p) for p in [updates, queries]}
        for name in methods:
            folder = output / name
            prefix = folder / 'result'
            if name == 'cap':
                command = [updates, prefix, epsilon.numerator, epsilon.denominator, args.k, math.ceil(2 * args.k / epsilon), 'CAP', 0, queries, '--query-mode', 'full', '--evidence-mode', 'delta', '--full-checkpoints', 64]
            elif name == 'udshp':
                command = [updates, prefix, '--epsilon', float(epsilon), '--queries', queries, '--check', 0, '--dump-assignments', 0, '--checkpoint', 1000, '--budget-seconds', 0]
            elif name == 'cq':
                command = [updates, queries, prefix, float(epsilon), args.k, 0, 0]
            else:
                raise ValueError('Unknown dynamic method ' + name)
            result = run(name, command, folder, args, hashes)
            if result['status'] == 'ok':
                stats_path = folder / ('result.summary.json' if name == 'udshp' else 'result.stats.json')
                stats = json.loads(stats_path.read_text())
                expected = sum(1 for line in queries.read_text().splitlines() if line.strip())
                if stats.get('status') != 'complete' or int(stats.get('queries', -1)) != expected:
                    result.update(status='invalid', validation_error='Incomplete baseline query schedule')
                result['query_count'] = int(stats.get('queries', -1))
                save(folder / 'record.json', result)
            results.append(result)
    if any(row['status'] != 'ok' for row in results):
        raise SystemExit(1)

if __name__ == '__main__':
    main()
