import json
from common import ROOT

DATASETS = json.loads((ROOT / 'benchmarks/datasets.json').read_text())

def job(dataset, **changes):
    source = DATASETS[dataset]
    result = dict(dataset=dataset, version='guarded', k=4, epsilon='0.01', batch=source['batch'], batches=20, mode=source['mode'], hold=source['hold'], check=source['check'], seed=1, threads=24, strategy='hybrid', q=64, env={}, limit=7200)
    result.update(changes)
    return result

def plan(suite):
    rows = []
    if suite == 'main':
        rows = [job(name) for name in DATASETS if name != 'friendster']
    elif suite == 'strict':
        for repeat in range(1, 4):
            for name in ['dblp', 'wiki', 'amazon', 'enron', 'google']:
                arms = ['maintained', 'warm', 'cold']
                for arm in arms[repeat - 1:] + arms[:repeat - 1]:
                    rows.append(job(name, version='trace', arm=arm, repeat=repeat, batch=100000 if name == 'wiki' else 1000, env={'KCDS_NOCHECK': '1'}, limit=14400))
    elif suite == 'triggered':
        for name in ['wiki', 'skitter', 'lj']:
            for arm in ['local', 'cold', 'warm']:
                rows.append(job(name, version='triggered', arm=arm, env={'KCDS_REPAIR_ARM': arm}, limit=14400))
    elif suite in ('ablations', 'work'):
        for repeat in range(1, 4):
            for name in ['wiki', 'so', 'skitter']:
                for arm, env in [('full', {}), ('noSAT', {'KCDS_SAT': '0'}), ('noABS', {'KCDS_QABS': '0'}), ('uniform', {'KCDS_UNIFORM': '1'}), ('noSplit', {'KCDS_NOSPLIT': '1'}), ('noFallback', {'KCDS_FALLBACK': '0'})]:
                    rows.append(job(name, arm=arm, repeat=repeat, version='profile' if suite == 'work' else 'guarded', env=dict(env) if suite == 'work' else dict(env, KCDS_NOCHECK='1')))
    elif suite == 'strategies':
        for name in ['dblp', 'skitter', 'wiki', 'so', 'lj']:
            for strategy in ['hybrid', 'newton', 'targeted']:
                rows.append(job(name, strategy=strategy))
    elif suite == 'adversarial':
        for name in ['dblp', 'skitter', 'lj']:
            for seed in [1, 2, 3]:
                for fallback in ['cold', 'warm']:
                    rows.append(job(name, mode='core', seed=seed, arm=fallback, env={'KCDS_WARMFB' if fallback == 'warm' else 'KCDS_COLDFB': '1'}))
        rows += [job(name, mode='plant') for name in ['dblp', 'skitter']]
    elif suite == 'sensitivity':
        for name in ['dblp', 'skitter', 'wiki', 'so']:
            rows += [job(name, k=k) for k in [3, 5, 6, 7]]
        for name in ['skitter', 'wiki']:
            rows += [job(name, epsilon=epsilon) for epsilon in ['0.005', '0.02', '0.05']]
            rows += [job(name, threads=threads) for threads in [1, 2, 4, 8, 16, 24]]
            rows += [job(name, batch=batch) for batch in ([100, 1000, 10000, 100000] if name == 'skitter' else [1000, 10000, 100000])]
        for name in ['dblp', 'skitter', 'lj']:
            rows += [job(name, seed=seed) for seed in [1, 2, 3]]
    elif suite == 'long':
        for name, count, check in [('wiki', 200, 50), ('so', 1000, 250), ('frwiki', 1000, 250)]:
            for rank in ['rerank', 'fixed']:
                rows.append(job(name, batches=count, check=check, arm=rank, env={} if rank == 'rerank' else {'KCDS_NORERANK': '1'}, limit=86400))
    elif suite == 'memory':
        for name in ['dblp', 'skitter', 'lj', 'orkut', 'wiki', 'so', 'friendster']:
            rows.append(job(name, env={'KCDS_NOCHECK': '1'}, arm='standard', limit=86400))
            if DATASETS[name]['mode'] == 'replay':
                rows.append(job(name, env={'KCDS_NOCHECK': '1', 'KCDS_COMPACT': '1'}, arm='compact', limit=86400))
    elif suite == 'orientation':
        for name in ['skitter', 'wiki']:
            for arm, env in [('rebuild', {}), ('reuse', {'KCDS_REUSE_ORIENTATION': '1'})]:
                rows.append(job(name, version='control', arm=arm, env=env))
    elif suite == 'policies':
        for arm in ['hybrid', 'targeted', 'newton-nosat']:
            rows.append(job('dblp', version='policy', mode='core', arm=arm, strategy=arm.removesuffix('-nosat'), fixed_edges=True, env={'KCDS_SAT': '0'} if arm.endswith('-nosat') else {}, check=1))
        for name in ['wiki', 'frwiki']:
            for epsilon in ['0.01', '0.005']:
                for box in [150, 300, 340]:
                    rows.append(job(name, version='policy', epsilon=epsilon, arm=f'box{box}', env={'KCDS_TEST_BOX': str(box), 'KCDS_TEST_PCG': '.0001'}, check=1))
        for q, pcg in [(32, '.0001'), (128, '.0001'), (64, '.001'), (64, '.00001')]:
            rows.append(job('wiki', version='policy', q=q, arm='pcg' + pcg, env={'KCDS_TEST_BOX': '300', 'KCDS_TEST_PCG': pcg}))
    elif suite == 'friend':
        rows = [job('friendster', version='friend', batches=10, check=1, hold='0', threads=8, env={'KCDS_COMPACT': '1', 'KCDS_NOCHECK': '1'}, limit=86400)]
    elif suite == 'guarded':
        rows = [job(name, version='guarded_audit', check=1) for name in ['dblp', 'skitter', 'lj', 'wiki', 'so', 'frwiki']]
        rows += [job('dblp', version='guarded_audit', mode=mode, check=1) for mode in ['core', 'plant']]
        rows += [job(name, version='guarded_audit', k=k, check=1) for name, k in [('amazon', 7), ('dblp', 5), ('wiki', 5)]]
        rows += [job(name, version='guarded_audit', epsilon='.005', check=1) for name in ['wiki', 'skitter']]
    elif suite == 'scaling':
        for name in ['lj', 'orkut']:
            rows += [job(name, sample=percentage, arm=f'v{percentage}') for percentage in [20, 40, 60, 80, 100]]
    else:
        raise ValueError('Unknown suite ' + suite)
    for index, row in enumerate(rows):
        row['name'] = f'{index + 1:03d}_{row["dataset"]}_k{row["k"]}_{row["mode"]}_{row.get("arm", row["strategy"])}_s{row["seed"]}_r{row.get("repeat", 1)}'
    return rows

SUITES = ['main', 'strict', 'triggered', 'ablations', 'work', 'strategies', 'adversarial', 'sensitivity', 'long', 'memory', 'orientation', 'policies', 'friend', 'scaling', 'guarded']
