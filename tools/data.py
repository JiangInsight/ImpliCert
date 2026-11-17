import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request
from common import ROOT, digest, save

def fetch(url, target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        temporary = target.with_suffix(target.suffix + '.part')
        with urllib.request.urlopen(url, timeout=90) as response, temporary.open('wb') as output:
            shutil.copyfileobj(response, output, 8 << 20)
        temporary.replace(target)
    return target

def convert(source, target, kind, native=False):
    source, target = Path(source), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if native and kind == 'static':
        with (gzip.open(source, 'rb') if source.suffix == '.gz' else source.open('rb')) as handle:
            with subprocess.Popen([str(ROOT / 'build/prep_big.exe'), str(target)], stdin=subprocess.PIPE) as process:
                shutil.copyfileobj(handle, process.stdin, 8 << 20)
                process.stdin.close()
                if process.wait() != 0:
                    raise RuntimeError('Native graph conversion failed')
    elif kind == 'events':
        with tarfile.open(source) as archive:
            members = [m for m in archive.getmembers() if m.isfile() and Path(m.name).name.startswith('out.')]
            if len(members) != 1:
                raise ValueError('Expected one KONECT out.* member')
            temporary = target.with_suffix('.input')
            with archive.extractfile(members[0]) as handle, temporary.open('wb') as output:
                shutil.copyfileobj(handle, output)
        subprocess.run([str(ROOT / 'build/prep_events.exe'), str(temporary), str(target)], check=True)
        temporary.unlink()
    else:
        import numpy as np
        import pandas as pd
        columns = ['u', 'v', 't'] if kind == 'temporal' else ['u', 'v']
        table = pd.read_csv(source, sep=r'\s+', header=None, names=columns, usecols=range(len(columns)), dtype=np.int64, comment='#')
        table = table[table.u != table.v]
        if kind == 'temporal':
            table = table.sort_values('t', kind='stable')
        labels, inverse = np.unique(np.concatenate([table.u.values, table.v.values]), return_inverse=True)
        size = len(table)
        u, v = inverse[:size].astype('<u4'), inverse[size:].astype('<u4')
        low, high = np.minimum(u, v), np.maximum(u, v)
        with target.open('wb') as output:
            if kind == 'temporal':
                np.array([len(labels), size], dtype='<u8').tofile(output)
                low.tofile(output)
                high.tofile(output)
                table.t.values.astype('<i8').tofile(output)
            else:
                pairs = np.unique(low.astype(np.uint64) * len(labels) + high)
                low, high = pairs // len(labels), pairs % len(labels)
                starts, ends = np.concatenate([low, high]), np.concatenate([high, low])
                order = np.lexsort((ends, starts))
                offsets = np.r_[0, np.cumsum(np.bincount(starts.astype(np.int64), minlength=len(labels)))].astype('<u8')
                np.array([len(labels), len(pairs)], dtype='<u8').tofile(output)
                offsets.tofile(output)
                ends[order].astype('<u4').tofile(output)
    save(target.with_suffix(target.suffix + '.json'), {'source': source.name, 'source_sha256': digest(source), 'format': kind, 'sha256': digest(target), 'bytes': target.stat().st_size})

def sample(source, target, fraction, seed):
    import numpy as np
    if not 0 < fraction <= 1:
        raise ValueError('Fraction must be in (0,1]')
    with Path(source).open('rb') as handle:
        n, m = map(int, np.fromfile(handle, '<u8', count=2))
        offset = np.fromfile(handle, '<u8', count=n + 1)
        adjacency = np.fromfile(handle, '<u4', count=2 * m)
    selected = np.flatnonzero(np.random.default_rng(seed).permutation(n) < int(round(n * fraction)))
    mapping = np.full(n, -1, dtype=np.int64)
    mapping[selected] = np.arange(len(selected))
    parts, lengths = [], []
    for vertex in selected:
        row = mapping[adjacency[offset[vertex]:offset[vertex + 1]]]
        row = row[row >= 0].astype('<u4')
        lengths.append(len(row))
        parts.append(row)
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('wb') as handle:
        np.array([len(selected), sum(lengths) // 2], dtype='<u8').tofile(handle)
        np.r_[0, np.cumsum(lengths)].astype('<u8').tofile(handle)
        np.concatenate(parts).tofile(handle)
    save(target.with_suffix(target.suffix + '.json'), {'source_sha256': digest(source), 'seed': seed, 'vertex_fraction': fraction, 'sha256': digest(target)})

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='action', required=True)
    obtain = sub.add_parser('fetch')
    obtain.add_argument('names', nargs='+')
    obtain.add_argument('--directory', type=Path, default=ROOT / 'datasets')
    obtain.add_argument('--source', type=Path)
    obtain.add_argument('--native', action='store_true')
    convert_parser = sub.add_parser('convert')
    convert_parser.add_argument('source', type=Path)
    convert_parser.add_argument('target', type=Path)
    convert_parser.add_argument('--kind', choices=['static', 'temporal', 'events'], default='static')
    convert_parser.add_argument('--native', action='store_true')
    subset = sub.add_parser('sample')
    subset.add_argument('source', type=Path)
    subset.add_argument('target', type=Path)
    subset.add_argument('--fraction', type=float, required=True)
    subset.add_argument('--seed', type=int, default=47)
    args = parser.parse_args()
    if args.action == 'convert':
        convert(args.source, args.target, args.kind, args.native)
    elif args.action == 'sample':
        sample(args.source, args.target, args.fraction, args.seed)
    else:
        manifest = json.loads((ROOT / 'benchmarks/datasets.json').read_text())
        if args.source and len(args.names) != 1:
            raise ValueError('--source requires exactly one dataset')
        for name in manifest if args.names == ['all'] else args.names:
            item = manifest[name]
            raw = args.source or fetch(item['url'], args.directory / 'raw' / item['archive'])
            target = args.directory / item['file']
            if not target.exists():
                convert(raw, target, item['kind'], args.native or name == 'friendster')
            expected = item.get('recorded_csr_sha256')
            if expected and digest(target) != expected:
                raise ValueError(f'{name}: conversion does not match the recorded experiment input')
            print(name, digest(target), flush=True)

if __name__ == '__main__':
    main()
