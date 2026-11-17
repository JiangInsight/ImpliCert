import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
from common import ROOT, digest, save
from data import fetch
from build import compile_msvc

SOURCES = {
    'ccas': ('https://github.com/forxenn/ccas.git', 'c1a8032ec4f415958e78e1d9c3a49207dbda9cf4'),
    'kcca': ('https://github.com/forxenn/KCCA.git', 'de420e1badc878989443bec75a95eca8d6e277b1'),
    'kclist': ('https://github.com/btsun/kClistpp.git', '2704e72c8fa26a4be1adaa4ee8128456811ade99'),
    'cap': ('https://github.com/chenjb1997/cap-artifact.git', '6bd103bf86e7c40b32f122efc6b8bb599b868e70')
}

def checkout(name):
    target = ROOT / 'external' / name
    url, revision = SOURCES[name]
    if not target.exists():
        subprocess.run(['git', 'clone', '--no-checkout', url, str(target)], check=True)
        subprocess.run(['git', '-C', str(target), 'checkout', '--detach', revision], check=True)
    current = subprocess.check_output(['git', '-C', str(target), 'rev-parse', 'HEAD'], text=True).strip()
    changed = subprocess.check_output(['git', '-C', str(target), 'status', '--porcelain'], text=True).strip()
    if current != revision or changed:
        raise RuntimeError(f'{target}: expected clean pinned revision {revision}; preserve local changes and use a fresh checkout')
    return target

def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Adapter anchor is not unique')
    return text.replace(old, new)

def stripped(text):
    pattern = r'R"([^ ()\\\t\r\n]{0,16})\([\s\S]*?\)\1"|"(?:\\[\s\S]|[^"\\])*"|\'(?:\\[\s\S]|[^\'\\])*\'|//[^\n]*|/\*[\s\S]*?\*/'
    def replace(match):
        token = match.group()
        return ('\n' * token.count('\n') if '\n' in token else ' ') if token.startswith(('//', '/*')) else token
    return re.sub(r'\n[ \t]*\n(?:[ \t]*\n)+', '\n\n', re.sub(pattern, replace, text)).strip() + '\n'

def adapted_kclist(source, output):
    text = stripped((source / 'kCDensestNoAtom.c').read_text())
    patch = json.loads((ROOT / 'benchmarks/kclist_adapter.json').read_text())
    if hashlib.sha256(text.encode()).hexdigest() != patch['stripped_sha256']:
        raise ValueError('Unexpected KClist++ source')
    lines = text.splitlines(keepends=True)
    for change in reversed(patch['changes']):
        if ''.join(lines[change['start']:change['end']]) != change['old']:
            raise ValueError('KClist++ adapter mismatch')
        lines[change['start']:change['end']] = [change['new']]
    text = ''.join(lines)
    if hashlib.sha256(text.encode()).hexdigest() != patch['result_sha256']:
        raise ValueError('KClist++ adapted source mismatch')
    target = output / 'kCDensestNoAtom.c'
    target.write_text(text)
    return target

def adapted_psctl(source, output):
    text = (source / 'sct.cpp').read_text()
    for callback in ('flag ? iter_kcca(k, r) : iter_psctl(k, r);', 'iter_pava(k, r, rr);'):
        text = replace_once(text, 'Psize -= sz;\n        ' + callback, callback + '\n        Psize -= sz;')
    start, end = text.index('void iter_pava('), text.index('void iter_psctl(')
    text = text[:start] + text[start:end].replace('r[pvertex[i].second] < r[hvertex[0].second]', 'pvertex[i] < hvertex[0]') + text[end:]
    text = text.replace('void init_c(){', 'void init_c(){\n    static bool ready = false;\n    if (ready) return;\n    ready = true;')
    text = text.replace('sz = sz_;', 'delete[] fIndex; delete[] vSet; sz = sz_;')
    for marker in ['void build_pava(', 'bool cmpT(']:
        before, after = text.split(marker, 1)
        end = before.rindex('}')
        text = before[:end] + 'delete[] Hnode; delete[] Pnode; delete[] id; delete[] deg;\n' + before[end:] + marker + after
    text = text.replace('return cnm[n][m];', 'if (n < 5000) return cnm[n][m];\n    long double result = 1;\n    for (unsigned j=1; j<=m; ++j) result = result * (n-j+1) / j;\n    return static_cast<double>(result);')
    target = output / 'psctl_sct.cpp'
    target.write_text(text)
    return target

def adapted_timers(source, output):
    text = (source / 'src/cap.cpp').read_text()
    old = ' auto batch_start=Clock::now(),last_progress=batch_start;\n'
    text = replace_once(text, old, ' std::ofstream query_times(prefix+".query_times.tsv");if(!query_times)throw std::runtime_error("Cannot open query timing output");query_times<<"step\\tlabel\\tupdate_ns\\tquery_ns\\n";\n' + old)
    text = replace_once(text, '   ++next_report;++queries_emitted;\n  }\n', '   ++next_report;++queries_emitted;\n   query_times<<q.step<<\'\\t\'<<q.label<<\'\\t\'<<update_ns<<\'\\t\'<<query_ns<<\'\\n\';\n  }\n')
    cap = output / 'cap.cpp'
    cap.write_text(text)
    text = (source / 'baselines/cq/run.cpp').read_text()
    old = 'std::ofstream progress(prefix+".progress.jsonl");\n'
    text = replace_once(text, old, old + '  std::ofstream query_times(prefix+".query_times.tsv");query_times.precision(17);query_times<<"step\\tlabel\\tupdate_s\\tquery_s\\n";\n')
    old = 'if(!output)throw std::runtime_error("gzip write failed");++qi;'
    text = replace_once(text, old, 'if(!output)throw std::runtime_error("gzip write failed");query_times<<step<<\'\\t\'<<asks[qi].second<<\'\\t\'<<update_s<<\'\\t\'<<query_s<<\'\\n\';++qi;')
    cq = output / 'cq.cpp'
    cq.write_text(text)
    return cap, cq

def gcc_build(cxx, name, sources, flags, env):
    output = ROOT / 'build'
    command = [str(cxx)] + flags + list(map(str, sources)) + ['-o', str(output / (name + '.exe')), '-static', '-Wl,--stack,268435456']
    with (output / (name + '-build.log')).open('w') as log:
        subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    save(output / (name + '-build.json'), {'command': command, 'compiler': subprocess.check_output([str(cxx), '--version'], env=env, text=True).splitlines()[0], 'sources': {str(p): digest(p) for p in sources}, 'binary_sha256': digest(output / (name + '.exe'))})
    print(name, 'built', flush=True)

def build_zlib(cxx, env):
    archive = fetch('https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz', ROOT / 'external/zlib-1.3.1-release.tar.gz')
    if digest(archive) != '9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23':
        raise ValueError('zlib archive checksum mismatch')
    source = ROOT / 'external/zlib-1.3.1'
    if not source.exists():
        with tarfile.open(archive) as packed:
            for item in packed.getmembers():
                if not (ROOT / 'external' / item.name).resolve().is_relative_to((ROOT / 'external').resolve()):
                    raise ValueError('Unexpected archive member')
            packed.extractall(ROOT / 'external', filter='data')
    output = ROOT / 'build/zlib'
    output.mkdir(exist_ok=True)
    gcc, ar = cxx.with_name('gcc.exe'), cxx.with_name('ar.exe')
    names = 'adler32 compress crc32 deflate gzclose gzlib gzread gzwrite infback inffast inflate inftrees trees uncompr zutil'.split()
    if not (output / 'libz.a').exists():
        for name in names:
            subprocess.run([str(gcc), '-O2', '-c', str(source / (name + '.c')), '-o', str(output / (name + '.o'))], check=True, env=env)
        subprocess.run([str(ar), 'rcs', str(output / 'libz.a')] + [str(output / (name + '.o')) for name in names], check=True, env=env)
        for header in ['zconf.h', 'zlib.h']:
            shutil.copyfile(source / header, output / header)
    return output

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('names', nargs='+', choices=['all', 'ccas', 'kcca', 'kclist', 'psctl', 'cap', 'udshp', 'cq'])
    parser.add_argument('--cxx', type=Path, default=shutil.which('g++'))
    parser.add_argument('--fetch-only', action='store_true')
    args = parser.parse_args()
    names = ['ccas', 'kcca', 'kclist', 'psctl', 'cap', 'udshp', 'cq'] if 'all' in args.names else args.names
    required = set('ccas' if name == 'psctl' else 'cap' if name in ('udshp', 'cq') else name for name in names)
    trees = {name: checkout(name) for name in sorted(required)}
    if args.fetch_only:
        return
    if any(name != 'kclist' for name in names) and (not args.cxx or not args.cxx.is_file()):
        raise ValueError('Supply --cxx PATH/TO/g++.exe (WinLibs GCC with C++20 and __int128)')
    cxx = args.cxx.resolve() if args.cxx else None
    env = os.environ.copy()
    if cxx:
        env['PATH'] = str(cxx.parent) + os.pathsep + env['PATH']
    output = ROOT / 'build/adapters'
    output.mkdir(parents=True, exist_ok=True)
    zlib = build_zlib(cxx, env) if any(name in names for name in ('cap', 'udshp', 'cq')) else None
    timer_sources = adapted_timers(trees['cap'], output) if any(name in names for name in ('cap', 'cq')) else None
    for name in names:
        source = trees['ccas' if name == 'psctl' else 'cap' if name in ('udshp', 'cq') else name]
        if name == 'ccas':
            gcc_build(cxx, name, [source / f for f in ['ccas.cpp', 'graph.cpp', 'pava.cpp', 'sct.cpp', 'heap.cpp']], ['-Ofast', '-mcmodel=medium', '-march=native', '-std=c++14', '-mavx', '-include', 'cassert'], env)
        elif name == 'kcca':
            gcc_build(cxx, name, [source / f for f in ['kcca.cpp', 'graph.cpp', 'sct.cpp', 'pava.cpp']], ['-O3'], env)
        elif name == 'kclist':
            compile_msvc(name, [adapted_kclist(source, output), source / 'BinaryHeap.c', source / 'GraphDupAdj.c'], [f'/I{source}', '/D_CRT_SECURE_NO_WARNINGS'], c=True)
        elif name == 'psctl':
            gcc_build(cxx, name, [ROOT / 'tools/psctl.cpp', adapted_psctl(source, output)] + [source / f for f in ['graph.cpp', 'pava.cpp', 'heap.cpp']], ['-std=c++17', '-O3', '-DNDEBUG', '-msse4.2', '-include', 'immintrin.h', '-include', 'cassert', '-fopenmp', '-I' + str(source)], env)
        else:
            includes = [source / 'src', source / 'baselines/cq', source / 'baselines/udshp/core', source / 'baselines/udshp/include', source / 'baselines/common', zlib]
            flags = ['-std=c++20', '-O3', '-DNDEBUG', '-DHB_PROFILE=0'] + ['-I' + str(p) for p in includes]
            if name == 'udshp':
                sources = [source / f for f in ['baselines/udshp/run.cpp', 'baselines/udshp/core/DynamicGraph-singleInstance.cpp', 'baselines/udshp/core/DynOpManager-singleInstance.cpp', 'baselines/udshp/include/EdgeManager.cpp']]
            else:
                sources = [timer_sources[0 if name == 'cap' else 1]]
            gcc_build(cxx, name, sources + [zlib / 'libz.a'], flags + ['-lpsapi'], env)

if __name__ == '__main__':
    main()
