import argparse
import os
from pathlib import Path
import shutil
import subprocess
from common import ROOT, digest, save

TARGETS = {
    'guarded': 'src/implicert.cpp', 'trace': 'src/trace.cpp',
    'export_checkpoint': 'src/export_checkpoint.cpp', 'verify_k4': 'src/verify_k4.cpp',
    'v7b': 'versions/v7b/main.cpp', 'v7c': 'versions/v7c/main.cpp',
    'triggered': 'versions/triggered/main.cpp', 'control': 'versions/control/main.cpp',
    'profile': 'versions/profile/main.cpp', 'policy': 'versions/policy/main.cpp',
    'friend': 'versions/friend/main.cpp', 'prep_big': 'tools/prep_big.cpp',
    'prep_events': 'tools/prep_events.cpp', 'maxclique': 'tools/maxclique.cpp',
    'memcap': 'tools/memcap.cpp', 'guarded_audit': 'versions/guarded_audit/main.cpp',
    'audit_micro': 'versions/guarded_audit/audit_micro.cpp',
    'test_absorber': 'tools/test_absorber.cpp', 'test_contract': 'tools/test_contract.cpp'
}

def compile_msvc(name, sources, extra=None, c=False):
    if os.name != 'nt':
        raise RuntimeError('The measured implementation requires x64 Windows and MSVC')
    output = ROOT / 'build'
    output.mkdir(exist_ok=True)
    setup = []
    if not os.environ.get('VSCMD_VER'):
        finder = Path(os.environ.get('ProgramFiles(x86)', 'C:/Program Files (x86)')) / 'Microsoft Visual Studio/Installer/vswhere.exe'
        installation = subprocess.check_output([str(finder), '-latest', '-products', '*', '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64', '-property', 'installationPath'], text=True).strip()
        if not installation:
            raise RuntimeError('Install Visual Studio C++ Build Tools with LLVM OpenMP')
        setup = [f'call "{installation}\\VC\\Auxiliary\\Build\\vcvars64.bat" >nul', 'if errorlevel 1 exit /b 1']
    command = ['cl', '/nologo', '/O2', '/TC' if c else '/std:c++20', '/openmp:llvm', '/EHsc']
    if name == 'verify_k4':
        command.append('/fp:strict')
    command += [f'/I{ROOT / "src"}'] + list(extra or []) + list(map(str, sources)) + [f'/Fe:{name}.exe', '/link', 'Psapi.lib', '/STACK:268435456']
    script = output / f'{name}-build.cmd'
    script.write_text('\n'.join(['@echo off', 'setlocal'] + setup + [subprocess.list2cmdline(command), 'exit /b %errorlevel%']) + '\n')
    with (output / f'{name}-build.log').open('w') as log:
        subprocess.run(['cmd', '/c', str(script)], cwd=output, stdout=log, stderr=subprocess.STDOUT, check=True)
    executable = output / f'{name}.exe'
    dependencies = set(map(Path, sources)) | set((ROOT / 'src').glob('*.hpp'))
    for source in sources:
        dependencies.update(Path(source).parent.glob('*.hpp'))
    save(output / f'{name}-build.json', {'compiler': 'MSVC x64', 'command': command, 'sources': {str(p.relative_to(ROOT)): digest(p) for p in sorted(dependencies)}, 'binary_sha256': digest(executable)})
    print(name, 'built', flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('targets', nargs='+', choices=['all'] + list(TARGETS))
    args = parser.parse_args()
    for name in TARGETS if 'all' in args.targets else args.targets:
        compile_msvc(name, [ROOT / TARGETS[name]])
    if (ROOT / 'build/guarded.exe').exists():
        shutil.copyfile(ROOT / 'build/guarded.exe', ROOT / 'build/implicert.exe')

if __name__ == '__main__':
    main()
