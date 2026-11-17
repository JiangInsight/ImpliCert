import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(8 << 20), b''):
            result.update(block)
    return result.hexdigest()

def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)

def clean_env():
    return {k: v for k, v in os.environ.items() if not k.startswith(('KCDS_', 'OMP_', 'KMP_', 'TRACE_', 'MEMCAP_'))}

def execute(command, output, environment=None, limit=3600, identity=None, cwd=None):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    record = output / 'record.json'
    identity = dict(identity or {}, command=list(map(str, command)), environment=environment or {}, limit_s=limit)
    if record.exists():
        previous = json.loads(record.read_text())
        if previous['identity'] != identity:
            raise ValueError(f'Existing run has a different identity: {output}')
        if previous['status'] != 'running':
            return previous
        raise ValueError(f'Unresolved previous run: {output}; use a new output directory')
    result = {'identity': identity, 'status': 'running', 'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    save(record, result)
    env = clean_env()
    env.update(environment or {})
    then = time.perf_counter()
    log = output / 'run.log'
    with log.open('x', encoding='utf-8') as handle:
        try:
            with subprocess.Popen(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT, env=env, cwd=cwd) as process:
                try:
                    code = process.wait(timeout=limit)
                    result.update(status='ok' if code == 0 else 'error', exit_code=code)
                except subprocess.TimeoutExpired:
                    if os.name == 'nt':
                        subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'], capture_output=True)
                    process.kill()
                    process.wait()
                    result.update(status='timeout', exit_code=process.returncode)
        except OSError as error:
            result.update(status='error', error=str(error))
    result.update(wall_s=time.perf_counter() - then, log_sha256=digest(log))
    save(record, result)
    print(output.name, result['status'], round(result['wall_s'], 3), flush=True)
    return result
