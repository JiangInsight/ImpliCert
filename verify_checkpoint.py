from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys

def main(native,checkpoint,epsilon):
    if epsilon<=0:raise ValueError('positive epsilon required')
    proc=subprocess.run([native,checkpoint],capture_output=True,text=True)
    if proc.returncode:
        raise ValueError(proc.stderr.strip() or 'native verifier failed')
    lines=[s for s in proc.stdout.splitlines() if s.startswith('checked ')]
    if len(lines)!=1:raise ValueError('missing verifier output')
    x=dict(t.split('=') for t in lines[0].split()[1:])
    upper=Fraction(float(x['safe']))
    raw=Fraction(float(x['raw']))
    lower=Fraction(int(x['count']),int(x['size']))
    passed=upper<=(1+epsilon)*lower
    print(json.dumps(dict(status='ok' if passed else 'bound_open',cliques=int(x['cliques']),
         witness_cliques=int(x['count']),lower=str(lower),upper=str(upper),raw=str(raw),
         upward_inflation=float(max(Fraction(),upper-raw)),relative_inflation=float(max(Fraction(),upper-raw)/raw) if raw else 0,
         rigorous=passed,epsilon=str(epsilon))))
    return 0 if passed else 3
if __name__=='__main__':
    try:sys.exit(main(Path(sys.argv[1]),Path(sys.argv[2]),Fraction(sys.argv[3])))
    except (ValueError,OSError) as e:print(e,file=sys.stderr);sys.exit(1)
