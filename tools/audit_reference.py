import csv
import json
from pathlib import Path
import re
import statistics
from common import ROOT, digest, save

def main():
    manifest = json.loads((ROOT / 'benchmarks/source_manifest.json').read_text())
    for target, identity in manifest.items():
        if digest(ROOT / target) != identity['artifact_sha256']:
            raise ValueError('Packaged source changed: ' + target)
    folder = ROOT / 'results/reference'
    output = []
    with (folder / 'version_map.csv').open() as handle:
        for row in csv.DictReader(handle):
            log = folder / 'main' / Path(row['frozen_log']).name
            if digest(log) != row['frozen_log_sha256']:
                raise ValueError('Historical log changed: ' + str(log))
            batches = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in log.read_text().splitlines() if line.startswith('batch=')]
            values = [float(batch['upd'].rstrip('s')) + float(batch['repair'].rstrip('s')) for batch in batches]
            output.append({'dataset': row['dataset_key'], 'historical_build': row['table_maintenance_and_warm_cold_build'], 'observed_batches': len(values), 'mean_maintenance_s': statistics.mean(values), 'max_maintenance_s': max(values), 'archived_log_sha256': digest(log), 'new_timing_run': False})
    save(folder / 'main-summary.json', output)
    print(len(manifest), 'source identities and', len(output), 'historical log hashes verified')

if __name__ == '__main__':
    main()
