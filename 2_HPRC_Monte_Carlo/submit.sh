#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
last=$(python - <<'PYTHON'
import csv, math
from pathlib import Path
with open('mc_accepted_runs.csv') as f:
    states=[r for r in csv.DictReader(f) if r['accepted'].lower() == 'true']
with open('mc_group2_samples.csv') as f:
    samples=list(csv.DictReader(f))
missing=[r['run_id'] for r in states if not Path(f"steady_state_fsed/mc_{int(r['run_id'])}_fsed_steady.dat").is_file()]
if missing: raise SystemExit(f'Missing steady states: {missing}')
if not states or not samples: raise SystemExit('No accepted states or samples')
print(math.ceil(len(states)*len(samples)/50)-1)
PYTHON
)
sbatch --array="0-${last}" main_job.sh
