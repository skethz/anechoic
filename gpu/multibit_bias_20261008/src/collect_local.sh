#!/usr/bin/env bash
# Package the bias study's results for the local mirror (run on gpu-host from the project root): summaries, batch files,
# checker/identity records, logs and per-trial CSVs; raw spin/trace files stay on gpu-host (sha256 in raw_manifest.sha256).
set -euo pipefail
cd /scratch/USER/anechoic_gpu_multibit_bias_20261008
OUT=tmp/collect; rm -rf $OUT; mkdir -p $OUT/results $OUT/remote_meta
find results -type f \( -name '*.spins.bin' -o -name '*.trace.bin' -o -name '*.trials.jsonl' \) | sort | xargs -r sha256sum > $OUT/results/raw_manifest.sha256
for ph in $(cd results && find . -mindepth 1 -maxdepth 2 -type d | sed 's#^\./##' | sort); do
  mkdir -p $OUT/results/$ph
  find results/$ph -maxdepth 1 -type f \( -name '*.summary.json' -o -name '*.batches.jsonl' -o -name '*.out' -o -name '*.jsonl' -o -name '*.json' -o -name '*.txt' -o -name '*.csv' -o -name '*.tsv' \) ! -name '*.trials.jsonl' -exec cp {} $OUT/results/$ph/ \;
  python3 - "$ph" "$OUT" <<'PY'
import sys, glob, json, gzip, os
ph, out = sys.argv[1], sys.argv[2]
fs = sorted(glob.glob(f'results/{ph}/*.trials.jsonl'))
if fs:
    keys = ['batch', 'chain', 'trial_id', 'cut', 'device_cut', 'sum_sh', 'sum_bs', 'energy', 'device_energy', 'field_mismatches', 'flips', 'success']
    with gzip.open(f"{out}/results/{ph.replace('/', '_')}_trials.csv.gz", 'wt') as g:
        g.write('run,' + ','.join(keys) + '\n')
        for f in fs:
            run = os.path.basename(f)[:-len('.trials.jsonl')]
            for l in open(f):
                d = json.loads(l)
                g.write(run + ',' + ','.join('' if k not in d else str(int(d[k])) if isinstance(d[k], bool) else str(d[k]) for k in keys) + '\n')
PY
done
cp results/tables.json results/tables.md results/GS_tables.json results/GS_tables.md $OUT/results/ 2>/dev/null || true
cp PROTOCOL*.md PROTOCOL*.sha256 $OUT/remote_meta/ 2>/dev/null || true
cp -r logs $OUT/remote_meta/logs
mkdir -p $OUT/remote_meta/build $OUT/remote_meta/data
cp build/BUILD.sha256 build/ptxas_sca_gpu_bias.log $OUT/remote_meta/build/
cp data/inst/SHA256SUMS $OUT/remote_meta/data/inst_SHA256SUMS; cp data/inst/*.json $OUT/remote_meta/data/ 2>/dev/null || true
cp data/gset/SHA256SUMS $OUT/remote_meta/data/gset_SHA256SUMS; cp data/gset/download_time_utc.txt $OUT/remote_meta/data/gset_download_time_utc.txt
[ -d data/hw_export ] && { mkdir -p $OUT/remote_meta/data/hw_export; cp data/hw_export/*.json data/hw_export/SHA256SUMS $OUT/remote_meta/data/hw_export/ 2>/dev/null || true; }
cp data/gset_study_results_SHA256SUMS $OUT/remote_meta/data/ 2>/dev/null || true
sha256sum data/v80_a3/* > $OUT/remote_meta/data/v80_a3_SHA256SUMS 2>/dev/null || true
cp environment.sh $OUT/remote_meta/
tar czf tmp/collect.tgz -C $OUT .
ls -la tmp/collect.tgz
