#!/usr/bin/env bash
# After the revision-2 chain (PROTOCOL_HW_MB Amendment 2): the frozen analyses on the r2 runs, the r1-r2 comparison and the
# negative-score exact-gate count; then pack the records (no PDI) as tmp/mbpack_board_mb2r2_e12_250mhz.tgz.
set -uo pipefail
cd /scratch/USER/sca_v80_20261003 && source scripts/environment.sh > /dev/null
R1=results/hwmb_board_mb2_e12_250mhz; R2=results/hwmb_board_mb2r2_e12_250mhz
python3 scripts/analyze_hw_mb.py $R2 results/hwv6_board_v64_e12_250mhz --predictions PROTOCOL_HW_MB_predictions.json > $R2/analysis.txt 2>&1
python3 scripts/analyze_hw_mb_a1.py ${R2}_a1 PROTOCOL_HW_MB_A1_predictions.json > ${R2}_a1/analysis.txt 2>&1
python3 scripts/compare_r1_r2.py $R1 $R2 rounds gset > $R2/compare.txt 2>&1
python3 scripts/compare_r1_r2.py ${R1}_a1 ${R2}_a1 gset > ${R2}_a1/compare.txt 2>&1
python3 - <<'PY' > ${R2}_a2/neg_gate.txt 2>&1
import json, glob
tot = ex = 0; cuts = set()
for f in sorted(glob.glob('results/hwmb_board_mb2r2_e12_250mhz_a2/verify/*.jsonl')):
    rows = [json.loads(l) for l in open(f)]
    n = len(rows); e = sum(r['exact_reference'] == 'PASS' for r in rows); tot += n; ex += e
    cuts |= {r['cut'] for r in rows}
    print(f, n, e, 'engines', len({r['engine'] for r in rows}), 'cuts', sorted({r['cut'] for r in rows})[:5],
          'scores_agree', all(r['scores_agree'] for r in rows))
print('negative-score exact gate', ex, '/', tot, 'PASS' if ex == tot and tot > 0 else 'FAIL')
PY
for f in $R2/analysis.txt ${R2}_a1/analysis.txt $R2/compare.txt ${R2}_a1/compare.txt ${R2}_a2/neg_gate.txt; do echo "== $f"; head -25 $f; done
pk=tmp/mbpack_board_mb2r2_e12_250mhz; rm -rf $pk; mkdir -p $pk/build
cp -r $R2 ${R2}_a1 ${R2}_a2 $pk/
(cd $pk/hwmb_board_mb2r2_e12_250mhz_a1 && tar czf gset_cohorts.tgz gset && sha256sum gset_cohorts.tgz > gset_cohorts.tgz.sha256 && rm -rf gset)
cp logs/bringup_mb2r2_board_mb2r2_e12_250mhz.log logs/pipeline_board_mb2r2_e12_250mhz.log logs/pipeline_board_mb2r2_e12_250mhz.env $pk/ 2>/dev/null
cp results/programming_board_mb2r2_e12_250mhz_*.log $pk/ 2>/dev/null; cp $(ls -t results/programming_board_v64_e12_250mhz_*.log | head -1) $pk/
for b in results/staged_image_backup_*; do mkdir -p $pk/$(basename $b); cp $b/SHA256SUMS $b/staged_image.json $pk/$(basename $b)/; done
for f in image_manifest.json timing_signoff.json kernel_clock.json timer_structure.json utilization.rpt route_status.rpt; do cp build/board_mb2r2_e12_250mhz/$f $pk/build/ 2>/dev/null; done
(cd tmp && tar czf mbpack_board_mb2r2_e12_250mhz.tgz mbpack_board_mb2r2_e12_250mhz) && ls -la tmp/mbpack_board_mb2r2_e12_250mhz.tgz
