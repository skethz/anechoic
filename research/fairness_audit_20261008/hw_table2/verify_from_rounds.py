"""Independent recomputation of the measured Table 2 numbers from the per-trial board records (rounds/*.jsonl), not from
hwm_summary.json: p = successes/trials; per round: success if any of the 12 trials succeeds, round time = max device time
over its 12 engines (round_device_ms), rounds after round 0 (tables load in round 0); primary = t_R if P_round = 1 else
t_R ln0.01/ln(1-P_round); secondary = t_R ln0.01/(12 ln(1-p))."""
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / 'fpga/v80_sca/results/hwv6_board_v64_e12_250mhz'
summ = json.loads((D / 'hwm_summary.json').read_text())['rounds']
out = {}
for f in sorted((D / 'rounds').glob('*.jsonl')):
    k = f.stem
    recs = [json.loads(x) for x in f.read_text().splitlines() if x.strip()]
    n = len(recs); succ = sum(r['success'] for r in recs); p = succ / n
    assert all(r['scores_agree'] for r in recs)
    rounds = defaultdict(list)
    for r in recs:
        rounds[r['round']].append(r)
    rs = sorted(rounds)
    P_round = sum(any(x['success'] for x in rounds[i]) for i in rs) / len(rs)
    t_all = [max(x['round_device_ms'] for x in rounds[i]) for i in rs]
    t_R = sum(t_all[1:]) / (len(t_all) - 1)
    prim = t_R if P_round >= 1 else t_R * math.log(0.01) / math.log1p(-P_round)
    sec = t_R * math.log(0.01) / (12 * math.log1p(-p))
    s = summ[k]
    out[k] = dict(p=p, P_round=P_round, t_R=t_R, primary=prim, secondary=sec, rounds=len(rs))
    flag = all(abs(a - b) <= 1e-6 * max(1, abs(b)) for a, b in ((p, s['p']), (P_round, s['P_round']), (prim, s['tts99_ms']), (sec, s['tts99_ind_ms'])))
    print(f"{k:3s} n={n} p={p:.4f} P_round={P_round:.4f} t_R={t_R:.5f} primary={prim:.5f} secondary={sec:.5f} "
          f"| summary primary={s['tts99_ms']:.5f} secondary={s['tts99_ind_ms']:.5f} t_round={s['t_round_ms']:.5f} {'MATCH' if flag else 'DIFF'}")
Path(__file__).with_name('verify_from_rounds.json').write_text(json.dumps(out, indent=1) + '\n')
