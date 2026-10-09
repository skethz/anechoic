"""Independent recomputation of Table 8a (both grid variants) from the PER-TRIAL board records of PROTOCOL_HW_MB
Amendment 3 (gset_cohorts_jsonl.tgz), not from a3_summary.json's per-cohort TTS: per cohort, rounds of 12 engines,
P_round = fraction of rounds with at least one success, t_R = mean round device time over rounds after round 0,
TTS99 = t_R if P_round = 1 else t_R ln0.01/ln(1 - P_round). Variant membership of each cohort from a3_summary.json.
Also checks every trial's scores_agree flag and that the recomputed TTS equals a3_summary.json's value."""
import json
import math
import tarfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
A3 = ROOT / 'fpga/v80_sca/results/multibit_20261007/board/board_mb2n_e12_250mhz/hwmb_board_mb2n_e12_250mhz_a3'
RULE = {'plain': 'SCA', 'tec': 'TEC', 'tecT': 'Onsager-kT', 'ons': 'Onsager-online'}


def wclass(g):
    if 1 <= g <= 5 or 22 <= g <= 26 or 43 <= g <= 47: return 'Random, +1'
    if 6 <= g <= 10 or 27 <= g <= 31: return 'Random, +-1'
    if 11 <= g <= 13 or 32 <= g <= 34: return 'Toroidal, +-1'
    if 14 <= g <= 17 or 35 <= g <= 38 or 51 <= g <= 54: return 'Planar-like, +1'
    return 'Planar-like, +-1'


def tts_from_records(recs):
    rounds = defaultdict(list)
    for r in recs:
        assert r['scores_agree'] and r['padding_zero']
        rounds[r['round']].append(r)
    rs = sorted(rounds)
    assert all(len(rounds[i]) == 12 for i in rs)
    P = sum(any(x['success'] for x in rounds[i]) for i in rs) / len(rs)
    t = [max(x['round_device_ms'] for x in rounds[i]) for i in rs]
    tR = sum(t[1:]) / (len(t) - 1)
    return (tR if P >= 1 else tR * math.log(0.01) / math.log1p(-P)), P, tR


def main():
    summ = json.loads((A3 / 'a3_summary.json').read_text())['cohorts']
    tts = {}
    with tarfile.open(A3 / 'gset_cohorts_jsonl.tgz') as tf:
        for m in tf.getmembers():
            if not m.name.endswith('.jsonl'):
                continue
            key = Path(m.name).stem
            recs = [json.loads(x) for x in tf.extractfile(m).read().decode().splitlines() if x.strip()]
            v, P, tR = tts_from_records(recs)
            assert abs(v - summ[key]['tts99_ms']) <= 1e-9 * max(1.0, v), (key, v, summ[key]['tts99_ms'])
            tts[key] = v
    print(len(tts), 'cohorts recomputed from per-trial records; all equal a3_summary.json')
    table = {}
    for variant in ('original', 'extended_amendment1'):
        per = defaultdict(dict)
        for key, c in summ.items():
            if variant in c['variants']:
                per[c['instance']][c['rule']] = tts[key]
        assert len(per) == 51
        rows = {}
        for cls in ('Random, +1', 'Random, +-1', 'Toroidal, +-1', 'Planar-like, +1', 'Planar-like, +-1', 'All 51'):
            inst = [g for g in per if cls == 'All 51' or wclass(int(g[1:])) == cls]
            row = {}
            for a, b in (('Onsager-kT', 'SCA'), ('Onsager-kT', 'TEC'), ('Onsager-online', 'SCA'), ('Onsager-online', 'TEC')):
                r = [per[g][b] / per[g][a] for g in inst]
                row[f'{a} vs {b}'] = (round(math.exp(sum(map(math.log, r)) / len(r)), 2), sum(x > 1 for x in r), len(r))
            row['geomean_ms'] = {k: round(math.exp(sum(math.log(per[g][k]) for g in inst) / len(inst)), 3)
                                 for k in ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')}
            rows[cls] = row
        table[variant] = rows
        print('==', variant)
        for cls, row in rows.items():
            print(f'{cls:17s}', {k: v for k, v in row.items()})
    Path(__file__).with_name('verify_board_variants.json').write_text(json.dumps(table, indent=1) + '\n')


if __name__ == '__main__':
    main()
