"""Table 8a from the board (PROTOCOL_HW_MB Amendment 3, a3_summary.json) for BOTH G-set grid variants.
No new data: each instance/rule has one cohort per variant ('original' and/or 'extended_amendment1' in its 'variants').
Speedup = geometric mean over instances of TTS99(reference)/TTS99(corrected rule), primary 12-engine estimator, with the
number of instances on which the corrected rule is faster (strictly lower TTS99)."""
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
A3 = ROOT / 'fpga/v80_sca/results/multibit_20261007/board/board_mb2n_e12_250mhz/hwmb_board_mb2n_e12_250mhz_a3/a3_summary.json'
RULES = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')


def wclass(g):
    if 1 <= g <= 5 or 22 <= g <= 26 or 43 <= g <= 47: return 'Random, +1'
    if 6 <= g <= 10 or 27 <= g <= 31: return 'Random, +-1'
    if 11 <= g <= 13 or 32 <= g <= 34: return 'Toroidal, +-1'
    if 14 <= g <= 17 or 35 <= g <= 38 or 51 <= g <= 54: return 'Planar-like, +1'
    if 18 <= g <= 21 or 39 <= g <= 42: return 'Planar-like, +-1'
    raise ValueError(g)


CLASSES = ('Random, +1', 'Random, +-1', 'Toroidal, +-1', 'Planar-like, +1', 'Planar-like, +-1')


def load():
    d = json.loads(A3.read_text())
    tts = {'original': defaultdict(dict), 'extended_amendment1': defaultdict(dict)}
    meta = {}
    for key, c in d['cohorts'].items():
        for v in c['variants']:
            assert c['rule'] not in tts[v][c['instance']], (key, v)
            tts[v][c['instance']][c['rule']] = (c['tts99_ms'], key, c['p'], c['S'])
        meta[c['instance']] = c['N']
    return tts


def table(tv):
    pairs = (('Onsager-kT', 'SCA'), ('Onsager-kT', 'TEC'), ('Onsager-online', 'SCA'), ('Onsager-online', 'TEC'), ('TEC', 'SCA'))
    out = {}
    for cls in CLASSES + ('All 51',):
        inst = [g for g in tv if cls == 'All 51' or wclass(int(g[1:])) == cls]
        row = {}
        for a, b in pairs:
            r = [tv[g][b][0] / tv[g][a][0] for g in inst]
            row[f'{a} vs {b}'] = dict(geomean=math.exp(sum(map(math.log, r)) / len(r)), wins=sum(x > 1 for x in r), n=len(r))
        for a in RULES:
            row[f'geomean {a} ms'] = math.exp(sum(math.log(tv[g][a][0]) for g in inst) / len(inst))
        out[cls] = row
    return out


def main():
    tts = load()
    res = {}
    for v in ('original', 'extended_amendment1'):
        tv = tts[v]
        assert len(tv) == 51 and all(len(tv[g]) == 4 for g in tv), v
        res[v] = table(tv)
        print(f'=== {v}')
        for cls, row in res[v].items():
            print(f"{cls:18s} " + '  '.join(f"{k.replace('Onsager-','')}: {x['geomean']:.2f} ({x['wins']}/{x['n']})"
                                             for k, x in row.items() if 'vs' in k) +
                  '  | gm ms: ' + ' '.join(f"{a}={row['geomean ' + a + ' ms']:.3f}" for a in RULES))
    # which cohorts differ between variants
    diff = {r: sum(tts['original'][g][r][1] != tts['extended_amendment1'][g][r][1] for g in tts['original']) for r in RULES}
    print('cohorts that differ between variants (per rule, of 51):', diff)
    res['cohorts_differing'] = diff
    Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).with_name('board_variants.json')).write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    main()
