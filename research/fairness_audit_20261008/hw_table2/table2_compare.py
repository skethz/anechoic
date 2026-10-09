"""Table 2 comparisons recomputed from the measured 12 x v6.4 board summary (no new data, no board action).
Paper labels: O1 = X5, O5 = X4 (re-selected Onsager-kT); baselines: P1-P4 and T1-T3 (the 'original' schedules of the
4 Oct protocol PROTOCOL_HW.md, from our joint sweep J) and B1-B4 (re-selected 6 Oct). Primary = measured 12-engine rounds;
secondary = per-trial p. Ratios = baseline TTS / O TTS (> 1: the correction is faster)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
d = json.loads((ROOT / 'fpga/v80_sca/results/hwv6_board_v64_e12_250mhz/hwm_summary.json').read_text())['rounds']
cfg = {}
for k in d:
    rec = json.loads(open(ROOT / f'fpga/v80_sca/results/hwv6_board_v64_e12_250mhz/rounds/{k}.jsonl').readline())
    cfg[k] = {x: rec[x] for x in ('steps', 't0', 'q', 'lambda', 'ramp', 'tec_jv', 'tecT_kappa')}
O1, O5 = 'X5', 'X4'
groups = {'original plain (P1-P4)': ['P1', 'P2', 'P3', 'P4'], 'original TEC (T1-T3)': ['T1', 'T2', 'T3'],
          're-selected plain (B1, B2)': ['B1', 'B2'], 're-selected TEC (B3, B4)': ['B3', 'B4']}
out = {'measured': {k: dict(cfg=cfg[k], p=v['p'], P_round=v['P_round'], t_round_ms=v['t_round_ms'], primary_ms=v['tts99_ms'],
                            secondary_ms=v['tts99_ind_ms']) for k, v in d.items()}, 'ratios': {}}
lines = []
for est, key, o in (('primary', 'tts99_ms', O1), ('secondary', 'tts99_ind_ms', O5)):
    for g, ks in groups.items():
        best = min(ks, key=lambda k: d[k][key])
        r = {k: d[k][key] / d[o][key] for k in ks}
        out['ratios'][f'{est}: {g}'] = dict(vs=o, best=best, best_ratio=r[best], per_config=r)
        lines.append(f"{est:9s} {o}({'O1' if o == 'X5' else 'O5'}) {d[o][key]:.4f} ms vs {g:28s}: best {best} {d[best][key]:.4f} ms -> "
                     f"{r[best]:.2f}x | " + ', '.join(f'{k} {d[k][key]:.4f} ({r[k]:.2f}x)' for k in ks))
# P3 = the paper's "original plain SCA" (Fig. 4 caption: q = 8, T0 = 30, S = 1560)
for est, key, o in (('primary', 'tts99_ms', O1), ('secondary', 'tts99_ind_ms', O5)):
    lines.append(f"{est}: P3 (q8,T30,S1560) {d['P3'][key]:.4f} ms -> {d['P3'][key] / d[o][key]:.2f}x")
t = '\n'.join(lines)
print(t)
(Path(__file__).with_name('table2_compare.json')).write_text(json.dumps(out, indent=1) + '\n')
(Path(__file__).with_name('table2_compare.txt')).write_text(t + '\n')
