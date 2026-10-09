"""Summarise a measure_power_mb.sh run (hwmon.txt + load_*.json) and compute energy-to-solution from measured TTS99.
(Copy of analyze_power_v6.py; TTS is looked up in mb_summary.json under rounds or gset.)
Usage: python3 analyze_power_mb.py <power_dir> <mb_summary.json of the same image> [--settle 10]
Phases: idle* = image loaded, engines idle; load_<key> = all engines busy back to back (host --power-run, rev 2).
The first `settle` seconds of every phase are dropped (sensor and thermal settling)."""
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path


def main():
    d = Path(sys.argv[1]); summ = json.loads(Path(sys.argv[2]).read_text())
    settle = float(sys.argv[sys.argv.index('--settle') + 1]) if '--settle' in sys.argv else 10.0
    rows = defaultdict(list); t_first = {}
    for line in (d / 'hwmon.txt').read_text().splitlines():
        ph, t, p_uw, vccint_ma, pex_ma = line.split()
        t = float(t); t_first.setdefault(ph, t)
        if t - t_first[ph] >= settle:
            rows[ph].append((float(p_uw) / 1e6, float(vccint_ma) / 1e3, float(pex_ma) / 1e3))
    phases = {}
    for ph, r in rows.items():
        P = [x[0] for x in r]; I = [x[1] for x in r]
        phases[ph] = dict(n=len(r), board_W_mean=st.mean(P), board_W_sd=st.pstdev(P), board_W_min=min(P), board_W_max=max(P),
                          vccint_A_mean=st.mean(I), vccint_W_mean=0.8 * st.mean(I))
    idle_keys = sorted(k for k in phases if k.startswith('idle'))
    order = [k for k in sorted(phases, key=lambda k: t_first[k])]
    out = dict(power_dir=str(d), settle_s=settle, phases=phases, loads={})
    for k in order:
        if not k.startswith('load_'):
            continue
        key = k[5:]; i = order.index(k)
        before = [x for x in order[:i] if x.startswith('idle')][-1]; after = [x for x in order[i + 1:] if x.startswith('idle')][0]
        idle = 0.5 * (phases[before]['board_W_mean'] + phases[after]['board_W_mean'])
        host = json.loads((d / f'load_{key}.json').read_text())
        P = phases[k]['board_W_mean']
        rec = dict(board_W=P, idle_W_bracket=idle, dynamic_W=P - idle, vccint_W=phases[k]['vccint_W_mean'],
                   vccint_dynamic_W=phases[k]['vccint_W_mean'] - 0.5 * (phases[before]['vccint_W_mean'] + phases[after]['vccint_W_mean']),
                   device_busy_fraction=host['device_busy_fraction'], host_over_device=host.get('host_over_device'),
                   trials=host['trials'], seconds=host['seconds'])
        rec['energy_per_trial_uJ_board'] = P * host['seconds'] / host['trials'] * 1e6
        tts = summ['rounds'].get(key, summ.get('gset', {}).get(key, {})).get('tts99_ms')
        if tts:
            rec['tts99_ms_primary'] = tts
            rec['ets99_mJ_board'] = P * tts            # W * ms = mJ
            rec['ets99_mJ_dynamic'] = (P - idle) * tts
        out['loads'][key] = rec
    print(json.dumps(out, indent=1))
    (d / 'power_summary.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
