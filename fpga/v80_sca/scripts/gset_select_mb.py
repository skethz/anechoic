#!/usr/bin/env python3
"""Contingency G-set configuration selection for the multi-bit board run (used only if research/gset_20261007/selected_configs.json
is not available when PROTOCOL_HW_MB.md is frozen). Reference model only (sca_ref.hpp run_trial(..., dense) via build/explore_mb),
fixed procedure, written before any multi-bit board data:
  1. pilot: grid below, seed 0x5EED, trial ids 0.., per instance and rule family;
  2. selection: minimum predicted 12-engine TTS99 with at least one round, t_R * max(1, r99(P_round)) (rule v2; see tts()),
     using the Wilson 95% LOWER bound of the pilot p, ties -> shorter S;
     cycle model c = max(22.6 S, 0.1424 F + 18.09 S) + 886 (fitted on the MB RTL simulations; F = flips per trial),
     round time t_R = (1.03 c + 600) / f with f = 250 MHz;
  3. validation: the selected configurations re-run with a fresh seed 0xBA11D (512 trials) -> unbiased p prediction for the board,
     whose cohort uses yet another seed (20261004, the protocol seed).
Usage: gset_select_mb.py jobs <outdir>          writes jobs.txt (pilot)
       gset_select_mb.py select <outdir>        reads pilot.jsonl, writes selected.json and validate_jobs.txt
       gset_select_mb.py report <outdir>        reads validate.jsonl, writes predictions.json"""
import json
import math
import sys
from pathlib import Path

ROOT = '/scratch/USER/sca_v80_20261003'
INST = {  # instance: (best known cut, target = ceil(0.99 best), grid)
    'G22': (13359, 13226, dict(t0=[5, 6, 8, 10], t1=[0.3, 0.5, 0.7], q=[5, 6, 7], S=[512, 1024, 2048, 4096], trials=48)),
    'G32': (1410, 1396, dict(t0=[1.5, 2, 3], t1=[0.3, 0.5, 0.7], q=[0.5, 1, 1.5], S=[1024, 2048, 4096], trials=32)),
}
RULES = {'G22': [('plain', 0, 0), ('ons', 0.5, 0), ('ons', 0.9, 0), ('ons', 1.3, 0), ('tecT', 0, 1.0), ('tecT', 0, 1.5), ('tecT', 0, 2.0)],
         'G32': [('plain', 0, 0), ('ons', 0.5, 0), ('ons', 0.9, 0), ('tecT', 0, 1.0), ('tecT', 0, 2.0)]}
F_MHZ, E = 250.0, 12


def wilson_lo(k, n, z=1.96):
    if n == 0: return 0.0
    ph = k / n; d = 1 + z * z / n; c = ph + z * z / (2 * n); h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return max(0.0, (c - h) / d)


def cycles(S, F):
    return max(22.6 * S, 0.1424 * F + 18.09 * S) + 886


def tts(p, S, F, floor=True):
    """12-engine TTS99 from the per-trial p. floor=True (selection rule v2): at least one round, t_R * max(1, r99(P_round)).
    Rule v1 used the continuous form without the floor, which drops below one round as P_round -> 1 (an artefact that made
    the Wilson lower bound look faster than the point estimate); it is kept only for the record (selected_v1.json)."""
    tR = (1.03 * cycles(S, F) + 600) / (F_MHZ * 1e3)   # ms
    if p <= 0: return math.inf
    Pr = 1 - (1 - p) ** E
    r = 1.0 if Pr >= 1 else math.log(0.01) / math.log1p(-Pr)
    return tR * (max(1.0, r) if floor else r)


def line(inst, trials, seed, off, t0, t1, S, q, lam, ramp, kappa, scale, name):
    return f"gset:{ROOT}/data/gset/{inst} {INST[inst][1]} {trials} {seed} {off} {t0} {t1} {S} {q} {lam} {ramp} 0 {kappa} {scale} {name}"


def main():
    mode, d = sys.argv[1], Path(sys.argv[2])
    if mode == 'jobs':
        L = []
        for inst, (_, _, g) in INST.items():
            for fam, lam, kap in RULES[inst]:
                name = f"{inst}_{fam}" + (f"{lam}" if fam == 'ons' else f"{kap}" if fam == 'tecT' else '')
                for t0 in g['t0']:
                    for t1 in g['t1']:
                        for q in g['q']:
                            for S in g['S']:
                                ramp = 0 if fam == 'plain' else 1
                                scale = '1' if fam == 'plain' else 'auto'
                                L.append(line(inst, g['trials'], '0x5EED', 0, t0, t1, S, q, lam, ramp, kap, scale, name))
        (d / 'jobs.txt').write_text('\n'.join(L) + '\n'); print(len(L), 'pilot jobs')
    elif mode == 'select':
        rows = [json.loads(l) for l in (d / 'pilot.jsonl').read_text().splitlines() if l.strip()]
        best = {}
        for r in rows:
            inst = r['name'].split('_')[0]; fam = r['name'].split('_')[1].rstrip('0123456789.')
            pl = wilson_lo(r['hits'], r['trials'])
            r['tts_lo_ms'] = tts(pl, r['S'], r['mean_flips']); r['tts_pt_ms'] = tts(r['p'], r['S'], r['mean_flips']); r['p_wilson_lo'] = pl
            r['tts_pt_continuous_ms'] = tts(r['p'], r['S'], r['mean_flips'], floor=False)
            k = (inst, fam)
            if k not in best or (r['tts_lo_ms'], r['S']) < (best[k]['tts_lo_ms'], best[k]['S']): best[k] = r
        sel = {f"{i}_{f}": r for (i, f), r in sorted(best.items())}
        (d / 'selected.json').write_text(json.dumps(sel, indent=1) + '\n')
        V = []
        for key, r in sel.items():
            inst = key.split('_')[0]
            for c in range(8):   # 8 chunks x 64 = 512 validation trials, fresh seed
                V.append(line(inst, 64, '0xBA11D', 64 * c, r['t0'], r['t1'], r['S'], r['q'], r['lambda'], r['ramp'], r['kappa'],
                              'auto' if r['scale'] != 1 else '1', key))
        (d / 'validate_jobs.txt').write_text('\n'.join(V) + '\n')
        for key, r in sel.items():
            print(f"{key:12s} t0 {r['t0']} t1 {r['t1']} q {r['q']} S {r['S']} lam {r['lambda']} kappa {r['kappa']} scale {r['scale']:.6g} "
                  f"pilot p {r['p']:.3f} ({r['hits']}/{r['trials']}) lo {r['p_wilson_lo']:.3f} TTS(lo) {r['tts_lo_ms']:.3f} ms (pt {r['tts_pt_ms']:.3f})")
    elif mode == 'report':
        sel = json.loads((d / 'selected.json').read_text())
        rows = [json.loads(l) for l in (d / 'validate.jsonl').read_text().splitlines() if l.strip()]
        out = {}
        for key, r in sel.items():
            vs = [v for v in rows if v['name'] == key]
            k = sum(v['hits'] for v in vs); n = sum(v['trials'] for v in vs); F = sum(v['mean_flips'] * v['trials'] for v in vs) / n
            p = k / n; lo = wilson_lo(k, n)
            z = 1.96; ph = p; dd = 1 + z * z / n; cc = ph + z * z / (2 * n); hh = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
            out[key] = dict(config=dict(t0=r['t0'], t1=r['t1'], q=r['q'], S=r['S'], lam=r['lambda'], ramp=r['ramp'], kappa=r['kappa'],
                                        corr_scale='auto' if r['scale'] != 1 else '1'),
                            target=INST[key.split('_')[0]][1], validation_hits=k, validation_trials=n, p=p,
                            p_wilson95=[lo, min(1.0, (cc + hh) / dd)], mean_flips=F, model_cycles=cycles(r['S'], F),
                            predicted_tts99_ms=tts(p, r['S'], F), predicted_round_ms=(1.03 * cycles(r['S'], F) + 600) / (F_MHZ * 1e3))
            print(f"{key:12s} p {p:.3f} ({k}/{n}) [{lo:.3f},{out[key]['p_wilson95'][1]:.3f}] flips {F:.0f} cycles {cycles(r['S'], F):.0f} "
                  f"TTS99 {out[key]['predicted_tts99_ms']:.4f} ms")
        (d / 'predictions.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
