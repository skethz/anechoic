"""EXPLORATORY additions (PROTOCOL.md section 6), run only after the pre-registered held-out stages.
  e1: S-scaling of the selected Onsager-online schedule beyond the grid (S = 32000), with and without the selected
      finishing phase; 4,096 runs, seeds [20261007, 20, chunk].
  e2: 'step-repeat' tables: the selected S > 4096 schedule compressed into 4,096 table rows, each row held for S/4096
      consecutive steps (what a per-row repeat count in the engine would execute), vs the smooth table on the same
      starts and RNG seeds; 8,192 runs each, seeds [20261007, 21, chunk].
  e3: the overall-selected Onsager-online parameters at the existing engine's maximum table length (S = 4096, T1 as
      selected) and with T1 = 5 plus the reference's selected finish compressed into 4,096 rows (S = 3996 + 100);
      8,192 runs each, seeds [20261007, 22, chunk]. Added after the held-out stages because the pre-registered
      existing-engine Onsager-online pick (S = 2000) failed on held-out.
Usage: python3 explore.py {e1|e2|e3}"""
import json
import sys
import time

import numpy as np

import sca_fast as F
import study as S

BASE = S.BASE


def repeat_tables(cfg, rows=4096):
    T, q, a, b = F.tables(cfg)
    n = int(cfg['S'])
    if n <= rows:
        return T, q, a, b
    t = np.arange(n)
    row = (t * rows) // n                       # row index used at step t
    first = np.searchsorted(row, np.arange(rows))  # first step of each row
    src = first[row]
    T2, q2, a2, b2 = T.copy(), q.copy(), a.copy(), b.copy()
    T2[:n], q2[:n], b2[:n] = T[src], q[src], b[src]
    # kcorr = lam_t / (2 T_{t-1}) is a row value too; keep a_0 = 0 (no s_prev at t = 0)
    a2[:n] = a[np.maximum(src, 1)]
    a2[0] = 0.0
    return T2, q2, a2, b2


def run_with_tables(cfg, tabs, B, seed, fin=None):
    s0, rs = F.make_batch(S.J, B, seed)
    st = F.State(S.J, None, s0, rs)
    F.advance(S.J2, st, tabs, 0, int(cfg['S']))
    out = dict(cut=F.cuts(st, S.SUMW), flips=st.flips.copy(), best=F.best_cuts(st, S.SUMW))
    out['greedy'], _ = S.greedy_cuts(st)
    if fin is not None:
        c2 = dict(cfg, **fin)
        x = S.copy_state(st)
        F.advance(S.J2, x, F.tables(c2), int(cfg['S']), int(cfg['S']) + fin['fin_L'])
        out['cut_fin'] = F.cuts(x, S.SUMW); out['flips_fin'] = x.flips.copy()
        out['greedy_fin'], _ = S.greedy_cuts(x)
    return out


def cat(parts):
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def e1():
    sel = S.load('a2_confirm')['selection']['onsager_overall']['cfg']
    cp = S.load('c_pilot')['configs']
    fin = cp[S.key(sel)]['selected'] if S.key(sel) in cp else None
    out = {}
    for Sx in (32000,):
        cfg = dict(sel, S=Sx); t = time.time()
        parts = [run_with_tables(cfg, F.tables(cfg), 1024, [BASE, 20, ch], fin) for ch in range(4)]
        a = cat(parts)
        r0 = S.summarize(cfg, a['cut'], a['flips'], Sx, a['best'], a['greedy'])
        res = dict(no_finish=r0)
        if fin is not None:
            c2 = dict(cfg, **fin)
            res['finish'] = S.summarize(c2, a['cut_fin'], a['flips_fin'], Sx + fin['fin_L'], None, a['greedy_fin'], dict(opt=fin))
        out[str(Sx)] = res
        S.log('E1', S.key(cfg), f"k337 {r0['k33337']}/{r0['n']} kg {r0['kgreedy33337']} TTS12 V80 {r0['tts']['V80_250MHz']['tts12_ms']:.1f}",
              ('finish k337 %d TTS12 %.1f' % (res['finish']['k33337'], res['finish']['tts']['V80_250MHz']['tts12_ms'])) if fin else '',
              f'({time.time() - t:.0f}s)')
        S.dump('explore_e1', dict(hashes=S.file_hashes(), base=sel, finish=fin, results=out))


def e2():
    sel = S.load('a2_confirm')['selection']
    cands = [v['cfg'] for k, v in sel.items() if k.endswith('_overall') and int(v['cfg']['S']) > 4096]
    out = {}
    for cfg in cands[:2]:
        t = time.time(); res = {}
        for name, tabs in (('smooth', F.tables(cfg)), ('repeat4096', repeat_tables(cfg))):
            parts = [run_with_tables(cfg, tabs, 2048, [BASE, 21, ch]) for ch in range(4)]
            a = cat(parts)
            res[name] = S.summarize(cfg, a['cut'], a['flips'], int(cfg['S']), a['best'], a['greedy'])
            S.log('E2', name, S.key(cfg), f"mean {res[name]['mean_cut']:.2f} k300 {res[name]['k33300']} k320 {res[name]['k33320']} "
                  f"k337 {res[name]['k33337']}/{res[name]['n']} kg337 {res[name]['kgreedy33337']}")
        out[S.key(cfg)] = res
        S.log('E2 done', f'({time.time() - t:.0f}s)')
        S.dump('explore_e2', dict(hashes=S.file_hashes(), results=out))


def e3():
    sel = S.load('a2_confirm')['selection']['onsager_overall']['cfg']
    cp = S.load('c_pilot')['configs']
    fin_ref = cp[S.key(S.REFERENCE)]['selected']
    cases = [('S4096_T1sel', dict(sel, S=4096), None),
             ('S3996_T1=5_plus_fin', dict(sel, S=4096 - fin_ref['fin_L'], T1=5.0), fin_ref)]
    out = {}
    for name, cfg, fin in cases:
        t = time.time()
        parts = [run_with_tables(cfg, F.tables(cfg), 2048, [BASE, 22, ch], fin) for ch in range(4)]
        a = cat(parts)
        res = dict(no_finish=S.summarize(cfg, a['cut'], a['flips'], int(cfg['S']), a['best'], a['greedy']))
        if fin is not None:
            c2 = dict(cfg, **fin)
            res['finish'] = S.summarize(c2, a['cut_fin'], a['flips_fin'], int(cfg['S']) + fin['fin_L'], None, a['greedy_fin'],
                                        dict(opt=fin))
        out[name] = res
        for k, r in res.items():
            S.log('E3', name, k, S.key(r['cfg']), f"k337 {r['k33337']}/{r['n']} kg {r['kgreedy33337']} "
                  f"TTS12 V80 {r['tts']['V80_250MHz']['tts12_ms']:.1f} [{r['tts']['V80_250MHz']['tts12_ms_lo']:.1f}, "
                  f"{r['tts']['V80_250MHz']['tts12_ms_hi']:.1f}] ms ({time.time() - t:.0f}s)")
        S.dump('explore_e3', dict(hashes=S.file_hashes(), results=out))


if __name__ == '__main__':
    {'e1': e1, 'e2': e2, 'e3': e3}[sys.argv[1]]()
