"""Track B pilot (exploratory, paired seeds): O1 Onsager SCA + extra echo terms -lam_x * r(t) * chat_k(t) * s(t-k), k=2..K.
chat_k(t): measured on O1 dynamics (multilag_chat.json), 21-step moving average; r(t) = the same end ramp as lambda_1."""
import json, sys, time
import numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m
import multilag as ML

J, sumw = m.load()
S, t0, q = 960, 12.0, 8.0
chat = np.array(json.load(open('multilag_chat.json'))['O1'])
ker = np.ones(21) / 21
sm = {k: np.convolve(np.pad(chat[:, k], 10, mode='edge'), ker, mode='valid') for k in range(2, ML.KMAX + 1)}
ramp = ML.lam_sched(1.0, True, S)


def hw_ms(fl): return (0.522 * fl + 40.2 * S + 1884) / 300e3


def tts(p, t): return float('inf') if p <= 0 else t * np.log(0.01) / np.log1p(-min(p, 1 - 1e-9))


def evaluate(cfgs, B, seed):
    res = []
    for lam1, K, lamx in cfgs:
        rng = np.random.default_rng(seed)          # paired: same starts and noise stream per config
        extra = None if K == 0 else {k: lamx * ramp * sm[k][:S] for k in range(2, K + 1)}
        t1 = time.time(); cut, flips, _ = ML.run(J, sumw, S, t0, q, lam1, True, B, rng, extra=extra)
        p = float((cut >= 33000).mean()); tr = hw_ms(flips.mean())
        res.append(dict(lam1=lam1, K=K, lamx=lamx, B=B, p=p, mean_cut=float(cut.mean()), flips=float(flips.mean()), t_run_ms=tr, tts_ms=tts(p, tr)))
        print(f"lam1 {lam1:.2f} K {K} lamx {lamx:.2f}: p={p:.3f} cut={cut.mean():.0f} flips={flips.mean():.0f} TTS={tts(p, tr):.3f} ms ({time.time()-t1:.0f}s)", flush=True)
    return res


grid = [(l1, 0, 0.0) for l1 in (0.9, 1.05)] + [(l1, K, lx) for l1 in (0.9, 1.05) for K in (2, 4, 6) for lx in (0.5, 1.0)]
print("PILOT (512 paired runs each)")
pilot = evaluate(grid, 512, 20261004)
best_extra = min((r for r in pilot if r['K'] > 0), key=lambda r: r['tts_ms'])
base = [r for r in pilot if r['K'] == 0]
best_base = min(base, key=lambda r: r['tts_ms'])
print("CONFIRM (2048 fresh paired runs): best base vs best multi-lag")
confirm = evaluate([(best_base['lam1'], 0, 0.0), (best_extra['lam1'], best_extra['K'], best_extra['lamx'])], 2048, 777)
json.dump(dict(pilot=pilot, confirm=confirm), open('multilag_test.json', 'w'), indent=1)
