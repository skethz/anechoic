"""PROTOCOL_PUB_V2.md: TEC validation of OUR kernel (pub_methods._tec_glauber_kernel, as used for Figure 2) against the
TEC paper's Fig. 3(b), with the targets and criteria the COP study pre-declared (PROTOCOL_ADDENDUM4.md V3):
K2000, J_v in {0, 30, 6, -6}, k_BT 100 -> 0.1 geometric over 3,000 cycles, 64 runs each; measure = first cycle at which
the 64-run mean cut reaches 31,670 (95% of 33,337). Paper (read from Fig. 3b): J_v = 0: 1,392; J_v = 30: 740 (+-10);
pass: each within +-15% and the ratio within +-20%; also Fig. 3(a) ordering (J_v = 6 faster, -6 slower than 0).
Seeds SeedSequence([20261008, 611, J_v index]) for the initial states and per-run seeds. Usage: python3 validate_tec.py [workers]"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pub_methods as P  # noqa: E402
from common import M  # noqa: E402

HERE = Path(__file__).resolve().parent
JVS = (0.0, 30.0, 6.0, -6.0)
S, B, CH = 3000, 64, 8


def job(arg):
    ji, c0 = arg
    J, sumw = M.m.load()
    s0, seeds = P.draw_sequential_inputs(np.random.SeedSequence([20261008, 611, ji]), B, len(J))
    s0, seeds = s0[c0:c0 + CH], seeds[c0:c0 + CH]
    Ts = P.geom_sched(100.0, 0.1, S)
    s = s0.astype(np.float64).copy(); E = np.zeros((S + 1, s.shape[0]))
    P._tec_glauber_kernel(J.astype(np.float64), s, Ts, float(JVS[ji]), np.asarray(seeds, np.int64), E)
    return ji, c0, (sumw - E) / 2.0


def main(workers):
    cuts = {ji: np.zeros((S + 1, B)) for ji in range(len(JVS))}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for ji, c0, c in ex.map(job, [(ji, c0) for ji in range(len(JVS)) for c0 in range(0, B, CH)]):
            cuts[ji][:, c0:c0 + CH] = c
    out = {}
    for ji, jv in enumerate(JVS):
        mean = cuts[ji].mean(1)
        first = int(np.argmax(mean >= 31670)) if (mean >= 31670).any() else None
        fp = [int(np.argmax(cuts[ji][:, b] >= 31670)) if (cuts[ji][:, b] >= 31670).any() else None for b in range(B)]
        fpv = [x for x in fp if x is not None]
        out[str(jv)] = dict(first_cycle_mean_cut_31670=first, median_first_passage=float(np.median(fpv)) if fpv else None,
                            final_mean_cut=float(mean[-1]), p33000_final=float((cuts[ji][-1] >= 33000).mean()),
                            mean_cut_trace_every_100=mean[::100].tolist())
        print(jv, out[str(jv)]['first_cycle_mean_cut_31670'], out[str(jv)]['median_first_passage'], out[str(jv)]['final_mean_cut'], flush=True)
    a, b = out['0.0']['first_cycle_mean_cut_31670'], out['30.0']['first_cycle_mean_cut_31670']
    out['verdict'] = dict(jv0=a, jv30=b, ratio=(a / b if a and b else None),
                          pass_jv0=bool(a and abs(a - 1392) <= 0.15 * 1392), pass_jv30=bool(b and abs(b - 740) <= 0.15 * 740),
                          pass_ratio=bool(a and b and abs(a / b - 1.88) <= 0.2 * 1.88))
    print(out['verdict'])
    (HERE / 'results_validation').mkdir(exist_ok=True)
    (HERE / 'results_validation' / 'tec_fig3b.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 32)
