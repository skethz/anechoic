"""Cross-platform reproduction fingerprint (laptop vs gpu-host): every method on small GPP/TSP/MCP cases with fixed seeds.
Prints the SHA-256 of the final spins and the mean normalized quality per case."""
import hashlib
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import hwmodel as HW  # noqa: E402
import problems as PB  # noqa: E402
import solvers as SV  # noqa: E402


def main():
    bkv = json.loads((HERE.parent / 'gset_20261007' / 'bkv.json').read_text())['instances']
    probs = [PB.gpp(14, 4, ref=1100, target=1111), PB.tsp('gr17'), PB.mcp(11, bkv['G11'])]
    out = {}
    for P in probs:
        a = P.sigma_T; st = P.stats()
        cfgs = [dict(family='plain', q=1.0 * a, T0=1.5 * a, tfin=0.3 * a, S=200),
                dict(family='tec', q=1.0 * a, jv=-0.2 * a, T0=1.5 * a, tfin=0.3 * a, S=200),
                dict(family='tecT', q=1.0 * a, kappa=0.01, ramp=True, T0=1.5 * a, tfin=0.3 * a, S=200),
                dict(family='onsager', q=1.0 * a, lam=0.05 * st['lam_factor'], ramp=True, T0=1.5 * a, tfin=0.3 * a, S=200),
                dict(family='apc', q_reset=2.0 * a, r_q=0.9, q_lim=0.5 * a, T0=1.5 * a, tfin=0.3 * a, S=200),
                dict(family='SA', T0=1.0 * a, T1=0.05 * a, S=100), dict(family='ReAIM', kset=(4, 8, 16, 32), T1=0.05, S=200),
                dict(family='aSB', dt=0.5, xi=1.0, S=200), dict(family='bSB', dt=0.75, xi=1.0, S=200),
                dict(family='dSB', dt=1.0, xi=1.0, S=200)]
        for cfg in cfgs:
            s, f, fin = SV.run_method(P, cfg, 16, np.random.SeedSequence([20261008, 77]))
            ev = P.evaluate(s)
            key = f"{P.name}/{cfg['family']}"
            out[key] = dict(sha=hashlib.sha256(np.ascontiguousarray(s, np.float32).tobytes()).hexdigest()[:16],
                            quality=float(ev['quality'].mean()), feasible=int(ev['feasible'].sum()))
            print(key, json.dumps(out[key]), flush=True)
        Q, _ = P.quantize(4)
        o = HW.run(Q, dict(t0=1.5 * Q.sigma_T, t1=0.3 * Q.sigma_T, S=200, q=1.0 * Q.sigma_T), 7, np.arange(4))
        out[f'{P.name}/hw_K4'] = dict(sha=hashlib.sha256(o['spins'].tobytes()).hexdigest()[:16], flips=o['flips'].tolist())
        print(f'{P.name}/hw_K4', json.dumps(out[f'{P.name}/hw_K4']), flush=True)
    tag = sys.argv[1] if len(sys.argv) > 1 else 'local'
    (HERE / f'repro_{tag}.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
