"""Bit-for-bit check of hwmodel.py (numba port) against the golden reference sca_ref_bias.hpp::run_trial_bias, through
hw_check.cpp (tables from mb_common.hpp::tables). Cases: quantized GPP, TSP and MCP instances, the four engine rules
(plain, TEC, Onsager-kT with ramp, Onsager-online with ramp), several trials each. Compares the tables, final spins,
flips, sum s*h, sum b*s and the full n_lin trace. Run on gpu-host after building hw_check (see hw_check.cpp)."""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import hwmodel as HW  # noqa: E402
import problems as PB  # noqa: E402

BIN = HERE / 'hw_check'


def unpack(hexs, n):
    words = [int(hexs[16 * w:16 * w + 16], 16) for w in range(len(hexs) // 16)]
    return np.array([1 if (words[i >> 6] >> (i & 63)) & 1 else -1 for i in range(n)], np.int64)


def cases():
    bkv = json.loads((HERE.parent / 'gset_20261007' / 'bkv.json').read_text())['instances']
    out = []
    g1, _ = PB.gpp(1, 4).quantize(4); out.append(g1)
    g14, _ = PB.gpp(14, 4).quantize(3); out.append(g14)
    t17, _ = PB.tsp('gr17').quantize(8); out.append(t17)
    t29, _ = PB.tsp('bayg29').quantize(4); out.append(t29)
    m1, _ = PB.mcp(1, bkv['G1']).quantize(2); out.append(m1)
    return out


def main():
    res = []; ok = True
    tmp = Path(tempfile.mkdtemp(dir=os.environ.get('TMPDIR')))
    for P in cases():
        J8 = P.jint8(); jf = tmp / f'{P.name}.jint8'; bf = tmp / f'{P.name}.bias'
        PB.write_jint8(jf, J8); PB.write_bias(bf, P.b.astype(np.int32))
        a = P.sigma_T; st = P.stats()
        rules = [dict(rule='plain', t0=1.5 * a, t1=0.3 * a, S=300, q=1.0 * a),
                 dict(rule='tec', t0=1.5 * a, t1=0.3 * a, S=300, q=1.0 * a, jv=-0.2 * a),
                 dict(rule='kT', t0=1.5 * a, t1=0.3 * a, S=300, q=1.0 * a, kappa=0.003 * st['kappa_factor'], ramp=True),
                 dict(rule='online', t0=1.5 * a, t1=0.3 * a, S=300, q=1.0 * a, lam=0.02 * st['lam_factor'], ramp=True)]
        for sch in rules:
            args = [str(BIN), str(jf), str(bf), repr(sch['t0']), repr(sch['t1']), str(sch['S']), repr(sch['q']),
                    repr(sch.get('lam', 0.0)), repr(sch.get('jv', 0.0)), repr(sch.get('kappa', 0.0)),
                    '1' if sch.get('ramp') else '0', '20261008', '3', '3']
            lines = subprocess.run(args, capture_output=True, text=True, check=True).stdout.strip().split('\n')
            ct = json.loads(lines[0])['tables']
            ft, qq, kc, kk = HW.tables(sch)
            tab_ok = (ct['fourT'] == ft.tolist() and ct['q'] == qq.tolist() and ct['kcorr'] == kc.tolist()
                      and ct['kconst'] == kk.tolist())
            out = HW.run(P, sch, 20261008, np.array([3, 4, 5]), trace=True)
            for k, ln in enumerate(lines[1:]):
                c = json.loads(ln)
                sp = unpack(c['spins'], P.N)
                same = (np.array_equal(sp, out['spins'][k].astype(np.int64)) and c['flips'] == int(out['flips'][k])
                        and c['sum_sh'] == int(out['sum_sh'][k]) and c['sum_bs'] == int(out['sum_bs'][k])
                        and sum(c['nlin']) == int(out['nlin_total'][k]))
                if k == 2:
                    same = same and c['nlin'] == out['nlin_trace_last'].tolist()
                # independent objective check: energy from the spins equals -(sum_sh + sum_bs)/2
                e = float(P.energy(sp[None, :].astype(np.float64))[0])
                e_ok = e == -(c['sum_sh'] + c['sum_bs']) / 2
                ok &= bool(same and tab_ok and e_ok)
                res.append(dict(instance=P.name, rule=sch['rule'], trial=c['trial'], tables_identical=bool(tab_ok),
                                run_identical=bool(same), energy_consistent=bool(e_ok), flips=c['flips']))
                print(f"{P.name:10s} {sch['rule']:6s} trial {c['trial']}: tables={tab_ok} run={same} energy={e_ok} "
                      f"flips={c['flips']}", flush=True)
    import platform
    (HERE / 'verify_hw.json').write_text(json.dumps(dict(all_identical=bool(ok), cases=res, node=platform.node(),
                                                         machine=platform.machine()), indent=1) + '\n')
    print('ALL IDENTICAL' if ok else 'MISMATCH')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
