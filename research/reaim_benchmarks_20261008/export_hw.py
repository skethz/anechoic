"""Hardware-ready export (PROTOCOL.md): for every GPP and TSP instance at K in {2, 4, 8}:
  hw_export/<inst>_K<K>.jint8      SCAJINT8: 8-byte magic, uint32 N, uint32 0, N*N int8 row-major (zero diagonal)
  hw_export/<inst>_K<K>.bias.bin   int32 little-endian, N entries
  hw_export/selected_configs.json  per rule: schedule in that K's units (mb_common.hpp Sched fields: t0, t1, S, q, lam,
                                   jv, kappa, ramp; scale = 1), range checks, precision-study result at that K, objective
                                   mapping and SHA-256 of every file.
Then (on gpu-host) python3 export_hw.py check: runs hw_check (the C++ golden reference) on the exported files for trials 0-1
of every (instance, K, rule) and compares with hwmodel bit for bit -> hw_export/export_check.json."""
import json
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import hwmodel as HW  # noqa: E402
import precision as PR  # noqa: E402
import problems as PB  # noqa: E402

OUT = HERE / 'hw_export'
KX = (2, 4, 8)
SEED = 20261008


def objective_mapping(P):
    if P.kind == 'gpp':
        return dict(problem='balanced min-cut bisection (Lucas 2014 Sec. 2.2, A = B = 1)',
                    spins='spin i <-> vertex i+1 of the G-set file research/gset_20261007/data/' + P.name,
                    decode='side(i) = sign(s_i)', feasible='sum_i s_i == 0 (400/400)',
                    objective='cut = number of G-set edges (i, j) with s_i != s_j, on the original unweighted graph',
                    quality='R / cut if feasible else 0', reference_R=P.meta['ref'], target='feasible and cut <= ' + str(P.meta['target']),
                    energy_note='E = -(sum_sh + sum_bs)/2 is the energy of the quantized matrix; evaluate the objective from the spins')
    n = P.meta['n']
    return dict(problem='TSP, one-hot (Lucas 2014 Sec. 7.1, B = 1, A = max distance = MQC)',
                spins=f'spin v*{n} + j <-> city v (TSPLIB node v+1) at tour position j (0..{n - 1}, cyclic)',
                decode='x_vj = (1 + s_vj)/2', feasible='x is a permutation matrix (every row and column sums to 1)',
                objective='tour length from the TSPLIB distance matrix (data/tsplib/' + P.name + '.tsp.gz), cities in position order',
                quality='L* / L if feasible else 0', optimum_Lstar=P.meta['ref'], target='feasible and L <= ' + str(P.meta['target']),
                energy_note='E = -(sum_sh + sum_bs)/2 is the energy of the quantized matrix; evaluate the objective from the spins')


def build():
    OUT.mkdir(exist_ok=True)
    sched = [s for s in json.loads((PR.OUT / 'schedules.json').read_text()) if s['problem'] in ('gpp', 'tsp')]
    res = {}
    for f in (PR.OUT / 'runs').glob('*.json'):
        d = json.loads(f.read_text())
        res[(d['instance'], d['rule'], d['K'])] = d
    out = dict(spec='fpga/v80_sca/MULTIBIT_SPEC.md (couplings and bias sections); tables: mb_common.hpp::tables with scale = 1',
               reference='fpga/v80_sca/src/sca_ref_bias.hpp::run_trial_bias(n, dense, bias, tables, seed, trial)',
               note=('Schedules: the held-out-selected configuration per rule at S_exp (precision/schedules.json) mapped to '
                     'K units: t0, t1, q, jv times alpha_K; kappa (kappa_eff) and lam (lambda_eff) recomputed from J^K. '
                     'precision_study gives the measured result of exactly this schedule and matrix in the bit-exact model '
                     '(1024 trials, seed in the entry, trials 0..1023).'),
               instances={})
    for s in sched:
        P = PR.problem(s['problem'], s['instance'])
        name = P.name
        ent = out['instances'].setdefault(name, dict(problem=s['problem'], N=P.N, objective=objective_mapping(P),
                                                     full_precision_stats=P.stats(), K={}))
        rel = PR.rel_params(s['problem'], s['cfg'], P)
        for K in KX:
            Q, alpha = P.quantize(K)
            kk = ent['K'].setdefault(str(K), None)
            if kk is None:
                jf = OUT / f'{name}_K{K}.jint8'; bf = OUT / f'{name}_K{K}.bias.bin'
                PB.write_jint8(jf, Q.jint8()); PB.write_bias(bf, Q.b.astype(np.int64))
                assert np.all(np.abs(Q.b) < 2 ** 31)
                qs = Q.stats()
                kk = dict(jint8=jf.name, jint8_sha256=PB.sha256(jf), bias=bf.name, bias_sha256=PB.sha256(bf), alpha=alpha,
                          maxabs_J=qs['maxabs_J'], distinct_absJ_levels=qs['distinct_absJ_levels'], absJ_levels=qs['absJ_levels'],
                          max_abs_b=qs['max_abs_b'], bits_b_signed=qs['bits_b_signed'], bits_field_signed=qs['bits_field_signed'],
                          lam_factor=qs['lam_factor'], kappa_factor=qs['kappa_factor'], sigma_T=Q.sigma_T, rules={})
                ent['K'][str(K)] = kk
            sch = PR.hw_schedule(rel, Q, s['S'])
            r = res.get((name, s['rule'], K))
            seed = 20261008 + 1000 * PR.PCODE[s['problem']] + 10 * (PR.icode(s['problem'], s['instance']) % 100) + PR.RULES.index(s['rule'])
            kk['rules'][s['rule']] = dict(family=rel['family'], S_exp=s['S'], S_exp_rule=s['why'], **sch,
                                          checks=HW.table_checks(sch, K), seed=seed,
                                          precision_study=None if r is None else dict(
                                              mean_quality=r['mean_quality'], p_feasible=r['p_feasible'], p_target=r['p_target'],
                                              trials=r['trials']),
                                          selected_cfg_float_model=s['cfg'])
    (OUT / 'selected_configs.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('export written:', len(out['instances']), 'instances')


def unpack(hexs, n):
    words = [int(hexs[16 * w:16 * w + 16], 16) for w in range(len(hexs) // 16)]
    return np.array([1 if (words[i >> 6] >> (i & 63)) & 1 else -1 for i in range(n)], np.int64)


def check():
    sc = json.loads((OUT / 'selected_configs.json').read_text())
    res = []; ok = True
    for name, ent in sc['instances'].items():
        P = PR.problem(ent['problem'], int(name[1:]) if ent['problem'] == 'gpp' else name)
        for K, kk in ent['K'].items():
            Q, _ = P.quantize(int(K))
            J8 = Q.jint8()
            # the exported file is the matrix the model runs
            raw = (OUT / kk['jint8']).read_bytes()
            same_file = raw[:8] == b'SCAJINT8' and np.array_equal(np.frombuffer(raw[16:], np.int8).reshape(P.N, P.N), J8)
            same_bias = np.array_equal(np.fromfile(OUT / kk['bias'], '<i4'), Q.b.astype(np.int32))
            for rule, rr in kk['rules'].items():
                sch = {k: rr[k] for k in ('t0', 't1', 'S', 'q', 'lam', 'jv', 'kappa', 'ramp')}
                if not rr['checks']['ref_int32_tables']:
                    res.append(dict(instance=name, K=K, rule=rule, skipped='tables exceed int32')); continue
                args = [str(HERE / 'hw_check'), str(OUT / kk['jint8']), str(OUT / kk['bias']), repr(sch['t0']), repr(sch['t1']),
                        str(sch['S']), repr(sch['q']), repr(sch['lam']), repr(sch['jv']), repr(sch['kappa']),
                        '1' if sch['ramp'] else '0', str(rr['seed']), '0', '2']
                lines = subprocess.run(args, capture_output=True, text=True, check=True).stdout.strip().split('\n')
                o = HW.run(Q, sch, rr['seed'], np.array([0, 1]), trace=True)
                ft, qq, kc, kc2 = HW.tables(sch); ct = json.loads(lines[0])['tables']
                tab_ok = ct['fourT'] == ft.tolist() and ct['q'] == qq.tolist() and ct['kcorr'] == kc.tolist() and ct['kconst'] == kc2.tolist()
                for k, ln in enumerate(lines[1:]):
                    c = json.loads(ln)
                    same = (np.array_equal(unpack(c['spins'], P.N), o['spins'][k].astype(np.int64)) and c['flips'] == int(o['flips'][k])
                            and c['sum_sh'] == int(o['sum_sh'][k]) and c['sum_bs'] == int(o['sum_bs'][k])
                            and sum(c['nlin']) == int(o['nlin_total'][k]) and (k < 1 or c['nlin'] == o['nlin_trace_last'].tolist()))
                    good = bool(same and tab_ok and same_file and same_bias); ok &= good
                    res.append(dict(instance=name, K=K, rule=rule, trial=c['trial'], identical=good, flips=c['flips']))
                    print(f'{name:7s} K={K} {rule:15s} trial {c["trial"]}: identical={good}', flush=True)
    (OUT / 'export_check.json').write_text(json.dumps(dict(all_identical=bool(ok), cases=res), indent=1) + '\n')
    print('ALL IDENTICAL' if ok else 'MISMATCH')


if __name__ == '__main__':
    build() if sys.argv[1] == 'build' else check()
