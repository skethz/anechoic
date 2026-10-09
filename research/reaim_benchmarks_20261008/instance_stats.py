"""Coupling-value statistics of every instance (full integer precision and after K-bit quantization):
distinct |J| levels, max |J|, signed bits for J, bias range and bits, field bound, density, scaling factors.
Writes data/instance_stats.json."""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import problems as PB  # noqa: E402

KS = (2, 3, 4, 6, 8)


def main():
    bkv = json.loads((HERE.parent / 'gset_20261007' / 'bkv.json').read_text())['instances']
    ref = json.loads((HERE / 'data' / 'gpp_reference.json').read_text())['instances']
    probs = [PB.mcp(g, bkv[f'G{g}']) for g in range(1, 21)]
    probs += [PB.gpp(g, 4, ref=ref[f'G{g}']['R'], target=ref[f'G{g}']['target']) for g in (1, 2, 3, 4, 5, 14, 15, 16, 17)]
    probs += [PB.tsp(t) for t in ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29')]
    out = {}
    for P in probs:
        st = dict(P.stats()); st['kind'] = P.kind
        st['penalty'] = {k: P.meta[k] for k in ('P', 'A') if k in P.meta}
        st['quantized'] = {}
        for K in KS:
            Q, alpha = P.quantize(K)
            qs = Q.stats()
            # coupling-information loss: distinct levels kept, nonzero couplings rounded to zero
            st['quantized'][K] = dict(alpha=alpha, maxabs_J=qs['maxabs_J'], distinct_absJ_levels=qs['distinct_absJ_levels'],
                                      absJ_levels=qs['absJ_levels'][:16], max_abs_b=qs['max_abs_b'],
                                      bits_b_signed=qs['bits_b_signed'], bits_field_signed=qs['bits_field_signed'],
                                      nnz_offdiag=qs['nnz_offdiag'], lost_couplings=st['nnz_offdiag'] - qs['nnz_offdiag'])
        out[P.name] = st
        q2 = st['quantized'][2]; q4 = st['quantized'][4]; q8 = st['quantized'][8]
        print(f"{P.name:7s} {P.kind} N={P.N:4d} dens={st['density']:.3f} |J| levels={st['distinct_absJ_levels']:4d} "
              f"max|J|={st['maxabs_J']:.0f} Jbits={st['bits_J_signed']:2d} max|b|={st['max_abs_b']:.0f} bbits={st['bits_b_signed']:2d} "
              f"| K2 levels {q2['distinct_absJ_levels']} lost {q2['lost_couplings']} | K4 levels {q4['distinct_absJ_levels']} "
              f"lost {q4['lost_couplings']} max|b| {q4['max_abs_b']:.0f} | K8 levels {q8['distinct_absJ_levels']} lost {q8['lost_couplings']} "
              f"bbits {q8['bits_b_signed']}")
    (HERE / 'data' / 'instance_stats.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
