#!/usr/bin/env python3
"""Independent evaluation of graph-partitioning (bisection) and TSP runs on the true, unquantized problem, from the final
spins written by sca_gpu_bias_tg (spins.bin: 32 x uint64 per trial, bit i = 1 means s_i = +1) and the per-trial records.
Conventions follow the export's objective mapping (research/reaim_benchmarks_20261008/export_hw.py, read only):
  GPP: spin i = vertex i+1 of the G-set file; feasible iff sum_i s_i = 0; value = cut on the original graph;
       quality = R / cut if feasible else 0; success iff feasible and cut <= target.
  TSP: spin v*n + j = city v (TSPLIB node v+1) at position j (cyclic); feasible iff x = (s > 0) is a permutation matrix;
       value = tour length (TSPLIB distances, cities in position order); quality = L*/L if feasible else 0;
       success iff feasible and L <= target.
TSPLIB distances: this file's own parser (EXPLICIT FULL_MATRIX / UPPER_ROW / LOWER_DIAG_ROW / UPPER_DIAG_ROW, GEO).

  eval_problem.py --spec OUT.json --kind gpp|tsp --name NAME --source FILE --ref R --target T --jint8 J.jint8 --bias B.bin
                  --alpha a --P P
        builds the instance spec (edges or W) and checks the exported quantized matrix and bias against a rebuild (mapping check)
  eval_problem.py SPEC.json PREFIX [PREFIX ...]   -> PREFIX.eval.jsonl, PREFIX.eval.json
  eval_problem.py --tsplib FILE.tsp [FILE.opt.tour]"""
import json, math, sys, os, gzip
import numpy as np

def _txt(path):
    raw = open(path, 'rb').read()
    return (gzip.decompress(raw) if path.endswith('.gz') else raw).decode('ascii')

def tsplib_matrix(path):
    t = _txt(path).replace(':', ' : ')
    hdr, sec = {}, None
    toks = t.split(); i = 0
    while i < len(toks):
        k = toks[i]
        if k in ('EDGE_WEIGHT_SECTION', 'NODE_COORD_SECTION'): sec = k; i += 1; break
        if i + 2 < len(toks) and toks[i + 1] == ':': hdr[k] = toks[i + 2]; i += 3
        else: i += 1
    n = int(hdr['DIMENSION']); wt = hdr.get('EDGE_WEIGHT_TYPE'); wf = hdr.get('EDGE_WEIGHT_FORMAT')
    W = np.zeros((n, n), np.int64)
    if sec == 'EDGE_WEIGHT_SECTION' and wt == 'EXPLICIT':
        vals = []
        while i < len(toks) and len(vals) < n * n:
            try: vals.append(int(float(toks[i])))
            except ValueError: break
            i += 1
        k = 0
        for a in range(n):
            rng = {'FULL_MATRIX': range(n), 'UPPER_ROW': range(a + 1, n), 'LOWER_DIAG_ROW': range(a + 1), 'UPPER_DIAG_ROW': range(a, n)}[wf]
            for b in rng:
                W[a, b] = vals[k]
                if wf != 'FULL_MATRIX': W[b, a] = vals[k]
                k += 1
    elif sec == 'NODE_COORD_SECTION' and wt == 'GEO':
        xy = []
        for _ in range(n): xy.append((float(toks[i + 1]), float(toks[i + 2]))); i += 3
        rad = lambda x: math.pi * (int(x) + 5.0 * (x - int(x)) / 3.0) / 180.0
        lat = [rad(a) for a, _ in xy]; lon = [rad(b) for _, b in xy]
        for a in range(n):
            for b in range(n):
                if a != b:
                    q1 = math.cos(lon[a] - lon[b]); q2 = math.cos(lat[a] - lat[b]); q3 = math.cos(lat[a] + lat[b])
                    W[a, b] = int(6378.388 * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)
    else:
        raise ValueError(f'unsupported TSPLIB file {path}: {wt} {wf}')
    np.fill_diagonal(W, 0)
    return n, W

def tour_file(path):
    t = _txt(path).split(); i = t.index('TOUR_SECTION') + 1; tour = []
    while t[i] != '-1': tour.append(int(t[i]) - 1); i += 1
    return tour

def read_jint8(path):
    b = open(path, 'rb').read(); assert b[:8] == b'SCAJINT8'
    n = int.from_bytes(b[8:12], 'little'); return np.frombuffer(b[16:16 + n * n], np.int8).reshape(n, n).astype(np.int64)

def rnd(x):  # round half away from zero (the export's quantizer)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)

def read_bias_raw(path, N):
    b = np.fromfile(path, dtype='<i4'); assert b.size == N, 'bias size'; return b.astype(np.int64)

def build_spec(a):
    """Instance spec for the evaluator, with a full mapping check: the exported quantized matrix and bias must equal
    round(alpha * J_full) and round(alpha * b_full) rebuilt here from the original graph / TSPLIB distances (P = penalty)."""
    J = read_jint8(a['jint8']); N = J.shape[0]; bias = read_bias_raw(a['bias'], N); alpha = float(a['alpha']); P = float(a['P'])
    if a['kind'] == 'gpp':
        lines = open(a['source']).read().split('\n'); n, m = (int(x) for x in lines[0].split())
        E = np.array([l.split() for l in lines[1:] if l.strip()], np.int64); assert len(E) == m and n == N
        ei, ej, w = E[:, 0] - 1, E[:, 1] - 1, E[:, 2]
        full = np.full((N, N), -P); full[ei, ej] += w; full[ej, ei] += w; np.fill_diagonal(full, 0)
        ok_J = bool(np.array_equal(rnd(alpha * full), J)); ok_b = bool((bias == 0).all())
        spec = dict(kind='gpp', name=a['name'], N=N, ei=ei.tolist(), ej=ej.tolist(), w=w.tolist(), ref=int(a['ref']), target=int(a['target']))
    else:
        n, W = tsplib_matrix(a['source']); assert n * n == N
        idx = lambda v, j: v * n + (j % n)
        full = np.zeros((N, N))
        for v in range(n):
            for j in range(n):
                r = idx(v, j)
                for k in range(n):
                    if k != j: full[r, idx(v, k)] = -P          # same city, other positions
                for u in range(n):
                    if u != v:
                        full[r, idx(u, j)] = -P                 # same position, other cities
                        full[r, idx(u, j + 1)] = -W[v, u]       # adjacent positions, other cities
                        full[r, idx(u, j - 1)] = -W[v, u]
        D = W.sum(1); bfull = np.array([-2.0 * P * (n - 2) - 2.0 * D[v] for v in range(n) for j in range(n)])
        ok_J = bool(np.array_equal(rnd(alpha * full), J)); ok_b = bool(np.array_equal(rnd(alpha * bfull).astype(np.int64), bias))
        spec = dict(kind='tsp', name=a['name'], N=N, n=n, W=W.tolist(), ref=int(a['ref']), target=int(a['target']))
    spec.update(mapping_check_J=ok_J, mapping_check_b=ok_b, alpha=alpha, P=P)
    json.dump(spec, open(a['spec'], 'w'))
    print(json.dumps(dict(name=a['name'], kind=a['kind'], N=N, mapping_check_J=ok_J, mapping_check_b=ok_b)))

def spins_of(path, ntr, N):
    raw = np.fromfile(path, dtype='<u8'); assert raw.size == 32 * ntr, 'spins.bin size'
    words = raw.reshape(ntr, 32)
    idx = np.arange(N)
    bits = (words[:, idx >> 6] >> (idx & 63).astype(np.uint64)) & np.uint64(1)
    return np.where(bits == 1, 1, -1).astype(np.int64)

def evaluate(spec, S):
    T = S.shape[0]
    if spec['kind'] == 'gpp':
        ei = np.array(spec['ei']); ej = np.array(spec['ej']); w = np.array(spec['w'])
        cut = ((S[:, ei] != S[:, ej]) * w[None, :]).sum(1); M = S.sum(1); feas = M == 0
        return dict(feasible=feas, value=np.where(feas, cut, -1), cut_any=cut, imbalance=M,
                    quality=np.where(feas, spec['ref'] / np.maximum(cut, 1), 0.0), success=feas & (cut <= spec['target']), opt=feas & (cut <= spec['ref']))
    n = spec['n']; W = np.array(spec['W']); x = S.reshape(T, n, n) > 0
    rows = (x.sum(2) == 1).all(1); cols = (x.sum(1) == 1).all(1); feas = rows & cols
    L = np.full(T, -1, np.int64)
    for r in np.nonzero(feas)[0]:
        tour = np.argmax(x[r], axis=0); L[r] = int(W[tour, np.roll(tour, -1)].sum())
    return dict(feasible=feas, value=L, rows_ok=rows, cols_ok=cols, quality=np.where(feas, spec['ref'] / np.maximum(L, 1), 0.0),
                success=feas & (L <= spec['target']), opt=feas & (L <= spec['ref']))

def main():
    if sys.argv[1] == '--tsplib':
        n, W = tsplib_matrix(sys.argv[2]); out = dict(n=n)
        if len(sys.argv) > 3:
            t = tour_file(sys.argv[3]); out['opt_tour_length'] = int(sum(W[t[k], t[(k + 1) % len(t)]] for k in range(len(t))))
        print(json.dumps(out)); return
    if sys.argv[1] == '--spec':
        a = {}; it = iter(sys.argv[1:])
        for k in it: a[k.lstrip('-')] = next(it)
        build_spec(a); return
    spec = json.load(open(sys.argv[1]))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import analyze as A
    for pre in sys.argv[2:]:
        tr = [json.loads(l) for l in open(pre + '.trials.jsonl')]
        bat = [json.loads(l) for l in open(pre + '.batches.jsonl')]
        S = spins_of(pre + '.spins.bin', len(tr), spec['N'])
        ev = evaluate(spec, S)
        batch = np.array([d['batch'] for d in tr])
        with open(pre + '.eval.jsonl', 'w') as f:
            for r, d in enumerate(tr):
                f.write(json.dumps(dict(batch=d['batch'], chain=d['chain'], trial_id=d['trial_id'], energy=d['energy'],
                                        **{k: (v[r].item() if hasattr(v[r], 'item') else v[r]) for k, v in ev.items()})) + '\n')
        B = bat[0]['chains']; nb = len(bat); t = sum(b['device_ms'] for b in bat) / nb
        succ = ev['success']; k = int(succ.sum()); n = len(tr)
        kb = int(sum(1 for b in range(nb) if succ[batch == b].any()))
        feas = ev['feasible']; nf = int(feas.sum())
        p_lo, p_hi = A.wilson(k, n); Pb_lo, Pb_hi = A.wilson(kb, nb); f_lo, f_hi = A.wilson(nf, n)
        summ = dict(prefix=os.path.basename(pre), B=B, batches=nb, trials=n, t_batch_ms=t, successes=k, p=k / n, p_wilson95=[p_lo, p_hi],
                    batches_with_success=kb, P_batch=kb / nb, P_batch_wilson95=[Pb_lo, Pb_hi], feasible=nf, feasibility=nf / n,
                    feasibility_wilson95=[f_lo, f_hi], mean_quality=float(ev['quality'].mean()),
                    mean_quality_feasible=float(ev['quality'][feas].mean()) if nf else None,
                    best_value=int(ev['value'][feas].min()) if nf else None, optimum_hits=int(ev['opt'].sum()),
                    tts_primary_ms=A.tts(kb / nb, t), tts_primary_ci_ms=[A.tts(Pb_hi, t), A.tts(Pb_lo, t)],
                    tts_secondary_ms=A.tts_sec(k / n, B, t), tts_secondary_ci_ms=[A.tts_sec(p_hi, B, t), A.tts_sec(p_lo, B, t)])
        json.dump(summ, open(pre + '.eval.json', 'w'), indent=1)
        print(json.dumps({k2: summ[k2] for k2 in ('prefix', 'B', 'trials', 'feasibility', 'p', 'P_batch', 'mean_quality', 'best_value', 'tts_primary_ms', 'tts_secondary_ms')}))

if __name__ == '__main__':
    main()
