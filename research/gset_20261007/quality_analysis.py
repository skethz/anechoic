"""Solution quality on G-set (post hoc analysis of the existing held-out runs; no new runs).
For every instance, method and step budget S, results/<variant>/<inst>/<method>_S<S>.json holds 256 held-out runs of
the configuration selected by pilot mean cut on separate seeds. Quality measures per instance:
  mean gap  = 100 * (BKV - mean final cut) / BKV
  best gap  = 100 * (BKV - best final cut of the 256 runs) / BKV
  BKV hit   = at least one of the 256 runs reaches the best-known cut
Aggregates: average over instances; 'wins' = instances on which the method has the highest held-out mean cut
among the compared methods (ties shared). Usage: python3 quality_analysis.py [budget|ext_budget]"""
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
VAR = sys.argv[1] if len(sys.argv) > 1 else 'budget'
RES = HERE / 'results' / VAR
BUDGETS = (250, 500, 1000, 2000, 4000)
ENGINE = ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online')
GROUPS = {'Random, +1': ('R800+', 'R2000+', 'R1000+'), 'Random, ±1': ('R800+-', 'R2000+-'),
          'Toroidal, ±1': ('T800+-', 'T2000+-'), 'Planar-like, +1': ('P800+', 'P2000+', 'P1000+'),
          'Planar-like, ±1': ('P800+-', 'P2000+-')}
CLS = {k: v['cls'] for k, v in json.loads((HERE / 'results_summary.json').read_text())['instances'].items()} \
    if isinstance(json.loads((HERE / 'results_summary.json').read_text())['instances'], dict) else {}


def load():
    data = defaultdict(dict)   # data[inst][(method, S)] = final dict
    for d in sorted(RES.iterdir()):
        for f in d.glob('*_S*.json'):
            m, s = f.stem.rsplit('_S', 1)
            data[d.name][(m, int(s))] = json.loads(f.read_text())['final']
    return data


def table(data, methods, insts, S):
    rows = {m: dict(mean_gap=[], best_gap=[], hits=0, wins=0.0) for m in methods}
    for i in insts:
        fin = {m: data[i].get((m, S)) for m in methods}
        if any(v is None for v in fin.values()):
            continue
        for m, v in fin.items():
            b = v['BKV']
            rows[m]['mean_gap'].append(100 * (b - v['mean_cut']) / b)
            rows[m]['best_gap'].append(100 * (b - v['max_cut']) / b)
            rows[m]['hits'] += int(v['k_bkv'] > 0)
        top = max(v['mean_cut'] for v in fin.values())
        winners = [m for m, v in fin.items() if v['mean_cut'] == top]
        for m in winners:
            rows[m]['wins'] += 1 / len(winners)
    n = len(rows[methods[0]]['mean_gap'])
    out = {}
    for m, r in rows.items():
        out[m] = dict(n=n, mean_gap=sum(r['mean_gap']) / n, best_gap=sum(r['best_gap']) / n, bkv_instances=r['hits'],
                      wins=r['wins'])
    return out


def main():
    data = load()
    insts = sorted(data, key=lambda s: int(s[1:]))
    subset = [i for i in insts if len({m for m, _ in data[i]}) >= 10]
    summary = {'variant': VAR, 'subset': subset, 'engine_all51': {}, 'all_methods_subset': {}, 'engine_by_group_S4000': {}}
    print(f'variant {VAR}: {len(insts)} instances; ten-method subset {subset}')
    for S in BUDGETS:
        summary['engine_all51'][S] = table(data, ENGINE, insts, S)
    allm = sorted({m for i in subset for m, _ in data[i]})
    for S in (BUDGETS if subset else ()):
        summary['all_methods_subset'][S] = table(data, allm, subset, S)
    if CLS:
        for g, cl in GROUPS.items():
            gi = [i for i in insts if CLS.get(i) in cl]
            summary['engine_by_group_S4000'][g] = table(data, ENGINE, gi, 4000)
    (HERE / f'quality_summary_{VAR}.json').write_text(json.dumps(summary, indent=1))

    def show(title, tab):
        print(f'\n{title}')
        print(f"{'method':16s} {'n':>3s} {'mean gap %':>10s} {'best gap %':>10s} {'BKV inst':>8s} {'wins':>6s}")
        for m, r in sorted(tab.items(), key=lambda kv: kv[1]['mean_gap']):
            print(f"{m:16s} {r['n']:3d} {r['mean_gap']:10.3f} {r['best_gap']:10.3f} {r['bkv_instances']:8d} {r['wins']:6.1f}")
    for S in BUDGETS:
        show(f'Four engine rules, all instances, S = {S}', summary['engine_all51'][S])
    for S in ((1000, 4000) if subset else ()):
        show(f'All methods, 12-instance subset, S = {S}', summary['all_methods_subset'][S])
    for g, tab in summary['engine_by_group_S4000'].items():
        show(f'Four engine rules, {g}, S = 4000', tab)


if __name__ == '__main__':
    main()
