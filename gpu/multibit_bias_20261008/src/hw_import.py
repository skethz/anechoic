#!/usr/bin/env python3
"""Phase B import (PROTOCOL.md Amendment 3). Reads data/hw_export/selected_configs.json (copy of the ReAIM study's export),
checks the SHA-256 of every exported file, and for K = 4:
  * converts each raw int32 bias to SCABIAS1 (values unchanged);
  * builds the evaluator spec from the original problem (G-set graph / TSPLIB distances) with a full mapping check (the
    exported matrix and bias must equal round(alpha * J_full), round(alpha * b_full));
  * maps every rule's schedule onto the engine flags and checks table identity with mb_common.hpp::tables (tables_equal);
  * records the engine's range plan, table placement and residency (60 chains per wave required);
  * writes data/hw_export/plan_b.tsv and data/hw_export/import_check.json. Any failed check raises.
  --power KIND prints "tag seed jint8 bias flags" of the KIND cell with the largest S * N (first in plan order on ties)."""
import hashlib, json, os, re, subprocess, sys, struct
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eval_problem as EP

X = 'data/hw_export'; BIN = 'build/sca_gpu_bias_tg2'; K = '4'
K4 = ['--kbits', '4', '--r', '2', '--preg', '2', '--cs', '8', '--groups', '1']

def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()

def flags_of(r):
    f = f"--t0 {r['t0']!r} --t1 {r['t1']!r} --steps {int(r['S'])} --q {r['q']!r} --lambda {r['lam']!r} --tec-jv {r['jv']!r} --tecT {r['kappa']!r}"
    return f + (' --ramp' if r['ramp'] else '')

def main():
    sc = json.load(open(f'{X}/selected_configs.json'))
    if len(sys.argv) > 2 and sys.argv[1] == '--power':
        rows = [l.rstrip('\n').split('\t') for l in open(f'{X}/plan_b.tsv')]
        rows = [r for r in rows if r[2] == sys.argv[2]]
        best = max(rows, key=lambda r: int(r[5]) * int(r[3]))
        print(best[0], best[6], best[7], best[8], best[10]); return
    checks = dict(files={}, instances={}, schedules=[]); plan = []
    for name, ent in sc['instances'].items():
        for Kx, kk in ent['K'].items():
            for key in ('jint8', 'bias'):
                p = f'{X}/{kk[key]}'; ok = os.path.exists(p) and sha(p) == kk[key + '_sha256']
                checks['files'][kk[key]] = ok
                if not ok: raise SystemExit(f'hash mismatch or missing file: {p}')
        kk = ent['K'][K]; kind = ent['problem']; N = ent['N']; obj = ent['objective']
        jf = f"{X}/{kk['jint8']}"; braw = f"{X}/{kk['bias']}"; bf = f"{X}/{name}_K{K}.bias"
        raw = open(braw, 'rb').read(); assert len(raw) == 4 * N
        vals = struct.unpack(f'<{N}i', raw); assert all(-32767 <= v <= 32767 for v in vals), 'bias beyond 16 bits'
        with open(bf, 'wb') as f: f.write(b'SCABIAS1'); f.write(struct.pack('<II', N, 0)); f.write(raw)
        target = int(re.search(r'<=\s*(-?\d+)', obj['target']).group(1))
        ref = obj['reference_R'] if kind == 'gpp' else obj['optimum_Lstar']
        src = f'data/gset/{name}' if kind == 'gpp' else f'data/tsplib/{name}.tsp'
        P = ent['full_precision_stats']['maxabs_J']
        spec = f'{X}/{name}_K{K}.spec.json'
        EP.build_spec(dict(spec=spec, kind=kind, name=name, source=src, ref=ref, target=target, jint8=jf, bias=braw, alpha=kk['alpha'], P=P))
        sp = json.load(open(spec))
        checks['instances'][name] = dict(kind=kind, N=N, K=K, alpha=kk['alpha'], P=P, ref=ref, target=target, maxabs_J=kk['maxabs_J'], max_abs_b=kk['max_abs_b'],
                                         mapping_check_J=sp['mapping_check_J'], mapping_check_b=sp['mapping_check_b'])
        if not (sp['mapping_check_J'] and sp['mapping_check_b']): raise SystemExit(f'mapping check failed for {name}')
        for rule, r in kk['rules'].items():
            fl = flags_of(r)
            te = subprocess.run(['build/tables_equal'], input=f"{r['t0']!r} {r['t1']!r} {int(r['S'])} {r['q']!r} {r['lam']!r} {r['jv']!r} {r['kappa']!r} {1 if r['ramp'] else 0}\n",
                                capture_output=True, text=True)
            same = te.returncode == 0 and '"different":0' in te.stdout
            info = json.loads(subprocess.run([BIN, '--info', '--jint8', jf, '--bias', bf, '--target-energy', '0', *K4, *fl.split()], capture_output=True, text=True, check=True).stdout)
            row = dict(instance=name, rule=rule, S=int(r['S']), seed=int(r['seed']), tables_identical=same, plain_ok=info['plain_ok'], wide=info['wide'],
                       tables_global=info['tables_global'], chains_per_wave=info['chains_per_wave'], export_checks=r.get('checks'),
                       precision_study=r.get('precision_study'), flags=fl)
            checks['schedules'].append(row)
            if not same: raise SystemExit(f'table identity failed: {name} {rule}')
            if info['chains_per_wave'] < 60: raise SystemExit(f'k4 cannot hold 60 chains at S = {r["S"]}: {name} {rule}')
            tag = f"{name}_{rule.replace('-', '')}"
            plan.append([tag, name, kind, str(N), rule, str(int(r['S'])), str(int(r['seed'])), jf, bf, spec, fl])
    json.dump(checks, open(f'{X}/import_check.json', 'w'), indent=1, default=str)
    with open(f'{X}/plan_b.tsv', 'w') as f:
        for p in plan: f.write('\t'.join(p) + '\n')
    print('import ok:', len(checks['instances']), 'instances,', len(plan), 'cells; all files, mappings, tables and residency checked')

if __name__ == '__main__':
    main()
