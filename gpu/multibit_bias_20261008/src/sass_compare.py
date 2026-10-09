#!/usr/bin/env python3
"""Static check for phase T: compare the SASS of the bias kernels (WIDE = false) with the multibit kernels of the same
template arguments. Reports instruction counts and an opcode-level diff (difflib) per kernel pair, plus the opcodes of
the inserted instructions. Usage (on gpu-host): python3 src/sass_compare.py BIAS_BINARY MULTIBIT_BINARY > results/T/sass_compare.json"""
import collections, difflib, json, re, subprocess, sys

def kernels(binary):
    txt = subprocess.run(['cuobjdump', '-sass', binary], capture_output=True, text=True, check=True).stdout
    out, name = {}, None
    for line in txt.splitlines():
        m = re.match(r'\s*Function : (\S+)', line)
        if m: name = m.group(1); out[name] = []; continue
        m = re.match(r'\s*/\*[0-9a-f]{4,}\*/\s+(.*?);', line)
        if m and name:
            ins = m.group(1).strip()
            ins = re.sub(r'^@!?U?P[T0-9]+\s+', '', ins)   # drop predicate guard
            out[name].append(ins.split()[0])
    return out

def key(mangled):
    m = re.search(r'mb_kernelILi(\d)ELi(\d)ELi(\d)ELi(\d+)ELi(\d)ELb([01])ELb([01])(?:ELb([01]))?E', mangled)
    if not m: return None
    K, R, P, CS, G, O, L, W = m.groups()
    return (f'K{K}R{R}P{P}G{G}' + ('ONS' if O == '1' else 'LIN' if L == '1' else 'CONST'), W)

bias, mb = kernels(sys.argv[1]), kernels(sys.argv[2])
B = {key(k): v for k, v in bias.items() if key(k) and key(k)[1] == '0'}
M = {key(k)[0]: v for k, v in mb.items() if key(k)}
res = {}
for (name, _), seq in sorted(B.items()):
    if name not in M: continue
    sm = difflib.SequenceMatcher(None, M[name], seq, autojunk=False)
    ins, dele = collections.Counter(), collections.Counter()
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag in ('replace', 'delete'): dele.update(M[name][i1:i2])
        if tag in ('replace', 'insert'): ins.update(seq[j1:j2])
    res[name] = dict(mb=len(M[name]), bias=len(seq), inserted=sum(ins.values()), deleted=sum(dele.values()),
                     inserted_ops=dict(ins.most_common(8)), deleted_ops=dict(dele.most_common(8)),
                     bmma_mb=M[name].count('BMMA.168256.AND.POPC'), bmma_bias=seq.count('BMMA.168256.AND.POPC'))
json.dump(res, sys.stdout, indent=1)
