#!/usr/bin/env python3
"""Markdown tables from an Amendment-3 summary (a3_summary.json): per instance the original-grid TTS99 of each rule (board p
in brackets) and the post-hoc best over rules and both grids; geometric means per rule and class. Reporting only.
Usage: a3_tables.py <a3_summary.json>"""
import json
import sys

RULES = [('SCA', 'plain'), ('TEC', 'tec'), ('Onsager-kT', 'tecT'), ('Onsager-online', 'ons')]


def fmt(x):
    return '–' if x is None else (f'{x:.3f}' if x < 10 else f'{x:.2f}' if x < 100 else f'{x:.1f}')


def main():
    S = json.load(open(sys.argv[1])); C = {k: c for k, c in S['cohorts'].items() if isinstance(c, dict)}
    insts = sorted({c['instance'] for c in C.values()}, key=lambda s: int(s[1:]))
    print('| Instance | N | Class | Target | SCA | TEC | Onsager-kT | Onsager-online | Best (key) | Study best E12 |')
    print('|---|---:|---|---:|---:|---:|---:|---:|---:|---:|')
    for i in insts:
        row = []
        cs = [c for c in C.values() if c['instance'] == i]
        c0 = cs[0]
        for rule, short in RULES:
            o = [c for c in cs if c['rule'] == rule and 'original' in c['variants']]
            row.append(f"{fmt(o[0]['tts99_ms'])} ({o[0]['p']:.2f})" if o else '–')
        b = S['best_per_instance'][i]
        sb = min(c['study_tts_E12_ms'] for c in cs)
        print(f"| {i} | {c0['N']} | {c0['class']} | {c0['target']} | " + ' | '.join(row) + f" | {fmt(b['tts99_ms'])} ({b['key']}) | {fmt(sb)} |")
    print()
    print('| Geometric mean TTS99 (ms) | random | toroidal | planar | all |')
    print('|---|---:|---:|---:|---:|')
    for rule in [r for r, _ in RULES] + ['best']:
        g = S['geomeans']
        print(f"| {rule}{' (post-hoc)' if rule == 'best' else ' (original grid)'} | " + ' | '.join(
            f"{fmt(g[f'{rule}/{c}']['geomean_tts99_ms'])} ({g[f'{rule}/{c}']['instances']})" for c in ('random', 'toroidal', 'planar', 'all')) + ' |')


if __name__ == '__main__':
    main()
