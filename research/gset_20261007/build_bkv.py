"""Best-known cuts for G1-G47, G51-G54, parsed from the text of two published tables (no value typed by hand).

Source A: U. Benlic, J.-K. Hao, "Breakout local search for the Max-Cut problem", Engineering Applications of Artificial
  Intelligence 26(3):1162-1173, 2013, doi:10.1016/j.engappai.2012.09.001. Read from the authors' accepted manuscript
  (dated 3 Sept 2012), https://leria-info.univ-angers.fr/~jinkao.hao/papers/BenlicHaoMaxCut2012.pdf, Table 2
  (columns f_prev = previous best known, f_best = BLS best of 20 runs).
Source B: F. Ma, J.-K. Hao, "A multiple search operator heuristic for the max-k-cut problem", Annals of Operations
  Research 248(1-2):365-403, 2017, doi:10.1007/s10479-016-2234-0. Read from arXiv:1510.09156v1, Table 5 (column f_pre =
  best known reported by any max-cut algorithm in the literature incl. parallel GES; column MOH = their best).

BKV = max(A.f_best, B.f_pre, B.MOH). A.f_prev is recorded but not used where it exceeds every later value (G23 only).
Target = ceil(0.99 * BKV), computed exactly as (99*BKV + 99) // 100.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IDS = list(range(1, 48)) + [51, 52, 53, 54]


def target(bkv):
    return (99 * bkv + 99) // 100   # exact ceil(0.99*bkv)


def parse_A():
    rows = {}
    for line in (ROOT / 'sources/BenlicHaoMaxCut2012.txt').read_text().split('\n'):
        m = re.match(r'^\s*G(\d+)\s+(\d+)\s+(\d+|–)\s+(\d+)\((\d+)\)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$', line)
        if m:
            rows[int(m.group(1))] = dict(N=int(m.group(2)), f_prev=None if m.group(3) == '–' else int(m.group(3)),
                                         f_best=int(m.group(4)), line=line.strip())
    return rows


def parse_B():
    txt = (ROOT / 'sources/MaHao_maxkcut_arXiv1510.09156.txt').read_text()
    t5 = txt[txt.index('Table 5: Comparative results of the proposed MOH'):]
    rows = {}
    for line in t5.split('\n'):
        m = re.match(r'^G(\d+)\s+(\d+)\s+(\d+)\s+(.*)$', line.strip())
        if m:
            g = int(m.group(1))
            if g in rows:
                break
            rest = m.group(4).split()
            rows[g] = dict(N=int(m.group(2)), f_pre=int(m.group(3)), MOH=int(rest[-1]), line=line.strip())
    return rows


def main():
    A, B = parse_A(), parse_B()
    out = dict(
        sources=dict(
            A=dict(citation='U. Benlic and J.-K. Hao, Breakout local search for the Max-Cut problem, Engineering Applications '
                            'of Artificial Intelligence 26(3):1162-1173, 2013. doi:10.1016/j.engappai.2012.09.001',
                   read_from='authors\' accepted manuscript (3 Sept 2012), '
                             'https://leria-info.univ-angers.fr/~jinkao.hao/papers/BenlicHaoMaxCut2012.pdf, Table 2',
                   local_copy='sources/BenlicHaoMaxCut2012.pdf'),
            B=dict(citation='F. Ma and J.-K. Hao, A multiple search operator heuristic for the max-k-cut problem, Annals of '
                            'Operations Research 248(1-2):365-403, 2017. doi:10.1007/s10479-016-2234-0',
                   read_from='arXiv:1510.09156v1 (30 Oct 2015), Table 5 (f_pre, MOH)',
                   local_copy='sources/MaHao_maxkcut_arXiv1510.09156.pdf')),
        rule='BKV = max(A.f_best, B.f_pre, B.MOH); target = ceil(0.99*BKV). Both tables were read from preprint versions; '
             'the published versions were not accessible (paywalled) and were not compared.',
        instances={})
    for g in IDS:
        a, b = A[g], B[g]
        assert a['N'] == b['N']
        bkv = max(a['f_best'], b['f_pre'], b['MOH'])
        src = [k for k, v in (('A.f_best', a['f_best']), ('B.f_pre', b['f_pre']), ('B.MOH', b['MOH'])) if v == bkv]
        flags = []
        if a['f_prev'] is not None and a['f_prev'] > bkv:
            flags.append(f"A.f_prev = {a['f_prev']} exceeds every other listed value (A.f_best {a['f_best']}, B.f_pre "
                         f"{b['f_pre']}, B.MOH {b['MOH']}); not used. Sensitivity target ceil(0.99*{a['f_prev']}) = "
                         f"{target(a['f_prev'])} is also reported.")
        if b['f_pre'] > a['f_best']:
            flags.append(f"improved after BLS (A.f_best {a['f_best']} < B.f_pre {b['f_pre']})")
        if b['MOH'] > b['f_pre']:
            flags.append(f"improved by MOH ({b['MOH']} > f_pre {b['f_pre']})")
        out['instances'][f'G{g}'] = dict(N=a['N'], A_f_prev=a['f_prev'], A_f_best=a['f_best'], B_f_pre=b['f_pre'], B_MOH=b['MOH'],
                                         BKV=bkv, BKV_from=src, target=target(bkv), flags=flags,
                                         A_line=a['line'], B_line=b['line'])
        print(f"G{g:<3d} N={a['N']:5d} A.prev={a['f_prev']} A.best={a['f_best']} B.pre={b['f_pre']} B.MOH={b['MOH']} -> "
              f"BKV {bkv} target {target(bkv)} {'; '.join(flags)}")
    (ROOT / 'bkv.json').write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')


if __name__ == '__main__':
    main()
