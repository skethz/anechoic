"""Addendum 3.1 analysis: analyze_a3.py with APC-SCA and ReAIM ASA taken from results_a3/final_bv (their papers' output
rules). Writes analysis_a3b.txt, results_summary_a3b.json, compact_table_a3b.md. Usage: python3 analyze_a3b.py"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyze_a3 as A  # noqa: E402

BV = ('APC-SCA', 'ReAIM ASA')
_orig = A.load_a3


def load_a3b():
    out = _orig()
    for prob in out:
        for t in out[prob]:
            for m in BV:
                for S in list(out[prob][t][m]):
                    f = HERE / 'results_a3' / 'final_bv' / prob / t / f"{m.replace(' ', '_')}_S{S}.json"
                    out[prob][t][m][S] = json.loads(f.read_text())
    return out


def main():
    A.load_a3 = load_a3b
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        A.main()
    # move the outputs written by analyze_a3.main to the a3b names, then restore the A3 (final-state) outputs
    for src, dst in (('analysis_a3.txt', 'analysis_a3b.txt'), ('results_summary_a3.json', 'results_summary_a3b.json'),
                     ('compact_table_a3.md', 'compact_table_a3b.md')):
        (HERE / dst).write_text((HERE / src).read_text().replace('# Addendum 3:', '# Addendum 3.1 (APC-SCA and ReAIM with their papers\' output rules):', 1))
    A.load_a3 = _orig
    with contextlib.redirect_stdout(io.StringIO()):
        A.main()
    print((HERE / 'analysis_a3b.txt').read_text())


if __name__ == '__main__':
    main()
