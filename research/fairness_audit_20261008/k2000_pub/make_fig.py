"""Regenerate Figure 2 from the PROTOCOL_PUB data with an UNMODIFIED copy of research/algorithm_compare_20261007/plot_alg.py
(it reads <its folder>/results/S*.json). Two variants:
  fig/     pre-registered PROTOCOL_PUB results (aSB: every batch flagged invalid by the frozen batch rule)
  fig_A1/  the same with aSB per run (PROTOCOL_PUB_A1.md: diverged runs count as failures; curve over non-diverged runs)
Outputs: alg_compare_published_settings{,_A1}.pdf/.png here; the A1 variant is also copied to
fpga27/images/alg_compare_published_settings.pdf (a NEW file; alg_compare.pdf is not touched)."""
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
S_LIST = (250, 500, 1000, 2000, 4000)


def main():
    a1 = json.loads((HERE / 'results_pub' / 'asb_perrun_A1.json').read_text())
    for d in ('fig', 'fig_A1'):
        (HERE / d / 'results').mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / 'research/algorithm_compare_20261007/plot_alg.py', HERE / d / 'plot_alg.py')
    for S in S_LIST:
        base = json.loads((HERE / 'fig' / 'results' / f'S{S}.json').read_text())   # written by analyze_pub.py
        x = a1[str(S)]
        f = dict(base['aSB']['final'])
        f.update(p33000=x['p33000'], p33000_wilson95=x['p33000_wilson95'], k33000=x['k33000'], mcs99=x['mcs99'],
                 mean_cut=x['mean_cut_nondiverged'], H_mean=x['H_mean_nondiverged'], H_p10=x['H_p10_nondiverged'],
                 H_p90=x['H_p90_nondiverged'], diverged=x['diverged'], amendment='PROTOCOL_PUB_A1 (per-run)')
        base['aSB'] = dict(final=f)
        (HERE / 'fig_A1' / 'results' / f'S{S}.json').write_text(json.dumps(base) + '\n')
    for d, out in (('fig', 'alg_compare_published_settings_prereg.pdf'), ('fig_A1', 'alg_compare_published_settings.pdf')):
        subprocess.run([sys.executable, str(HERE / d / 'plot_alg.py'), str(HERE / out)], check=True)
    shutil.copy(HERE / 'alg_compare_published_settings.pdf', ROOT / 'fpga27/images/alg_compare_published_settings.pdf')
    shutil.copy(HERE / 'alg_compare_published_settings.png', ROOT / 'fpga27/images/alg_compare_published_settings.png')
    print('written; fpga27/images/alg_compare.pdf untouched')


if __name__ == '__main__':
    main()
