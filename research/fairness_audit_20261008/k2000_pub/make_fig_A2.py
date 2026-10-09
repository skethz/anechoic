"""Figure 2 from results_pub_A2 (COP-aligned readings; assemble_A2.py) with an UNMODIFIED copy of plot_alg.py.
Output: alg_compare_published_settings_A2.pdf/.png here and the same as a NEW file
fpga27/images/alg_compare_published_settings_A2.pdf (alg_compare.pdf is not touched)."""
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    d = HERE / 'fig_A2'
    (d / 'results').mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / 'research/algorithm_compare_20261007/plot_alg.py', d / 'plot_alg.py')
    for S in (250, 500, 1000, 2000, 4000):
        shutil.copy(HERE / 'results_pub_A2' / f'S{S}.json', d / 'results' / f'S{S}.json')
    out = HERE / 'alg_compare_published_settings_A2.pdf'
    subprocess.run([sys.executable, str(d / 'plot_alg.py'), str(out)], check=True)
    shutil.copy(out, ROOT / 'fpga27/images/alg_compare_published_settings_A2.pdf')
    shutil.copy(out.with_suffix('.png'), ROOT / 'fpga27/images/alg_compare_published_settings_A2.png')
    print('written; fpga27/images/alg_compare.pdf untouched')


if __name__ == '__main__':
    main()
