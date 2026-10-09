"""After gset_pub.py and gset_pub_sigma.py: write the board pre-registration of the published-setting baselines for
Table 8a (PROTOCOL_GSET_PUB.md last section; Amendment A1) and STOP. Lines follow fpga/v80_sca/PROTOCOL_HW_MB_A3_gset.txt
(key instance target trials host-arguments; seed 20261004, trial ids 0..trials-1, one trial per launch, tables
resident). Two readings: 'pub' = the papers' K2000 values unchanged; 'pubs' = the COP study's sigma transfer.
Nothing here touches the board."""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRIALS = 6156
VARIANTS = (('results_gset_pub', {'SCA-pub': 'plainpub', 'TEC-pub': 'tecpub'}, 'K2000 values unchanged'),
            ('results_gset_pub_sigma', {'SCA-pub-sigma': 'plainpubs', 'TEC-pub-sigma': 'tecpubs'}, 'sigma transfer (COP rule)'))


def args(cfg, tfin):
    a = f"--q {cfg['q']:.12g} --t0 {cfg['T0']:.12g} --t1 {tfin:.12g} --steps {cfg['S']}"
    if cfg['family'] == 'tec':
        return a + f" --tec-jv {cfg['jv']:.12g} --corr-scale 1"
    return a + ' --lambda 0 --corr-scale 1'


def main():
    lines = ['# PROTOCOL_GSET_PUB board pre-registration (frozen, no board data): key instance target trials host-arguments;',
             '# seed 20261004, trial ids 0..trials-1, one trial per launch, tables resident. Published-setting baselines for',
             '# Table 8a: plainpub/tecpub = STATICA q 4, T 40/30/50 -> 5 (+ TEC J_v +30), K2000 values unchanged;',
             '# plainpubs/tecpubs = the same transferred by alpha = sigma/sigma_K2000 (COP rule). Not run by this audit.']
    pred = {}
    for folder, rules, label in VARIANTS:
        rows = {int(k): v for k, v in json.loads((HERE / folder / 'per_instance.json').read_text()).items()}
        for g in sorted(rows):
            for rule, tag in rules.items():
                r = rows[g][rule]
                tfin = r.get('tfin', 5.0)
                key = f'G{g}_{tag}'
                lines.append(f"{key} G{g} {r['target']} {TRIALS} {args(r['cfg'], tfin)}")
                pred[key] = dict(instance=f'G{g}', rule=rule, reading=label, cfg=r['cfg'], tfin=tfin, model_p=r['p'],
                                 model_tts12_ms=r['tts12_ms'] if math.isfinite(r['tts12_ms']) else None,
                                 model_tts12_secondary_ms=r['tts12_secondary_ms'] if math.isfinite(r['tts12_secondary_ms']) else None,
                                 model_t_round12_ms=r['t_round12_ms'], model_mean_cut=r['mean_cut'])
    (HERE / 'PROTOCOL_GSET_PUB_board.txt').write_text('\n'.join(lines) + '\n')
    (HERE / 'PROTOCOL_GSET_PUB_predictions.json').write_text(json.dumps(pred, indent=1) + '\n')
    n_none = sum(1 for v in pred.values() if v['model_tts12_ms'] is None)
    print(len(pred), 'cohorts;', n_none, 'predicted never to reach the target in the model')


if __name__ == '__main__':
    main()
