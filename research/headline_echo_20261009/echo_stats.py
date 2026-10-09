"""Echo statistics for the headline figure: per-step flips and flips that undo the same spin's flip one step earlier,
from the engine's float model (same update as flips_profile_v2.run), 1,024 runs per configuration on K2000.
Used for Figure 1(c) of the paper; writes echo_stats.json next to this file."""
import json, sys
from pathlib import Path
import numpy as np
R = Path(str(__import__('pathlib').Path(__file__).resolve().parent / '..'))
for d in ('theory_ideas_20261003', 'reaim_reproduction_20261003', 'ablation_20261005', 'algorithm_compare_20261007'):
    sys.path.insert(0, str(R / d))
import methods as M
import abl as _abl
_abl.J_RESULTS = R / 'theory_ideas_20261003' / 'J_results.json'
CFGS = {
 'B5':       dict(family='plain', q=4.0, T0=30.0, S=560),                                   # STATICA's short schedule
 'O1_off':   dict(family='plain', q=8.0, T0=15.0, S=280),                                   # O1's schedule, correction off
 'O1':       dict(family='tecT', q=8.0, kappa=1.75, ramp=True, T0=15.0, S=280),             # Anechoic, Onsager-kT
}
def run(J, cfg, B, seed):
    rng = np.random.default_rng(seed)
    s = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
    fam = cfg['family']; S = cfg['S']; T = M.m.sched(cfg['T0'], S)
    h = s @ J; s_prev = None; T_prev = None; prev_flip = np.zeros(s.shape, bool)
    q = np.full(s.shape, float(cfg.get('q', 0.0)), np.float32)
    flips = np.zeros((S, B)); back = np.zeros((S, B))
    for t, Tt in enumerate(T):
        field = h
        if s_prev is not None and fam == 'tecT':
            field = h - cfg['kappa'] * T_prev * M.abl.ramp_factor(t, S, cfg['ramp']) * s_prev
        z = s * field + q
        flip = np.clip(z / (4 * Tt) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev = s.copy(); T_prev = Tt
        s += d; h += d @ J
        flips[t] = flip.sum(1); back[t] = (flip & prev_flip).sum(1); prev_flip = flip
    return flips, back, s, h
J, sumw = M.m.load()
out = {}
for i, (k, cfg) in enumerate(CFGS.items()):
    flips, back, s, h = run(J, cfg, 1024, 97000 + i)
    cuts = (sumw + 0.5 * np.einsum('bi,bi->b', s, h)) / 2
    out[k] = dict(cfg=cfg, echo_fraction=float(back.sum() / flips.sum()), flips_per_trial=float(flips.sum(0).mean()),
                  p_33000=float((cuts >= 33000).mean()), mean_cut=float(cuts.mean()),
                  flips_step=flips.mean(1).tolist(), back_step=back.mean(1).tolist())
    print('%-7s echo %.3f  flips/trial %8.0f  p(cut>=33000) %.3f  mean cut %.0f' % (k, out[k]['echo_fraction'], out[k]['flips_per_trial'], out[k]['p_33000'], out[k]['mean_cut']))
Path(__file__).with_name('echo_stats.json').write_text(json.dumps(out) + '\n')
