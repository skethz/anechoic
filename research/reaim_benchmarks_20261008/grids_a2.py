"""Addendum 2 (PROTOCOL_ADDENDUM2.md): grid families and axes of all ten methods, edge counting, and the automatic
edge-extension rule. One procedure for every method:

  * Grid families: 'gpp', 'tsp' (relative units of grids.py) and the three Max-Cut weight classes of G1-G20
    'mcp:P+' (G1-G5), 'mcp:Q+' (G14-G17), 'mcp:M+-' (G6-G13, G18-G20) (relative units of research/gset_20261007/run_gset.py).
  * A1 grid = the frozen grid of the earlier run (grids.py for GPP/TSP, run_gset.py for Max-Cut). The A2 base adds
    (i) ReAIM's own options: F in {max, min} (ReAIM Step 4) and its published temperature schedule (ReAIM Table II:
        GPP T 1 -> 0.01, TSP T 0.5 -> 0.1; the Max-Cut schedule 1 -> 0.1 is already on its grid);
    (ii) APC-SCA's T_fin axis {0.15, 0.6} on GPP/TSP, the values every other SCA-family rule already has there.
  * An axis = a numeric parameter with >= 2 values (categorical F is not an axis; fixed parameters are not axes).
  * Edge count of an axis end = fraction of the family's selections (instances x 5 budgets) whose value equals that end.
    A selection whose pilot means are identical at every grid point (an all-tied grid, e.g. all infeasible) carries no
    direction and counts at both ends of every axis.
  * Trigger: an end with fraction > 0.20. Extension = one new value beyond that end with the existing spacing:
      positive or sign-definite axes: geometric, v_new = v_end^2 / v_neighbour;
      axes whose end and neighbour differ in sign: linear, v_new = 2 v_end - v_neighbour;
      r_q (APC decay factor): geometric in (1 - r_q);
      ReAIM k set: halve (low end, floor 1, stops at (1,1,1,1)) or double (high end, cap N) every flip cap.
  * Constraints (drop points): T1 <= T0 (SA, ReAIM), T_fin <= T0 (SCA family), q_lim < q_reset (APC).
  * Cap: a (method, family) grid may hold at most A1 size + CAP_EXTRA points (CAP_EXTRA = 224 for every method and family);
    every A2 addition counts, including ReAIM's F = min copies and published-schedule points and APC's T_fin values.
    Within a round, triggers are applied in decreasing order of their fraction (ties: axis order, low end first); an
    extension that would exceed the cap is skipped.
  * Up to three extension rounds; round 1 is triggered by the A1 selections, rounds 2 and 3 by the A2 selections of the
    previous round. Residual edges after round 3 are reported.
"""
import itertools
import json
import math
import sys
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE.parent
for _d in ('gset_20261007', 'ablation_20261005', 'algorithm_compare_20261007', 'reaim_reproduction_20261003',
           'theory_ideas_20261003', 'statica_reproduction_20261003'):
    if str(RES / _d) not in sys.path:
        sys.path.insert(0, str(RES / _d))
sys.path.insert(0, str(HERE))

THRESH = 0.20
MAX_ROUNDS = 3
CAP_EXTRA = 224
METHODS = ('SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'Onsager-kT', 'Onsager-online')
MCP_INST = tuple(range(1, 21))
GPP_INST = (1, 2, 3, 4, 5, 14, 15, 16, 17)
TSP_INST = ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29')
INSTANCES = {'mcp': MCP_INST, 'gpp': GPP_INST, 'tsp': TSP_INST}
S_LIST = {'mcp': (250, 500, 1000, 2000, 4000), 'gpp': (256, 512, 1024, 2048, 4096), 'tsp': (512, 1024, 2048, 4096, 8192)}
REAIM_PUBLISHED = {'gpp': (1.0, 0.01), 'tsp': (0.5, 0.1)}       # (T_init, T_final), ReAIM Table II
APC_TFIN_A2 = (0.15, 0.6)


def mcp_class(g):
    return 'P+' if 1 <= g <= 5 else ('Q+' if 14 <= g <= 17 else 'M+-')


def family_of(prob, inst):
    return f'mcp:{mcp_class(inst)}' if prob == 'mcp' else prob


def families():
    return ['mcp:P+', 'mcp:Q+', 'mcp:M+-', 'gpp', 'tsp']


def family_instances(fam):
    if fam == 'gpp':
        return [('gpp', g) for g in GPP_INST]
    if fam == 'tsp':
        return [('tsp', t) for t in TSP_INST]
    c = fam.split(':')[1]
    return [('mcp', g) for g in MCP_INST if mcp_class(g) == c]


# ------------------------------------------------------------------------------------------------ A1 axes
def _ax(values, kind='num'):
    return dict(values=list(values), kind=kind)


def a1_axes(fam, method):
    """(axes, fixed): the frozen A1 grid as ordered axes (product order = A1 point order) and fixed relative params."""
    if fam in ('gpp', 'tsp'):
        import grids
        g = grids.G[fam]; sh = g['shared']
        if method == 'SA':
            return OrderedDict(T0=_ax(g['SA']['T0']), T1=_ax(g['SA']['T1'])), {}
        if method == 'SCA':
            return OrderedDict(q=_ax(g['SCA']['q']), T0=_ax(g['SCA']['T0']), tfin=_ax(g['SCA']['tfin'])), {}
        if method == 'TEC':
            return OrderedDict(jv=_ax(g['TEC']['jv']), q=_ax(sh['q']), T0=_ax(sh['T0']), tfin=_ax(sh['tfin'])), {}
        if method == 'Onsager-kT':
            return OrderedDict(kappa=_ax(g['kT']['kappa']), q=_ax(sh['q']), T0=_ax(sh['T0']), tfin=_ax(sh['tfin'])), {}
        if method == 'Onsager-online':
            return OrderedDict(lam=_ax(g['online']['lam']), q=_ax(sh['q']), T0=_ax(sh['T0']), tfin=_ax(sh['tfin'])), {}
        if method == 'APC-SCA':
            c = g['APC']
            return (OrderedDict(q_reset=_ax(c['q_reset']), r_q=_ax(c['r_q'], 'rq'), q_lim=_ax(c['q_lim']), T0=_ax(c['T0'])),
                    dict(tfin=c['tfin']))
        if method == 'ReAIM ASA':
            return OrderedDict(kset=_ax(grids.REAIM_KSETS, 'kset'), T1=_ax(grids.REAIM_T1)), dict(T0=1.0, F='max')
        if method in ('aSB', 'bSB', 'dSB'):
            return OrderedDict(dt=_ax(g[method]['dt']), xi=_ax(g[method]['xi'])), {}
    else:
        import run_gset as RG  # only constants are used
        c = RG.SHARED[fam.split(':')[1]]
        if method == 'SA':
            return OrderedDict(T0=_ax(RG.SA_T0), T1=_ax(RG.SA_T1)), {}
        if method == 'SCA':
            return OrderedDict(q=_ax(c['plain_q']), T0=_ax(c['plain_T0'])), dict(tfin=c['tfin'])
        if method == 'TEC':
            return OrderedDict(jv=_ax(c['jv']), q=_ax(c['q']), T0=_ax(c['T0'])), dict(tfin=c['tfin'])
        if method == 'Onsager-kT':
            return OrderedDict(kappa=_ax(RG.KAPPA), q=_ax(c['q']), T0=_ax(c['T0'])), dict(tfin=c['tfin'])
        if method == 'Onsager-online':
            return OrderedDict(lam=_ax(RG.LAMBDA), q=_ax(c['q']), T0=_ax(c['T0'])), dict(tfin=c['tfin'])
        if method == 'APC-SCA':
            return (OrderedDict(q_reset=_ax(c['q_reset']), r_q=_ax((0.9, 0.97), 'rq'), q_lim=_ax(c['q_lim']), T0=_ax(c['T0'])),
                    dict(tfin=c['tfin']))
        if method == 'ReAIM ASA':
            return OrderedDict(kset=_ax(RG.REAIM_K, 'kset'), T1=_ax(RG.REAIM_T1)), dict(T0=1.0, F='max')
        if method == 'aSB':
            return OrderedDict(dt=_ax((0.25, 0.5, 0.9, 1.25)), xi=_ax((0.5, 0.75, 1.0, 1.5))), {}
        if method in ('bSB', 'dSB'):
            return OrderedDict(dt=_ax((0.5, 0.75, 1.0, 1.25)), xi=_ax((0.5, 0.75, 1.0, 2.0))), {}
    raise ValueError((fam, method))


# ------------------------------------------------------------------------------------------------ grid state
class Grid:
    """The A2 grid of one (method, family): ordered axes, fixed params, categorical option F (ReAIM), published
    schedule points (ReAIM), the ordered point list (A1 points first, then additions in the order added), and the
    extension log."""

    def __init__(self, fam, method):
        self.fam, self.method = fam, method
        self.axes, self.fixed = a1_axes(fam, method)
        self.a1_keyset = set()
        self.a1_keyset = {self._key(p) for p in self._a1_points_with_options()}
        self.a1_size = len(self.a1_keyset)
        self.options = {}
        self.published = None
        if method == 'ReAIM ASA':
            self.options['F'] = ['max', 'min']
            if fam in REAIM_PUBLISHED:
                self.published = REAIM_PUBLISHED[fam]
        if method == 'APC-SCA' and fam in ('gpp', 'tsp'):
            self.axes['tfin'] = _ax(APC_TFIN_A2)
            del self.fixed['tfin']
        self.log = []
        self.order = []          # keys in point order
        self._sync()

    # points are dicts of relative values; the key is a hashable tuple
    def _key(self, p):
        return tuple((k, (tuple(v) if isinstance(v, (list, tuple)) else (round(v, 12) if isinstance(v, float) else v)))
                     for k, v in sorted(p.items()))

    def _valid(self, p):
        T0 = p.get('T0', self.fixed.get('T0'))
        if self.method == 'SA' and p['T1'] > p['T0'] * (1 + 1e-12):
            return False
        if self.method == 'ReAIM ASA' and p['T1'] >= T0 * (1 - 1e-12):
            return False
        tf = p.get('tfin', self.fixed.get('tfin'))
        if self.method in ('SCA', 'TEC', 'Onsager-kT', 'Onsager-online', 'APC-SCA') and tf is not None and tf > p['T0'] * (1 + 1e-12):
            return False
        if self.method == 'APC-SCA' and p['q_lim'] >= p['q_reset']:
            return False
        return True

    def _product(self, axes=None, options=None, filt=True):
        """Product points; constraints apply to points that are not in the A1 grid (all A1 points are kept)."""
        axes = self.axes if axes is None else axes
        names = list(axes)
        out = []
        for vals in itertools.product(*[axes[n]['values'] for n in names]):
            base = dict(zip(names, [list(v) if isinstance(v, (list, tuple)) else v for v in vals]))
            for opt in itertools.product(*[[(k, v) for v in vs] for k, vs in (options or {}).items()]):
                p = dict(base); p.update(dict(opt))
                if not filt or self._key(p) in self.a1_keyset or self._valid(p):
                    out.append(p)
        return out

    def _all_points(self, axes=None):
        pts = self._product(axes, self.options)
        if self.published is not None:
            axes = self.axes if axes is None else axes
            T0, T1 = self.published
            for ks in axes['kset']['values']:
                for F in self.options['F']:
                    p = dict(kset=list(ks), T1=T1, F=F, T0=T0)
                    if self._key(p) in self.a1_keyset or self._valid(p):
                        pts.append(p)
        return pts

    def _sync(self):
        """Rebuild the ordered point list: existing order kept, new points appended in product order."""
        have = set(self.order)
        pts = {}
        # A1 points first (with the A1 option values), in A1 order
        for p in self._all_points():
            k = self._key(p); pts[k] = p
        a1_first = []
        for p in self._a1_points_with_options():
            k = self._key(p)
            if k in pts and k not in have and k not in a1_first:
                a1_first.append(k)
        for k in a1_first:
            self.order.append(k); have.add(k)
        for k in pts:
            if k not in have:
                self.order.append(k); have.add(k)
        self.points = {k: pts[k] for k in self.order if k in pts}
        self.order = [k for k in self.order if k in pts]

    def _a1_points_with_options(self):
        a1_axes_, fixed = a1_axes(self.fam, self.method)
        out = []
        for p in self._product(a1_axes_, None, filt=False):
            q = dict(p)
            if self.method == 'ReAIM ASA':
                q['F'] = 'max'
            if self.method == 'APC-SCA' and self.fam in ('gpp', 'tsp'):
                q['tfin'] = fixed['tfin']
            out.append(q)
        return out

    def size(self):
        return len(self.order)

    def cap(self):
        return self.a1_size + CAP_EXTRA

    def point_list(self):
        return [self.points[k] for k in self.order]

    def index_of(self, p):
        return self.order.index(self._key(p))

    # ------------------------------------------------------------------ extension
    @staticmethod
    def _new_value(ax, end, N=None):
        v = ax['values']; kind = ax['kind']
        if kind == 'kset':
            if end == 'min':
                lo = v[0]
                if all(k == 1 for k in lo):
                    return None
                return tuple(max(1, int(math.floor(k / 2))) if k > 1 else 1 for k in lo)
            hi = v[-1]
            new = tuple(int(min(N or 10 ** 9, 2 * k)) for k in hi)
            return None if new == tuple(hi) else new
        s = sorted(v)
        if kind == 'rq':
            u = sorted(1 - x for x in s)            # u = 1 - r_q
            if end == 'min':                        # smaller r_q = larger u
                un = u[-1] ** 2 / u[-2]
                return None if un >= 1 else 1 - un
            un = u[0] ** 2 / u[1]
            return 1 - un
        a, b = (s[0], s[1]) if end == 'min' else (s[-1], s[-2])
        if a * b > 0:
            return a * a / b
        return 2 * a - b

    def axis_end_value(self, name, end):
        ax = self.axes[name]
        if ax['kind'] == 'kset':
            return tuple(ax['values'][0] if end == 'min' else ax['values'][-1])
        return min(ax['values']) if end == 'min' else max(ax['values'])

    def extend(self, triggers, round_no, N=None):
        """triggers: list of (fraction, axis, end). Applies them in priority order under the cap; returns applied list."""
        applied = []
        order_axes = list(self.axes)
        for frac, name, end in sorted(triggers, key=lambda t: (-t[0], order_axes.index(t[1]), 0 if t[2] == 'min' else 1)):
            ax = self.axes[name]
            nv = self._new_value(ax, end, N)
            if nv is None:
                self.log.append(dict(round=round_no, axis=name, end=end, frac=frac, result='at bound'))
                continue
            trial = OrderedDict((k, dict(values=list(a['values']), kind=a['kind'])) for k, a in self.axes.items())
            if ax['kind'] == 'kset':
                trial[name]['values'] = ([list(nv)] + trial[name]['values']) if end == 'min' else (trial[name]['values'] + [list(nv)])
            else:
                vals = trial[name]['values'] + [nv]
                trial[name]['values'] = sorted(vals)
            keys = set(self.order) | {self._key(p) for p in self._all_points(trial)}
            if len(keys) > self.cap():
                self.log.append(dict(round=round_no, axis=name, end=end, frac=frac, value=nv, result='skipped (cap)'))
                continue
            self.axes = trial
            self._sync()
            applied.append((name, end, nv))
            self.log.append(dict(round=round_no, axis=name, end=end, frac=frac, value=nv, result='added', size=self.size()))
        return applied

    def to_json(self):
        return dict(fam=self.fam, method=self.method, a1_size=self.a1_size, cap=self.cap(), size=self.size(),
                    axes={k: v for k, v in self.axes.items()}, fixed=self.fixed, options=self.options,
                    published=self.published, log=self.log, points=self.point_list())


# ------------------------------------------------------------------------------------------------ edges
def _eq(a, b):
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        return tuple(a) == tuple(b)
    return abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def edge_counts(grid, selections):
    """selections: list of (point dict, tied flag). Returns {axis: dict(n, min, max, frac_min, frac_max)} over numeric axes
    with >= 2 values. Tied selections count at both ends."""
    out = OrderedDict()
    n = len(selections)
    for name, ax in grid.axes.items():
        if len(ax['values']) < 2:
            continue
        lo = grid.axis_end_value(name, 'min'); hi = grid.axis_end_value(name, 'max')
        cmin = cmax = 0
        for p, tied in selections:
            if tied:
                cmin += 1; cmax += 1; continue
            v = p.get(name)
            if v is None:
                continue
            cmin += _eq(v, lo); cmax += _eq(v, hi)
        out[name] = dict(n=n, min=cmin, max=cmax, frac_min=cmin / n if n else 0.0, frac_max=cmax / n if n else 0.0)
    return out


def triggers_from(counts):
    t = []
    for name, c in counts.items():
        if c['frac_min'] > THRESH:
            t.append((c['frac_min'], name, 'min'))
        if c['frac_max'] > THRESH:
            t.append((c['frac_max'], name, 'max'))
    return t


# ------------------------------------------------------------------------------------------------ point -> cfg
def cfg_of(prob, P, method, p, S):
    """Absolute configuration of a relative point on Problem P (grids.py / run_gset.py conventions)."""
    a = P.sigma_T; st = P.stats(); N = P.N
    if prob == 'mcp':
        c_lamf = (2.0 * P.meta['m'] / N) / N / (1999.0 / 2000.0)          # engine.Graph.lam_factor, exactly
        kf = 1.0
    else:
        c_lamf = st['lam_factor']; kf = st['kappa_factor']
    tf = p.get('tfin')
    if method == 'SA':
        return dict(family='SA', T0=p['T0'] * a, T1=p['T1'] * a, S=S, rel=dict(T0=p['T0'], T1=p['T1']))
    if method == 'SCA':
        return dict(family='plain', q=p['q'] * a, T0=p['T0'] * a, tfin=tf * a, S=S, rel=dict(q=p['q'], T0=p['T0'], tfin=tf))
    if method == 'TEC':
        return dict(family='tec', jv=p['jv'] * a, q=p['q'] * a, T0=p['T0'] * a, tfin=tf * a, S=S,
                    rel=dict(jv=p['jv'], q=p['q'], T0=p['T0'], tfin=tf))
    if method == 'Onsager-kT':
        return dict(family='tecT', kappa=p['kappa'] * kf, q=p['q'] * a, T0=p['T0'] * a, ramp=True, tfin=tf * a, S=S,
                    rel=dict(kappa=p['kappa'], kappa_factor=kf, q=p['q'], T0=p['T0'], tfin=tf))
    if method == 'Onsager-online':
        return dict(family='onsager', lam=p['lam'] * c_lamf, q=p['q'] * a, T0=p['T0'] * a, ramp=True, tfin=tf * a, S=S,
                    rel=dict(lam_k2000=p['lam'], lam_factor=c_lamf, q=p['q'], T0=p['T0'], tfin=tf))
    if method == 'APC-SCA':
        return dict(family='apc', q_reset=p['q_reset'] * a, r_q=p['r_q'], q_lim=p['q_lim'] * a, T0=p['T0'] * a,
                    tfin=tf * a, S=S, rel=dict(q_reset=p['q_reset'], r_q=p['r_q'], q_lim=p['q_lim'], T0=p['T0'], tfin=tf))
    if method == 'ReAIM ASA':
        ks = p['kset']
        if prob == 'mcp':
            kk = tuple(int(min(N, max(1, round(k * N / 2000)))) for k in ks)
            rel = dict(kset_n2000=tuple(ks), T1=p['T1'])
        else:
            kk = tuple(int(min(N, k)) for k in ks)
            rel = dict(kset=tuple(ks), T1=p['T1'])
        rel.update(F=p['F'], T0=p['T0'])
        return dict(family='ReAIM', kset=kk, T1=p['T1'], T0=p['T0'], F=p['F'], S=S, rel=rel)
    if method in ('aSB', 'bSB', 'dSB'):
        return dict(family=method, dt=p['dt'], xi=p['xi'], S=S)
    raise ValueError(method)


def point_from_cfg(prob, method, cfg):
    """Relative point of an A1 selected configuration (for the round-1 trigger)."""
    r = cfg.get('rel', {})
    if method == 'SA':
        return dict(T0=r['T0'], T1=r['T1'])
    if method == 'SCA':
        return dict(q=r['q'], T0=r['T0'], tfin=r.get('tfin'))
    if method == 'TEC':
        return dict(jv=r['jv'], q=r['q'], T0=r['T0'], tfin=r.get('tfin'))
    if method == 'Onsager-kT':
        return dict(kappa=r['kappa'], q=r['q'], T0=r['T0'], tfin=r.get('tfin'))
    if method == 'Onsager-online':
        return dict(lam=r['lam_k2000'], q=r['q'], T0=r['T0'], tfin=r.get('tfin'))
    if method == 'APC-SCA':
        return dict(q_reset=r['q_reset'], r_q=r['r_q'], q_lim=r['q_lim'], T0=r['T0'], tfin=r.get('tfin'))
    if method == 'ReAIM ASA':
        ks = r.get('kset_n2000', r.get('kset'))
        return dict(kset=list(ks), T1=r['T1'], F=r.get('F', 'max'), T0=r.get('T0', 1.0))
    return dict(dt=cfg['dt'], xi=cfg['xi'])
