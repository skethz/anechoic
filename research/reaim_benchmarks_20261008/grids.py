"""Tuning grids of PROTOCOL.md: 32 points for every method on GPP and TSP. Field-valued parameters (T0, T_fin, T1 of SA,
q, J_v, q_reset, q_lim) are in units of the instance's sigma_T (GPP: sqrt(mean graph degree); TSP: A = max distance).
kappa is in units of kappa_factor = mean_i sum_j J_ij^2 / sigma_T^2 and lambda in units of
lam_factor = (mean_i sum_j J_ij^2 / N) / (1999/2000) (K2000 units), i.e. kappa_eff = kappa * kappa_factor,
lambda_eff = lambda * lam_factor. SB parameters are dimensionless (Goto normalisation with sigma_J of the full J).
ReAIM's temperature is relative (T: 1 -> T1) and its k values are absolute flip caps.
Placed by the exploratory synthetic-instance calibration (calib.py, calib/ and logs/calib_*.log); see PROTOCOL.md."""
import itertools

P_ = itertools.product
RAMP = True
REAIM_KSETS = ((1, 1, 2, 2), (1, 2, 3, 4), (1, 2, 4, 8), (2, 4, 8, 16))
REAIM_T1 = (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125)

G = {
    'gpp': dict(
        SA=dict(T0=(0.25, 0.5, 1.0, 2.0), T1=(0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.48, 0.64)),
        SCA=dict(q=(4.0, 8.0, 16.0, 32.0), T0=(0.5, 1.0, 2.0, 4.0), tfin=(0.15, 0.6)),
        shared=dict(q=(8.0, 32.0), T0=(1.0, 4.0), tfin=(0.15, 0.6)),
        TEC=dict(jv=(-16.0, -8.0, -4.0, -1.0)),
        kT=dict(kappa=(2.0 ** -10, 2.0 ** -8, 2.0 ** -6, 2.0 ** -4)),
        online=dict(lam=(2.0 ** -8, 2.0 ** -5, 2.0 ** -2, 2.0)),
        APC=dict(q_reset=(64.0, 128.0, 256.0, 512.0), r_q=(0.97, 0.99), q_lim=(2.0, 8.0), T0=(1.0, 4.0), tfin=0.15),
        bSB=dict(dt=(0.03125, 0.0625, 0.125, 0.25), xi=(6.0, 8.0, 11.0, 16.0, 22.0, 32.0, 45.0, 64.0)),
        dSB=dict(dt=(0.03125, 0.0625, 0.125, 0.25), xi=(0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)),
        aSB=dict(dt=(0.03125, 0.0625, 0.125, 0.25), xi=(1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0, 32.0)),
    ),
    'tsp': dict(
        SA=dict(T0=(0.5, 1.0, 2.0, 4.0), T1=(0.04, 0.08, 0.16, 0.24, 0.32, 0.48, 0.64, 0.96)),
        SCA=dict(q=(8.0, 12.0, 16.0, 24.0), T0=(8.0, 16.0, 32.0, 64.0), tfin=(0.15, 0.6)),
        shared=dict(q=(12.0, 16.0), T0=(16.0, 32.0), tfin=(0.15, 0.6)),
        TEC=dict(jv=(-16.0, -8.0, -4.0, -2.0)),
        kT=dict(kappa=(2.0 ** -10, 2.0 ** -8, 2.0 ** -6, 2.0 ** -4)),
        online=dict(lam=(2.0 ** -8, 2.0 ** -6, 2.0 ** -4, 2.0 ** -2)),
        APC=dict(q_reset=(32.0, 64.0, 256.0, 512.0), r_q=(0.97, 0.99), q_lim=(0.125, 0.5), T0=(1.0, 4.0), tfin=0.15),
        bSB=dict(dt=(0.125, 0.25, 0.5, 0.75), xi=(0.5, 0.71, 1.0, 1.41, 2.0, 2.83, 4.0, 8.0)),
        dSB=dict(dt=(0.0625, 0.125, 0.25, 0.5), xi=(1.0, 1.41, 2.0, 2.83, 4.0, 5.66, 8.0, 16.0)),
        aSB=dict(dt=(0.03125, 0.0625, 0.125, 0.25), xi=(0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)),
    ),
}


def grid(prob, method, P, S):
    g = G[prob]; a = P.sigma_T; st = P.stats(); N = P.N; out = []
    sh = g['shared']
    if method == 'SA':
        for t0, t1 in P_(g['SA']['T0'], g['SA']['T1']):
            out.append(dict(family='SA', T0=t0 * a, T1=t1 * a, S=S, rel=dict(T0=t0, T1=t1)))
    elif method == 'SCA':
        c = g['SCA']
        for q, t0, tf in P_(c['q'], c['T0'], c['tfin']):
            out.append(dict(family='plain', q=q * a, T0=t0 * a, tfin=tf * a, S=S, rel=dict(q=q, T0=t0, tfin=tf)))
    elif method == 'TEC':
        for j, q, t0, tf in P_(g['TEC']['jv'], sh['q'], sh['T0'], sh['tfin']):
            out.append(dict(family='tec', jv=j * a, q=q * a, T0=t0 * a, tfin=tf * a, S=S,
                            rel=dict(jv=j, q=q, T0=t0, tfin=tf)))
    elif method == 'Onsager-kT':
        for k, q, t0, tf in P_(g['kT']['kappa'], sh['q'], sh['T0'], sh['tfin']):
            out.append(dict(family='tecT', kappa=k * st['kappa_factor'], q=q * a, T0=t0 * a, ramp=RAMP, tfin=tf * a,
                            S=S, rel=dict(kappa=k, kappa_factor=st['kappa_factor'], q=q, T0=t0, tfin=tf)))
    elif method == 'Onsager-online':
        for lam, q, t0, tf in P_(g['online']['lam'], sh['q'], sh['T0'], sh['tfin']):
            out.append(dict(family='onsager', lam=lam * st['lam_factor'], q=q * a, T0=t0 * a, ramp=RAMP, tfin=tf * a,
                            S=S, rel=dict(lam_k2000=lam, lam_factor=st['lam_factor'], q=q, T0=t0, tfin=tf)))
    elif method == 'APC-SCA':
        c = g['APC']
        for x, r, y, t0 in P_(c['q_reset'], c['r_q'], c['q_lim'], c['T0']):
            out.append(dict(family='apc', q_reset=x * a, r_q=r, q_lim=y * a, T0=t0 * a, tfin=c['tfin'] * a, S=S,
                            rel=dict(q_reset=x, r_q=r, q_lim=y, T0=t0, tfin=c['tfin'])))
    elif method == 'ReAIM ASA':
        for ks, t1 in P_(REAIM_KSETS, REAIM_T1):
            out.append(dict(family='ReAIM', kset=tuple(int(min(N, k)) for k in ks), T1=t1, S=S, rel=dict(kset=ks, T1=t1)))
    elif method in ('aSB', 'bSB', 'dSB'):
        c = g[method]
        for d, x in P_(c['dt'], c['xi']):
            out.append(dict(family=method, dt=d, xi=x, S=S))
    assert len(out) == 32, (prob, method, len(out))
    return out
