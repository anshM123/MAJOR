"""
dc_sens (DIRCV.md §3.2): linear response of the fitted r (F1, F1B) to a relative change (1 + eps) of the STATE record at one node,
through the full pipeline (a_h refit, bases, de-shift, fit).  r_dir - r_cv = sum_nu S_nu . d_nu (d_nu = scalar MC log deviation).
"""
import dc_guard, os, sys
import numpy as np
import z3_nodes as Z3, fw_gamma as FG, pw_gamma as PG
import dc_pipe as DP, dc_fit as F

J, D = 100, 1100.0
N4 = [13.3, 16.3, 20.8, 25.0]
EPS = {}
_tr = Z3.to_record


def to_record(r, DD, JJ, fac=1.0, sel=None):
    e = EPS.get((float(DD), round(float(r.get('nu', -1)), 6)), 0.0)
    return _tr(r, DD, JJ, fac*(1.0 + e), sel)


Z3.to_record = to_record


def rfit(est, tags, kf, nus, model):
    FG.EST = est; PG.TAGS[:] = list(tags); PG.KEEPFW = kf; PG.STATE[J] = D
    g = FG.row(J, D); m = g['m']; xs = g['xs']
    o = dict(J=J, D=D, m=m, nus=xs*m)
    out = {}
    for h in ('++', '+-'):
        img = FG.CF.curv_image('C2', h, J, m, xs, FG.Q.A.env, FG.Q.A.SIG)
        y = np.asarray(g[('C2', h)], float) - img
        o[('C2', h)] = dict(y=y, yj=y[None, :].repeat(2, 0))
        out[h] = F.fit(o, h, nus, model)['r']
    return out


def sens(est, tags, kf, nus, model, h=1e-4):
    EPS.clear(); r0 = rfit(est, tags, kf, nus, model)
    S = {}
    for nu in nus:
        row = {}
        for comp, e in (('re', h), ('im', 1j*h)):
            EPS.clear(); EPS[(D, round(nu, 6))] = e
            r1 = rfit(est, tags, kf, nus, model)
            row[comp] = {hh: (r1[hh] - r0[hh])/h for hh in r0}
        S[nu] = row
    EPS.clear()
    return r0, S


if __name__ == '__main__':
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    import pickle
    w = pickle.load(open(os.path.join(dc_guard.SPIN, 'cache', 'pw_pts_J100_D1100_w.pkl'), 'rb'))
    import fw_scalar as FS
    FS.install(J, D)
    dev = {}
    for nu in N4:
        sw, rf = PG.to_z2(w[(D, nu)], 'EE'), PG.to_z2(w[(Z3.DREF, nu)], 'EE')
        dev[nu] = Z3.scalar_check(J, nu, {D: sw, Z3.DREF: rf}, D)['d']
    res = {}
    for model, kf, nus in (('F1', False, N4), ('F1B', True, [13.3, 16.3, 16.7, 20.8, 25.0, 33.4])):
        r0, S = sens('dir', ['_w'], kf, nus, model)
        res[model] = S
        pr(f"# {model} (dir, _w, nodes {nus}): r0 = {r0};  S_nu = dr/d(Re eps), dr/d(Im eps)  [c2 per unit relative change]")
        for nu in nus:
            pr(f"   nu={nu:5.1f}: " + "  ".join(f"{hh}: {S[nu]['re'][hh]:+.3e}, {S[nu]['im'][hh]:+.3e}" for hh in ('++', '+-')))
        # predicted r_dir - r_cv from the scalar deviations d_nu (cv multiplies the record by exp(-d))
        for hh in ('++', '+-'):
            pred = sum(S[nu]['re'][hh]*dev[nu].real + S[nu]['im'][hh]*dev[nu].imag for nu in N4 if nu in S)
            pr(f"   predicted r_dir - r_cv ({hh}) from the scalar deviations d (pw nodes) = {pred:+.0f};  "
               f"per node: " + ", ".join(f"{nu}: {S[nu]['re'][hh]*dev[nu].real + S[nu]['im'][hh]*dev[nu].imag:+.0f}" for nu in N4 if nu in S))
    pr("   scalar deviations d (waves, seed 4242): " + ", ".join(f"{nu}: {dev[nu].real:+.5f}{dev[nu].imag:+.5f}i" for nu in N4))
    pickle.dump(dict(S=res, dev=dev), open(os.path.join(dc_guard.CACHE, 'dc_sens.pkl'), 'wb'))
    open(os.path.join(dc_guard.CACHE, 'dc_sens.txt'), 'w').write('\n'.join(lines) + '\n')
