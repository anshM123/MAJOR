"""
k = 2 (the C+_2 moment) under the even kernels S_ff, K* (tensor weights, t18) and the leading modes, both helicities: diagnosis and a
cleaned extraction (SPINNING_TTTT_EXT.md sect. 7).  Existing caches only.
Diagnosis tools: (i) the noise floor of the cached modes, sigma(node) = median over modes of |Im G_q| (the modes used are real by symmetry;
     their imaginary parts are pure numerical noise), smoothed over neighbouring t-nodes; (ii) the k-weight (rho rhobar)^{(k-4)/2} =
     S^{2k-8} u^{k-4}, which enhances the tiny-t corner (S -> 0) by S^-4 and the light-cone wedge (u -> 0) by u^-2 at k = 2 relative to k = 4.
Cleaned extraction (no fit to the prediction):
  1. per u-row, the first t-node i0 with three consecutive accepted nodes whose |Re d| > KAPPA x (propagated noise): nodes before i0 are
     replaced by a power law in t fitted to nodes i0..i0+2 (exponent clipped to [0, 40])  ["row tail"];
  2. per t-column with t >= T_W, the dropped light-cone wedge u < u_first (radial series not converged) is filled with H(u) = d/F extrapolated
     as a power law in u from the first three accepted nodes, times the exact kernel F(xi(u))  ["wedge"];
  3. R = (clean accepted + row tail + wedge)/prediction.  Variants (fit nodes shifted by one, two-node fits, KAPPA = 3 and 30, constant
     extrapolation) give the model spread; constant tails and the looser production-level mask (cut 32) are added to form the envelope.
Output: out_t19_k2.txt
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from t13_analysis_l2 import Modes, design, Tpred, CLASSES
from predict4 import weyl
from moments4 import kernel_np
from moments import frame_data
import t5_analysis_l1 as T5
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

KS = {'lead': None,
      'Sff': {(0, 0): 1/72, (2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/32, (-1, 1): -1/32},
      'K*': {(2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/64, (-1, 1): -1/64},
      'A': {(2, -2): 1/(4j), (-2, 2): -1/(4j)}}
LEAD = {('same', 'same'): (2, -2, -2, 2), ('opp', 'opp'): (2, -2, 2, -2)}
FILES = ['cache/modes_l2_D45.861045_J20_d.npz', 'cache/modes_l2_D80.484300_J36_d.npz', 'cache/modes_l2_D115.107554_J52_d.npz']
T_W = 3.0
CUT = 40.0

_X = {}
def crit(mo):
    """X - 10 log X of e2e/moments.series_converged at every grid point (the production mask is crit > 25)"""
    if id(mo) not in _X:
        g = mo.g; N = int(np.load(mo.fn)['N']) if hasattr(mo, 'fn') else 8000
        fd = frame_data(g['U'], np.minimum(g['S'], 0.999))
        X = -N*(np.log(fd['rs']) + fd['chi']) - mo.J*fd['chi']
        with np.errstate(divide='ignore', invalid='ignore'):
            c = X - 10.0*np.log(np.maximum(X, 1.0))
        c[X <= 8] = -np.inf
        _X[id(mo)] = c
    return _X[id(mo)]

_NOISE = {}
def noise(mo, wave):
    key = (id(mo), wave)
    if key not in _NOISE:
        im = np.array([np.abs(v.imag) for v in mo.modes[wave].values()])
        im = im[im.reshape(len(im), -1).max(1) > 0]
        s = np.median(im, axis=0)
        sm = s.copy(); sm[:, 1:] = np.maximum(sm[:, 1:], s[:, :-1]); sm[:, :-1] = np.maximum(sm[:, :-1], s[:, 1:])
        _NOISE[key] = sm
    return _NOISE[key]

_SET = {}
def setup(mo, k, nu, wave, K):
    key = (id(mo), k, nu, wave, id(K))
    if key not in _SET: _SET[key] = _setup(mo, k, nu, wave, K)
    return _SET[key]

def _setup(mo, k, nu, wave, K):
    g = mo.g; U, S = g['U'], g['S']; xi = -np.log(g['us'])
    lst = [(LEAD[wave], 1.0)] if K is None else [((a, b, c, d), w1*w2) for (a, b), w1 in K.items() for (c, d), w2 in K.items()]
    mus = {(q[0] + q[1], q[2] + q[3]) for q, w in lst}; assert mus == {(0, 0)}
    F = kernel_np(nu, [(0, 0)], xi)[(0, 0)]
    with np.errstate(divide='ignore', invalid='ignore'):
        logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(S**2 - (S*U)**2) + 0.5*(k - 4)*np.log(S**4*U**2) + mo.logscale
    logw[~g['ok']] = -np.inf
    sh = np.max(np.where(mo.conv, logw, -np.inf))
    wgt = np.where(np.isfinite(logw), np.exp(logw - sh), 0.0)
    GK = sum(w*mo.modes[wave][q] for q, w in lst).real
    A = design(mo.Delta, mo.J, k, nu); Tp = Tpred(wave[0], wave[1], mo.Delta, mo.J, 8.0)
    p = sum(w*sum(A[q][c]*Tp[c] for c in CLASSES) for q, w in lst)
    f = (-1)**mo.J*np.exp(sh - np.log(abs(p)))/(p/abs(p))
    H = wgt*GK                                              # density without the kernel F and the normalisation (real)
    Fc = (F*f)                                              # complex factor per u-row; density d = Re(Fc H)
    d = (Fc[:, None]*H).real
    sig = wgt*np.abs(Fc[:, None])*noise(mo, wave)*np.sqrt(sum(w*w for q, w in lst))
    return d, sig, H, Fc

def smooth_start(x, y, good, lo=0, need=3):
    """first index i >= lo with need consecutive good nodes of one sign whose successive log-slopes agree (|dp| < 1 + 0.3 max|p|)"""
    n = len(x)
    for i in range(lo, n - need + 1):
        idx = np.arange(i, i + need)
        if not np.all(good[idx]): continue
        yy = y[idx]
        if not (np.all(yy > 0) or np.all(yy < 0)): continue
        ps = np.diff(np.log(np.abs(yy)))/np.diff(np.log(x[idx]))
        if np.all(np.abs(np.diff(ps)) < 1.0 + 0.3*np.max(np.abs(ps))): return i
    return None

def powfit(x, y, idx, xfill, clip):
    yy = y[idx]
    pp = np.polyfit(np.log(x[idx]), np.log(np.abs(yy)), 1)[0] if len(idx) > 1 else 0.0
    pp = min(max(pp, clip[0]), clip[1])
    return yy[0]*(xfill/x[idx[0]])**pp, pp

def cleaned(mo, k, nu, wave, K, kappa=10.0, shift=0, nfit=3, const=False, bound=False, tw=None, cut=None):
    tw = T_W if tw is None else tw; cut = CUT if cut is None else cut
    g = mo.g; W = g['W']; ts, us = g['ts'], g['us']; ok = g['ok']
    d, sig, H, Fc = setup(mo, k, nu, wave, K)
    conv = mo.conv & ok & (crit(mo) > cut)
    raw = float((d*W*mo.conv*ok).sum())
    good = conv & (np.abs(d) > kappa*sig)
    Hf = np.where(conv, H, np.nan)                         # working copy of H (wedge-filled)
    twedge = 0.0; tdisc = 0.0
    # light-cone wedge, columns t >= tw: discard the region before the smooth start, extrapolate H in u
    for it in np.where(ts >= tw)[0]:
        i0 = smooth_start(us, H[:, it], good[:, it])
        if i0 is None or i0 == 0: continue
        fill = np.arange(i0); fill = fill[ok[fill, it]]
        if len(fill) == 0: continue
        idx = np.arange(i0 + shift, min(i0 + shift + nfit, len(us)))
        if const or bound: hv = np.full(len(fill), H[i0, it])
        else: hv, _ = powfit(us, H[:, it], idx, us[fill], (-1.49, 30.0))
        dv = np.abs(Fc[fill]*hv) if bound else (Fc[fill]*hv).real
        tdisc += np.sum(np.where(conv[fill, it], d[fill, it], 0.0)*W[fill, it])
        twedge += np.sum(dv*W[fill, it]); Hf[fill, it] = hv
    dF = (Fc[:, None]*np.nan_to_num(Hf)).real
    filled = ~np.isnan(Hf)
    goodF = filled & ((np.abs(dF) > kappa*sig) | (~conv & filled))
    tot = 0.0; trow = 0.0
    for iu in range(len(us)):
        i0 = smooth_start(ts, dF[iu], goodF[iu])
        if i0 is None: continue
        keep = np.arange(i0, len(ts)); keep = keep[filled[iu, keep]]
        tot += np.sum(dF[iu, keep]*W[iu, keep])
        fill = np.arange(i0); fill = fill[ok[iu, fill]]
        if len(fill) == 0: continue
        tdisc += np.sum(np.where(conv[iu, fill] & ~(ts[fill] >= tw), d[iu, fill], 0.0)*W[iu, fill])
        idx = np.arange(i0 + shift, min(i0 + shift + nfit, len(ts)))
        if const or bound: yv = np.full(len(fill), abs(dF[iu, i0]) if bound else dF[iu, i0])
        else: yv, _ = powfit(ts, dF[iu], idx, ts[fill], (0.0, 40.0))
        trow += np.sum(yv*W[iu, fill])
    return dict(raw=raw, R=tot + trow, row=trow, wedge=twedge, disc=tdisc)

VARIANTS = {'nominal': {}, 'shift1': dict(shift=1), 'nfit2': dict(nfit=2), 'kappa3': dict(kappa=3.0), 'kappa30': dict(kappa=30.0),
            'tw5': dict(tw=5.0), 'cut50': dict(cut=50.0)}
ENVELOPE = {'const': dict(const=True), 'cut32': dict(cut=32.0)}     # constant tails (ignores the measured slopes); looser mask

def analyse(mos, ms, k, nu, wave, kn):
    res = {v: [cleaned(mo, k, nu, wave, KS[kn], **kw) for mo in mos] for v, kw in list(VARIANTS.items()) + list(ENVELOPE.items())}
    Rn = [r['R'] for r in res['nominal']]
    spread = [max(abs(res[v][i]['R'] - Rn[i]) for v in VARIANTS) for i in range(len(mos))]
    env = [(min(res[v][i]['R'] for v in res), max(res[v][i]['R'] for v in res)) for i in range(len(mos))]
    return res, env, Rn, spread

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    mos = [Modes(os.path.join(HERE, f)) for f in FILES]
    for mo, f in zip(mos, FILES): mo.fn = os.path.join(HERE, f)
    ms = [float(weyl(mo.Delta, mo.J)[2]) for mo in mos]
    pr("# k = 2 (C+_2) under S_ff, K* and the leading modes: cleaned extraction (fourpt/tttt/t19_k2.py)")
    pr(f"# beta = 1, m = " + ", ".join(f"{m:.1f}" for m in ms) + f"; T_W = {T_W}; variants {list(VARIANTS)}")
    for mo in mos:
        c = crit(mo); agree = np.mean((c > 25.0)[mo.g['ok']] == mo.conv[mo.g['ok']])
        pr(f"  mask check m = {float(weyl(mo.Delta, mo.J)[2]):.1f}: crit > 25 reproduces the production mask on {100*agree:.2f} % of the grid; nodes kept at cut 25/32/40/50: "
           + "/".join(str(int(((c > cc) & mo.g['ok']).sum())) for cc in (25, 32, 40, 50)))
    pr("\n# (a) noise floor sigma (median |Im G_q|, mode units) against the leading-mode signal, at the smallest accepted t of three rows (m = 100.7)")
    mo = mos[-1]; g = mo.g
    for wave in (('same', 'same'), ('opp', 'opp')):
        s = noise(mo, wave); Gl = mo.modes[wave][LEAD[wave]].real
        for iu in (40, 55, 65):
            acc = np.where(mo.conv[iu])[0][:6]
            pr(f"  {wave[0]} u={g['us'][iu]:.3f}: " + " ".join(f"t{g['ts'][i]:.1f}: sig {s[iu, i]:.1e} lead {Gl[iu, i]:+.1e} iso {mo.modes[wave][(0,0,0,0)][iu, i].real:+.1e}" for i in acc))
    pr("\n# (b) raw ratio versus a hard tiny-t cut t >= t_c (t_c = 0, 2, 3, 4, 5), k = 2, at m = 39 | 70 | 101")
    for nu in (0.0, 3.0):
        for wave in (('same', 'same'), ('opp', 'opp')):
            for kn in KS:
                cells = []
                for mo in mos:
                    d, sig, H, F = setup(mo, 2, nu, wave, KS[kn])
                    cells.append("/".join(f"{(d*mo.g['W']*(mo.conv & (mo.g['ts'][None, :] >= tc))).sum():.3f}" for tc in (0, 2, 3, 4, 5)))
                pr(f"  nu={nu} {wave[0]} {kn}: " + " | ".join(cells))
    ks = (2, 4, 6) if which == 'all' else (2,)
    pr("\n# (c) cleaned ratios R (nominal) [model spread: max deviation over the variants] {envelope: min..max including constant tails and the looser mask}; parts: row tail, wedge; m^2(R-1); c1 of the fit c1/m + c2/m^2 + c4/m^4")
    for k in ks:
        for nu in (0.0, 1.0, 3.0):
            for wave in (('same', 'same'), ('opp', 'opp')):
                for kn in KS:
                    if k > 2 and kn == 'A': continue
                    res, env, Rn, sp = analyse(mos, ms, k, nu, wave, kn)
                    f = T5.fits(ms, Rn)
                    pr(f"  k={k} nu={nu} {wave[0]} {kn}: R " + "/".join(f"{R:.4f}" for R in Rn) + " [" + "/".join(f"{s:.4f}" for s in sp) + "] {"
                       + "/".join(f"{a:.3f}..{b:.3f}" for a, b in env) + "}; raw " + "/".join(f"{r['raw']:.3f}" for r in res['nominal'])
                       + "; row " + "/".join(f"{r['row']:+.4f}" for r in res['nominal']) + "; wedge " + "/".join(f"{r['wedge']:+.4f}" for r in res['nominal'])
                       + "; discarded " + "/".join(f"{r['disc']:+.3f}" for r in res['nominal'])
                       + "; m2d " + "/".join(f"{m*m*(R-1):+.0f}" for m, R in zip(ms, Rn)) + f"; c1 {f['A'][0]:+.2f}")
    pr("\n# (d) beta = 0.4, m = 72.0 (J = 14), one mass, k = 2: cleaned R [model spread] {envelope}; raw")
    mo4 = Modes(os.path.join(HERE, 'cache', 'modes_l2_D74.964099_J14_d.npz')); mo4.fn = os.path.join(HERE, 'cache', 'modes_l2_D74.964099_J14_d.npz')
    for nu in (0.0, 1.0, 3.0):
        for wave in (('same', 'same'), ('opp', 'opp')):
            cells = []
            for kn in KS:
                res, env, Rn, sp = analyse([mo4], [72.0], 2, nu, wave, kn)
                cells.append(f"{kn} {Rn[0]:.4f} [{sp[0]:.4f}] {{{env[0][0]:.3f}..{env[0][1]:.3f}}} raw {res['nominal'][0]['raw']:.3f}")
            pr(f"  k=2 nu={nu} {wave[0]}: " + " | ".join(cells))
    open(os.path.join(HERE, 'out_t19_k2.txt' if which == 'all' else 'cache/t19_k2only.txt'), 'w').write("\n".join(out) + "\n")
