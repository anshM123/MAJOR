"""
gm_cert8 (GAP_MORSE.md §4.4): Taylor-model data valid for a whole xi-INTERVAL [xi0 - rx, xi0 + rx] (inside the complex disc
|delta| <= rx, delta = xi - xi0 the 8th Taylor-model variable of gm_ball8), in the form consumed by gm_chart (Lemma M), uniform in delta:
  * crystal z*(xi) exact as a Taylor model in delta (a = xi/(1 + sqrt(1 - xi^2))); phase coefficient i xi/2 a Taylor model; the Morse
    matrix L is the float one at xi0 (so the Hessian is -(1 + E'(delta)), E' = O(delta), carried in EF).
  * v-degree-k Bombieri norms sup_delta [P_k(.; delta)] <= sum_j [p_{k,j}] rx^j; the Taylor-model remainder R (sup <= rem on the product
    domain, vanishing to TOTAL order D+1) is split at fixed delta into its v-homogeneous parts R_k: k = 0, 1 are dropped exactly (the
    constant is irrelevant; P_1 + R_1 = 0 since the crystal is exact for every xi), k = 2..6 go into the chart polynomial with
    [R_k]_B <= sqrt(C(k+6,6)) rem/rho^k ([p]_B^2 = dim H_k x mean of |p|^2 on the sphere), and R~ = R_{>=7} stays in the exponent:
    |R~| <= 8 rem (|v|/rho)^7.
  * amplitude: sup bound sum_k (sum_j [h_{k,j}] rx^j) t^k + remh (flat remainder), divided by lam_min = min over the disc of
    env(xi)/env(xi0) (so the bound is relative to the envelope at each xi).
Usage: GC_DEG=8 python3 gm_cert8.py J xi0 rx rho  -> cache/gm_int_J<J>.json
"""
import os, sys, time, json, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gm_guard  # noqa: F401
import numpy as np
import gm_ball8 as B
import gm_scalar as GS
from flint import acb, arb
sys.path.insert(0, os.path.join(HERE, '..', 'gapc3'))
import gc_node as N

DREF = 1002.0


def up(z):
    return float(abs(z).upper()) if isinstance(z, acb) else float(z.upper())


def lap_part(c, k):
    """coefficients of Lap_v(p) for the part of v-degree k of p (8-variable monomials; Laplacian in the first 7)."""
    idx = {m: i for i, m in enumerate(B.MONS)}
    out = B.zeros(B.NM)
    for a, m in enumerate(B.MONS):
        if sum(m[:7]) != k:
            continue
        for i in range(7):
            if m[i] >= 2:
                mm = list(m); mm[i] -= 2
                out[idx[tuple(mm)]] = out[idx[tuple(mm)]] + c[a]*(m[i]*(m[i] - 1))
    return out


def cert(J, xi0, rx, rho, Delta=DREF):
    t0 = time.time()
    B.BAD.clear(); B.WMAX[0] = 0.0; B.RHO[0] = rho; B.RHOX[0] = rx
    D = B.DEG
    nu = xi0*Delta
    zs, a, L = N.morse(nu, Delta)
    LP = L*np.sqrt(Delta)
    xi = B.const(acb(arb(nu)/arb(Delta))) + B.var(7)
    one = B.const(1.0)
    sq = B.powr(one - B.mul(xi, xi), 0.5)
    at = B.mul(xi, B.recip(one + sq))
    I = acb(0, 1)
    zc = [B.const(0.0), -I*at, I*at, I*at, -I*at, B.const(0.0), B.const(0.0)]
    zv = B.linsub_ball(zc, LP)
    K = GS.TMB(B)
    snap = {}

    def mark():
        snap['bad'] = len(B.BAD); snap['q'] = B.WMAX[0]
    r = GS.build(zv, nu, Delta, J, K, at[0], xi=xi, mark=mark)
    NM = B.NM
    dPhi = r['dPhi']
    rem = up(dPhi[NM])
    kv = np.array([sum(m[:7]) for m in B.MONS])
    # gradient: every coefficient of v-degree 1 must enclose 0
    lin = [dPhi[i] for i in range(NM) if kv[i] == 1]
    grad_ok = all(z.contains(0) for z in lin)
    gradmax = max(up(z) for z in lin)
    # E'(delta): v-degree-2 part + |v|^2/2  ->  [ -v.E'v/2 ]_B = ||E'||_F / 2
    p2 = dPhi[:NM].copy()
    for i in range(7):
        e = [0]*8; e[i] = 2
        k = B.MONS.index(tuple(e))
        p2[k] = p2[k] + acb(0.5)
    p2[kv != 2] = B.ZERO
    rx_ = arb(rx); r_ = arb(rho)
    C = lambda k: arb(math.comb(k + 6, 6)).sqrt()
    vb2 = B.vbounds(p2[None])[0]
    EF = 2*vb2[2] + 2*C(2)*arb(rem)/r_**2
    trE = arb(7).sqrt()*EF
    vb = B.vbounds(dPhi[:NM][None])[0]
    bPhi = [0.0, 0.0, 1.3228756555322962] + [up(vb[k] + (C(k)*arb(rem)/r_**k if k <= 6 else 0)) for k in range(3, D + 1)]
    bTr = [0.0, 0.0, 0.0]
    for k in range(3, D + 1):
        lp = lap_part(dPhi[:NM], k)/acb(k)
        v = B.vbounds(lp[None])[0][k - 2]
        if k <= 6:
            v = v + (k - 1)*arb(7).sqrt()*C(k)*arb(rem)/r_**k
        bTr.append(up(v))
    # amplitude and envelope ratio
    if J > 0:
        h = r['hJ']
        envt = B.mul(r['grel'], B.powr(r['bJ'], J/2.0))       # env(xi)/env(xi0) on the v = 0 slice (modulus)
    else:
        h = r['grel']; envt = r['grel']
    vbh = B.vbounds(h[:NM][None])[0]
    e0 = envt[:NM].copy(); e0[kv != 0] = B.ZERO
    bge = B.bigraded(e0[None])[0]
    lam = abs(envt[0]) - sum(bge[0, j]*rx_**j for j in range(1, D + 1)) - abs(envt[NM])
    lam_lo = float(lam.lower())
    out = dict(J=J, xi0=xi0, xi=xi0, rx=rx, rho=rho, deg=D, q=B.WMAX[0], bad=len(B.BAD), qPhi=snap['q'], badPhi=snap['bad'],
               grad_ok=bool(grad_ok), gradmax=gradmax, EF=up(EF), trE=up(trE), bPhi=bPhi, remPhi=8*rem, ordR=7, bTr=bTr,
               bh=[up(vbh[k])/lam_lo if lam_lo > 0 else float('inf') for k in range(D + 1)],
               remh=up(h[NM])/lam_lo if lam_lo > 0 else float('inf'), remh_flat=True, lam_min=lam_lo, rem_raw=rem,
               t=time.time() - t0)
    return out


if __name__ == '__main__':
    J, xi0, rx = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    rhos = [float(a) for a in sys.argv[4].split(',')]
    fn = os.path.join(HERE, 'cache', 'gm_int_J%d.json' % J)
    for rho in rhos:
        key = '%.4f_%.4f_%.4f_%d' % (xi0, rx, rho, B.DEG)
        db = json.load(open(fn)) if os.path.exists(fn) else {}
        if key in db:
            continue
        try:
            o = cert(J, xi0, rx, rho)
        except Exception as ex:
            o = dict(J=J, xi0=xi0, rx=rx, rho=rho, deg=B.DEG, bad=-1, err=repr(ex)[:300])
        db = json.load(open(fn)) if os.path.exists(fn) else {}
        db[key] = o
        json.dump(db, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
        print(json.dumps({k: o.get(k) for k in ('J', 'xi0', 'rx', 'rho', 'bad', 'badPhi', 'grad_ok', 'EF', 'remPhi', 'remh', 'lam_min',
                                                'bPhi', 'bh', 't', 'err')}), flush=True)
