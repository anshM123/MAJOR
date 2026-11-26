"""
ct_check.py [DV]: exact nu-contour ("Regge") decomposition of the ideal-weight action at finite twist.
  L(G_A^psi) = (2/Lam) int_0^Lam psi2(nu/Lam) K(nu) dnu,  K(nu) = L(F_{3,nu}) (exact terminating ladder),
  = (2/Lam) { int_0^delta psi K dnu + 2 Re[ int_{delta}^{delta+iY} - int_{Lam}^{Lam+iY} + int_{delta+iY}^{Lam+iY} ] psi kappa dnu },
kappa = finite-mass outgoing kernel (kplus.Kplus); infinite twist: kappa_inf = outgoing (1-z)^{-3}P+ at c.
Prints, per piece: Re value (fm, inf), (2/Lam) int|kappa_inf psi|, (2/Lam) int|(kappa_fm - kappa_inf) psi|, max|kappa_fm/kappa_inf - 1|.
Usage: python3 ct_check.py J a b Y [ntop] [prec]
"""
import sys, time, json, numpy as np
from an_common import *
from kplus import *
J, a, b = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]); Y = float(sys.argv[4])
NT = int(sys.argv[5]) if len(sys.argv) > 5 else 120
ctx.prec = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
st = state(J, a, b); U = ladder(st)
beta = float((2*(st['c'].sqrt()).atanh()).mid()); d = float(st['d'].mid())
LAM = arb(LAM2).sqrt(); Lf = float(LAM.mid())
p2, p3 = nf02()
P2 = [acb(arb(x.numerator)/x.denominator) for x in p2]
def psi(p):
    v = acb(0)
    for co in reversed(P2): v = v*p + co
    return v
k = 2/Lf
print(f'state J={J} a={a} b={b}: c={float(st["c"].mid()):.6f} beta_c={beta:.4f} L=-log d={-np.log(d):.4f}  Y={Y}', flush=True)
t0 = time.time()
res = {}
fails = [0]
def piece(name, fn, a_, b_, nseg, nn):
    edges = np.linspace(a_, b_, nseg+1)
    Sf = acb(0); Si = acb(0); Ad = 0.0; Ai = 0.0; mr = 0.0
    x, w = np.polynomial.legendre.leggauss(nn)
    for e0, e1 in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(x, w):
            t = (e0+e1)/2 + (e1-e0)/2*xi; ww = (e1-e0)/2*wi
            nu, jac, wt = fn(t)
            try:
                Kf, _ = Kplus(nu, st, U, jmax=40000)
            except RuntimeError:
                fails[0] += 1; continue          # j-series unusable at this node (small c near the turning point): node skipped
            Ki, _ = Kplus(nu, st, U, inf=True)
            vf = wt*Kf*jac; vi = wt*Ki*jac
            Sf += vf*ww; Si += vi*ww
            Ad += float(abs(vf - vi).mid())*ww; Ai += float(abs(vi).mid())*ww
            r = float((abs(Kf - Ki)/abs(Ki)).mid()); mr = max(mr, r)
    if fails[0]: print(f'  [{name}: {fails[0]} nodes skipped (j-series failed)]'); fails[0] = 0
    out = dict(fm=float((2*k*Sf.real).mid()), inf=float((2*k*Si.real).mid()), absinf=2*k*Ai, absdiff=2*k*Ad, maxrel=mr)
    print(f'  {name:7s}: Re fm {out["fm"]:+.8e}  inf {out["inf"]:+.8e}  |inf| {out["absinf"]:.3e}  |fm-inf| {out["absdiff"]:.3e}  max|m-1| {mr:.3e}  [{time.time()-t0:.0f}s]', flush=True)
    res[name] = out
    return Sf, Si
delta = 0.25
# bottom: real axis [0, delta] with the full kernel K (no split): K = 2 Re kappa for real nu
piece('bottom', lambda t: (acb(t), 1, psi(acb(t)/LAM)), 1e-12, delta, 1, 24)
piece('left', lambda y: (acb(delta, y), I, psi(acb(delta, y)/LAM)), 0, Y, max(2, int(2*Y)), 12)
piece('right', lambda y: (LAM + I*y, I, psi((LAM + I*y)/LAM)), 0, Y, max(2, int(Y)), 10)
piece('top', lambda t: (acb(t, Y), 1, psi(acb(t, Y)/LAM)), delta, Lf, NT, 10)
tf = res['bottom']['fm'] + res['left']['fm'] - res['right']['fm'] + res['top']['fm']
ti = res['bottom']['inf'] + res['left']['inf'] - res['right']['inf'] + res['top']['inf']
print(f' TOTAL fm {tf:.10e}  inf {ti:.10e}  rel diff {(tf-ti)/ti:.3e}  sum|diff| {sum(r["absdiff"] for r in res.values()):.3e}  sum|inf| {sum(r["absinf"] for r in res.values()):.3e}')
res.update(dict(state=[J, a, b], beta=beta, Y=Y, total_fm=tf, total_inf=ti))
json.dump(res, open(f'{CACHE}/ct_{J}_{a}_{b}_{Y}.json', 'w'))
