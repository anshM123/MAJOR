"""
fullpos3/fp3_grouped_dv.py (FULLPOS3.md §5) [DV]: size of the GROUPED carrier bound for the oscillatory local components at a point (c, h).
  psi2 endpoint component G_e(c e^t) = sum_{+-} e^{+-i lam t} B_+-(t),  lam = Lam sqrt(c)/(1-c),
  |B_+-| <= 2|1-z|^{-3} (sqrt2/pi)(1/2) J_e(Re xi) e^{Lam |Im phi|},  |Im phi| <= (R^2/2) max |xi''|,  xi'' = sqrt z (1+z)/(2(1-z)^2).
  Grouped local term:  |(Q_{M1}(theta) - 1) G_e (c)| <= sum_{+-} S_B(R) sum_i |q_i(+-i lam)| i! R^{-i},
  q_i = Taylor coefficients of Q(k) = sum_{m=1}^{M1} h^m P_m(k) at k = +-i lam.   Compared with the per-order bound and with F_cert.
Usage: fp3_grouped_dv.py beta a0 [M1]
"""
import sys, os, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, acb, arb_poly, acb_poly, ctx
import fp3_ht as H
import fp3_htdata as D

def shifted_coeffs(P, k0):
    """Taylor coefficients of polynomial P (list of arb) at k0 (acb), Horner with polynomials."""
    res = acb_poly([0])
    lin = acb_poly([k0, 1])
    for co in reversed(P):
        res = res*lin + acb_poly([acb(co)])
    return [res[i] for i in range(res.degree() + 1)]

def run(beta, a0, M1=70, b=0):
    t0 = time.time()
    s = 2*b + 1
    ctx.prec = 2000
    c = (arb(beta)/2).tanh()**2
    nmin, h = H.state_params(c, a0, s); h = arb(h.mid())
    lam = D.LAM*c.sqrt()/(1 - c)
    P = H.Pm_exact(c, s, M1)
    Qc = [arb(0)]*(3*M1 + 1)
    for m in range(1, M1 + 1):
        for j in range(P[m].degree() + 1):
            Qc[j] += h**m*P[m][j]
    sh = shifted_coeffs(Qc, acb(0, lam))
    # per-order version for comparison
    perord = []
    for m in range(1, M1 + 1):
        co = [h**m*P[m][j] for j in range(P[m].degree() + 1)]
        perord.append(shifted_coeffs(co, acb(0, lam)))
    Jtab = H.jt(2)
    out = []
    for R in [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6]:
        if R > 0.9*float(-(c.log()).mid()):
            continue
        # geometry on the circle |t| = R: min Re xi, max |1-z|^{-1}, max |xi''|
        bmin = None; om = None; x2 = None
        for i in range(64):
            ph = arb(i + 0.5)/64*2*arb.pi()
            t = acb(R*ph.cos(), R*ph.sin()); z = c*t.exp()
            xi = ((1 + z)/(1 - z)).acosh()
            bb = float(xi.real.mid()); oo = float(abs(1 - z).mid()); xx = float(abs(z.sqrt()*(1 + z)/(2*(1 - z)**2)).mid())
            bmin = bb if bmin is None else min(bmin, bb); om = oo if om is None else min(om, oo); x2 = xx if x2 is None else max(x2, xx)
        Je = float(H.jval(2, arb(bmin))[2].mid())
        SB = 2*om**(-3)*(2**0.5/math.pi)*0.5*Je*math.exp(float(D.LAM.mid())*R*R/2*x2)
        RA = arb(R)
        grp = float(sum((abs(sh[i])*arb.fac_ui(i)/RA**i for i in range(len(sh))), arb(0)).mid())
        po = float(sum((sum((abs(v[i])*arb.fac_ui(i)/RA**i for i in range(len(v))), arb(0)) for v in perord), arb(0)).mid())
        out.append((R, SB, grp, po, 2*SB*grp, 2*SB*po))
    best = min(out, key=lambda r: r[4]); bestp = min(out, key=lambda r: r[5])
    fc = [F for (x, y, F) in D.fcert_cells() if x <= Fr(beta).limit_denominator(1000) <= y]
    F = min(fc)
    w3 = float(((1 - c)**3).mid())
    print(f'beta={beta} a0={a0} c={float(c.mid()):.4f} lam={float(lam.mid()):.1f} h={float(h.mid()):.3e} M1={M1}: F_cert={F:.3e}', flush=True)
    for r in out:
        print(f'   R={r[0]:.2f}: S_B<={r[1]:.2e}  sum_i|q_i|i!R^-i grouped={r[2]:.2e} per-order={r[3]:.2e}  -> (1-c)^3*bound grouped={w3*r[4]/F:.2e} F, per-order={w3*r[5]/F:.2e} F')
    print(f'   BEST grouped (1-c)^3 bound / F_cert = {w3*best[4]/F:.3e} (R={best[0]}),  per-order = {w3*bestp[5]/F:.3e}  [{time.time()-t0:.0f}s]', flush=True)

if __name__ == '__main__':
    beta = float(sys.argv[1]); a0 = int(sys.argv[2]); M1 = int(sys.argv[3]) if len(sys.argv) > 3 else 70
    run(beta, a0, M1)
