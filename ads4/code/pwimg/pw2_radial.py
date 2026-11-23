"""
pwimg/pw2_radial.py (PW_IMAGES.md §3): is the zero of c(mu) at mu = +-i/2 (Regge spin l = 0, -1) a property of the
polarisation/angular structure (robust) or of the specific radial weights (accidental)?

Leading-Regge heavy image by the exact split-representation method of H7 §23.1 / EE_FUNCTIONAL §2 (h7i14/v4_hclosed.py,
spinfunc/u2_hclosed_struct.py; read-only imports): h_-+ = 2 pi K_A K_B R_-+, R rational.  Here the radial weight of side A (pair 2,4)
and/or side B (pair 1,3) is shifted, |x|^c -> |x|^{c + sA}, |y|^c -> |y|^{c + sB} (the reference K is shifted with it, so the
nu-only normalisation cancels in EE/scalar ratios).  For each shift: c_s(mu) := R_-^EE(mu - i)/R_-^sc(mu - i) and the consistency
R_+^EE(nu)/R_+^sc(nu) = c_s(nu - i); zeros of c_s.
Output: cache/pw2_radial.txt
"""
import os, sys, io, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14')); sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
import numpy as np
import sympy as sp
import v3_regge_kernel as V        # read-only
import v4_hclosed as H4            # read-only
import u2_hclosed_struct as U2     # read-only
import s1_offdiag_kernel as S1     # read-only

d = V._dl; I = sp.I; nu = H4.nu_s


def side(monos, side_name, shift, spin_e):
    z = (1, 0, 1)
    s2 = 1/sp.sqrt(2)
    em = (0, s2, -I*s2); ep = (0, s2, I*s2)
    if side_name == 'A':
        e, f = em, ep; kref, mref = 1 - d, -(1 - d)
    else:
        e, f = ep, em; kref, mref = d, -d
    cargs = (H4._md(z, e), H4._md(z, f), H4._md(e, f))
    cref = 3 + kref + shift; lref = -mref
    V3 = [0, 0, 0]; cache = {}
    for (coef, kf, mf, i, j) in monos:
        k = sp.sympify(kf); m = sp.sympify(mf)
        c = 3 + k + shift; lam = -m; kap = (2*lam - c - 3)/2
        if side_name == 'A':
            cv = sp.sympify(coef).subs({V._a2: cargs[0], V._a4: cargs[1], V._e24: cargs[2]})
        else:
            cv = sp.sympify(coef).subs({V._A1: cargs[0], V._A3: cargs[1], V._E13: cargs[2]})
        n1 = sp.simplify(c/2 - cref/2); n2 = sp.simplify((c + 3)/2 - lam - (cref + 3)/2 + lref)
        dmu = sp.simplify((c + 3 - lam) - (cref + 3 - lref))
        Kr = sp.Integer(2)**dmu*H4._gratio(1 + cref/2, n1)*H4._gratio((cref + 3)/2 - lref, n2)
        for mu in range(3):
            key = (i, j, mu)
            if key not in cache:
                cache[key] = H4._dpat(e, f, z, i, j, mu)
            V3[mu] += cv*Kr*cache[key].subs({H4._lam: lam, H4._kap: kap})
    V3 = [sp.factor(sp.simplify(v)) for v in V3]
    assert sp.simplify(V3[1]) == 0
    c2 = V3[2]; c1 = sp.factor(sp.simplify(-V3[0] - c2))
    return c1, c2


def Rpm(Ac, Bc):
    (c1, c2), (d1, d2) = Ac, Bc
    l = d - 1
    Rp = c1*d1*(l + 1)/(2*l + 1) + c1*d2
    Rm = c1*d1*l/(2*l + 1) + c2*d1
    sub = {d: sp.Rational(1, 2) + I*nu}
    return sp.factor(sp.simplify(Rm.subs(sub))), sp.factor(sp.simplify(Rp.subs(sub)))


def main():
    out = io.StringIO()
    def pr(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); out.write(s + '\n')
    pr(__doc__)
    A_EE = U2.monos_struct_sym('EE')
    A_sc = H4._monos_sym(False, 'A'); B_sp = H4._monos_sym(True, 'B'); B_sc = H4._monos_sym(False, 'B')
    t0 = time.time()
    for (sA, sB) in ((0, 0), (1, 1), (2, 0), (0, 2), (2, 2), (sp.Rational(1, 2), sp.Rational(-1, 2)), (3, 3)):
        RmE, RpE = Rpm(side(A_EE, 'A', sA, True), side(B_sp, 'B', sB, True))
        RmS, RpS = Rpm(side(A_sc, 'A', sA, False), side(B_sc, 'B', sB, False))
        cm = sp.factor(sp.simplify(RmE/RmS))                 # = c_s(nu + i)
        cp = sp.factor(sp.simplify(RpE/RpS))                 # should be c_s(nu - i)
        cs = sp.factor(sp.simplify(cm.subs(nu, nu - I)))     # c_s(mu)
        cons = sp.simplify(cp - cs.subs(nu, nu - I))
        num, den = sp.fraction(sp.together(cs))
        rz = sorted([complex(r) for r in sp.Poly(sp.expand(num), nu).nroots()], key=lambda z: (z.imag, z.real)) if sp.degree(num, nu) > 0 else []
        rp = sorted([complex(r) for r in sp.Poly(sp.expand(den), nu).nroots()], key=lambda z: (z.imag, z.real)) if sp.degree(den, nu) > 0 else []
        fm = lambda r: ', '.join(f'{z.real:+.4f}{z.imag:+.4f}i' for z in r)
        pr(f'shift (sA, sB) = ({sA}, {sB}):  c_s(mu) = {cs}\n    consistency R_+^EE/R_+^sc - c_s(nu - i) = {cons}\n'
           f'    zeros of c_s: {fm(rz)}\n    poles of c_s: {fm(rp)}   [{time.time() - t0:.0f}s]')
    open(os.path.join(HERE, 'cache', 'pw2_radial.txt'), 'w').write(out.getvalue())


if __name__ == '__main__':
    main()
