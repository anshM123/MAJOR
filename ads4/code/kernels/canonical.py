"""
Canonical structures S^{++}, S^{+-}, S^{odd} (even J >= 4), S^{oddJ} (odd J >= 5), and J = 0, 2; comparison of the
pipeline amplitudes with closed_forms.py; export of the structure tables to data/structures_J*.json.
Normalisations (all coefficient vectors are polynomials in s = (Delta+J)/2 without common factor):
   S^{+-}  : no amplitude for m <= 0;  entry H14^2 H34^2 V4^{J-4}        has leading term  8 s^4
   S^{++}  : no amplitude for m >= 4;  entry V1 V3 V4^J H13              has leading term  8 s^4
   S^{odd} : (even J >= 4)             entry eps14 X14 V1 V3 V4^{J-2} H34 equals J(J-1)(J-2)
   S^{oddJ}: (odd J >= 5)              entry eps13 X13 V1 V3 V4^J        has leading term  8 s^3
   J = 0   : even: entry V1^2 V3^2 = 1 (rational in s);  odd: (beta_1, beta_2) = (2s+1, 2s-3)
   J = 2   : even: entry V1^2 V3^2 V4^2 = 2(s-5)(s-4)(7s-2); odd: entry eps14 X14 V1 V3^2 V4 = 8s+4
Run: python3 canonical.py
"""
import os, sys, json, pickle
import sympy as sp, mpmath as mp
import kernels_lib as kl, closed_forms as cf, pipeline as pl
mp.mp.dps = 30
S = kl.S

def load(J, par):
    return pickle.load(open(os.path.join(kl.CACHE, f'kern_J{J}_{par}.pkl'), 'rb'))

def lead(c):
    return sp.Poly(c, S).LC()

def canonical(J):
    """returns dict name -> (parity, monos, coeff vector, K, {m: P_m})"""
    out = {}
    if J % 2 == 0:
        monos, basis, K, res = load(J, 'even')
        if J >= 4:
            (r0, r1) = res
            a0, a1, b0, b1 = r0[0], r1[0], r0[4], r1[4]
            g = sp.gcd(a0, a1); gb = sp.gcd(b0, b1)
            A0, A1, B0, B1 = sp.cancel(a0/g), sp.cancel(a1/g), sp.cancel(b0/gb), sp.cancel(b1/gb)
            for name, (c0, c1), ref, target in (('pm', (A1, -A0), (0, 0, J - 4, 0, 2, 2, ''), 8),
                                                ('pp', (B1, -B0), (1, 1, J, 1, 0, 0, ''), 8)):
                raw = [sp.expand(c0*x + c1*y) for x, y in zip(*basis)]
                vec = kl._clean(raw)
                i = monos.index(ref)
                f = sp.Rational(target)/lead(vec[i])
                vec = [sp.factor(f*v) for v in vec]
                j = next(q for q, v in enumerate(vec) if v != 0)
                fac = sp.cancel(vec[j]/raw[j])
                P = {m: sp.factor(fac*(c0*r0.get(m, 0) + c1*r1.get(m, 0))) for m in range(-J, J + 1)}
                out[name] = ('even', monos, vec, K, P)
        else:
            vec, P = basis[0], res[0]
            if J == 0:
                f = 1/vec[monos.index((2, 2, 0, 0, 0, 0, ''))]
                vec = [sp.factor(f*v) for v in vec]; P = {m: sp.factor(f*p) for m, p in P.items()}
            out['pp'] = ('even', monos, vec, K, P)
        monos, basis, K, res = load(J, 'odd')
        vec, P = basis[0], res[0]
        if J >= 4:
            f = sp.Integer(J*(J - 1)*(J - 2))/vec[monos.index((1, 1, J - 2, 0, 0, 1, '14'))]
            vec = [sp.factor(f*v) for v in vec]; P = {m: sp.factor(f*p) for m, p in P.items()}
        out['odd'] = ('odd', monos, vec, K, P)
    elif J >= 5:
        monos, basis, K, res = load(J, 'odd')
        vec, P = basis[0], res[0]
        f = sp.Integer(8)/lead(vec[monos.index((1, 1, J, 0, 0, 0, '13'))])
        vec = [sp.factor(f*v) for v in vec]; P = {m: sp.factor(f*p) for m, p in P.items()}
        out['oddJ'] = ('odd', monos, vec, K, P)
    return out

def M_num(parity, K, P, J, Delta, m):
    sv = (mp.mpf(Delta) + J)/2
    val = sp.N(sp.sympify(P.get(m, 0)).subs(S, sp.Float(str(sv), 30)), 30)
    return kl.G_prefactor(parity, K, sv)*mp.mpc(str(sp.re(val)), str(sp.im(val)))

if __name__ == '__main__':
    os.makedirs('data', exist_ok=True)
    Js = [int(a) for a in sys.argv[1:]] or [0, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    for J in Js:
        try:
            can = canonical(J)
        except FileNotFoundError:
            print(f"J={J}: kernels not computed"); continue
        dump = {}
        for name, (par, monos, vec, K, P) in can.items():
            dump[name] = dict(parity=par, K=list(K), monomials=[list(m) for m in monos], coefficients=[str(v) for v in vec],
                              P={str(m): str(p) for m, p in P.items() if p != 0})
            errs = []
            for Delta in (mp.mpf(J) + mp.mpf('1.37'), mp.mpf(J) + mp.mpf('4.9'), mp.mpf(2*J) + mp.mpf('10.23')):
                if J >= 4 and name in ('pm', 'pp', 'odd', 'oddJ'):
                    m0 = 4 if name in ('pm', 'oddJ') else 0
                    closed = {'pm': cf.M4_pm, 'pp': cf.M0_pp, 'odd': cf.M0_odd, 'oddJ': cf.M4_oddJ}[name](Delta, J)
                    num = M_num(par, K, P, J, Delta, m0)
                    errs.append(abs(num/closed - 1))
                elif J == 0:
                    num = M_num(par, K, P, J, Delta, 0)
                    k2 = cf.K2(num, Delta, 0, 0)
                    closed = cf.K0_even_scalar2(Delta) if name == 'pp' else cf.K0_odd_scalar2(Delta)
                    errs.append(abs(k2/closed - 1))
            if errs:
                print(f"J={J:2d} {name:4s}: closed form vs pipeline, max rel. deviation = {mp.nstr(max(errs), 3)}")
        json.dump(dump, open(f'data/structures_J{J}.json', 'w'), indent=1)
