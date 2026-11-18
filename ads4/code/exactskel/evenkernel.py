"""
Even (2,4) kernels for the even moments C+_{2,4,6} (SPINNING_TTTT_EXT.md sect. 6, EXACT_SKELETON_NEARGAP.md sect. 9):
   S_ff = D1^2 + D1 D2            (project choice, SUPERCONVERGENT_UPLIFT.md item 6)
   K*   = D1^2 + D1 D2/2 - D2^2/2 = (D1 + D2)(D1 - D2/2)   (flip-free and without the isotropic circle mode)
Contraction weights on the circle modes in the project normalisation (superconv/s1_responses.py: D1^2 = 4[(2,-2)+(-2,2)] + (2/9)(0,0)):
   D1D2 = -(1/2)[(1,-1)+(-1,1)],  D2^2 = (4/9)(0,0);  S_ff: (+-2) 4, (+-1) -1/2, (0,0) 2/9;  K*: (+-2) 4, (+-1) -1/4, (0,0) 0.
Exact one-pair responses A(h1,h3) = sum_q w_q M^q_{h1h3}(nu; a) from e2e/laplace.M_table(2) (a = k + 5), and
   N_k^K(nu) = pres[(1,3) j'=0 celestial] / pres[K]      (k = 4, 6; the C+_k -> Psi_k/N_k^K replacement)
   M_K(nu; 7)                                              (k = 2: C+_2 keeps its moment, normalised by M_K(nu;7)^2)
Closed forms checked here (sympy, exact):  pres[K*](nu;a) = (4nu^2+(4a-7)^2)/(4nu^2+(2a-5)^2) x const,
   N_4^K* = 8(16nu^4+216nu^2+263)/((4nu^2+81)(4nu^2+841)),  N_6^K* = 8(16nu^4+536nu^2+2023)/((4nu^2+169)(4nu^2+1369)).
Graviton pole of C+_2 (I9): the Regge density enters a flip-free kernel only through its preserved pair responses, which the
normalisation by M_K(nu;7)^2 divides out, so the normalised residue factor is 1 for every flip-free K (no IP needed for the flip part).
Usage: python3 evenkernel.py   -> cache/evenkernel.txt.   Importable: Nk(k, nu, kernel), Mpres(nu, a, kernel) (numpy, closed forms).
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))

KERNELS = ('Sff', 'Kstar')

def Mpres(nu, a, kernel):
    """preserved one-pair response (project normalisation, up to the common factor), closed form"""
    x = 4*np.asarray(nu, dtype=float)**2
    if kernel == 'Sff':
        return (x + 8*(a - 2)**2 + 1)/(2*(x + (2*a - 5)**2))
    if kernel == 'Kstar':
        return (x + (4*a - 7)**2)/(4*(x + (2*a - 5)**2))
    raise ValueError(kernel)

def Nk(k, nu, kernel):
    """N_k^K(nu) for k = 4, 6 (ratio of the (1,3) j'=0 response to the kernel's response)"""
    x = 4*np.asarray(nu, dtype=float)**2
    c = {4: (216, 263, 81, 393, 841), 6: (536, 2023, 169, 649, 1369)}[k]
    num = 16*np.asarray(nu, dtype=float)**4 + c[0]*np.asarray(nu, dtype=float)**2 + c[1]
    if kernel == 'Sff':
        return 4*num/((x + c[2])*(x + c[3]))
    if kernel == 'Kstar':
        return 8*num/((x + c[2])*(x + c[4]))
    raise ValueError(kernel)

def residue_factor_C2(nu, kernel):
    """normalised graviton-pole residue factor of C+_2 relative to its value for S_ff: preserved Einstein density x M_K^2 / M_K^2"""
    return np.ones_like(np.asarray(nu, dtype=float))

if __name__ == "__main__":
    import sympy as sp
    sys.path.insert(0, os.path.join(HERE, '..', 'superconv')); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
    cwd = os.getcwd(); os.chdir(os.path.join(HERE, '..', 'e2e'))
    from laplace import M_table, nu, a, g1
    from s1_responses import w_celestial
    os.chdir(cwd)
    T = M_table(2)
    W = {'D1sq': {(2, -2): 4, (-2, 2): 4, (0, 0): sp.Rational(2, 9)},
         'D1D2': {(1, -1): -sp.Rational(1, 2), (-1, 1): -sp.Rational(1, 2)},
         'D2sq': {(0, 0): sp.Rational(4, 9)}}
    def comb(**c):
        w = {}
        for nme, cf in c.items():
            for q, v in W[nme].items():
                w[q] = w.get(q, 0) + cf*v
        return w
    KW = {'Sff': comb(D1sq=1, D1D2=1), 'Kstar': comb(D1sq=1, D1D2=sp.Rational(1, 2), D2sq=-sp.Rational(1, 2))}
    def resp(w, h1, h3, aa):
        tot = 0
        for q, ww in w.items():
            R, S, mu = T[(q[0], q[1], h1, h3)]
            tot += ww*(R + g1*S)
        return sp.simplify(sp.together(sp.expand(tot.subs(a, aa))))
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr("# evenkernel: exact one-pair responses of S_ff and K* (project normalisation) and the normalisations N_k^K (exactskel/evenkernel.py)")
    for nmx in ('D1sq', 'D1D2', 'D2sq'):
        aa = 7
        fe = sp.factor(sp.simplify((resp(W[nmx], 1, 1, aa) + resp(W[nmx], -1, -1, aa))/2))
        pr(f"  Lemma L check a=7: flip[{nmx}] = {fe}")
    for k in (2, 4, 6):
        aa = k + 5
        pres = {}; flip = {}
        for kn, w in list(KW.items()) + [('SC', w_celestial(0))]:
            pe = sp.factor(sp.simplify((resp(w, 1, -1, aa) + resp(w, -1, 1, aa))/2))
            fe = sp.factor(sp.simplify((resp(w, 1, 1, aa) + resp(w, -1, -1, aa))/2))
            pres[kn], flip[kn] = pe, fe
            if kn != 'SC':
                xs = np.linspace(0, 60, 601)
                fnum = sp.lambdify(nu, pe, 'numpy')
                ok = np.allclose(np.real(fnum(xs + 0j)), Mpres(xs, aa, kn)*complex(fnum(0j)).real/Mpres(0.0, aa, kn), rtol=1e-12)
                pr(f"  k={k} a={aa} {kn}: pres = {pe};  flip = {fe};  closed form (up to a constant) {'OK' if ok else 'MISMATCH'};  min pres on [0,60]: {np.min(np.real(fnum(xs + 0j))):.4g}")
        if k in (4, 6):
            for kn in KW:
                Nsym = sp.factor(sp.simplify(pres['SC']/pres[kn]))
                f = sp.lambdify(nu, Nsym, 'numpy'); xs = np.linspace(0, 60, 601)
                ok = np.allclose(np.real(f(xs + 0j)), Nk(k, xs, kn), rtol=1e-10)
                pr(f"  N_{k}^{kn} = {Nsym}   [closed form in Nk(): {'OK' if ok else 'MISMATCH'}];  N(0) = {float(np.real(f(0j))):.5f}, N(inf) -> {float(np.real(f(1e6+0j))):.4f}")
        else:
            r = sp.factor(sp.simplify(pres['Kstar']/pres['Sff']))
            pr(f"  k=2: C+_2 keeps its moment (normalised by M_K(nu;7)^2); M_K*/M_Sff = {r};  (M_K*/M_Sff)^2 at nu = 0, 3, 10, inf: "
               + ", ".join(f"{float(sp.N(r.subs(nu, v)))**2:.4f}" for v in (0, 3, 10)) + f", {float(sp.limit(r, nu, sp.oo))**2:.4f}")
            pr(f"       j'=0 replacement (reference only): N_2^K* = {sp.factor(sp.simplify(pres['SC']/pres['Kstar']))}")
    pr("  I9: normalised graviton-pole residue factor of C+_2 under K* relative to S_ff = 1 identically (flip-free; preserved responses divided out)")
    open(os.path.join(HERE, 'cache', 'evenkernel.txt'), 'w').write("\n".join(out) + "\n")
