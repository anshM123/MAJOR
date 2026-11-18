"""
Referee check (b), (c): even-kernel circle weights n_|q| = (1/36, 1/16, 1), flip-freeness of K* and S_ff, closed forms N_4^{K*}, N_6^{K*}.

Independent pieces:
 * M-table entries are recomputed from laplace.M_coeff (NOT the pickled cache), l = 2.
 * Pointwise circle Fourier coefficients of D1^2, D1 D2, D2^2 are computed here by FFT of cos^2, -cos, 1 on the circle
   (D1 = cos(th_out - th_in), D2 = -1), and the tensor weights are n_|q| * w~_q with n to be TESTED, not assumed:
   we fit n from the flip responses (36:16:1) and then test the 18 Lemma-L-type preserved ratios and the closed forms.
 * The j'=0 celestial weights are re-derived here: w_q = (-1)^{q3} (delta/2)_n/(1-delta/2)_n, n = |q1+q3|/2, delta = 1/2 + i nu
   (Fourier coefficients of |sin psi|^{-delta}); we check this Fourier statement numerically too.
Output: cache/rf1_kernels.txt
"""
import os, sys, time
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
E2E = os.path.join(HERE, '..', 'e2e')
sys.path.insert(0, E2E)
cwd = os.getcwd(); os.chdir(E2E)
import laplace
os.chdir(cwd)
nu, a, g1 = laplace.nu, laplace.a, laplace.g1
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

t0 = time.time()
T = {}
for q1 in range(-2, 3):
    for q3 in range(-2, 3):
        for h1 in (1, -1):
            for h3 in (1, -1):
                T[(q1, q3, h1, h3)] = laplace.M_coeff(2, q1, q3, h1, h3)
pr(f"# rf1: M_coeff(l=2) recomputed from source in {time.time()-t0:.0f} s")

def g1num(nuv, aa):
    from mpmath import gamma, mpf, mpc
    ga = lambda A: gamma((1 + A - 1.5 - 1j*nuv)/2)*gamma((1 + A - 1.5 + 1j*nuv)/2)
    return complex(ga(aa - 1)/ga(aa))

def resp_num(w, h1, h3, nuv, aa):
    tot = 0j
    G = g1num(nuv, aa)
    for q, ww in w.items():
        R, S, mu = T[(q[0], q[1], h1, h3)]
        Rv = complex(sp.N(R.subs({nu: nuv, a: aa}))) if R != 0 else 0
        Sv = complex(sp.N(S.subs({nu: nuv, a: aa}))) if S != 0 else 0
        tot += complex(ww)*(Rv + G*Sv)
    return tot

# 1. flip responses per circle mode (q,-q): ratio should be C(4,2+q)^2 = 36 : 16 : 1
pr("## 1. one-pair flip responses per circle mode (q,-q) relative to (2,-2)")
for aa in (7, 9, 11):
    for nuv in (0.0, 1.1, 3.0, 10.0):
        f = {}
        for q in (0, 1, 2):
            wq = {(q, -q): 1}
            f[q] = (resp_num(wq, 1, 1, nuv, aa) + resp_num(wq, -1, -1, nuv, aa))/2
        pr(f"  a={aa} nu={nuv:5.1f}: flip(0,0)/flip(2,-2) = {f[0]/f[2]:.10f}   flip(1,-1)/flip(2,-2) = {f[1]/f[2]:.10f}")

# 2. pointwise Fourier coefficients on the circle (relative angle psi = th_out - th_in)
Npts = 64
psi = 2*np.pi*np.arange(Npts)/Npts
def fourier(fvals):
    c = np.fft.fft(fvals)/Npts
    return {q: c[q % Npts] for q in range(-2, 3)}
wt = {'D1sq': fourier(np.cos(psi)**2), 'D1D2': fourier(-np.cos(psi)), 'D2sq': fourier(np.ones(Npts))}
pr("## 2. pointwise circle coefficients w~_q (mode (q,-q) <-> e^{i q psi}):")
for k, v in wt.items():
    pr("  " + k + ": " + ", ".join(f"q={q}: {v[q].real:+.4f}" for q in range(-2, 3)))
n = {0: sp.Rational(1, 36), 1: sp.Rational(1, 16), 2: sp.Integer(1)}
def tensor_w(name, nn=n):
    return {(q, -q): nn[abs(q)]*sp.nsimplify(round(wt[name][q].real, 12)) for q in range(-2, 3) if abs(wt[name][q]) > 1e-12}
W = {k: tensor_w(k) for k in wt}
pr("  tensor weights with n=(1/36,1/16,1): " + "; ".join(f"{k}: {W[k]}" for k in W))

def comb(**c):
    w = {}
    for nm, cf in c.items():
        for q, v in W[nm].items():
            w[q] = w.get(q, 0) + sp.nsimplify(cf)*v
    return w
KW = {'D1sq': W['D1sq'], 'D1D2': W['D1D2'], 'D2sq': W['D2sq'],
      'Sff': comb(D1sq=1, D1D2=1), 'Kstar': comb(D1sq=1, D1D2=sp.Rational(1, 2), D2sq=-sp.Rational(1, 2))}

pr("## 3. Lemma-L flip ratios and flip-freeness (numerical, several nu, a)")
for aa in (7, 9, 11):
    for nuv in (0.0, 1.1, 3.0, 7.5):
        fl = {k: (resp_num(w, 1, 1, nuv, aa) + resp_num(w, -1, -1, nuv, aa))/2 for k, w in KW.items()}
        pr(f"  a={aa} nu={nuv}: flip D1D2/D1sq = {fl['D1D2']/fl['D1sq']:.8f}, D2sq/D1sq = {fl['D2sq']/fl['D1sq']:.8f};"
           f" |flip Sff|/|flip D1sq| = {abs(fl['Sff'])/abs(fl['D1sq']):.2e}, |flip K*|/|flip D1sq| = {abs(fl['Kstar'])/abs(fl['D1sq']):.2e}")

pr("## 4. preserved responses: Lemma-L ratios quoted in SPINNING_TTTT_EXT 6.2 (nu=0, a=7: D1D2/D1sq = -0.270947, D2sq/D1sq = +0.129489)")
for aa in (7, 9, 11):
    for nuv in (0.0, 1.1, 3.0):
        p = {k: (resp_num(w, 1, -1, nuv, aa) + resp_num(w, -1, 1, nuv, aa))/2 for k, w in KW.items()}
        pr(f"  a={aa} nu={nuv}: D1D2/D1sq = {(p['D1D2']/p['D1sq']).real:+.6f}  D2sq/D1sq = {(p['D2sq']/p['D1sq']).real:+.6f}"
           f"  K*/Sff = {(p['Kstar']/p['Sff']).real:.6f}  (closed form (4nu^2+(4a-7)^2)/(2(4nu^2+8(a-2)^2+1)) = "
           f"{(4*nuv**2+(4*aa-7)**2)/(2*(4*nuv**2+8*(aa-2)**2+1)):.6f})")

pr("## 5. alternative normalisation n = 1/C(4,2+q) (natural SO(3) spinor inner product) for comparison: does it pass the flip test?")
n1 = {0: sp.Rational(1, 6), 1: sp.Rational(1, 4), 2: sp.Integer(1)}
W1 = {k: tensor_w(k, n1) for k in wt}
for aa in (7,):
    for nuv in (0.0, 3.0):
        fl = {k: (resp_num(w, 1, 1, nuv, aa) + resp_num(w, -1, -1, nuv, aa))/2 for k, w in W1.items()}
        pr(f"  n=1/C: a={aa} nu={nuv}: flip D1D2/D1sq = {fl['D1D2']/fl['D1sq']:.5f}, D2sq/D1sq = {fl['D2sq']/fl['D1sq']:.5f} (Lemma L: -1, +1)")

pr("## 6. N_k^{K*} = pres[(1,3) j'=0 celestial] / pres[K*] (symbolic) vs the closed forms")
delta = sp.Rational(1, 2) + sp.I*nu
def w_cel0():
    w = {}
    for q1 in range(-2, 3):
        for q3 in range(-2, 3):
            M = q1 + q3
            if M % 2: continue
            nn = abs(M)//2
            w[(q1, q3)] = sp.Integer(-1)**q3*sp.rf(delta/2, nn)/sp.rf(1 - delta/2, nn)
    return w
# check the Fourier statement of the celestial weights: |sin psi|^{-delta} coefficients of e^{2 i n psi} ratio to n=0
import mpmath as mp
for nuv in (0.0, 2.0):
    d = mp.mpf(0.5) + 1j*nuv
    c0 = mp.quad(lambda x: abs(mp.sin(x))**(-d), [0, mp.pi/2, mp.pi])
    for nn in (1, 2):
        cn = mp.quad(lambda x: abs(mp.sin(x))**(-d)*mp.cos(2*nn*x), [0, mp.pi/2, mp.pi])
        pred = mp.rf(d/2, nn)/mp.rf(1 - d/2, nn)
        pr(f"  |sin|^-delta Fourier: nu={nuv} n={nn}: ratio c_n/c_0 = {complex(cn/c0):.8f}  vs (-1)^n (d/2)_n/(1-d/2)_n = {complex((-1)**nn*pred):.8f}")
def resp_sym(w, h1, h3, aa):
    tot = 0
    for q, ww in w.items():
        R, S, mu = T[(q[0], q[1], h1, h3)]
        tot += ww*(R + g1*S)
    return sp.together(sp.expand(tot.subs(a, aa)))
for k in (4, 6):
    aa = k + 5
    pc = sp.simplify((resp_sym(w_cel0(), 1, -1, aa) + resp_sym(w_cel0(), -1, 1, aa))/2)
    pk = sp.simplify((resp_sym(KW['Kstar'], 1, -1, aa) + resp_sym(KW['Kstar'], -1, 1, aa))/2)
    fk = sp.simplify((resp_sym(KW['Kstar'], 1, 1, aa) + resp_sym(KW['Kstar'], -1, -1, aa))/2)
    Nsym = sp.factor(sp.simplify(pc/pk))
    claim = {4: 8*(16*nu**4 + 216*nu**2 + 263)/((4*nu**2 + 81)*(4*nu**2 + 841)),
             6: 8*(16*nu**4 + 536*nu**2 + 2023)/((4*nu**2 + 169)*(4*nu**2 + 1369))}[k]
    diff = sp.simplify(Nsym - claim)
    pr(f"  k={k}: flip[K*] = {sp.simplify(fk)};  N_{k}^K* = {Nsym};  N - claim = {diff};  g1-free: {not Nsym.has(g1)}")
    f = sp.lambdify(nu, Nsym, 'numpy'); xs = np.linspace(0, 200, 4001)
    vals = np.real(f(xs + 0j))
    pr(f"     min on [0,200] = {vals.min():.5f} at nu = {xs[vals.argmin()]:.2f}; N(0) = {vals[0]:.5f}; N(200) = {vals[-1]:.4f}")
open(os.path.join(HERE, 'cache', 'rf1_kernels.txt'), 'w').write("\n".join(out) + "\n")
