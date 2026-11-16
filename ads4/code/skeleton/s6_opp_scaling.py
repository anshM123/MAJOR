"""
s6 (follow-up 2): scaling of the opposite-helicity non-uniformity d(x, m) (s4 definition, c = 0, i.e. relative to scalar skeleton x constant),
diagnostic fits [NE, never derived], and an INFERRED projection to the graviton +- near-gap states of UNIFORM_REMAINDER.md.

(1a) d m / x at fixed x vs m.
(1b) one-parameter helicity-shift fit: d ~ c g(x, m),  g = (d[c=0] - d[c=8])/8  (the exact shape of eta -> eta - c/m^2 on the preserved part,
     linearised); c fitted per case over x <= X.
(1c) polynomial ansaetze over all cases at fixed (beta, k, kernel):  A: d = c1 x/m + c2 x^2/m + c3 x^2/m^2;  B: d = c1 x/(m sinh beta);
     C: d = c1 x^2 (flat-level, m-independent).  Reported: coefficients and rms residual / rms d.
(2)  graviton projection (c_g = 2 J_z^2 = 32 for +-; inference): envelope of the curved helicity shift on near-gap states, relative to the
     bulk-point envelope E_k:  dA = max_x |Z(eta - c/m^2) - Z(eta)|/|Z(eta)| * 2/E_k(x)    (all of it; Z = Harish-Chandra amplitude, |P| for nu beta < 1.5)
     dB = its curvature part only: max_x nu c |1/sinh(beta) - 1/beta| / m^2 * 2/E_k(x)       (the flat part being in F_flat via d~_44)
     with m = Delta_g m_flat (AdS units), sinh(beta/2) = (2J+1)/(2m), nu = x m, x in (0, 1), k = 2 and 3.
Writes cache/out_s6.txt.
"""
import os, pickle
from common import CACHE
import numpy as np, mpmath as mp
mp.mp.dps = 20
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

TAILMAX = 1e-4
allres = pickle.load(open(os.path.join(CACHE, 's4_allres.pkl'), 'rb'))
rows = sorted(((m, b, res) for (res, m, b) in allres.values()), key=lambda t: (round(t[1], 1), t[0]))

def pts(res, k, kern, h, Df, c):
    return [p for p in res[(k, kern, h, Df, c)]['pts'] if p[2] < TAILMAX and p[0] > 0.01]

# ---------------- (1a)
pr("# (1a) opp, Df = 2: d*m/x  (d = m^2 d / m^2) at x = 0.1, 0.2, 0.3   [odd k=6 | odd k=4 | even k=6];  envelope ref c/sinh(beta) with c = 8")
for (m, b, res) in rows:
    cells = []
    for (k, kern) in ((6, 'odd'), (4, 'odd'), (6, 'even')):
        P = pts(res, k, kern, 'opp', 2.0, 0)
        v = []
        for x0 in (0.1, 0.2, 0.3):
            c = [p for p in P if abs(p[0] - x0) < 0.026]
            v.append(f"{min(c, key=lambda p: abs(p[0]-x0))[1]/m/x0:7.1f}" if c else "      -")
        cells.append("".join(v))
    pr(f"  beta={b:.2f} m={m:6.1f} | " + " | ".join(cells) + f" | 8/sinh(beta) = {8/np.sinh(b):.1f}")

# ---------------- (1b)
pr("\n# (1b) helicity-shift fit d = c g (g = exact linearised shape of eta -> eta - c/m^2), per case; c over x <= 0.2 / x <= 0.3; residual rms/(rms d)")
cfits = []
for (m, b, res) in rows:
    line = []
    for (k, kern) in ((6, 'odd'), (4, 'odd'), (6, 'even'), (4, 'even')):
        for Df in (2.0, 2.5):
            P0 = {round(p[0], 4): p[1] for p in pts(res, k, kern, 'opp', Df, 0)}
            P8 = {round(p[0], 4): p[1] for p in pts(res, k, kern, 'opp', Df, 8)}
            for X in (0.2, 0.3):
                xs = [x for x in P0 if x in P8 and x <= X + 1e-4]
                if len(xs) < 2:
                    continue
                d0 = np.array([P0[x] for x in xs]); g = (d0 - np.array([P8[x] for x in xs]))/8
                c = float(d0 @ g/(g @ g)); res_ = np.sqrt(np.mean((d0 - c*g)**2))/np.sqrt(np.mean(d0**2))
                cfits.append(dict(m=m, b=b, k=k, kern=kern, Df=Df, X=X, c=c, r=res_))
                if Df == 2.0:
                    line.append(f"{kern[0]}{k} X={X}: c={c:5.2f} r={res_:.2f}")
    pr(f"  beta={b:.2f} m={m:6.1f}: " + ";  ".join(line))
for b0 in (0.4, 1.0):
    for X in (0.2, 0.3):
        cs = [f['c'] for f in cfits if abs(f['b'] - b0) < 0.05 and f['X'] == X]
        cs_big = [f['c'] for f in cfits if abs(f['b'] - b0) < 0.05 and f['X'] == X and f['m'] >= 20]
        if cs:
            pr(f"  summary beta={b0} X<={X}: c in [{min(cs):.2f}, {max(cs):.2f}], median {np.median(cs):.2f} (all ansatz variants k, kernel, Df, m);"
               f" m >= 20: [{min(cs_big):.2f}, {max(cs_big):.2f}] median {np.median(cs_big):.2f}" if cs_big else "")

# ---------------- (1c)
pr("\n# (1c) global ansaetze per (beta, k, kernel), Df = 2, x <= 0.3, all m with data:  coefficients; rms residual / rms d")
for b0 in (0.4, 1.0):
    for (k, kern) in ((6, 'odd'), (4, 'odd'), (6, 'even'), (4, 'even')):
        X_, Y_ = [], []
        for (m, b, res) in rows:
            if abs(b - b0) > 0.05:
                continue
            for (x, m2d, t) in pts(res, k, kern, 'opp', 2.0, 0):
                if x <= 0.3 + 1e-4:
                    X_.append((x, m, b)); Y_.append(m2d/m**2)
        if len(Y_) < 4:
            continue
        Y = np.array(Y_); x = np.array([t[0] for t in X_]); mm = np.array([t[1] for t in X_]); bb = np.array([t[2] for t in X_])
        rmsY = np.sqrt(np.mean(Y**2))
        def fit(cols):
            A = np.stack(cols, axis=1); c, *_ = np.linalg.lstsq(A, Y, rcond=None)
            return c, np.sqrt(np.mean((Y - A @ c)**2))/rmsY
        cA, rA = fit([x/mm, x*x/mm, x*x/mm**2]); cB, rB = fit([x/(mm*np.sinh(bb))]); cC, rC = fit([x*x])
        pr(f"  beta={b0} {kern} k={k} ({len(Y)} pts): A c1={cA[0]:7.2f} c2={cA[1]:8.2f} c3={cA[2]:8.1f} r={rA:.2f} | B c1={cB[0]:6.2f} r={rB:.2f} | C c1={cC[0]:6.3f} r={rC:.2f}")

# ---------------- (2)+(3)
def Z(nu, eta):
    nu = mp.mpf(nu); eta = mp.mpf(eta)
    return 2*mp.gamma(1j*nu)/(mp.sqrt(mp.pi)*mp.gamma(0.5 + 1j*nu))*(2*eta)**(-0.5 + 1j*nu)*mp.hyp2f1(0.25 - 0.5j*nu, 0.75 - 0.5j*nu, 1 - 1j*nu, 1/eta**2)

def Pn(nu, eta):
    return mp.re(mp.hyp2f1(0.5 - 1j*nu, 0.5 + 1j*nu, 1, (1 - mp.mpf(eta))/2))

def Ek(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2)

def dproj(Dg, mfl, J, c, k):
    m = Dg*mfl
    beta = 2*np.arcsinh((2*J + 1)/(2*m)); eta = np.cosh(beta)
    best = (0, 0, 0)
    for x in np.linspace(0.02, 0.98, 49):
        nu = x*m
        if nu*beta >= 1.5:
            z0, z1 = Z(nu, eta), Z(nu, eta - c/m**2)
            rel = float(abs(z1 - z0)/abs(z0))
        else:
            p0 = Pn(nu, eta); rel = float(abs(Pn(nu, eta - c/m**2) - p0)/abs(p0))
        dA = rel*2/Ek(k, x)
        dB = nu*c*abs(1/np.sinh(beta) - 1/beta)/m**2*2/Ek(k, x)
        if dA > best[0]:
            best = (dA, best[1], x)
        best = (best[0], max(best[1], dB), best[2])
    return best, beta

KAPPA = {20: (0.58, 0.62), 50: (0.66, 0.53), 100: (0.44, 0.41), 300: (0.32, 0.32), 1000: (0.29, 0.29)}   # UNIFORM_REMAINDER.md, M02 +-, [Lam,1.5Lam), Df=3 | Df=2
pr("\n# (2) INFERENCE: graviton +- near-gap states, c_g = 32 (= 2 J_z^2, J_z = 4; = 4 x the current value 8).  m = Delta_g m_flat, J = 4..6, x in (0,1).")
pr("#     dA = all of the helicity shift (if F_flat carried none of it); dB = curvature part only (if F_flat's d~_44 carries the flat part).  max over J, m_flat in [1, 1.5]")
for Dg in (20, 50, 100, 300, 1000):
    for k in (2, 3):
        dA, dB, loc = 0, 0, None
        for J in (4, 5, 6):
            for mfl in (1.0, 1.1, 1.25, 1.5):
                (a, bb, x), beta = dproj(Dg, mfl, J, 32.0, k)
                if a > dA:
                    dA, loc = a, (J, mfl, round(x, 2), round(beta, 3))
                dB = max(dB, bb)
        kc = KAPPA[Dg]
        pr(f"  Dg={Dg:5d} k={k}: dA_max = {dA:.3f} at (J, m_flat, x, beta) = {loc};  dB_max = {dB:.2e};  kappa_crit(+-, M02, [1,1.5)Lam) = {kc[0]} (Df=3) / {kc[1]} (Df=2)"
           f"  -> A: {'below' if dA < min(kc) else 'ABOVE'}; B: {'below' if dB < min(kc) else 'ABOVE'}")
pr("\n# (2b) dA (k = 2, c = 32) vs m_flat at J = 4, 6 and at the largest J with (2J+1) ~ 0.3 m (b-hat = (2J+1)/m_flat ~ fixed)")
for Dg in (20, 100, 1000):
    s = []
    for mfl in (1.0, 1.5, 2.0, 3.0, 5.0):
        v = [dproj(Dg, mfl, J, 32.0, 2)[0][0] for J in (4, 6)]
        s.append(f"m_fl={mfl}: J4 {v[0]:.2f} J6 {v[1]:.2f}")
    pr(f"  Dg={Dg}: " + ";  ".join(s))
open(os.path.join(CACHE, 'out_s6.txt'), 'w').write("\n".join(out) + "\n")
