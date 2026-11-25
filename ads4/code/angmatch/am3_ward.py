"""ANGULAR_MATCH item 3: the Ward step C^(0) = C^(1) used in ANGULAR_ROUTE §13 to remove the delta'' contact.

(a) Per-harmonic form of the mixed commutator from the LADDER_ALLJ §4 identity
        sum_O [ Gbar(n1).F(n2) - Fbar(n2).G(n1) - (i/2)( Fbar(n1).DF(n2) + (D Fbar(n2)).F(n1) ) ],
        G = (i/6) W F,   D = 1 + 2 cos(phi) + sin(phi) d/dtheta2,   phi = theta2 - theta1,
    for a generic real amplitude vector and generic real symmetric W:  result = -i * C(phi) with
        C(phi) = sum_k |a_k|^2 [cos k phi + (1+k/2) cos(k+1)phi + (1-k/2) cos(k-1)phi] + (1/3) sum_k a_k (W a)_k cos k phi.
    This fixes the relative normalisation (the 1/3) shared by the D-term and the first-moment term.
(b) The T family is annihilated by W at (Delta, J) = (3, 2): unnormalised §7 coefficients, so d_T = 0 (no first-moment term
    at harmonics 0 and 4 of the T family).
(c) Exponential Fourier coefficients C^(n) of the state-summed C(phi) on a random spectrum (T family + heavy b = 0, 4 ladders):
    C^(0) = C^(1) = C^(-1) = |a_{T,0}|^2 exactly; and the q-pairing reproduces the endpoint rule.
Output: cache/am3_ward.txt
"""
import os, time
import numpy as np
import sympy as sp

out = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'am3_ward.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()
t0 = time.time()

# ---------------- (a) ----------------
th1, th2 = sp.symbols('theta1 theta2', real=True)
h = 2
ms = [-1, 0, 1, 2, 3, 4]                       # O-helicities in a small generic example (harmonic k = h - m)
nm = len(ms)
a = sp.symbols('a0:%d' % nm, real=True)
w = {}
for i in range(nm):
    for j in range(i, nm):
        w[(i, j)] = w[(j, i)] = sp.Symbol('w%d%d' % (i, j), real=True)
Wm = sp.Matrix(nm, nm, lambda i, j: w[(i, j)])
av = sp.Matrix(a)
Wa = Wm*av
F = lambda th: [a[i]*sp.exp(sp.I*(h-ms[i])*th) for i in range(nm)]
G = lambda th: [sp.I/6*Wa[i]*sp.exp(sp.I*(h-ms[i])*th) for i in range(nm)]
phi = th2 - th1
def Dop(expr):   # D acting on a function of theta2 (with theta1 fixed)
    return expr + 2*sp.cos(phi)*expr + sp.sin(phi)*sp.diff(expr, th2)
F1, F2, G1 = F(th1), F(th2), G(th1)
tot = 0
for i in range(nm):
    tot += sp.conjugate(G1[i])*F2[i] - sp.conjugate(F2[i])*G1[i] \
        - sp.I/2*(sp.conjugate(F1[i])*Dop(F2[i]) + sp.conjugate(Dop(F2[i]))*F1[i])
p = sp.Symbol('p', real=True)
Cclaim = 0
for i in range(nm):
    kk = h - ms[i]
    Cclaim += a[i]**2*(sp.cos(kk*p) + (1 + sp.Rational(kk, 2))*sp.cos((kk+1)*p) + (1 - sp.Rational(kk, 2))*sp.cos((kk-1)*p)) \
        + sp.Rational(1, 3)*a[i]*Wa[i]*sp.cos(kk*p)
z1, z2 = sp.symbols('z1 z2', nonzero=True)       # z_j = exp(i theta_j): the identity becomes a Laurent-polynomial identity
expr = (tot + sp.I*Cclaim.subs(p, phi)).rewrite(sp.exp)
expr = sp.expand(sp.powsimp(sp.expand(expr)))
expr = expr.subs({sp.exp(sp.I*th1): z1, sp.exp(sp.I*th2): z2})
expr = sp.expand(expr.replace(lambda e: e.func == sp.exp, lambda e: sp.exp(sp.expand(e.args[0]))))
expr = expr.subs({th1: sp.log(z1)/sp.I, th2: sp.log(z2)/sp.I})
diff = sp.simplify(sp.expand(sp.powsimp(sp.expand(expr), force=True)))
P('(a) LADDER_ALLJ §4 identity  ==  -i C(phi)  with the §12 per-harmonic form (coefficient 1/3 on a(Wa)); exact residual:', diff)
# numeric spot check at random point
vals = {s: np.random.default_rng(1).normal() for s in list(a) + list(set(w.values()))}
vals.update({th1: 0.37, th2: 1.91})
P('    numeric residual at a random point:', complex(sp.N((tot + sp.I*Cclaim.subs(p, phi)).subs(vals))))
ctrl = Cclaim - sum(sp.Rational(1, 6)*a[i]*Wa[i]*sp.cos((h-ms[i])*p) for i in range(nm))     # 1/3 -> 1/6
P('    control (1/3 replaced by 1/6): residual', complex(sp.N((tot + sp.I*ctrl.subs(p, phi)).subs(vals))), '(must be nonzero)')

# ---------------- (b) ----------------
Dl, Jl, m = sp.symbols('Delta J m')
T_ = (Dl-1)*(Dl-2); C_ = Jl*(Jl+1)
diag = T_ + C_ - 2*m**2
up = -(Dl-m-2)*(Jl+m+1)        # coef of e_m in W e_{m+1}
dn = -(Dl+m-2)*(Jl-m+1)        # coef of e_m in W e_{m-1}
sub = {Dl: 3, Jl: 2}
P('(b) at (Delta,J) = (3,2): diag(m=+-2) =', diag.subs(sub).subs(m, 2), diag.subs(sub).subs(m, -2),
  '; coef e_1 in W e_2 =', up.subs(sub).subs(m, 1), '; coef e_{-1} in W e_{-2} =', dn.subs(sub).subs(m, -1),
  '; coef e_3 in W e_2 =', dn.subs(sub).subs(m, 3), '(J=2: e_3 absent anyway); coef e_{-3} in W e_{-2} =', up.subs(sub).subs(m, -3))
P('    => W e_{+-2} = 0: the T family carries no first-moment term (d_T = 0) at harmonic 0 (m = h) nor at harmonic 4 (m = -h)')
P('    normalised off-diagonals b_m = sqrt((C-m(m+1))(T-m(m+1))) at (3,2): b_1 =', sp.sqrt(((C_-2)*(T_-2)).subs(sub)),
  ', b_-2 =', sp.sqrt(((C_-2)*(T_-2)).subs(sub)), '(m(m+1) = 2 for m = 1, -2)')

# ---------------- (c) ----------------
def ladder(Jv, Dv, b):
    T = (Dv-1)*(Dv-2); C = Jv*(Jv+1); n = Jv-b+1; av = np.zeros(n); av[0] = 1.0
    for kk in range(n-1):
        s = b+kk; av[kk+1] = av[kk]*(kk+4)/(kk+1)*np.sqrt((C-s*(s+1))/(T-s*(s+1)))
    s = np.arange(b, Jv+1, dtype=float); A = np.diag(T+C-2*s**2)
    for i in range(n-1): A[i, i+1] = A[i+1, i] = -np.sqrt((C-s[i]*(s[i]+1))*(T-s[i]*(s[i]+1)))
    return av, A@av, T
rng = np.random.default_rng(17)
worst = 0; worst_pair = 0
for trial in range(200):
    Nmax = 80
    Ch = {}                                           # exponential Fourier coefficients
    def add(n, v):
        Ch[n] = Ch.get(n, 0.0) + v
    def addcos(n, coef):                              # coef * cos(n phi) -> exponential coefficients
        if n == 0: add(0, coef)
        else: add(n, 0.5*coef); add(-n, 0.5*coef)
    def harm(kk, w2, dk):                             # one signed harmonic kk with weight |a|^2 = w2, first-moment a(Wa) = dk
        addcos(kk, w2 + dk/3); addcos(kk+1, (1+kk/2)*w2); addcos(kk-1, (1-kk/2)*w2)
    # T family: |a_{T,0}|^2 = 1 at harmonic 0, |c|^2 at harmonic 4, no first moment (part (b))
    c2 = rng.random()
    harm(0, 1.0, 0.0); harm(4, c2, 0.0)
    heavy_pred = 0.0
    for st in range(int(rng.integers(1, 6))):
        b = int(rng.choice([0, 4])); Jv = int(rng.integers(b + 2, 40)); Dv = Jv + 2 + 10**rng.uniform(-1, 2.5)
        av, Aav, T = ladder(Jv, Dv, b); wt = rng.random()*1e-3/av.max()**2
        sgn = 1 if b == 0 else -1                     # b = 0 ladder: harmonic +(2+k); b = 4: -(2+k)
        for kk in range(len(av)):
            harm(sgn*(2+kk), wt*av[kk]**2, wt*av[kk]*Aav[kk])
        heavy_pred += wt*((T - b*(b-1))/3 + 1.5)*av[0]**2
    worst = max(worst, abs(Ch.get(0, 0) - 1), abs(Ch.get(1, 0) - 1), abs(Ch.get(-1, 0) - 1))
    # q-pairing with the full (untruncated, finite spectrum) coefficients: sum_n q_n C^(n) = -3 + heavy endpoint sum
    q = lambda n: -3.0 if n == 0 else (0.0 if abs(n) == 1 else 6.0/((abs(n)-1)*abs(n)*(abs(n)+1)))
    pair = sum(q(n)*v for n, v in Ch.items())
    worst_pair = max(worst_pair, abs(pair - (-3 + heavy_pred)))
P('(c) random spectra (200): max |C^(n) - 1| over n = 0, +1, -1 :', worst)
P('    q-pairing sum_n q_n C^(n) vs -3 + sum_heavy[(T-b(b-1))/3 + 3/2]|a_0|^2 : max abs diff', worst_pair)
P('    (the finite-spectrum pairing is exact: the sum rule itself then reads 0 = -3 + P, i.e. P = 3 when D = 0)')
P('(d) contact: an even distribution supported at phi = 0 with C^(n) = o(n^4) is c0 delta + c2 delta\'\';'
  ' C^(n) = (c0 - c2 n^2)/(2 pi); C^(0) = C^(1) forces c2 = 0. [algebra]')
P('done in %.1f s' % (time.time()-t0))
