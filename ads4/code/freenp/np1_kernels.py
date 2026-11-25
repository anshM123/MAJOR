"""
np1_kernels.py -- null-plane light-ray operators of the d = 3 free scalar (improved, any xi) and free Dirac fermion,
as one-body Fock operators, from the covariant stress-tensor matrix elements.  Exact sympy.

Conventions (match ANGULAR_ROUTE §17-18):
  ds^2 = -du dv + dy^2,  u = t - x,  v = t + x;   g_uv = -1/2, g_yy = 1;  metric (-,+,+) in (t,x,y).
  Massless momentum p = (p^0, p^x, q).  Light-front variables:
      K = -p_v = (p^0 - p^x)/2 >= 0   (conjugate to v; the 'Galilean mass' is 2K),
      h = -p_u = (p^0 + p^x)/2 = q^2/(4K) >= 0   (conjugate to u; this is H on one particle).
  Plane wave e^{i p.x} = e^{-i h u - i K v + i q y}.  Heisenberg: O(u) = e^{iHu} O e^{-iHu}, H >= 0.
  Light-ray operators on u = 0:  calE(y) = int dv T_vv,  calN(y) = int dv T_vy,  A(f) = int f calE,  N(g) = int g calN.
  One-particle states |K,q> with <K,q|K',q'> = (2pi)^2 2K delta(K-K') delta(q-q')  (covariant normalisation:
  d^2p/(2p^0) = dK dq/(2K)).  <p1|T_mn(x)|p2> = M_mn(p1,p2) e^{i(p2-p1).x}; int dv gives 2 pi delta(K1-K2).

Output: cache/np1_kernels.txt
"""
import sympy as sp, os
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

K, K1, K2 = sp.symbols('K K1 K2', positive=True); q1, q2, xi = sp.symbols('q1 q2 xi', real=True)
k, l, Q = sp.symbols('k l Q', real=True)
I = sp.I

def lf_vec(Kv, qv):
    """contravariant (t,x,y) components of the massless momentum with LF data (K,q)."""
    h = qv**2/(4*Kv)
    return sp.Matrix([h + Kv, h - Kv, qv])
eta = sp.diag(-1, 1, 1)
# coordinate map (u,v,y) -> (t,x,y): t=(u+v)/2, x=(v-u)/2.  Covector components: w_v = (w_t + w_x)/2, w_u = (w_t - w_x)/2.
def low(p):  # covariant (t,x,y)
    return eta*p
def comp_v(pl):  # covariant v-component of a covector given covariant (t,x,y) comps
    return (pl[0] + pl[1])/2
def comp_u(pl):
    return (pl[0] - pl[1])/2
def comp_y(pl):
    return pl[2]
def dot(a, b):
    return (a.T*eta*b)[0]

p1, p2 = lf_vec(K1, q1), lf_vec(K2, q2)
assert sp.simplify(dot(p1, p1)) == 0
P("== LF dictionary: p_v = -K, p_u = -q^2/(4K), p_y = q :",
  sp.simplify(comp_v(low(p1)) + K1) == 0, sp.simplify(comp_u(low(p1)) + q1**2/(4*K1)) == 0)

# ------------------------------------------------------------------ scalar
P("\n=== 1. free scalar, T_mn = d_m phi d_n phi - 1/2 g (dphi)^2 - xi (d_m d_n - g box) phi^2 ===")
# one-particle (a^dag a) kernel: M_mn = p1_m p2_n + p2_m p1_n - g_mn p1.p2 + 2 xi (D_m D_n - g_mn D^2),  D = p2 - p1
# for the v v and v y components g_vv = g_vy = 0.
p1v, p2v = comp_v(low(p1)), comp_v(low(p2)); p1y, p2y = q1, q2
Dv, Dy = p2v - p1v, p2y - p1y
Mvv = 2*p1v*p2v + 2*xi*Dv*Dv
Mvy = p1v*p2y + p2v*p1y + 2*xi*Dv*Dy
onshell = {K2: K1}           # int dv -> delta(K1 - K2)
Evv = sp.simplify(Mvv.subs(onshell)); Evy = sp.simplify(Mvy.subs(onshell))
P("  <p1|T_vv|p2> on delta(K1-K2):", Evv, "   (xi drops: improvement term is int dv d_v(...))")
P("  <p1|T_vy|p2> on delta(K1-K2):", Evy)
assert sp.simplify(Evv - 2*K1**2) == 0 and sp.simplify(Evy + K1*(q1 + q2)) == 0
# pair-creation kernel <0|T|p1 p2>: M = -(p1 p2 + p2 p1) + g p1.p2 - 2 xi (S S - g S^2), S = p1+p2; int dv -> delta(K1+K2)
P("  pair kernels <0|T_vv|p1p2> = -2 p1v p2v - 2 xi S_v^2,  <0|T_vy|p1p2> = -(p1v q2 + p2v q1) - 2 xi S_v S_y,")
P("  supported (after int dv) on K1 = K2 = 0, i.e. on exactly-forward collinear pairs; see section 3 for their vanishing.")

# ------------------------------------------------------------------ Dirac fermion
P("\n=== 2. free Dirac fermion, T_mn = (i/4)[psibar g_m d_n psi - d_n psibar g_m psi + (m<->n)] ===")
g0 = sp.Matrix([[0, 1], [-1, 0]]); g1 = sp.Matrix([[0, 1], [1, 0]]); g2 = sp.Matrix([[1, 0], [0, -1]])
gam = [g0, g1, g2]
for a in range(3):
    for b in range(3):
        assert gam[a]*gam[b] + gam[b]*gam[a] == 2*eta[a, b]*sp.eye(2)
beta = -g0                  # psibar = psi^dag beta;  beta gamma^mu Hermitian, ubar gamma^0 u > 0
for a in range(3):
    assert (beta*gam[a]).H == beta*gam[a]
def u(Kv, qv):
    return sp.Matrix([sp.sqrt(2*Kv), qv/sp.sqrt(2*Kv)])
def slash(p):
    pl = low(p); return sum((gam[a]*pl[a] for a in range(3)), sp.zeros(2))
assert sp.simplify(slash(p1)*u(K1, q1)) == sp.zeros(2, 1)
cur = [sp.simplify((u(K1, q1).H*beta*gam[a]*u(K1, q1))[0]) for a in range(3)]
assert all(sp.simplify(cur[a] - 2*p1[a]) == 0 for a in range(3))
P("  u(K,q) = (sqrt(2K), q/sqrt(2K)):  p-slash u = 0,  ubar gamma^mu u = 2 p^mu   [checked]")
gu = g0 - g1                                     # gamma^u  (u = t - x)
gv_low = -gu/2                                   # gamma_v = g_vu gamma^u
gy_low = g2
def bil(G, Ka, qa, Kb, qb):
    return sp.simplify((u(Ka, qa).H*beta*G*u(Kb, qb))[0])
Pv = comp_v(low(p1)) + comp_v(low(p2)); Py = q1 + q2
FMvv = sp.Rational(1, 4)*2*bil(gv_low, K1, q1, K2, q2)*Pv
FMvy = sp.Rational(1, 4)*(bil(gv_low, K1, q1, K2, q2)*Py + bil(gy_low, K1, q1, K2, q2)*Pv)
FEvv = sp.simplify(FMvv.subs(onshell)); FEvy = sp.simplify(FMvy.subs(onshell))
P("  <p1|T_vv|p2> on delta(K1-K2):", FEvv)
P("  <p1|T_vy|p2> on delta(K1-K2):", FEvy)
P("  ubar(p1) gamma_v u(p2) =", bil(gv_low, K1, q1, K2, q2), "  (good component only, q-independent)")
P("  ubar(p1) gamma_y u(p2) =", bil(gy_low, K1, q1, K2, q2), "  (good x bad: the constrained component enters calN)")
assert sp.simplify(FEvv - Evv.subs(xi, 0)) == 0 and sp.simplify(FEvy - Evy) == 0
P("  => fermion one-body kernels are IDENTICAL to the scalar ones (no spin term: transverse space is 1-dim,")
P("     there is no eps_AB d_B (spin density) structure in d = 3).  Antiparticles: same, by charge conjugation.")

# ------------------------------------------------------------------ 3. pair terms vanish on collinear-forward pairs
P("\n=== 3. pair (vacuum-creation) terms: delta(K1+K2) lives on exactly-forward pairs ===")
r1, r2, t1, t2, s = sp.symbols('r1 r2 t1 t2 s', positive=True)
def cov(r, t):
    return sp.Matrix([r, r*sp.cos(t), r*sp.sin(t)])
c1, c2 = cov(r1, s*t1), cov(r2, s*t2)          # angles -> 0 as s -> 0
def Kc(p): return (p[0] - p[1])/2
def qc(p): return p[2]
S_v = -(Kc(c1) + Kc(c2)); S_y = qc(c1) + qc(c2)
c1v, c2v = -Kc(c1), -Kc(c2)
sc_vv = -2*c1v*c2v - 2*xi*S_v**2
sc_vy = -(c1v*qc(c2) + c2v*qc(c1)) - 2*xi*S_v*S_y
def vbar_u(G, pa, pb):   # v(p) = u(p) for massless (same Dirac eq.); v-bar(pa) G u(pb)
    ua = sp.Matrix([sp.sqrt(2*Kc(pa)), qc(pa)/sp.sqrt(2*Kc(pa))]); ub = sp.Matrix([sp.sqrt(2*Kc(pb)), qc(pb)/sp.sqrt(2*Kc(pb))])
    return (ua.T*beta*G*ub)[0]
Dm_v = comp_v(low(c1)) - comp_v(low(c2)); Dm_y = qc(c1) - qc(c2)
fe_vv = -sp.Rational(1, 4)*2*vbar_u(gv_low, c2, c1)*Dm_v
fe_vy = -sp.Rational(1, 4)*(vbar_u(gv_low, c2, c1)*Dm_y + vbar_u(gy_low, c2, c1)*Dm_v)
K_sum = sp.series(sp.simplify(Kc(c1) + Kc(c2)), s, 0, 3).removeO()
P("  K1 + K2 =", sp.factor(K_sum), "+ O(s^4): positive-definite quadratic in the two angles (s = overall angle scale)")
for name, M in (('scalar T_vv', sc_vv), ('scalar T_vy', sc_vy), ('Dirac  T_vv', fe_vv), ('Dirac  T_vy', fe_vy)):
    ser = sp.series(sp.simplify(M), s, 0, 6).removeO()
    order = min(sp.Poly(sp.expand(ser), s).monoms())[0] if ser != 0 else None
    P(f"  {name}: pair matrix element = O(s^{order})  [leading: {sp.factor(sp.expand(ser).coeff(s, order))}]")
P("  int d^2p1 d^2p2/(p1^0 p2^0) delta(K1+K2) M psi  ~  int r dr delta(c r^2) r^k = 0 for k >= 1, so on any state whose")
P("  wavefunction is bounded near the forward-collinear set: calE|0> = calN|0> = 0 (both theories, any xi).")
P("  Hence <0|[N(g),A(f)]|0> = 0: no Schwinger (c-number) term.")

# ------------------------------------------------------------------ 4. one-body operators and the algebra
P("\n=== 4. one-body operators (normalised modes b = a/(2pi sqrt(2K)), b~(K,y) = (2pi)^-1/2 int dq e^{iqy} b) ===")
P("  calE(y) = int dK K b~^dag b~(K,y);  calN(y) = (i/2) int dK [b~^dag d_y b~ - (d_y b~^dag) b~];  H = int dK dq q^2/(4K) b^dag b")
P("  i.e. on L^2(dK dy):  a(f) = K f(y),  n(g) = -(1/2){g(y), p},  h = p^2/(4K)   (p = -i d_y)")
P("  = 1+1 dim Galilean QM (mass 2K): calE = mass density, -calN = momentum density, h = kinetic energy.")
# kernel check with plane-wave profiles f = e^{iky}, g = e^{ily}: <q1|mult(e^{iky})|q2> = delta(q1 - q2 - k)
a_ker = lambda qa, qb: K            # coefficient of delta(qa - qb - k) for a(e^{iky}); same K on both sides
n_ker = lambda qa, qb: -(qa + qb)/2  # coefficient of delta(qa - qb - l) for n(e^{ily})
h_ = lambda qq: qq**2/(4*K)
# [n(g), a(f)] at (q1, q1 - k - l):
comm_na = n_ker(Q, Q - l)*a_ker(Q - l, Q - l - k) - a_ker(Q, Q - k)*n_ker(Q - k, Q - k - l)
target = I*(I*k)*K                    # i a(g f'),  g f' = i k e^{i(k+l)y}
P("  [n(g), a(f)] - i a(g f') (plane waves) =", sp.simplify(comm_na - target))
comm_ha = (h_(Q) - h_(Q - k))*a_ker(Q, Q - k)
target2 = (I/2)*n_ker(Q, Q - k)*(I*k)  # (i/2) n(f'),  f' = i k e^{iky}
P("  [h, a(f)] - (i/2) n(f')           =", sp.simplify(comm_ha - target2))
# double commutator with f = A e^{iky} + B e^{-iky} (real profile) is covered by bilinearity; direct check, two modes:
k2 = sp.symbols('k2', real=True)
def comm_ha_ker(qa, kk):   # [h, a(e^{i kk y})] coefficient at (qa, qa-kk)
    return (h_(qa) - h_(qa - kk))*K
dc = (K*comm_ha_ker(Q - k, k2) - comm_ha_ker(Q, k2)*K)   # [a(e^{iky}), [h, a(e^{ik2 y})]] at (Q, Q-k-k2)
target3 = sp.Rational(1, 2)*K*(I*k)*(I*k2)                # (1/2) a(f1' f2') symmetrised below
P("  [a(f1),[h,a(f2)]] + (f1<->f2) - a(f1' f2') =", sp.simplify(dc + dc.subs({k: k2, k2: k}, simultaneous=True) - 2*target3))
P("  => [N(g),A(f)] = i A(g f'),  [H,A(f)] = (i/2) N(f'),  [A(f),[H,A(f)]] = (1/2) A(f'^2): EXACT, as identities of")
P("     integral kernels for all (Q, k, l); the residual R(g,f) of the author's §18 is identically 0 (no delta''' etc.).")

# ------------------------------------------------------------------ 5. Cordova-Shao local form
y_, z_ = sp.symbols('y z', real=True)
P("\n=== 5. local form ===")
P("  From [N(g),A(f)] = i int g f' calE:  (1/i)[calN(y), calE(z)] = -delta(y-z) d_z calE(z) + delta'(y-z) calE(z)")
P("  (delta' = derivative w.r.t. its argument), i.e. Cordova-Shao's canonical form with no extra contact term.")
# check the smeared form of the local kernel against a test case with explicit functions
# int dy dz g(y) f(z)[-delta(y-z) E'(z) + delta'(y-z) E(z)] = int dz [ -g f E' + (d/dy g->) ...]:
# int dy g(y) delta'(y-z) = -g'(z)  -> = int dz f [-g E' - g' E] = -int f (g E)' = int f' g E.   (integration by parts)
P("  smeared: int dz f(z)[-g E' - g' E] = -int f (gE)' = int g f' E   [by parts]  -> i A(g f') ✓")

P("\n=== 6. the Galilean f-sum rule reading ===")
P("  [A(f),[H,A(f)]] = (1/2) A(f'^2) is the longitudinal f-sum (Thomas-Reiche-Kuhn) rule of the null-plane Galilean")
P("  system with particle masses m_i = 2K_i and rho_f = sum_i K_i f(y_i) = A(f): [rho_f,[H,rho_f]] = sum_i K_i^2 f'(y_i)^2/m_i")
P("  = (1/2) sum_i K_i f'(y_i)^2 = (1/2) A(f'^2).  Its only input is j = p/m (momentum density = mass current),")
P("  which is exactly [N(g),A(f)] = iA(gf').  Interacting Galilean systems keep it iff interactions are density-dependent.")

os.makedirs(os.path.join(os.path.dirname(__file__), 'cache'), exist_ok=True)
open(os.path.join(os.path.dirname(__file__), 'cache', 'np1_kernels.txt'), 'w').write('\n'.join(out) + '\n')
