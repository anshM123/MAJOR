# r3_endtoend.py -- end-to-end numerical test of the realisation (note sec 1, 2.1) and of the light-side normalisations
# (note sec 3) at a small mass parameter M (default 3), where the partial-fraction coefficients can be handled at
# high precision.  Everything is built from the note's definitions and the paper's formulas:
#   f(x)=(1-z)^3F(z), F=z^3(1-z)^M prod_j(1-rho_j z)^-1, rho_j=j/(100 m^2), m=M+6           (note 2.1)
#   E(nu)=(Lam/2) pi nu tanh(pi nu)/cosh(pi nu) sum_j C_j P_{-1/2+i nu}(A_j)                  (note 2.1)
#   chi = B q,  B=Phi(z+i)/R_-(z),  q = 1/4 int W(v) sech(pi(z-v)/2) dv,  W=E/Phi           (note 1)
# Tests:
#  T0  f = sum C_j/(x+A_j)  and  (2/Lam) int_0^oo E P(x) dnu = f(x)
#  T1  F_k = (2/Lam) int E g_{k,3}       (Taylor coefficients, k=0..4)
#  T2  -(2/Lam^7) int E Q = (100F_2+36F_3)/Lam^6
#  T3  I = int_1^oo f Q_1 dx (z-form) = (2/Lam) int E/(nu^2+9/4)
#  T4  leading-Regge image of chi:  (2/Lam) int_0^oo chi [R_- P_{-1/2+i(nu+i)}(eta)+R_+ P_{-1/2+i(nu-i)}(eta)] = f(eta)
#  T5  raw light integrals of chi:  int chi*P_K = int E Q ;  int chi*Ein = int E D     (contour shift)
#  T6  paper's moment route (fp-moments, fp-single): P_k = D_k^2 L_k from Z_k = int g z_{k+2} equals F_k
#  T7  decay chi(nu) ~ nu^-1 e^{-pi nu}
import sys, math
import mpmath as mp
M = int(sys.argv[1]) if len(sys.argv) > 1 else 3
mp.mp.dps = 60
m = M+6
Lam = mp.mpf(1)
rho = [mp.mpf(j)/(100*m*m) for j in range(1, m+1)]
A = [(1+r)/(1-r) for r in rho]
K0 = mp.mpf(2)**(M+3)/mp.fprod([1-r for r in rho])
C = []
for j in range(m):
    C.append(K0*(-A[j]-1)**3/mp.fprod([A[i]-A[j] for i in range(m) if i != j]))
print('M=%d m=%d  max|C_j|=%s' % (M, m, mp.nstr(max(abs(c) for c in C), 5)))
fx = lambda x: mp.mpf(2)**(M+3)*(x-1)**3/mp.fprod([(1-r)*x+(1+r) for r in rho])
print('T0a f(2) direct=%s  partial fractions=%s' % (mp.nstr(fx(mp.mpf(2)), 20), mp.nstr(sum(C[j]/(2+A[j]) for j in range(m)), 20)))
def Pleg(nu, x):  # P_{-1/2+i nu}(x), x>=1, via 2F1(1/2-i nu,1/2+i nu;1;(1-x)/2)
    return mp.re(mp.hyp2f1(0.5-1j*nu, 0.5+1j*nu, 1, (1-x)/2))
def E(nu):
    if nu == 0: return mp.mpf(0)
    s = sum(C[j]*Pleg(nu, A[j]) for j in range(m))
    return Lam/2*mp.pi*nu*mp.tanh(mp.pi*nu)/mp.cosh(mp.pi*nu)*s
# Gauss-Legendre panels
mp.mp.dps = 60
GLn = 16
def gl_nodes(n):
    # nodes/weights on [-1,1]
    xs = []; ws = []
    for i in range(1, n+1):
        x = mp.cos(mp.pi*(i-0.25)/(n+0.5))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for kk in range(2, n+1):
                p0, p1 = p1, ((2*kk-1)*x*p1-(kk-1)*p0)/kk
            dp = n*(x*p1-p0)/(x*x-1)
            dx = p1/dp; x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps+5): break
        xs.append(x); ws.append(2/((1-x*x)*dp*dp))
    return xs, ws
gx, gw = gl_nodes(GLn)
def panel_nodes(a, b, h):
    N = []; Wt = []
    t = mp.mpf(a)
    while t < b-1e-12:
        u = min(t+h, mp.mpf(b))
        for x, w in zip(gx, gw):
            N.append((u-t)/2*x+(u+t)/2); Wt.append((u-t)/2*w)
        t = u
    return N, Wt
VMAX = 70; NUMAX = 26
nodesW, wtsW = panel_nodes(0, VMAX, mp.mpf('0.5'))
print('evaluating E on %d nodes ...' % len(nodesW)); sys.stdout.flush()
Ev = [E(v) for v in nodesW]
mp.mp.dps = 30
Ev = [mp.mpf(e) for e in Ev]
inte = lambda vals: mp.fsum(w*v for w, v in zip(wtsW, vals))
# T0b
for x in [mp.mpf('1.5'), mp.mpf(4)]:
    Px = [Pleg(v, x) for v in nodesW]
    print('T0b (2/Lam) int E P(%s) = %s   f = %s' % (x, mp.nstr(2/Lam*inte([e*p for e, p in zip(Ev, Px)]), 16), mp.nstr(fx(x), 16)))
# Taylor coefficients of F
def taylorF(N):
    R = [mp.mpf(1)]+[mp.mpf(0)]*N
    for r in rho:
        for i in range(1, N+1): R[i] += r*R[i-1]
    B = [(-1)**l*mp.binomial(M, l) if l <= M else 0 for l in range(N+1)]
    Am = [mp.fsum(B[i]*R[k-i] for i in range(k+1)) for k in range(N+1)]
    return [Am[k-3] if k >= 3 else mp.mpf(0) for k in range(N+1)]
Fk = taylorF(6)
def gk3(k, nu):
    a = mp.mpc(0.5, -nu)
    tot = 0
    for j in range(k+1):
        tot += mp.rf(a, j)**2/mp.factorial(j)**2*mp.rf(3-a, k-j)/mp.factorial(k-j)
    return mp.re(tot)
for k in range(5):
    val = 2/Lam*inte([e*gk3(k, v) for e, v in zip(Ev, nodesW)])
    print('T1  k=%d  (2/Lam) int E g_k = %s   F_k = %s' % (k, mp.nstr(val, 15), mp.nstr(Fk[k], 15)))
Qp = lambda v: v**6-mp.mpf(245)/4*v**4+mp.mpf(9619)/16*v**2-mp.mpf(51575)/64
Kv = -2/Lam**7*inte([e*Qp(v) for e, v in zip(Ev, nodesW)])
print('T2  -(2/Lam^7) int E Q = %s   (100F2+36F3)/Lam^6 = %s' % (mp.nstr(Kv, 15), mp.nstr((100*Fk[2]+36*Fk[3])/Lam**6, 15)))
EDint = inte([e/(v**2+mp.mpf(9)/4) for e, v in zip(Ev, nodesW)])
Fz = lambda zz: zz**3*(1-zz)**M/mp.fprod([1-r*zz for r in rho])
Iz = mp.quad(lambda zz: Fz(zz)*((1+zz)*(-mp.log(zz))-2*(1-zz)), [0, 0.1, 0.5, 1])
print('T3  I(z-form) = %s   (2/Lam) int E D = %s' % (mp.nstr(Iz, 15), mp.nstr(2/Lam*EDint, 15)))
# ---- chi ----
Pp = lambda w: (w**2+0.25)**2*(w**2+2.25)*(w**2+6.25)*(w**2+20.25)
Phi = lambda w: Pp(w)*mp.cosh(mp.pi*w/2)*mp.gamma(1.25+1j*w/2)*mp.gamma(1.25-1j*w/2)/mp.cosh(mp.pi*w)
Rm = lambda w: -2*(2*w-1j)*(2*w+1j)**2*(2*w+3j)*(2*w+5j)*(2*w+11j)/(w*(2*w-7j)*(2*w-3j)*(2*w+7j))
Rp = lambda w: Rm(-w)
Bf = lambda w: mp.re(Phi(w+1j)/Rm(w))
Wv = [e/mp.re(Phi(v)) for e, v in zip(Ev, nodesW)]
nodesN, wtsN = panel_nodes(0, NUMAX, mp.mpf('0.5'))
def q(x):
    return mp.fsum(w*wv*(mp.sech(mp.pi*(x-v)/2)+mp.sech(mp.pi*(x+v)/2)) for w, wv, v in zip(wtsW, Wv, nodesW))/4
print('building chi on %d nodes ...' % len(nodesN)); sys.stdout.flush()
chi = [Bf(v)*q(v) for v in nodesN]
intN = lambda vals: mp.fsum(w*v for w, v in zip(wtsN, vals))
for eta in [mp.mpf('1.5'), mp.mpf(5)]:
    img = 2/Lam*intN([c*mp.re(Rm(v)*mp.legenp(-0.5+1j*(v+1j), 0, eta, type=3)+Rp(v)*mp.legenp(-0.5+1j*(v-1j), 0, eta, type=3))
                      for c, v in zip(chi, nodesN)])
    print('T4  LR image of chi at eta=%s : %s   f(eta)=%s' % (eta, mp.nstr(img, 14), mp.nstr(fx(eta), 14)))
N10 = lambda v: 1024*v**10-195328*v**8+8229504*v**6-90780512*v**4+251603732*v**2-114141375
PK = lambda v: -(4*v**2+1)*N10(v)/(8*(4*v**2+9)*(4*v**2+49))
Ein = lambda v: -32*(16*v**4-440*v**2+273)/((4*v**2+9)*(4*v**2+49))
print('T5  int chi P_K = %s   int E Q = %s' % (mp.nstr(intN([c*PK(v) for c, v in zip(chi, nodesN)]), 14), mp.nstr(inte([e*Qp(v) for e, v in zip(Ev, nodesW)]), 14)))
print('T5  int chi Ein = %s   int E D = %s' % (mp.nstr(intN([c*Ein(v) for c, v in zip(chi, nodesN)]), 14), mp.nstr(EDint, 14)))
# T6 paper moment route
def fm_list(nu, mmax):
    s = nu**2+mp.mpf(1)/4
    f = [mp.mpf(0), mp.mpf(1)]
    for mm in range(1, mmax):
        f.append((-(s-2*mm*mm)*f[mm]-(mm-1)**2*f[mm-1])/(mm+1)**2)
    return f
KM = 5
Z = [0]*KM
gfac = lambda v: -mp.mpf(9216)/Lam*(4*v**2+1)/((4*v**2+9)*(4*v**2+49))
zs = []
for c, v in zip(chi, nodesN):
    f = fm_list(v, KM+5)
    zs.append([12*f[mm]-(mm-1)*f[mm-1]+(mm+1)*f[mm+1] for mm in range(KM+3)])
for kk in range(KM):
    Z[kk] = intN([gfac(v)*c*zz[kk+2] for c, v, zz in zip(chi, nodesN, zs)])
Lk = [Z[0]/3]+[(Z[kk]-mp.mpf(kk)/(kk+3)*Z[kk-1])/(kk+3) for kk in range(1, KM)]
for kk in range(KM):
    print('T6  k=%d  D_k^2 L_k (paper moments of chi) = %s   F_k = %s' % (kk, mp.nstr(mp.binomial(kk+3, 3)**2*Lk[kk], 12), mp.nstr(Fk[kk], 12)))
for v in [5, 10, 15, 20, 25]:
    v = mp.mpf(v)
    print('T7  chi(%s)*nu*e^{pi nu} = %s' % (v, mp.nstr(Bf(v)*q(v)*v*mp.exp(mp.pi*v), 8)))
