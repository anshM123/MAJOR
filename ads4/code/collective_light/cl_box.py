"""cl_box: block-wise versus collective light side of a scalar CHMRS functional on an exactly crossing-symmetric toy.

Functional (CHMRS arXiv:2106.10274 eqs. (BMellinDefinition), (CfromBMellin); d = 3, k = 2; overall constant dropped):
    C_nu[M] = (1/|Gamma(c)|^4) int dz/(2 pi i) Gamma(z)^2 Gamma(c-z) Gamma(cbar-z) Bhat(2z),   c = (A + i nu)/2, A = 2Df + 1/2,
    Bhat(T)  = -4 sum_{poles S'} Res M(S', T + 2Df - S') / ((2Df - S')^2 (S' - T))            (fixed U' = 2Df - T).
Heavy normalisation is built in (1/gamma_a(nu)^2): a heavy state's image is O(1) for nu below its mass.  Check: C_nu[1] = 1/Gamma(A).

Toy amplitude (crossing symmetric, meromorphic, Regge spin 0, double zeros at the subtraction points S, T, U = 2Df):
    M = f(S) f(T) f(U),   f(X) = (X - 2Df)^2 / ((X - sL)(X - sH)),   S + T + U = 4 Df.
So Bhat(T) = Bhat_L(T) + Bhat_H(T) = 0 identically (no arc, no subtraction residues).  sL is a light twist, sH a heavy one.
A pole at twist s with residue function R(T') contributes (s- and t-channel together)
    Bhat_s[R](T) = 4 R(T + 2Df - s)/(2Df - s) * [ 1/(s - T)^2 - 1/((2Df - s)(s - T)) ].
'collective light' = C_nu[Bhat_L] with the exact residue function R_L;  identically = -C_nu[Bhat_H] ('heavy side').
'block-wise light'  = sum_J C_nu[Bhat_sL[c_J Q_J]],  R_L = sum_J c_J Q_J the spin (Mack n = 0) decomposition.  The n = 0 Mack
polynomials are orthogonal w.r.t. Gamma(T/2)^2 Gamma((s-T)/2)^2 on Re T = s/2 (cl_mackorth.py), so Q_J is built by Stieltjes from
that weight; R_L is even under T' <-> U (Bose), so only even J occur.  The functional contour is the same line Re T = sL/2.
Descendant poles of the twist-sL primaries are not tracked (pole-term bookkeeping; see README sec. 2).

Everything in Arb (python-flint).  Output: out/box_<tag>.txt and out/box_<tag>.json
usage: python3 -B cl_box.py [tag]   (parameters below)
"""
import os, sys, json, time
sys.dont_write_bytecode = True
import numpy as np
from flint import arb, acb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
TAG = sys.argv[1] if len(sys.argv) > 1 else 'main'
P = dict(main=dict(Df=1.0, sL=2.5, sH=40.0, h=0.05, Y=320.0, N=160, prec=700,
                   nus=[2, 5, 10, 20, 30, 40, 50, 60], Lams=[10, 20, 30, 40], np_=24),
         quick=dict(Df=1.0, sL=2.5, sH=40.0, h=0.08, Y=200.0, N=60, prec=400, nus=[2, 10, 20, 30], Lams=[10, 20], np_=12),
         sH20=dict(Df=1.0, sL=2.5, sH=20.0, h=0.05, Y=320.0, N=160, prec=700,
                   nus=[2, 5, 10, 20, 30, 40], Lams=[10, 20, 30], np_=24))[TAG]
ctx.prec = P['prec']
Df, sL, sH = arb(P['Df']), arb(P['sL']), arb(P['sH'])
A = 2*Df + arb(0.5)
lines = []
def log(s=''):
    print(s, flush=True); lines.append(s)

def f(X):
    return (X - 2*Df)**2/((X - sL)*(X - sH))
def RL(Tp):                      # Res_{S=sL} M(S,T'),  U = 4Df - sL - T'
    return (sL - 2*Df)**2/(sL - sH)*f(Tp)*f(4*Df - sL - Tp)
def RH(Tp):
    return (sH - 2*Df)**2/(sH - sL)*f(Tp)*f(4*Df - sH - Tp)
def Bfac(s, T):                  # Bhat_s[R](T) = R(T + 2Df - s) * Bfac(s, T)
    return 4/(2*Df - s)*(1/(s - T)**2 - 1/((2*Df - s)*(s - T)))

# ---------------- nodes on Re T = sL/2 :  T = sL/2 + i y,  z = T/2
h = arb(P['h']); M_ = int(round(P['Y']/P['h']))
ys = [h*k for k in range(-M_, M_ + 1)]
Ts = [acb(sL/2, y) for y in ys]
zs = [T/2 for T in Ts]
t0 = time.time()
gz2 = [z.gamma()**2 for z in zs]
wts = [h*(g*g.conjugate()).real for g in gz2]                 # |Gamma(z)|^4 = Gamma(T/2)^2 Gamma((sL-T)/2)^2 on the line
RLn = [RL(T + 2*Df - sL) for T in Ts]
BL = [Bfac(sL, T) for T in Ts]
FL = [r*b for r, b in zip(RLn, BL)]                          # collective light Bhat_L(T)
FH = [RH(T + 2*Df - sH)*Bfac(sH, T) for T in Ts]              # heavy Bhat_H(T)
log('# cl_box tag=%s  Df=%s sL=%s sH=%s  nodes=%d h=%s Y=%s N=%d prec=%d' % (TAG, P['Df'], P['sL'], P['sH'], len(ys), P['h'], P['Y'], P['N'], P['prec']))
mx = max(abs(float((a + b).abs_upper())) for a, b in zip(FL[::97], FH[::97]))
log('per-T sum rule check: max_T |Bhat_L + Bhat_H| on the line = %.2e  (max |Bhat_L| = %.2e)'
    % (mx, max(float(a.abs_upper()) for a in FL[::97])))

# ---------------- pass 1: Stieltjes (monic, symmetric weight => alpha_k = 0), spin coefficients c_J
N = P['N']
pm, p0 = [arb(0)]*len(ys), [arb(1)]*len(ys)
nrm_prev = None
beta, cJ, nrm = [], [], []
RLr = [r.real for r in RLn]
log('imag part of R_L on the line (should be 0): %.1e' % max(abs(float(r.imag.abs_upper())) for r in RLn[::50]))
for k in range(N + 1):
    nk = sum(w*p*p for w, p in zip(wts, p0))
    ck = sum(w*r*p for w, r, p in zip(wts, RLr, p0))/nk
    b = arb(0) if k == 0 else nk/nrm_prev
    beta.append(b); cJ.append(ck); nrm.append(nk)
    pn = [y*p - (b*q if k > 0 else 0) for y, p, q in zip(ys, p0, pm)]
    # next: p_{k+1} = y p_k - beta_{k+1} p_{k-1}; beta_{k+1} needs ||p_{k+1}||, computed at the next step
    nrm_prev = nk
    if k < N:
        # recompute properly: p_{k+1} = y p_k - beta_k' p_{k-1} with beta_k' = ||p_k||^2/||p_{k-1}||^2
        bk = arb(0) if k == 0 else nk/nrm[k - 1]
        pn = [y*p - bk*q for y, p, q in zip(ys, p0, pm)]
        pm, p0 = p0, pn
betas = [arb(0)] + [nrm[k]/nrm[k - 1] for k in range(1, N + 1)]
log('pass 1 (Stieltjes) done in %.0f s' % (time.time() - t0))
# Q_J(T) := p_J(y), monic in y.  Size of the spin-J piece of R_L at the point T' = 2Df (T = sL, y = -i sL/2):
yst = acb(0, -sL/2)
def pvals_at(y, N):
    out = [acb(1)]; pmm, pp = acb(0), acb(1)
    for k in range(N):
        pn = y*pp - betas[k]*pmm
        out.append(pn); pmm, pp = pp, pn
    return out
pst = pvals_at(yst, N)
log('\nspin data: c_J ||Q_J|| (L2 size), c_J Q_J(T\'=2Df) (the coefficient of each state\'s t-image e^{pi nu/2} term)')
part = acb(0)
for J in range(0, N + 1, 2):
    part += cJ[J]*pst[J]
    if J <= 12 or J % 20 == 0:
        log('  J=%3d  c_J||Q_J|| = %+.3e   c_J Q_J(2Df) = %+.3e   partial sum = %+.3e' %
            (J, float((cJ[J]*nrm[J].sqrt()).mid()), float(cJ[J].mid()*pst[J].real.mid()), float(part.real.mid())))
log('  exact R_L(2Df) = %s (double zero)' % RL(2*Df))

# ---------------- pass 2: images per nu
def kernel(nu):
    c = acb(A/2, arb(nu)/2)
    nrm4 = (c.gamma()*c.conjugate().gamma())**2              # |Gamma(c)|^4
    return [h/(4*arb.pi())*g*(c - z).gamma()*(c.conjugate() - z).gamma()/nrm4 for g, z in zip(gz2, zs)]

def images(nu):
    K = kernel(nu)
    one = sum(K, acb(0))
    CL = sum((k*v for k, v in zip(K, FL)), acb(0))
    CH = sum((k*v for k, v in zip(K, FH)), acb(0))
    G = [k*b for k, b in zip(K, BL)]
    I = []
    pmm, pp = [arb(0)]*len(ys), [arb(1)]*len(ys)
    for k in range(N + 1):
        I.append(cJ[k]*sum((g*p for g, p in zip(G, pp)), acb(0)))
        if k < N:
            pn = [y*p - betas[k]*q for y, p, q in zip(ys, pp, pmm)]
            pmm, pp = pp, pn
    return one, CL, CH, I

def fl(x):
    return float(x.real.mid())

res = {}
log('\nper-nu images.  Barnes check: C_nu[1]*Gamma(A) = 1.  CL = collective light, CH = heavy side, CL+CH = 0 exactly.')
log('  SB_N = block-wise light sum over J <= N;  max|I_J| = largest single-state image;  abs = sum_J |I_J| (light-data-only bound)')
for nu in P['nus']:
    t1 = time.time()
    one, CL, CH, I = images(nu)
    S = np.cumsum([fl(x) for x in I])
    absb = sum(abs(fl(x)) for x in I)
    big = max(abs(fl(x)) for x in I)
    jbig = int(np.argmax([abs(fl(x)) for x in I]))
    res[str(nu)] = dict(one=fl(one*A.gamma()), CL=fl(CL), CH=fl(CH), imCL=float(CL.imag.mid()), I=[fl(x) for x in I],
                        rad=float(CL.rad()))
    log('nu=%3d  Barnes %.12f  CL=%+.6e (rad %.1e)  CH=%+.6e  CL+CH=%.1e | max|I_J|=%.3e (J=%d)  abs=%.3e  SB_40=%+.3e  SB_80=%+.3e  SB_%d=%+.6e  [%.0fs]'
        % (nu, fl(one*A.gamma()), fl(CL), float(CL.rad()), fl(CH), fl(CL + CH), big, jbig, absb, S[min(40, N)], S[min(80, N)], N, S[-1], time.time() - t1))

log('\nconvergence of the block-wise sum to the collective value (SB_N - CL)/max|I_J|, even N:')
for nu in P['nus']:
    r = res[str(nu)]; S = np.cumsum(r['I']); big = max(abs(x) for x in r['I'])
    log('  nu=%3d ' % nu + ' '.join('N=%d:%.1e' % (n, (S[n] - r['CL'])/big) for n in (10, 20, 40, 60, 80, 120, 160) if n <= N))

# ---------------- smeared functional  Omega_Lam = int_0^1 dp psi2(p) C_{Lam p},  psi2 = (1-p)^3 p (65 + 155 p + 47 p^2)
xg, wg = np.polynomial.legendre.leggauss(P['np_'])
pg = 0.5*(xg + 1); wg = 0.5*wg
psi2 = lambda p: (1 - p)**3*p*(65 + 155*p + 47*p*p)
log('\nsmeared Omega_Lam = int dp psi2(p) C_{Lam p}  (Gauss %d nodes in p)' % P['np_'])
sm = {}
for Lam in P['Lams']:
    t1 = time.time()
    acc = dict(CL=0.0, CH=0.0, I=np.zeros(N + 1), absI=np.zeros(N + 1))
    for p, w in zip(pg, wg):
        one, CL, CH, I = images(Lam*p)
        acc['CL'] += w*psi2(p)*fl(CL); acc['CH'] += w*psi2(p)*fl(CH)
        acc['I'] += w*psi2(p)*np.array([fl(x) for x in I])
        acc['absI'] += w*psi2(p)*np.abs(np.array([fl(x) for x in I]))
    S = np.cumsum(acc['I'])
    sm[str(Lam)] = dict(CL=acc['CL'], CH=acc['CH'], I=list(acc['I']), absI=list(acc['absI']))
    log('Lam=%3d  Omega[light, collective]=%+.6e  Omega[heavy]=%+.6e  max|state|=%.3e  sum|state|=%.3e  block-wise SB_N=%+.6e  e^{pi Lam/2}=%.3e  [%.0fs]'
        % (Lam, acc['CL'], acc['CH'], np.max(np.abs(acc['I'])), np.sum(np.abs(acc['I'])), S[-1], np.exp(np.pi*Lam/2), time.time() - t1))

json.dump(dict(params=P, cJ=[float(c.mid()) for c in cJ], cQ2Df=[float((c*q.real).mid()) for c, q in zip(cJ, pst)],
               nu=res, smeared=sm), open(os.path.join(HERE, 'out', 'box_%s.json' % TAG), 'w'))
open(os.path.join(HERE, 'out', 'box_%s.txt' % TAG), 'w').write('\n'.join(lines) + '\n')
log('total %.0f s' % (time.time() - t0))
