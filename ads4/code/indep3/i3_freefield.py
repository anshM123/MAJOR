"""INDEP_CHECK3, claim 1 step (ii) and the end-to-end zero (Route B), written from scratch (C1_FINAL's code not read).

Free real conformally coupled scalar (B) and free massless Dirac fermion (F) in d = 3, mostly-plus metric.
T|0> is a two-particle state; the T family at momentum P (timelike) is spanned by psi_h = conj(eps_h(P)) . t(k1, k2),
eps_h = e_h e_h with e_+-(P) = Lambda_boost(P) (0, 1, +-i)/sqrt2.

Light-ray operators at future null infinity, direction n (nbar = (1, n)), one-body:
   E(n)     : sum_i  w(omega) a^dag a              (particle i along n)
   E1(n)    : sum_i  w(omega) (-i/2)(a^dag d_omega a - d_omega a^dag a)     [= int du u T at I+]
Matrix elements between T-family wave packets Psi = sum_h F_h(P) psi_h / sqrt(n(P)) are computed DIRECTLY from the two-particle
wavefunction (d/d omega at fixed spectator momentum, numerically, 4th order, mp.dps = 30), i.e. WITHOUT assuming the
decomposition E1 = -i H nbar.d + ... + M of C1_FINAL.  Phase-space bookkeeping (derived here):
   <Psi'|E(n)|Psi>  = int d^3P  sum_i [omega/(-P.nbar)] Psi'^* Psi |_{k_i = omega nbar}
   <Psi'|Psi>       = int d^3P  sum F'^* F,   n(P) = int dphi1 |psi|^2/(-P.nbar1)   (common constant dropped)
Checks: tracelessness/conservation of t, Ward identity int dphi nbar^mu H(P,phi) = P^mu, HM endpoints (|flip| = diag, B/F flips
opposite), M := E1 matrix on constant F vs the frame connection, and the pointwise O(t4) T-block term under the EE (2,4)
contraction (eta-contraction of the polarisation index pairs) with controls.
"""
import mpmath as mp
import random, sys, json

mp.mp.dps = 30
ETA = [-1, 1, 1]
I = mp.mpc(0, 1)

def dot(a, b): return -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def low(v): return [-v[0], v[1], v[2]]
def nbar(phi): return [mp.mpf(1), mp.cos(phi), mp.sin(phi)]

def boost(P):
    m = mp.sqrt(-dot(P, P)); g = P[0]/m; b = [P[1]/P[0], P[2]/P[0]]; b2 = b[0]**2 + b[1]**2
    L = [[g, g*b[0], g*b[1]], [g*b[0], 0, 0], [g*b[1], 0, 0]]
    for i in range(2):
        for j in range(2):
            L[1+i][1+j] = (1 if i == j else 0) + ((g - 1)*b[i]*b[j]/b2 if b2 > 0 else 0)
    return L

def frames(P):
    L = boost(P); s = 1/mp.sqrt(2)
    ep = [s*(L[mu][1] + I*L[mu][2]) for mu in range(3)]
    em = [s*(L[mu][1] - I*L[mu][2]) for mu in range(3)]
    return {+2: ep, -2: em}

def eps_up(P):
    fr = frames(P)
    return {h: [[fr[h][a]*fr[h][b] for b in range(3)] for a in range(3)] for h in (2, -2)}

def eps_low(P):
    e = eps_up(P)
    return {h: [[ETA[a]*ETA[b]*e[h][a][b] for b in range(3)] for a in range(3)] for h in (2, -2)}

# ---- stress-tensor two-particle matrix elements (lower indices) ----
def tB(k1, k2):
    P = [k1[i] + k2[i] for i in range(3)]; l1, l2, lP = low(k1), low(k2), low(P)
    k12 = dot(k1, k2); P2 = dot(P, P)
    return [[l1[a]*l2[b] + l2[a]*l1[b] - (ETA[a] if a == b else 0)*k12
             - (lP[a]*lP[b] - (ETA[a] if a == b else 0)*P2)/4 for b in range(3)] for a in range(3)]

G0 = mp.matrix([[0, 1], [-1, 0]]); G1 = mp.matrix([[0, 1], [1, 0]]); G2 = mp.matrix([[1, 0], [0, -1]])
GUP = [G0, G1, G2]; GLO = [-G0, G1, G2]
def spinor(k):   # k null, future; u = sqrt(omega) (sin(phi/2), cos(phi/2)); kslash u = 0
    w = k[0]; phi = mp.atan2(k[2], k[1])
    return mp.matrix([mp.sqrt(w)*mp.sin(phi/2), mp.sqrt(w)*mp.cos(phi/2)])
def tF(k1, k2):
    u1, v2 = spinor(k1), spinor(k2)
    ub = (u1.T * G0)   # ubar = u^dag gamma^0 (u real); overall constant irrelevant
    bil = [(ub * GLO[a] * v2)[0] for a in range(3)]
    d = low([k1[i] - k2[i] for i in range(3)])
    return [[(bil[a]*d[b] + bil[b]*d[a])/2 for b in range(3)] for a in range(3)]
TT = {'B': tB, 'F': tF}

def psi(th, k1, k2, el=None):
    P = [k1[i] + k2[i] for i in range(3)]
    e = eps_up(P) if el is None else el
    t = TT[th](k1, k2)
    return {h: sum(mp.conj(e[h][a][b])*t[a][b] for a in range(3) for b in range(3)) for h in (2, -2)}

def config(P, phi, which):   # particle `which` along nbar(phi)
    nb = nbar(phi); w = dot(P, P)/(2*dot(P, nb))
    k = [w*nb[i] for i in range(3)]; o = [P[i] - k[i] for i in range(3)]
    return (k, o, w, nb) if which == 1 else (o, k, w, nb)

NCACHE = {}
def norm_rest(th, N=64):     # n(P) at rest, m = 1; integrand is a trig polynomial at rest -> trapezoid exact
    if th in NCACHE: return NCACHE[th]
    P = [mp.mpf(1), mp.mpf(0), mp.mpf(0)]
    acc = {(a, b): mp.mpc(0) for a in (2, -2) for b in (2, -2)}
    for j in range(N):
        phi = 2*mp.pi*j/N; k1, k2, w, nb = config(P, phi, 1)
        ps = psi(th, k1, k2)
        for a in (2, -2):
            for b in (2, -2):
                acc[(a, b)] += mp.conj(ps[a])*ps[b]/(-dot(P, nb))*2*mp.pi/N
    NCACHE[th] = acc
    return acc
def nP(th, P):
    return norm_rest(th)[(2, 2)].real * (-dot(P, P))**mp.mpf(1.5)

def Hmat(th, P, phi):
    """normalised E(n) matrix H[h', h] at momentum P, direction phi."""
    nn = nP(th, P); H = {(a, b): mp.mpc(0) for a in (2, -2) for b in (2, -2)}
    for which in (1, 2):
        k1, k2, w, nb = config(P, phi, which); ps = psi(th, k1, k2)
        for a in (2, -2):
            for b in (2, -2):
                H[(a, b)] += w*mp.conj(ps[a])*ps[b]/(-dot(P, nb))/nn
    return H

def struct(fn, S, *args):
    """S = 'E' (Einstein, t4 = 0): (B+F)/2 ; 'W' (t4 part per unit t4/4): (B-F)/2 ; or 'B'/'F'."""
    if S in ('B', 'F'): return fn(S, *args)
    A, B = fn('B', *args), fn('F', *args)
    sg = 1 if S == 'E' else -1
    if isinstance(A, dict): return {k: (A[k] + sg*B[k])/2 for k in A}
    return (A + sg*B)/2

# ---------- direct first-moment density ----------
HSTEP = mp.mpf('1e-7')
def Psi_val(th, F, k1, k2):
    P = [k1[i] + k2[i] for i in range(3)]
    Fv = F(P); ps = psi(th, k1, k2); s = 1/mp.sqrt(nP(th, P))
    return sum(Fv[h]*ps[h] for h in (2, -2))*s
def D1(th, Fbra, Fket, P, phi):
    """density of <Psi_bra|E1(n)|Psi_ket> at P: sum_i [w/(-P.nbar)] (-i/2)(Psi_b^* dPsi_k - dPsi_b^* Psi_k), d = d/d omega_i."""
    tot = mp.mpc(0)
    for which in (1, 2):
        k1, k2, w, nb = config(P, phi, which)
        def at(dw):
            if which == 1: kk1 = [k1[i] + dw*nb[i] for i in range(3)]; kk2 = k2
            else: kk2 = [k2[i] + dw*nb[i] for i in range(3)]; kk1 = k1
            return Psi_val(th, Fbra, kk1, kk2), Psi_val(th, Fket, kk1, kk2)
        h = HSTEP
        vals = {s: at(s*h) for s in (-2, -1, 0, 1, 2)}
        def deriv(idx):
            return (vals[-2][idx] - 8*vals[-1][idx] + 8*vals[1][idx] - vals[2][idx])/(12*h)
        b0, k0 = vals[0]; db, dk = deriv(0), deriv(1)
        tot += w/(-dot(P, nb))*(-I/2)*(mp.conj(b0)*dk - mp.conj(db)*k0)
    return tot
def D0(th, Fbra, Fket, P, phi):
    tot = mp.mpc(0)
    for which in (1, 2):
        k1, k2, w, nb = config(P, phi, which)
        tot += w/(-dot(P, nb))*mp.conj(Psi_val(th, Fbra, k1, k2))*Psi_val(th, Fket, k1, k2)
    return tot

out = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s); sys.stdout.flush(); out.append(s)

def randP(rng, boosted=True):
    m = mp.mpf(rng.uniform(0.7, 1.6))
    if not boosted: return [m, mp.mpf(0), mp.mpf(0)]
    y = rng.uniform(0.1, 1.2); a = rng.uniform(0, 6.283)
    return [m*mp.cosh(y), m*mp.sinh(y)*mp.cos(a), m*mp.sinh(y)*mp.sin(a)]

def main():
    rng = random.Random(20261125)
    # ---- (1) t is conserved and traceless ----
    worst = 0
    for th in ('B', 'F'):
        for _ in range(4):
            P = randP(rng); k1, k2, _, _ = config(P, rng.uniform(0, 6.28), 1)
            t = TT[th](k1, k2); up = [[ETA[a]*ETA[b]*t[a][b] for b in range(3)] for a in range(3)]
            tr = sum(ETA[a]*t[a][a] for a in range(3))
            cons = max(abs(sum(P[a]*t[a][b] for a in range(3))) for b in range(3))
            sym = max(abs(t[a][b] - t[b][a]) for a in range(3) for b in range(3))
            worst = max(worst, abs(tr), cons, sym)
    log(f"(1) t_B, t_F: max |trace|, |P.t|, |t - t^T| = {mp.nstr(worst, 3)}")
    # ---- (2) Lorentz invariance of n(P) (boosted, brute-force quadrature) ----
    for th in ('B', 'F'):
        P = randP(rng); N = 400; acc = mp.mpf(0)
        for j in range(N):
            phi = 2*mp.pi*j/N; k1, k2, w, nb = config(P, phi, 1); ps = psi(th, k1, k2)
            acc += abs(ps[2])**2/(-dot(P, nb))*2*mp.pi/N
        off = norm_rest(th)[(2, -2)]
        log(f"(2) {th}: n(P) boosted quadrature / (n_rest m^3) - 1 = {mp.nstr(acc/nP(th, P) - 1, 3)};"
            f" rest off-diagonal n_(+,-) = {mp.nstr(abs(off), 3)}; n_(++)-n_(--) = {mp.nstr(abs(norm_rest(th)[(2,2)]-norm_rest(th)[(-2,-2)]), 3)}")
    # ---- (3) Ward identity and HM endpoints ----
    res = {}
    for th in ('B', 'F'):
        wmax = 0; Ps = [randP(rng, boosted=(j > 0)) for j in range(4)]
        for P in Ps:
            N = 256 if P[1] == 0 else 400
            acc = [[mp.mpc(0)]*3 for _ in range(3)]  # acc[hidx][mu]
            for j in range(N):
                phi = 2*mp.pi*j/N; H = Hmat(th, P, phi); nb = nbar(phi)
                for mu in range(3):
                    acc[0][mu] += nb[mu]*H[(2, 2)]*2*mp.pi/N
                    acc[1][mu] += nb[mu]*H[(-2, -2)]*2*mp.pi/N
                    acc[2][mu] += nb[mu]*H[(2, -2)]*2*mp.pi/N
            for mu in range(3):
                wmax = max(wmax, abs(acc[0][mu] - P[mu]), abs(acc[1][mu] - P[mu]), abs(acc[2][mu]))
        res[th] = wmax
    log(f"(3) Ward: max |int dphi nbar^mu H_hh - P^mu|, |int nbar^mu H_flip|: B {mp.nstr(res['B'],3)}, F {mp.nstr(res['F'],3)}")
    dd = fs = rat = 0; phases = []
    for _ in range(6):
        P = randP(rng); phi = rng.uniform(0, 6.28)
        HB, HF = Hmat('B', P, phi), Hmat('F', P, phi)
        dd = max(dd, abs(HB[(2, 2)] - HF[(2, 2)]), abs(HB[(-2, -2)] - HF[(-2, -2)]), abs(HB[(2,2)] - HB[(-2,-2)]))
        fs = max(fs, abs(HB[(2, -2)] + HF[(2, -2)]))
        rat = max(rat, abs(abs(HB[(2, -2)])/HB[(2, 2)].real - 1))
    Pr = [mp.mpf(1), 0, 0]
    for phi in (mp.mpf('0.37'), mp.mpf('1.1')):
        phases.append((float(phi), float(mp.arg(Hmat('B', Pr, phi)[(2, -2)]))))
    log(f"(3) HM: max|H_B,diag - H_F,diag| (and ++ vs --) = {mp.nstr(dd,3)}; max|H_B,flip + H_F,flip| = {mp.nstr(fs,3)};"
        f" max | |H_B,flip|/H_diag - 1 | = {mp.nstr(rat,3)}; rest-frame arg H_B(+2,-2)(phi): {phases}")
    # ---- (4) first-moment local matrix M (constant F) vs frame connection ----
    mx = {'flip': 0, 'BF': 0, 'conn': 0, 'herm': 0}
    for _ in range(5):
        P = randP(rng); phi = rng.uniform(0, 6.28)
        # frame connection along nbar: Gamma_h = -i eps_h^* . d eps_h   (eta contraction)
        nb = nbar(phi); h = HSTEP
        def epsl(Q): return eps_low(Q)
        eu = eps_up(P)
        der = {}
        for hh in (2, -2):
            vals = {s: eps_low([P[i] + s*h*nb[i] for i in range(3)])[hh] for s in (-2, -1, 1, 2)}
            der[hh] = [[(vals[-2][a][b] - 8*vals[-1][a][b] + 8*vals[1][a][b] - vals[2][a][b])/(12*h) for b in range(3)] for a in range(3)]
        Gam = {hh: (-I*sum(mp.conj(eu[hh][a][b])*der[hh][a][b] for a in range(3) for b in range(3))) for hh in (2, -2)}
        M = {}
        for th in ('B', 'F'):
            H = Hmat(th, P, phi)
            for a in (2, -2):
                for b in (2, -2):
                    Fa = (lambda Q, a=a: {a: mp.mpf(1), -a: mp.mpf(0)})
                    Fb = (lambda Q, b=b: {b: mp.mpf(1), -b: mp.mpf(0)})
                    M[(th, a, b)] = D1(th, Fa, Fb, P, phi)
            # my convention psi_h = eps_h^* . t  =>  predicted M = -(1/2){Gamma, H}
            for a in (2, -2):
                for b in (2, -2):
                    pred = -(Gam[a] + Gam[b])/2*H[(a, b)]
                    mx['conn'] = max(mx['conn'], abs(M[(th, a, b)] - pred))
            mx['flip'] = max(mx['flip'], abs(M[(th, 2, -2)]), abs(M[(th, -2, 2)]))
            mx['herm'] = max(mx['herm'], abs(M[(th, 2, -2)] - mp.conj(M[(th, -2, 2)])))
        mx['BF'] = max(mx['BF'], abs(M[('B', 2, 2)] - M[('F', 2, 2)]), abs(M[('B', -2, -2)] - M[('F', -2, -2)]))
        mx['gsum'] = max(mx.get('gsum', 0), abs(Gam[2] + Gam[-2]))
    log(f"(4) M (E1 on constant F): max|M_flip| = {mp.nstr(mx['flip'],3)}; max|M_B,diag - M_F,diag| = {mp.nstr(mx['BF'],3)};"
        f" max|M - (-1/2){{Gamma,H}}| = {mp.nstr(mx['conn'],3)}; |Gamma_+ + Gamma_-| = {mp.nstr(mx['gsum'],3)}")
    # ---- (5) end-to-end: pointwise O(t4) T-block term under EE ----
    def packet(rng):
        cs = [(mp.mpc(rng.gauss(0, 1), rng.gauss(0, 1)), [mp.mpf(rng.uniform(0.8, 2.0)), mp.mpf(rng.uniform(-.8, .8)), mp.mpf(rng.uniform(-.8, .8))],
               [mp.mpf(rng.gauss(0, 1)) for _ in range(3)]) for _ in range(3)]
        def g(P):
            return sum(c*mp.exp(-sum((P[i]-c0[i])**2 for i in range(3)) + I*sum(kk[i]*P[i] for i in range(3))) for c, c0, kk in cs)
        return g
    R = [1, 1, -1]
    def comp_F(g, mu, nu, refl=False):
        def F(P):
            el = eps_low(P); gv = g(P)
            if refl: return {h: R[mu]*R[nu]*el[h][mu][nu]*gv for h in (2, -2)}
            return {h: el[h][mu][nu]*gv for h in (2, -2)}
        return F
    def applyH(S, Fin, phi):
        def F(P):
            H = struct(Hmat, S, P, phi); v = Fin(P)
            return {a: sum(H[(a, b)]*v[b] for b in (2, -2)) for a in (2, -2)}
        return F
    def D1S(S, Fb, Fk, P, phi):
        if S in ('B', 'F'): return D1(S, Fb, Fk, P, phi)
        A, B = D1('B', Fb, Fk, P, phi), D1('F', Fb, Fk, P, phi)
        return (A + B)/2 if S == 'E' else (A - B)/2
    def term(kind, S3, S1, P, ph1, ph3, g2, g4, refl=False):
        """eta-contracted sum over polarisation component pairs (mu nu) of the pointwise density in P."""
        tot = mp.mpc(0)
        for mu in range(3):
            for nu in range(3):
                wgt = ETA[mu]*ETA[nu]
                F2 = comp_F(g2, mu, nu); F4 = comp_F(g4, mu, nu, refl)
                if kind == 'E3E1':
                    H3 = struct(Hmat, S3, P, ph3); H1 = struct(Hmat, S1, P, ph1); a4 = F4(P); a2 = F2(P)
                    v = {a: sum(H1[(a, b)]*a2[b] for b in (2, -2)) for a in (2, -2)}
                    v = {a: sum(H3[(a, b)]*v[b] for b in (2, -2)) for a in (2, -2)}
                    tot += wgt*sum(mp.conj(a4[a])*v[a] for a in (2, -2))
                elif kind == 'E3E11':   # <Psi4| E3 P_T E1^(1)(n1) |Psi2>
                    tot += wgt*D1S(S1, applyH(S3, F4, ph3), F2, P, ph1)
                elif kind == 'E31E1':   # <Psi4| E3^(1)(n3) P_T E1 |Psi2>
                    tot += wgt*D1S(S3, F4, applyH(S1, F2, ph1), P, ph3)
        return tot
    rows = []
    for ipt in range(int(sys.argv[1]) if len(sys.argv) > 1 else 6):
        P = randP(rng); ph1, ph3 = mp.mpf(rng.uniform(0, 6.28)), mp.mpf(rng.uniform(0, 6.28))
        g2, g4 = packet(rng), packet(rng)
        for kind in ('E3E1', 'E3E11', 'E31E1'):
            EEv = term(kind, 'E', 'E', P, ph1, ph3, g2, g4)
            EW = term(kind, 'E', 'W', P, ph1, ph3, g2, g4); WE = term(kind, 'W', 'E', P, ph1, ph3, g2, g4)
            WW = term(kind, 'W', 'W', P, ph1, ph3, g2, g4)
            EWr = term(kind, 'E', 'W', P, ph1, ph3, g2, g4, refl=True); WEr = term(kind, 'W', 'E', P, ph1, ph3, g2, g4, refl=True)
            row = dict(pt=ipt, kind=kind, EE=float(abs(EEv)), X_rel=float(abs(EW + WE)/abs(EEv)),
                       EW_rel=float(abs(EW)/abs(EEv)), WE_rel=float(abs(WE)/abs(EEv)),
                       WW_rel=float(abs(WW)/abs(EEv)), refl_rel=float(abs(EWr + WEr)/abs(EEv)))
            rows.append(row)
            log(f"(5) pt {ipt} {kind:6s}: |I[E,E]| = {row['EE']:.3e}; |X|/|I[E,E]| = {row['X_rel']:.2e} "
                f"(EW {row['EW_rel']:.1e}, WE {row['WE_rel']:.1e}); controls: W x W {row['WW_rel']:.3f}, reflected (2,4) X {row['refl_rel']:.3f}")
    json.dump(rows, open('/home/user/MAJOR/ads4/code/indep3/cache/i3_freefield_rows.json', 'w'), indent=1)
    log(f"(5) SUMMARY: max |X|/|I[E,E]| = {max(r['X_rel'] for r in rows):.2e}; max EW, WE separately = "
        f"{max(max(r['EW_rel'], r['WE_rel']) for r in rows):.2e}; WW range {min(r['WW_rel'] for r in rows):.3f}-{max(r['WW_rel'] for r in rows):.3f};"
        f" reflected-control range {min(r['refl_rel'] for r in rows):.3f}-{max(r['refl_rel'] for r in rows):.3f}")
    open('/home/user/MAJOR/ads4/code/indep3/cache/i3_freefield.txt', 'w').write("\n".join(out) + "\n")

if __name__ == '__main__':
    main()
