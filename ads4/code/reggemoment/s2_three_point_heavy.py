"""
Stage 2(a): heavy limit (Delta, l -> infinity, l/Delta = tanh(beta/2) fixed) of spinning three-point functions in
the symmetric frame of CHMRS (3.23):
      (x1, x2, x3, x4; x5) = (-y, x/x^2, y, -x/x^2; e),   x = X/m, y = Y/m,   n = z5 = (1, e_z),
      p = (Delta, l e_z),  pbar = (Delta, -l e_z),  phat, pbarhat in H^2, cosh(beta) = -phat.pbarhat.

Test 1  <J phi O_l>   (conserved current at x1 = -y):  conserved parity-even combination  c1 V1 V5^l + c2 H15 V5^{l-1}
        and the parity-odd structure.  "Vertex vector" t^mu := d/dz1_mu [structure / |V5|^l].
Test 2  <T phi O_l>   (stress tensor at x1): three even monomials, one conserved combination; quadratic form in z1.
Test 3  the same for the Lorentzian shadow operator (Delta, l) -> (3 - Delta, -1 - l) at (x3, x4; x5): right vertex.
Test 4  <T T O_l>: conservation at both points, parity even; helicity matrix of the two conserved structures.
Prediction (derived in REGGE_MOMENT_TTTT.md): the vertex is, up to O(1/m), a transverse tensor at the point of H^2
associated with the Regge pair (pbarhat for operators 1,3; phat for 2,4):  helicity +-l_1 only.
"""
import sys
import mpmath as mp
from threept import *

mp.mp.dps = 50
EZ = [0, mp.mpf(1), 0]; EY = [0, 0, mp.mpf(1)]; E0 = [mp.mpf(1), 0, 0]
N5 = [mp.mpf(1), mp.mpf(1), 0]

def frame(Delta, l, X, Y, m=None):
    m = mp.sqrt(Delta**2 - l**2) if m is None else m
    x = [c/m for c in X]; y = [c/m for c in Y]
    x2 = [c/dot(x, x) for c in x]                  # = -x/(-x^2): large, past
    x4 = [-c/dot(x, x) for c in x]
    return dict(x1=[-c for c in y], x2=x2, x3=y, x4=x4, x5=E0, x=x, y=y, m=m)

def hyper(Delta, l):
    m = mp.sqrt(Delta**2 - l**2)
    ph = [Delta/m, l/m, 0]; pbh = [Delta/m, -l/m, 0]
    npar_bar = [l/m, -Delta/m, 0]                  # unit tangent at pbarhat in the (t,z) plane (pointing away from phat)
    npar = [l/m, Delta/m, 0]                       # unit tangent at phat (pointing away from pbarhat)
    return ph, pbh, npar, npar_bar

def conserved_coeffs_J(dims, spins, kind='even', seed=1):
    """coefficients of the conserved <J phi O> (l1=1) or <T phi O> (l1=2) structures, generic configuration."""
    import random
    rng = random.Random(seed)
    l1 = spins[0]
    rows = []
    monos = [(0, k, 0) for k in range(l1+1)] if kind == 'even' else list(range(l1))
    for trial in range(3):
        x1 = [mp.mpf(rng.randint(-9, 9))/7 for _ in range(3)]
        x2 = [mp.mpf(rng.randint(-9, 9))/5 + 3 for _ in range(3)]
        x5 = [mp.mpf(rng.randint(-9, 9))/11 - 2 for _ in range(3)]
        z5 = [mp.mpf(5), 3, 4] if trial % 2 == 0 else [mp.mpf(13), -5, 12]
        base = log_scalar(blocks(x1, E3[0], x2, E3[0], x5, z5), dims, spins)
        cols = []
        for mono in monos:
            def S(xx, zz, mono=mono):
                B = blocks(xx, zz, x2, E3[0], x5, z5)
                pre = mp.exp(log_scalar(B, dims, spins) - base)
                if kind == 'even':
                    return pre*poly_even(mono, B, spins)
                return pre*poly_odd(mono, B, eps15(xx, zz, x2, x5, z5), spins)
            cols.append(divergence_at_1(S, x1, l1))
        for r in range(len(cols[0])):
            rows.append([cols[a][r] for a in range(len(monos))])
    sv, V = nullspace(rows, len(monos))
    return monos, sv, V

def run_J_or_T(l1, Delta, l, D2, X, Y, shadow=False, verbose=True):
    D1 = mp.mpf(l1 + 1)
    D5, l5 = (3 - Delta, -1 - l) if shadow else (Delta, l)
    dims, spins = (D1, D2, D5), (l1, 0, l5)
    fr = frame(Delta, l, X, Y)
    xa, xb = (fr['x3'], fr['x4']) if shadow else (fr['x1'], fr['x2'])
    ph, pbh, npar, npar_bar = hyper(Delta, l)
    res = {}
    for kind in ('even', 'odd'):
        monos, sv, V = conserved_coeffs_J(dims, spins, kind)
        nnull = sum(1 for s in sv if s < sv[0]*mp.mpf(10)**(-20)) + (len(monos) - len(sv))
        cvec = [V[len(monos)-1, a] for a in range(len(monos))]       # right-singular vector of the smallest singular value
        def poly(z):
            B = blocks(xa, z, xb, E3[0], fr['x5'], N5)
            if kind == 'even':
                return sum(cvec[a]*poly_even(monos[a], B, spins) for a in range(len(monos)))
            return sum(cvec[a]*poly_odd(monos[a], B, eps15(xa, z, xb, fr['x5'], N5), spins) for a in range(len(monos)))
        res[kind] = (tensor_at_1(poly, l1), sv, cvec)
    return res, (ph, pbh, npar, npar_bar)

def report_vector(tag, t, ph, pbh, npar_bar):
    nrm = mp.sqrt(abs(dot(t, t)))
    comps = [dot(t, pbh)/nrm, dot(t, npar_bar)/nrm, dot(t, EY)/nrm, dot(t, ph)/nrm]
    print(f"   {tag}: t.pbarhat = {mp.nstr(comps[0], 6):>12s}  t.n_par = {mp.nstr(comps[1], 8):>12s}  t.e_perp = {mp.nstr(comps[2], 8):>12s}   [t.phat = {mp.nstr(comps[3], 5)}]")

def report_form(tag, s, pbh, npar_bar):
    # components of the quadratic form in the frame (pbarhat, n_par, e_perp)
    fr3 = [pbh, npar_bar, EY]
    def comp(u, w):
        ul = [ETA[i]*u[i] for i in range(3)]; wl = [ETA[i]*w[i] for i in range(3)]
        return sum(s[a][b]*ul[a]*wl[b] for a in range(3) for b in range(3))
    # polynomial Q(z) = s_{ab} z^a z^b ; component along frame vectors f,g:  Q = sum q_fg (z.f)(z.g) with q_fg = s(f^,g^) and dual basis
    dual = [[-c for c in pbh], npar_bar, EY]        # dual vectors: (z.pbh) coefficient etc.  z = -(z.pbh) pbh + (z.n) n + (z.e) e
    q = [[sum(s[a][b]*dual[i][a]*dual[j][b] for a in range(3) for b in range(3)) for j in range(3)] for i in range(3)]
    # modulo z^2 = -(z.p)^2 + (z.n)^2 + (z.e)^2:  transverse-traceless  <=>  q_pp : q_nn : q_ee = (c - k) : (c + k)... we print the
    # invariant content: after removing a multiple of z^2 so that q_pp -> 0:
    k = q[0][0]                                     # add k*z^2:  q_pp -> 0, q_nn -> q_nn + k, q_ee -> q_ee + k
    qnn, qee = q[1][1] + k, q[2][2] + k
    sc = max(abs(qnn), abs(qee), abs(q[1][2]))
    print(f"   {tag}: mod z^2:  q_nn = {mp.nstr(qnn/sc, 8)}, q_ee = {mp.nstr(qee/sc, 8)}, q_ne = {mp.nstr(q[1][2]/sc, 8)};  non-transverse: q_pn = {mp.nstr(q[0][1]/sc, 5)}, q_pe = {mp.nstr(q[0][2]/sc, 5)}")

if __name__ == "__main__":
    X = [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5')]; Y = [mp.mpf('0.9'), mp.mpf('-0.3'), mp.mpf('0.2')]
    D2 = mp.mpf('1.7')
    ratio = mp.mpf(2)/5                                        # l/Delta  (beta = 2 atanh(0.4) = 0.847)
    print("# x = X/m, y = Y/m with X =", [mp.nstr(c, 3) for c in X], " Y =", [mp.nstr(c, 3) for c in Y], "; l/Delta = 0.4; Delta_phi(2) = 1.7")
    print("Test 1: <J phi O_l>, left vertex (operator 1 at x1 = -y).  Expect: even -> n_par(pbarhat), odd -> e_perp, t.pbarhat = O(1/m)")
    for Delta in (mp.mpf(50), mp.mpf(500), mp.mpf(5000)):
        res, (ph, pbh, npar, npar_bar) = run_J_or_T(1, Delta, ratio*Delta, D2, X, Y)
        for kind in ('even', 'odd'):
            t, sv, cvec = res[kind]
            tu = [ETA[i]*t[i] for i in range(3)]               # raise the index
            report_vector(f"Delta={mp.nstr(Delta, 5):>7s} {kind:4s}", tu, ph, pbh, npar_bar)
        print(f"      conserved even combination (c1 : c2) = {mp.nstr(res['even'][2][0]/res['even'][2][1], 10)};  singular values {[mp.nstr(s, 3) for s in res['even'][1]]}")
    print("Test 2: <T phi O_l>, left vertex.  Expect (mod z^2): even -> (z.n_par)^2 - (z.e_perp)^2, odd -> (z.n_par)(z.e_perp), i.e. helicity +-2 only")
    for Delta in (mp.mpf(50), mp.mpf(500), mp.mpf(5000)):
        res, (ph, pbh, npar, npar_bar) = run_J_or_T(2, Delta, ratio*Delta, D2, X, Y)
        for kind in ('even', 'odd'):
            s, sv, cvec = res[kind]
            report_form(f"Delta={mp.nstr(Delta, 5):>7s} {kind:4s}", s, pbh, npar_bar)
        print(f"      singular values (even) {[mp.nstr(s, 3) for s in res['even'][1]]}  (odd) {[mp.nstr(s, 3) for s in res['odd'][1]]}")
    print("Test 3: right vertex with the Lorentzian shadow operator (3-Delta, -1-l), operator 3 at x3 = +y")
    for l1 in (1, 2):
        for Delta in (mp.mpf(50), mp.mpf(5000)):
            res, (ph, pbh, npar, npar_bar) = run_J_or_T(l1, Delta, ratio*Delta, D2, X, Y, shadow=True)
            for kind in ('even', 'odd'):
                t, sv, cvec = res[kind]
                if l1 == 1:
                    report_vector(f"J  Delta={mp.nstr(Delta, 5):>7s} {kind:4s}", [ETA[i]*t[i] for i in range(3)], ph, pbh, npar_bar)
                else:
                    report_form(f"T  Delta={mp.nstr(Delta, 5):>7s} {kind:4s}", t, pbh, npar_bar)
