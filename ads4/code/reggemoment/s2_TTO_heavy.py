"""
Stage 2(a), Test 4: <T T O_l> in d = 3, parity even, conservation at both stress tensors (no Bose symmetry imposed),
at numerical (Delta, l) in the Regge regime l/Delta fixed.

 * 14 monomials V1^{m1} V2^{m2} V5^{m5} H12^{n12} H15^{n15} H25^{n25};  one linear relation among them in d = 3.
 * conserved structures: null space of the divergence equations at points 1 and 2  (expected dimension 2, modulo the relation).
 * vertex in the symmetric frame: tau(z1, z2') = structure/(scalar factor), z2 = I(x2) z2' (operator 2 sits at x2 = x/x^2 -> infinity).
   Frames at pbarhat (operator 1) and phat (operator 2), parallel along the geodesic of H^2:
        n_par(pbarhat) = (-l, Delta, 0)/m,  n_par(phat) = (l, Delta, 0)/m,  e_perp = e_y,  eps_pm = (n_par +- i e_perp)/sqrt2.
   Null polarisations z(th) = phat-direction + cos(th) n_par + sin(th) e_perp at each point; Fourier analysis in (th1, th2):
   a transverse (helicity +-2, +-2) vertex has only the harmonics e^{+-2i th1} e^{+-2i th2}.
"""
import mpmath as mp, random, sys
from threept import *
from s2_three_point_heavy import frame, hyper, EY, E0, N5

mp.mp.dps = 60
MONOS = monomials(2, 2)
NULLS = [(5, 3, 4), (5, 4, 3), (13, 5, 12), (13, 12, 5), (17, 8, 15), (17, 15, 8), (5, -3, 4), (25, 7, 24), (29, 20, 21)]

def poly(mono, x1, z1, x2, z2, x5, z5, spins):
    return poly_even(mono, blocks(x1, z1, x2, z2, x5, z5), spins)

def rand_pt(rng, shift):
    return [mp.mpf(rng.randint(-9, 9))/7 + shift for _ in range(3)]

def conserved_space(Delta, l, seed=5):
    dims, spins = (mp.mpf(3), mp.mpf(3), Delta), (2, 2, l)
    rng = random.Random(seed)
    rows = []
    for trial in range(8):
        x1, x2, x5 = rand_pt(rng, 0), rand_pt(rng, 4), rand_pt(rng, -3)
        z5 = [mp.mpf(5), 3, 4] if trial % 2 else [mp.mpf(13), -5, 12]
        zo = [mp.mpf(c) for c in rng.choice(NULLS)]                  # polarisation of the other stress tensor: must be null
        base = log_scalar(blocks(x1, E3[0], x2, E3[0], x5, z5), dims, spins)
        # conservation at 1
        cols = []
        for mono in MONOS:
            S = lambda xx, zz, mono=mono: mp.exp(log_scalar(blocks(xx, zz, x2, zo, x5, z5), dims, spins) - base)*poly(mono, xx, zz, x2, zo, x5, z5, spins)
            cols.append(divergence_at_1(S, x1, 2))
        for r in range(3):
            rows.append([cols[a][r] for a in range(len(MONOS))])
        # conservation at 2
        cols = []
        for mono in MONOS:
            S = lambda xx, zz, mono=mono: mp.exp(log_scalar(blocks(x1, zo, xx, zz, x5, z5), dims, spins) - base)*poly(mono, x1, zo, xx, zz, x5, z5, spins)
            cols.append(divergence_at_1(S, x2, 2))
        for r in range(3):
            rows.append([cols[a][r] for a in range(len(MONOS))])
    # rescale rows
    rows = [[v/max(abs(u) for u in r) for v in r] for r in rows]
    sv, V = nullspace(rows, len(MONOS))
    # linear relation among monomials (null polarisations!)
    ev = []
    for _ in range(24):
        x1, x2, x5 = rand_pt(rng, 0), rand_pt(rng, 4), rand_pt(rng, -3)
        z1, z2, z5 = ([mp.mpf(c) for c in rng.choice(NULLS)] for _ in range(3))
        r = [poly(mono, x1, z1, x2, z2, x5, z5, spins) for mono in MONOS]
        ev.append([v/max(abs(u) for u in r) for v in r])
    sv2, V2 = nullspace(ev, len(MONOS))
    return sv, V, sv2, V2

def vertex(coeffs, Delta, l, X, Y):
    spins = (2, 2, l)
    fr = frame(Delta, l, X, Y)
    x2 = fr['x2']; X2 = dot(x2, x2)
    def inv(z):                                    # z2 = I(x2) z2'
        s = 2*dot(x2, z)/X2
        return [z[i] - s*x2[i] for i in range(3)]
    def tau(z1, z2p):
        z2 = inv(z2p)
        B = blocks(fr['x1'], z1, x2, z2, fr['x5'], N5)
        return sum(coeffs[a]*poly_even(MONOS[a], B, spins) for a in range(len(MONOS)))/X2**2     # 1/X2^2: weight of the inversion for spin 2
    return tau

def analyse(Delta, l, X, Y, basis):
    m = mp.sqrt(Delta**2 - l**2)
    ph = [Delta/m, l/m, 0]; pbh = [Delta/m, -l/m, 0]
    n1 = [-l/m, Delta/m, 0]; n2 = [l/m, Delta/m, 0]
    nth = 8
    out = []
    for coeffs in basis:
        tau = vertex(coeffs, Delta, l, X, Y)
        vals = {}
        for i in range(nth):
            for j in range(nth):
                t1, t2 = 2*mp.pi*i/nth, 2*mp.pi*j/nth
                z1 = [pbh[k] + mp.cos(t1)*n1[k] + mp.sin(t1)*EY[k] for k in range(3)]
                z2 = [ph[k] + mp.cos(t2)*n2[k] + mp.sin(t2)*EY[k] for k in range(3)]
                vals[(i, j)] = tau(z1, z2)
        four = {}
        for a in range(-3, 5):
            for b in range(-3, 5):
                four[(a, b)] = sum(vals[(i, j)]*mp.exp(-2j*mp.pi*(a*i + b*j)/nth) for i in range(nth) for j in range(nth))/nth**2
        out.append(four)
    return out

if __name__ == "__main__":
    X = [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5')]; Y = [mp.mpf('0.9'), mp.mpf('-0.3'), mp.mpf('0.2')]
    ratio = mp.mpf(2)/5
    print(f"# <T T O_l>, parity even; l/Delta = 0.4; {len(MONOS)} monomials")
    for Delta in (mp.mpf(60), mp.mpf(600), mp.mpf(6000)):
        l = ratio*Delta
        sv, V, sv2, V2 = conserved_space(Delta, l)
        print(f"Delta = {mp.nstr(Delta, 5)}, l = {mp.nstr(l, 5)}:  smallest singular values of the conservation system: {[mp.nstr(s, 3) for s in sv[-5:]]}")
        print(f"     smallest singular values of the monomial evaluation matrix (linear relations): {[mp.nstr(s, 3) for s in sv2[-3:]]}")
        n = len(MONOS)
        null = [[V[n-1-i, a] for a in range(n)] for i in range(3)]
        rel = [V2[n-1, a] for a in range(n)]
        # remove the component along the trivial relation, then orthonormalise two conserved structures
        def proj(vv):
            c = sum(vv[a]*rel[a] for a in range(n))/sum(r*r for r in rel)
            return [vv[a] - c*rel[a] for a in range(n)]
        null = [proj(vv) for vv in null]
        # Gram-Schmidt: two independent vectors
        basis = []
        for vv in null:
            for b in basis:
                c = sum(vv[a]*b[a] for a in range(n))
                vv = [vv[a] - c*b[a] for a in range(n)]
            nr = mp.sqrt(sum(x*x for x in vv))
            if nr > mp.mpf(10)**(-15):
                basis.append([x/nr for x in vv])
        print(f"     dimension of the conserved space modulo the relation: {len(basis)}")
        fours = analyse(Delta, l, X, Y, basis)
        for ib, four in enumerate(fours):
            tot = mp.sqrt(sum(abs(v)**2 for v in four.values()))
            hel = {k: four[k]/tot for k in ((2, 2), (-2, -2), (2, -2), (-2, 2))}
            rest = mp.sqrt(sum(abs(v)**2 for k, v in four.items() if k not in hel))/tot
            print(f"     structure {ib}: harmonics (2,2) {mp.nstr(hel[(2,2)], 6)}  (-2,-2) {mp.nstr(hel[(-2,-2)], 6)}  (2,-2) {mp.nstr(hel[(2,-2)], 6)}  (-2,2) {mp.nstr(hel[(-2,2)], 6)};  all other harmonics / total = {mp.nstr(rest, 4)}")
        # combinations supported on a single helicity orbit
        A = mp.matrix([[fours[0][(2, 2)], fours[1][(2, 2)]], [fours[0][(2, -2)], fours[1][(2, -2)]]])
        print(f"     det of [(2,2),(2,-2)] x [structures] / norms = {mp.nstr(mp.det(A)/(mp.norm(A[0, :])*mp.norm(A[1, :])), 6)}  (non-zero: the two conserved structures span both helicity orbits)")
