r"""
hk2_k1spin (HK_RESIDUAL.md §2): K1' (uniform Regge-box bound, LEMMAS_CLOSE §2) on a SPINNING example: the exact free-scalar <TTTT> in
d = 3 (conformally coupled scalar, N copies; connected O(N) + disconnected O(N^2)), Wick engine of resultq/og_free.py (read-only import:
DFUN link derivatives, VERTEX, MATCHES, correlator).

Frame (Caron-Huot / KKSDZ, as lem2_k1prime): light-cone coordinates x^+- = x^1 +- t in the (t, x^1) plane, x^2 transverse,
    x1 = (+1, +1), x2 = -x1, x4 = (rho, rhob), x3 = -x4;   Regge sheet: 0 < rho < 1 < rhob;  rho_b := 1/rhob;  box B_c = {rho, rho_b <= c}.
Rindler reflection Theta: (t, x^1, x^2) -> (-t, -x^1, x^2) (x2 = Theta x1, x3 = Theta x4), polarisation z -> R z*, R = diag(-1, -1, 1).
Only the pairs (14) and (23) are timelike.  Wightman orderings <Theta[P] Q>, P, Q in {O1 O4, O4 O1}, fix the i0 of those two links:
    <O3 O2 O1 O4> (+,+) and <O2 O3 O4 O1> (-,-): Regge sheet (rhob -+ i0);   <O2 O3 O1 O4> (+,-), <O3 O2 O4 O1> (-,+): first sheet,
with x14^2 -> |x14^2| e^{i pi s14}, x23^2 -> |x23^2| e^{i pi s23}.
Cauchy-Schwarz (Rindler positivity):  |<Theta[P] Q>|^2 <= <Theta[P] P> <Theta[Q] Q>.  For given polarisations (z1, z2, z3, z4):
Q carries (z1, z4), P carries (R z2, R z3) at (x1, x4) (real null polarisations).  Normalisation: each Wightman value is divided by
sqrt(N_P N_Q), N = identity (s-channel, (12)(34) disconnected) part of the positive comparison correlator; then
    |G| <= sqrt(G_P G_Q) <= M := (G_P + G_Q)/2,     G_P := <Theta[P]P>/N_P >= 1,   G_Q := <Theta[Q]Q>/N_Q >= 1.
Polarisations "normalised in each pair's rest frame" (KKSDZ; EESTEP §7.3): z1, z2 = (1, cos a, sin a) in the lab (pair (12) at rest);
z3, z4 = Lambda_y (1, cos a, sin a), Lambda_y the boost along x^1 taking e_1 to x4/|x4|, y = (1/2) log(rhob/rho).
Checks:
  (s) scalar control O = phi^2 (same engine, VERTEX -> 1): the first-sheet orderings equal the all-spacelike value at (rho, rho_b),
      i.e. the comparison correlators are first-sheet correlators at the same cross-ratios (as for lem2's scalars);
  (i) Cauchy-Schwarz at every box point and polarisation set (implementation / phase check; it is a theorem);
  (ii) B(c) := sup over the box and polarisations of max(G_P, G_Q) (both orderings): finite, attained in the interior; the edge limits
      M -> 1 (identity) as rho_b -> 0 or rho -> 0; sup |G_Regge| and sup |dDisc G| / 2 against B(c), uniformly in eta;
  (iii) lab-frame polarisations for 3, 4 (not boosted with their pair): is the bound sensitive to the polarisation frame?
Output: cache/hk2_k1spin.txt.   ~ 2-4 min, one core.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
import numpy as np
import og_free as OG

OUT = []
T0 = time.time()
R = np.array([-1.0, -1.0, 1.0])
CLASSES = ('s', 't', 'u', '1234', '1243', '1324')


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def dot(a, b):
    return -a[..., 0]*b[..., 0] + a[..., 1]*b[..., 1] + a[..., 2]*b[..., 2]


def lc(xp, xm):
    """light-cone (x^+, x^-) -> (t, x^1, x^2)."""
    xp = np.asarray(xp, float); xm = np.asarray(xm, float)
    return np.stack([(xp - xm)/2, (xp + xm)/2, 0*xp], -1)


def nvec(a):
    return np.array([1.0, np.cos(a), np.sin(a)])


def boost(y):
    """boost along x^1 with rapidity y: e_1 -> (sinh y, cosh y, 0) (a spacelike unit vector at rapidity y)."""
    y = np.asarray(y, float)
    L = np.zeros(y.shape + (3, 3))
    L[..., 0, 0] = np.cosh(y); L[..., 0, 1] = np.sinh(y); L[..., 1, 0] = np.sinh(y); L[..., 1, 1] = np.cosh(y); L[..., 2, 2] = 1
    return L


def wightman(rho, rhob, pols, s14, s23, N=1.0, vertex=None):
    """<T T T T> (free scalar, N copies) at x1 = (1,1), x4 = (rho, rhob), x2 = -x1, x3 = -x4 with polarisations pols = {i: (...,3)};
    s14, s23 in {0, +1, -1}: phase of the (14), (23) bases (0: principal, the base must be > 0).  Returns (full, identity part s)."""
    x1 = lc(np.ones_like(rho), np.ones_like(rho)); x4 = lc(rho, rhob)
    X = {1: x1, 4: x4}
    X[2] = x1.copy(); X[2][..., 0] *= -1; X[2][..., 1] *= -1
    X[3] = x4.copy(); X[3][..., 0] *= -1; X[3][..., 1] *= -1
    L = {}
    for (i, j) in OG.PAIRS:
        Xij = X[i] - X[j]
        Y = dot(Xij, Xij)
        s = {(1, 4): s14, (2, 3): s23}.get((i, j), 0)
        if s == 0:
            assert np.all(Y > 0), (i, j)
            r = Y**-0.5 + 0j
        else:
            r = np.abs(Y)**-0.5*np.exp(-0.5j*np.pi*s)
        A = dot(pols[i], Xij); B = dot(pols[j], Xij); ab = dot(pols[i], pols[j])
        L[(i, j)] = {mn: OG.DFUN[mn](r, Y + 0j, A, B, ab, 1.0) for mn in OG.DFUN}
    old = OG.VERTEX
    if vertex is not None:
        OG.VERTEX = vertex
    try:
        C = OG.correlator(L, CLASSES)
    finally:
        OG.VERTEX = old
    disc = C['s'] + C['t'] + C['u']; conn = C['1234'] + C['1243'] + C['1324']
    return N**2*disc + N*conn, N**2*C['s']


def pol_sets(angles, rho, rhob, frame='pair'):
    a1, a2, a3, a4 = angles
    y = -0.5*np.log(rhob/rho)               # x4/|x4| = (-sinh y', cosh y', 0), y' = log(rhob/rho)/2
    L = boost(y) if frame == 'pair' else boost(0*y)
    z1 = np.broadcast_to(nvec(a1), rho.shape + (3,)); z2 = np.broadcast_to(nvec(a2), rho.shape + (3,))
    z3 = np.einsum('...ij,j->...i', L, nvec(a3)); z4 = np.einsum('...ij,j->...i', L, nvec(a4))
    return {1: z1, 2: z2, 3: z3, 4: z4}


def evaluate(rho, rb, angles, N=1.0, frame='pair', vertex=None):
    """all four orderings + the comparison correlators, normalised.  rho, rb arrays (rhob = 1/rb)."""
    rhob = 1/rb
    p = pol_sets(angles, rho, rhob, frame)
    # Q = (z1, z4) at (x1, x4); P = (R z2, R z3) at (x1, x4)
    pQ = {1: p[1], 4: p[4], 2: p[1]*R, 3: p[4]*R}
    pP = {1: p[2]*R, 4: p[3]*R, 2: p[2], 3: p[3]}
    out = {}
    for tag, pp in (('Q', pQ), ('P', pP)):
        for (s14, s23) in ((+1, -1), (-1, +1)):
            W, Nid = wightman(rho, rhob, pp, s14, s23, N, vertex)
            out[(tag, s14, s23)] = (W, Nid)
    NQ = out[('Q', 1, -1)][1]; NP = out[('P', 1, -1)][1]
    norm = np.sqrt(NP.real*NQ.real)
    G = {}
    for (s14, s23) in ((+1, +1), (-1, -1), (+1, -1), (-1, +1)):
        W, _ = wightman(rho, rhob, p, s14, s23, N, vertex)
        G[(s14, s23)] = W/norm
    comp = {k: v[0]/v[1] for k, v in out.items()}       # G_P, G_Q per first-sheet ordering
    return G, comp, NP, NQ, out


def scalar_control():
    pr('(s) scalar control O = phi^2 (vertex 1): first-sheet orderings vs the all-spacelike value at (rho, rho_b)')
    vtx = [(1.0, (0, 0))]
    g = np.array([1e-4, 1e-2, 0.2, 0.6, 0.9])
    Rr, Bb = np.meshgrid(g, g, indexing='ij'); Rr = Rr.ravel(); Bb = Bb.ravel()
    ang = (0.3, 1.1, 2.0, 4.0)
    G, comp, NP, NQ, _ = evaluate(Rr, Bb, ang, N=1.0, vertex=vtx)
    p0 = {i: np.broadcast_to(nvec(0.0), Rr.shape + (3,)) for i in (1, 2, 3, 4)}
    Wsp, Nsp = wightman(Rr, Bb, p0, 0, 0, 1.0, vtx)            # x4 = (rho, rho_b): all spacelike, first sheet
    Gsp = Wsp/Nsp
    d1 = np.max(np.abs(G[(1, -1)]/Gsp - 1)); d2 = np.max(np.abs(G[(-1, 1)]/Gsp - 1))
    dc = max(np.max(np.abs(comp[(t, s, -s)]/Gsp - 1)) for t in 'PQ' for s in (1, -1))
    pr(f'    max |G_(+-)/G_E(rho, rho_b) - 1| = {d1:.1e}, |G_(-+)/G_E - 1| = {d2:.1e}, comparison correlators {dc:.1e}')
    # Regge value for phi^2 (Delta = 1): G_R = 1 + 2 u + 2 u/v e^{-2 pi i s} ... check modulus bound only
    rg = np.abs(G[(1, 1)])
    pr(f'    |G_Regge| <= G_E(rho, rho_b) pointwise: {bool(np.all(rg <= Gsp.real*(1 + 1e-12)))}  (max ratio {np.max(rg/Gsp.real):.4f})')


def box(c, n, lo=1e-6):
    g = np.geomspace(lo, c, n)
    Rr, Bb = np.meshgrid(g, g, indexing='ij')
    return g, Rr.ravel(), Bb.ravel()


def main_scan(c, N, nang=60, n=61, frame='pair', label=''):
    rng = np.random.default_rng(11)
    angs = [tuple(rng.uniform(0, 2*np.pi, 4)) for _ in range(nang)] + [(0, 0, 0, 0), (0, np.pi, 0, np.pi), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)]
    g, Rr, Bb = box(c, n)
    cs_worst = 0.0; Bc = 0.0; argB = None; GR = 0.0; DD = 0.0; nneg = 0; minN = np.inf; edge = 0.0
    Ms = []
    pairs = {(1, 1): (('Q', 1, -1), ('P', -1, 1)), (-1, -1): (('Q', -1, 1), ('P', 1, -1)),
             (1, -1): (('Q', 1, -1), ('P', 1, -1)), (-1, 1): (('Q', -1, 1), ('P', -1, 1))}
    for ang in angs:
        G, comp, NP, NQ, raw = evaluate(Rr, Bb, ang, N=N, frame=frame)
        minN = min(minN, np.min(NP.real/np.abs(NP)), np.min(NQ.real/np.abs(NQ)))
        Mall = np.zeros_like(Rr)
        for k, (a, b) in pairs.items():
            ca, cb = comp[a], comp[b]
            nneg += int(np.sum(ca.real < 1 - 1e-9)) + int(np.sum(cb.real < 1 - 1e-9))
            cs_worst = max(cs_worst, np.max(np.abs(G[k])**2/(ca.real*cb.real)))
            Mall = np.maximum(Mall, 0.5*(ca.real + cb.real))
        if np.max(Mall) > Bc:
            Bc = np.max(Mall); argB = (Rr[np.argmax(Mall)], Bb[np.argmax(Mall)])
        GR = max(GR, np.max(np.abs(G[(1, 1)])), np.max(np.abs(G[(-1, -1)])))
        dd = 0.5*(G[(1, -1)] + G[(-1, 1)]) - 0.5*(G[(1, 1)] + G[(-1, -1)])
        DD = max(DD, np.max(np.abs(dd)))
        Ms.append(Mall)
    Ms = np.array(Ms)
    Mmax = Ms.max(0)
    # K1' edge statement: |M - 1| <= (B(c) - 1) (min(rho, rho_b)/c)^{tau_min/2}, tau_min = 1 (free scalar: phi^2 and the twist-1 currents)
    edge = np.max((Mmax - 1)/((Bc - 1)*(np.minimum(Rr, Bb)/c)**0.5))
    Mg = Mmax.reshape(n, n)
    sl = np.polyfit(np.log(g[:8]), np.log(Mg[-1, :8] - 1), 1)[0]
    pr(f'    {label} c = {c}, N = {N}, {len(angs)} polarisation sets, {n}x{n} log grid rho, rho_b in [1e-6, c]:')
    pr(f'       Cauchy-Schwarz max |G|^2/(G_P G_Q) = {cs_worst:.12f} (<= 1; = 1 when P ~ Q);  comparison values < 1: {nneg};  '
       f'min N/|N| = {minN:.3f} (identity parts > 0)')
    pr(f'       B(c) := sup max(M) = {Bc:.6g} at (rho, rho_b) = ({argB[0]:.3g}, {argB[1]:.3g});  sup |G_Regge| = {GR:.6g} (/B = {GR/Bc:.4f});  '
       f'sup |dDisc G| = {DD:.6g} (/2B = {DD/(2*Bc):.4f})')
    pr(f'       edge: max over box, pols of (M - 1)/((B - 1)(min(rho,rho_b)/c)^(1/2)) = {edge:.4f} (<= 1: K1\' edge row with tau_min = 1);'
       f'  slope dlog(M-1)/dlog rho_b at rho = c, rho_b <= {g[7]:.0e}: {sl:.4f};  M at the corner (1e-6, 1e-6) = {Mg[0, 0]:.8f}')
    return Bc, GR, DD, edge


def eta_profile(c, N):
    """M and |G_R| along the light-cone edge direction at fixed max(rho, rho_b) = c: rho_b/rho -> 0 (eta -> infinity)."""
    rng = np.random.default_rng(11)
    angs = [tuple(rng.uniform(0, 2*np.pi, 4)) for _ in range(12)]
    rb = np.geomspace(1e-8, c, 9); rho = np.full_like(rb, c)   # fixed rho = c: eta -> infinity along the light-cone edge
    eta = (rho + rb)/(2*np.sqrt(rho*rb))
    Mx = np.zeros_like(rb); Gx = np.zeros_like(rb)
    for ang in angs:
        G, comp, *_ = evaluate(rho, rb, ang, N=N)
        M = 0.5*(comp[('Q', 1, -1)].real + comp[('P', -1, 1)].real)
        Mx = np.maximum(Mx, M); Gx = np.maximum(Gx, np.abs(G[(1, 1)]))
    pr(f'    eta profile at rho = {c} (sup over 12 pol sets): eta = ' + ' '.join(f'{e:.1e}' for e in eta))
    pr('       M      = ' + ' '.join(f'{m:.4f}' for m in Mx))
    pr('       |G_R|  = ' + ' '.join(f'{m:.4f}' for m in Gx))


if __name__ == '__main__':
    pr('# hk2_k1spin: K1\' on the free-scalar <TTTT> (spin 2, Delta = 3), Rindler/Cauchy-Schwarz with generic polarisations')
    scalar_control()
    pr('\n(i)-(ii) pair-frame polarisations')
    res = {}
    for c in (0.5, 0.9):
        for N in (1.0, 10.0):
            res[(c, N)] = main_scan(c, N)
    for c in (0.5, 0.9):
        eta_profile(c, 1.0)
    pr('\n(iii) control: lab-frame polarisations for 3, 4 (not boosted with the pair)')
    main_scan(0.5, 1.0, nang=20, frame='lab', label='[lab frame]')
    pr(f'\n[{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'hk2_k1spin.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
