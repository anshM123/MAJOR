"""
g1 (M2): the CLS-D scattering crystal in D = 4 (d = 3), explicitly, for photons and gravitons.
CLS-D (arXiv:2311.04271) shock frame, eqs. (shock_frame_vecs_inCMframe), (shock_momenta_inbasisvec), with D = 4 (m = 1, nu = x):
   sbar = (s, 1, 0, i x), s_ = (s, -1, 0, i x), v1 = (i x, 0, 0, s), v2 = (0, 0, 1, 0),  s = sqrt(1 - x^2)   [components (t, x1, x2, x3)]
   p1 = p sbar + q1, p2 = p s_ - q3, p3 = -p sbar + q3, p4 = -p s_ - q1,  p = s/2,  q1 = -(i x/2) v1 + (x/2) v2,  q3 = (i x/2) v1 + (x/2) v2.
Checks / results:
 (1) kinematics are real CM kinematics, cos(theta_13) = 1 - 2x^2; shock polarisation sbar of 1 and 3 is a pure helicity state
     (self-dual field strength), coefficient |alpha| = sqrt2 x, the same for every helicity of 2, 4 (common factor);
 (2) the CLS-D kernel vectors V_2, W_2, V_4, W_4 (eq. kernel_VWvec_definition) on normalised helicity states of 2 and 4;
     photon kernels of CLS-D (eq. JJO_kernels_flatspace at D = 4: delta_{perp;D-4} = 0, j' = 2 absent): helicity tables;
 (3) H24 ~ VV + WW at the crystal (eq. H24_kernel_toflat) is nonzero only on equal all-incoming helicities of 2, 4 = physical helicity FLIP
     (consistent with tr(F2 F4) = 0 for opposite labels) -- an independent check of the helicity bookkeeping;
 (4) flip-free (physically helicity-preserving) kernels: photons: unique, VV - WW (= CLS-D K^(1)_{j'=0} at D = 4); gravitons: the j' = 0
     structures H20^{2-n}H40^{2-n}H24^n -> {P^2, PH, H^2} (P = VV, H = VV + WW); flip-free subspace and its channel ratio;
 (5) with CLS-D's polarisation-sum pairing (block sees the conjugate of the kernel's vector for 2 and 4), the channel ratio
     (opposite/same helicity of 2 relative to 1) of every flip-free kernel is (1-x^2)^{-2} (photons), (1-x^2)^{-4} (gravitons): the STRIPPED
     dictionary (d-tilde, reading B-tilde), not the x-independent ratio of the B' (reading-B) dictionary.
Usage: python3 g1_crystal4d.py -> cache/g1_crystal4d.txt
"""
import os, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(s=""):
    print(s, flush=True); out.append(str(s))

ETA = np.diag([-1.0, 1, 1, 1])
def dot(a, b): return a @ ETA @ b
LC = np.zeros((4, 4, 4, 4))
for perm in itertools.permutations(range(4)):
    sgn = np.linalg.det(np.eye(4)[list(perm)]); LC[perm] = sgn

def duality(p, e):
    """+1 if F = p^e is self-dual (*F = +iF), -1 if anti-self-dual, in the index-lowered convention used consistently."""
    F = np.outer(p, e) - np.outer(e, p)            # upper indices
    Fl = ETA @ F @ ETA                              # lower
    dual = 0.5*np.einsum('abcd,cd->ab', LC, F)      # *F_{ab} = 1/2 eps_{abcd} F^{cd}
    for s in (+1, -1):
        if np.allclose(dual, s*1j*Fl, atol=1e-10):
            return s
    return 0

def crystal(x):
    s = np.sqrt(1 - x*x)
    sbar = np.array([s, 1, 0, 1j*x]); s_ = np.array([s, -1, 0, 1j*x])
    v1 = np.array([1j*x, 0, 0, s]); v2 = np.array([0, 0, 1.0, 0])
    p = s/2; q1 = -(1j*x/2)*v1 + (x/2)*v2; q3 = (1j*x/2)*v1 + (x/2)*v2
    P = [p*sbar + q1, p*s_ - q3, -p*sbar + q3, -p*s_ - q1]
    m, nu = 1.0, x
    V2 = m**2*q1 + (m**2 - nu**2)*q3; V4 = (m**2 - nu**2)*q1 + m**2*q3
    W2 = m**2*q1 - (m**2 - nu**2)*q3; W4 = (m**2 - nu**2)*q1 - m**2*q3
    return s, sbar, P, (V2, W2, V4, W4)

def hel_vectors(pvec):
    """normalised transverse helicity vectors (real basis e_a, e_b, e_a x e_b = n) for the spatial direction of pvec (physical motion)."""
    n = np.real(pvec[1:]); n = n/np.linalg.norm(n)
    eb = np.array([0, 0, 1.0]); ea = np.cross(eb, n); ea = ea/np.linalg.norm(ea)   # ea x eb = n
    return {h: np.concatenate([[0], (ea + 1j*h*eb)/np.sqrt(2)]) for h in (+1, -1)}

pr("# g1: CLS-D scattering crystal in D = 4: helicity structure of the CFT -> flat dictionary")
for x in (0.1, 0.3, 0.5, 0.7, 0.9):
    s, sbar, P, (V2, W2, V4, W4) = crystal(x)
    p1, p2, p3, p4 = P
    assert all(abs(dot(pp, pp)) < 1e-12 for pp in P) and np.allclose(sum(P), 0)
    real = all(np.allclose(np.imag(pp), 0) for pp in P)
    c13 = np.real(p1[1:] @ (-p3[1:]))/(np.linalg.norm(np.real(p1[1:]))*np.linalg.norm(np.real(p3[1:])))
    # shock polarisations: pure helicity?  coefficient of sbar on the helicity basis of p1 (incoming) and of -p3 (outgoing 3)
    d1, d3 = duality(p1, sbar), duality(p3, sbar)
    H1 = hel_vectors(p1)
    g1 = sbar - (sbar[0]/p1[0])*p1                    # remove time component by gauge p1
    a1 = {h: np.vdot(H1[h][1:], g1[1:]) for h in (+1, -1)}
    # probe helicity states: physical directions: 2 incoming along p2, 4 outgoing along -p4 (polarisation entering amplitude: eps_h^* )
    H2 = hel_vectors(p2); H4 = hel_vectors(-p4)
    E2 = {h: H2[h] for h in (+1, -1)}                    # incoming: eps_h
    E4 = {h: np.conj(H4[h]) for h in (+1, -1)}           # outgoing: eps_h^*
    lab2 = {h: duality(p2, E2[h]) for h in (+1, -1)}     # all-incoming duality labels
    lab4 = {h: duality(p4, E4[h]) for h in (+1, -1)}
    lab1 = d1
    pr(f"x={x}: real CM kinematics {real}; cos th13 = {c13:.12f} (1-2x^2 = {1-2*x*x:.12f}); sbar: duality label (1, 3) = ({d1:+d}, {d3:+d}); "
       f"|helicity coefficients| of sbar on p1: {abs(a1[+1]):.3e}, {abs(a1[-1]):.6f} (sqrt2 x = {np.sqrt(2)*x:.6f})")
    pr(f"      all-incoming duality labels: particle 2 phys h=+1 -> {lab2[1]:+d}, h=-1 -> {lab2[-1]:+d}; particle 4 phys h=+1 -> {lab4[1]:+d}, h=-1 -> {lab4[-1]:+d}")
    # photon structures on (h2, h4) physical helicities
    def tab(A2, A4):
        return {(h2, h4): dot(A2, E2[h2])*dot(A4, E4[h4]) for h2 in (+1, -1) for h4 in (+1, -1)}
    VV = tab(V2, V4); WW = tab(W2, W4); WV = tab(W2, V4); VW = tab(V2, W4)
    H = {k: VV[k] + WW[k] for k in VV}; Pv = VV
    # CLS-D eq. (sumrule_bulk_polarizationmapall): e2, e4 are POLARISATION-SUMMED, sum_h (K.eps_h)(G.eps_h^*): the block (partial wave) is
    # evaluated with the CONJUGATE of the vector the kernel sees.  So kernel label h for particle 2 (incoming) = amplitude helicity -h;
    # for particle 4 likewise (both flip, so preserved/flip is unchanged; the channel relative to particle 1 -- whose s-bar enters the
    # amplitude directly -- is swapped).  "same" channel: amplitude helicity of 2 = helicity of 1 (both incoming; J_z = 0, d^J_00).
    h1phys = lab1
    hsame = -1 if lab2[1] == h1phys else +1               # kernel label whose conjugate has the helicity of 1
    def show(name, T):
        return f"{name}: " + "  ".join(f"({h2:+d},{h4:+d}) {T[(h2, h4)]:+.4e}" for h2 in (+1, -1) for h4 in (+1, -1))
    pr("      photon " + show("VV", VV)); pr("      photon " + show("WW", WW)); pr("      photon " + show("H=VV+WW", H))
    # preserved = (h2, h4) with h4 = h2 (physical helicity conserved along the 2 -> 4 line)
    pres = [(h, h) for h in (+1, -1)]; flip = [(h, -h) for h in (+1, -1)]
    pr(f"      H on preserved: {[f'{abs(H[k]):.1e}' for k in pres]}  on flips: {[f'{abs(H[k]):.3e}' for k in flip]}   (H24 is pure flip)")
    K = {k: VV[k] - WW[k] for k in VV}
    rat_ph = abs(K[(-hsame, -hsame)])/abs(K[(hsame, hsame)])
    pr(f"      photon flip-free VV-WW: preserved same {abs(K[(hsame, hsame)]):.4e}, opp {abs(K[(-hsame, -hsame)]):.4e}; flips "
       f"{max(abs(K[k]) for k in flip):.1e}; opp/same = {rat_ph:.6f}; (1-x^2)^2 = {(1-x*x)**2:.6f}, (1-x^2)^-2 = {(1-x*x)**-2:.6f}")
    # gravitons: structures P^2, PH, H^2 on graviton helicity states eps (x) eps -> products of photon weights
    G = {'P2': {k: Pv[k]**2 for k in Pv}, 'PH': {k: Pv[k]*H[k] for k in Pv}, 'H2': {k: H[k]**2 for k in Pv}}
    M = np.array([[G[n][k] for n in ('P2', 'PH', 'H2')] for k in flip])    # flip conditions
    ns = np.linalg.svd(M)[2].conj().T[:, np.linalg.matrix_rank(M, tol=1e-12):]
    pr(f"      graviton j'=0 span {{P^2, PH, H^2}}: flip-free subspace dimension {ns.shape[1]}")
    for i in range(ns.shape[1]):
        c = ns[:, i]
        w = {k: sum(c[j]*G[n][k] for j, n in enumerate(('P2', 'PH', 'H2'))) for k in Pv}
        if max(abs(w[k]) for k in pres) < 1e-14:
            pr(f"        basis vector {i}: zero on preserved too (trivial)"); continue
        rg = abs(w[(-hsame, -hsame)])/abs(w[(hsame, hsame)])
        pr(f"        basis vector {i}: opp/same = {rg:.6f}; (1-x^2)^4 = {(1-x*x)**4:.6f}, (1-x^2)^-4 = {(1-x*x)**-4:.6f}")
pr("# conclusion: every flip-free kernel has an x-dependent channel ratio (1-x^2)^{+-2} (photons) / (1-x^2)^{+-4} (gravitons);")
pr("#   the B' (unstripped, reading-B) dictionary needs ratio 1; the stripped B^(1) dictionary needs (1-x^2)^{-4} (+- enhanced, d-tilde).")
pr("#   (without the conjugate pairing the ratio would be (1-x^2)^{+4}: neither dictionary; the pairing is CLS-D's eq. sumrule_bulk_polarizationmapall)")
open(os.path.join(HERE, 'cache', 'g1_crystal4d.txt'), 'w').write("\n".join(out) + "\n")
