"""
g2 (M2b): structure of the CLS-D (2,4) dual structures at the D = 4 crystal (uses g1's kinematics and CLS-D's polarisation-sum pairing).
Photon-level crystal images (CLS-D eqs. H20H40_kernel_toflat, H24_kernel_toflat, JJO_kernels_flatspace at D = 4):
   P := V2 V4 (H20 H40),  H := V2 V4 + W2 W4 (H24, up to normalisation),  O := W2 V4 - V2 W4 (j' = 1; the symmetric W2V4 + V2W4 is pure flip), O1 := W2 V4 alone (CLS-D K_{j'=1}).
Checks at x = 0.05 ... 0.95:
 (a) H_flip = 2 P_flip and H_pres = 0  =>  the graviton j'=0 flip-free condition on a P^2 + b PH + c H^2 is a + 2b + 4c = 0, x-INDEPENDENT;
     every flip-free combination has preserved weight a P_pres^2 (universal shape; normalisation a);
 (b) the photon structure Q := H - 2P (= W2W4 - V2V4, CLS-D K^(1)_{j'=0} at D = 4) is exactly flip-free; Q^2 (a = 4, b = -4, c = 1) is a
     graviton j' = 0 structure flip-free at the crystal for all x;
 (c) odd sector: O and O*Q at the crystal: flip content, channel sign and ratio (expect antisymmetric, (1-x^2)^{-4} for gravitons).
Usage: python3 g2_structures.py -> cache/g2_structures.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location('g1', os.path.join(HERE, 'g1_crystal4d.py'))
src = open(os.path.join(HERE, 'g1_crystal4d.py')).read().split('pr("# g1:')[0]   # helpers only
G1 = {'__file__': os.path.join(HERE, 'g1_crystal4d.py')}; exec(compile(src, 'g1_helpers', 'exec'), G1)
crystal, hel_vectors, dot, duality = G1['crystal'], G1['hel_vectors'], G1['dot'], G1['duality']
out = []
def pr(s=""):
    print(s, flush=True); out.append(str(s))
pr("# g2: CLS-D (2,4) structures at the D = 4 crystal")
for x in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95):
    s, sbar, Pm, (V2, W2, V4, W4) = crystal(x)
    p1, p2, p3, p4 = Pm
    H2 = hel_vectors(p2); H4 = hel_vectors(-p4)
    E2 = {h: H2[h] for h in (1, -1)}; E4 = {h: np.conj(H4[h]) for h in (1, -1)}
    lab1 = duality(p1, sbar); lab2 = {h: duality(p2, E2[h]) for h in (1, -1)}
    hsame = -1 if lab2[1] == lab1 else +1           # CLS-D conjugate pairing (see g1)
    T = lambda A2, A4: {(a, b): dot(A2, E2[a])*dot(A4, E4[b]) for a in (1, -1) for b in (1, -1)}
    VV, WW, WV, VW = T(V2, V4), T(W2, W4), T(W2, V4), T(V2, W4)
    P = VV; H = {k: VV[k] + WW[k] for k in VV}; Q = {k: H[k] - 2*P[k] for k in VV}; O = {k: WV[k] - VW[k] for k in VV}; O1 = WV
    pres = [(hsame, hsame), (-hsame, -hsame)]; flip = [(1, -1), (-1, 1)]
    r = lambda A: max(abs(A[k]) for k in flip)/max(abs(A[k]) for k in pres)
    pr(f"x={x}: H_flip/P_flip = {H[flip[0]]/P[flip[0]]:.12f}, {H[flip[1]]/P[flip[1]]:.12f}; |H_pres| = {max(abs(H[k]) for k in pres):.1e}; "
       f"P_flip^2/(P_same P_opp) = {P[flip[0]]**2/(P[pres[0]]*P[pres[1]]):.12f}")
    Q2 = {k: Q[k]**2 for k in Q}
    pr(f"     photon Q = H-2P: flip/pres = {r(Q):.1e}; graviton Q^2: flip/pres = {r(Q2):.1e}, opp/same = {abs(Q2[pres[1]]/Q2[pres[0]]):.8f}"
       f" ((1-x^2)^-4 = {(1-x*x)**-4:.8f})")
    pr(f"     photon O1 = WV: same {O1[pres[0]]:+.4e} opp {O1[pres[1]]:+.4e} flips {O1[flip[0]]:+.2e} {O1[flip[1]]:+.2e}")
    pr(f"     photon O = WV-VW: same {O[pres[0]]:+.4e} opp {O[pres[1]]:+.4e} flips {O[flip[0]]:+.1e} {O[flip[1]]:+.1e};  "
       f"opp/same = {O[pres[1]]/O[pres[0]]:+.8f} (-(1-x^2)^-2 = {-(1-x*x)**-2:+.8f})")
    OQ = {k: O[k]*Q[k] for k in O}
    pr(f"     graviton O*Q: flip/pres = {r(OQ):.1e}; opp/same = {OQ[pres[1]]/OQ[pres[0]]:+.8f} (-(1-x^2)^-4 = {-(1-x*x)**-4:+.8f})")
open(os.path.join(HERE, 'cache', 'g2_structures.txt'), 'w').write("\n".join(out) + "\n")
