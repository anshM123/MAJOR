"""
Exact MFT data of [J J]_{N,J} (two conserved currents) in the exact helicity basis of fourpt/t1b (hhat_opp, hhat_same),
two-particle Fock space.  Checks: polynomiality of the structures at double-twist Delta; Gram matrix diagonal in the helicity basis
(CHL: MFT is helicity-diagonal) [exact]; P_h = 2^J / G_hh.
Output: out_t1_mft_l1.txt, cache/mft_l1.pkl
"""
import sys, os, pickle, time
import sympy as sp
from mft2 import structure_inf, gram
from infframe import Delta, J, CACHE
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
m2 = (Delta - sp.Rational(3, 2))**2 - (J + sp.Rational(1, 2))**2
v1 = {'AC': (Delta**2 - 3*Delta - J**2 - J)*(Delta**2 - 3*Delta - J**2 - J + 4), 'AE': -(Delta + J)*(Delta - J - 1)*(Delta**2 - 3*Delta + J**2 + J),
      'BC': -(Delta + J)*(Delta - J - 1)*(Delta**2 - 3*Delta + J**2 + J), 'BE': (Delta**2 - 3*Delta + J**2 + J)**2, 'H': 0}
v2 = {'AC': 2, 'AE': 0, 'BC': 0, 'BE': 0, 'H': Delta**2 - 3*Delta + J**2 + J}
f_opp = J*(J + 1)*(m2 + 2); f_same = -(Delta - 1)*(Delta - 2)*(m2 - 2)
hopp = {w: v1[w] + f_opp*v2[w] for w in v1}; hsame = {w: v1[w] + f_same*v2[w] for w in v1}
if __name__ == "__main__":
    pr("# exact MFT data of [JJ]_{N,J} in the exact helicity basis (fourpt/tttt/t1_mft_l1.py); Delta = 4 + N, n = (2 + N - J)/2")
    res = {}
    for N in range(0, 7):
        for Jv in range(N + 2, -1, -2):
            t = time.time()
            try:
                Fo = structure_inf(1, hopp, N, Jv); Fs = structure_inf(1, hsame, N, Jv)
            except AssertionError as e:
                pr(f"  N={N} J={Jv}: not polynomial ({e})"); continue
            structs = []; names = []
            for nm, F in (('opp', Fo), ('same', Fs)):
                if F != 0: structs.append(F); names.append(nm)
            G = gram(1, structs)
            if len(names) == 2 and sp.simplify(G.det()) == 0:
                # degenerate (leading twist n = 0): the two helicity structures are proportional; one primary
                ratio = sp.nsimplify(G[0, 1]/G[1, 1])
                P = {'opp': 2**Jv/G[0, 0], 'same': 2**Jv/G[1, 1]}
                res[(N, Jv)] = (names, G)
                pr(f"  N={N} J={Jv} n={sp.Rational(2 + N - Jv, 2)}: hhat_opp = {ratio} hhat_same (degenerate, one primary);  P[hhat_opp] = {P['opp']}, P[hhat_same] = {P['same']}   ({time.time()-t:.0f} s)")
                continue
            P = {nm: (2**Jv/G[i, i] if G[i, i] != 0 else None) for i, nm in enumerate(names)}
            off = G[0, 1] if len(names) == 2 else 0
            res[(N, Jv)] = (names, G)
            pr(f"  N={N} J={Jv} n={sp.Rational(2 + N - Jv, 2)}: structures {names}; G_offdiag = {off};  P = {P}   ({time.time()-t:.0f} s)")
    pickle.dump(res, open(os.path.join(CACHE, 'mft_l1.pkl'), 'wb'))
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t1_mft_l1.txt'), 'w').write("\n".join(out) + "\n")
