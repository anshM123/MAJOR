"""
Analysis of the l = 2 heavy vertices (t12_vertex_l2 caches): with T X_h(Delta) = lambda_h X_h(Delta~) (lambda_h exact,
lambda_h = 2^{2J+1} P(Delta~)/(A_h P(Delta)), P = ((2-Delta+J)/2)_4 ((Delta+J-7)/2)_4  [from kappa_h A_h = 2^{3J+3/2} and C_2]),
the leading integrand is  T^{hL hR} = g_hL(Delta) lambda_hR g'_hR,  g = left vertex of X_h, g' = right vertex of X_h(Delta~).
Printed per ray and scale s (m ~ s m_0):  E_h := m^8 T^{hh}/lambda_h (Weyl-even part; -> 8 <=> 1/4 per helicity quadruple, see
t15), the same for the geometric mean of the mixed components, and the asymmetry T^{so}/T^{os} (leading order 1; the O(1/m) term of the
integrand is antisymmetric under the exchange of the vertices, as for l = 1).
Representative components: same (h1,h2,h3,h4) = (+,-,-,+), opp = (+,+,-,-) (preserved-preserved), mixed (+,-,+,+)/(+,+,+,-).
"""
import os, sys, pickle
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def P(D, J):
    return mp.rf((2 - D + J)/2, 4)*mp.rf((D + J - 7)/2, 4)

def A(h, D, J):
    return mp.rf(-J, 4)/mp.rf(J + 1, 4) if h == 'opp' else mp.rf(D - 1, 4)/mp.rf(2 - D, 4)

def lam(h, D, J):
    return 2**(2*J + 1)*P(3 - D, -1 - J)/(A(h, D, J)*P(D, J))

REP = {('same', 'same'): (1, -1, -1, 1), ('opp', 'opp'): (1, 1, -1, -1), ('same', 'opp'): (1, -1, 1, 1), ('opp', 'same'): (1, 1, 1, -1)}

if __name__ == "__main__":
    for ray in ('b0.4', 'b0.9'):
        fn = os.path.join(HERE, 'cache', 't12_res_%s.pkl' % ray)
        if not os.path.exists(fn): continue
        data = pickle.load(open(fn, 'rb'))
        pr(f"# ray {ray}")
        for s, (Dv, Jv, m, R) in sorted(data.items()):
            D = mp.mpf(Dv.p)/Dv.q; J = mp.mpf(Jv.p)/Jv.q
            T = {c: R[c][0][REP[c]] for c in REP}
            ls_, lo_ = lam('same', D, J), lam('opp', D, J)
            Es = T[('same', 'same')]*m**8/ls_; Eo = T[('opp', 'opp')]*m**8/lo_
            Eso = mp.sqrt(T[('same', 'opp')]*T[('opp', 'same')]/(ls_*lo_))*m**8
            asym = T[('same', 'opp')]/T[('opp', 'same')]
            pr(f"  s={s} Delta={mp.nstr(D, 8)} J={mp.nstr(J, 6)} m={mp.nstr(m, 7)}: m^8 T^ss/lambda_s = {mp.nstr(Es, 8)}, m^8 T^oo/lambda_o = {mp.nstr(Eo, 8)}, "
               f"m^8 sqrt(T^so T^os/(lambda_s lambda_o)) = {mp.nstr(Eso, 8)};  T^so/T^os = {mp.nstr(asym, 6)} (lambda_o/lambda_s = {mp.nstr(lo_/ls_, 6)})")
    open(os.path.join(HERE, 'out_t14_vertex_fit_l2.txt'), 'w').write("\n".join(out) + "\n")
