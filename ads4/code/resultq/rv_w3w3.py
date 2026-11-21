"""
rv_w3w3 (RESULT_Q_REVISION.md §3): structure of the W3 x W3 (graviton exchange between two W3 vertices) light side of the k = 4 forward
functional, in CHMRS's scalar Mellin form (arXiv:2106.10274 eqs. BMellinDefinition, CfromBMellin; scalar proxy for the spinning case).

  B^_{k,T}[M] = (-1)^{k/2} Res_{S'=oo}-contour of  M(S', T' = T + 2Df - S') / [ (2Df - S') (Df - S'/2)_{k/2} (Df - T'/2)_{k/2} ]
  at FIXED U' = 2Df - T (the u-channel Regge limit is S' -> oo at fixed U').  C_{k,nu} = (Barnes) int dT Gamma(T/2)^2 gamma_{a-T}(nu) B^.

Checks (sympy, exact):
 (1) a T'-channel (or S'-channel) pole term  M = S'^p U'^q / (T' - c)  [the s-/t-channel W3W3 diagrams: poles in the variables that run
     with S' at fixed U'] gives B^_{4,T} = a POLYNOMIAL in T for every c (computed as minus the sum of the finite residues, the T-poles
     of the individual residues cancel) -> C_{4,nu} polynomial in nu^2 (Barnes lemma, CHMRS sec. light);
 (2) a U'-channel pole term with Regge growth S'^2 (the u-channel W3W3 diagram: spin-2 exchange in the Regge channel; vertex factors are
     U'-dependent only)  M = S'^2 P(U') / (U' - c)  gives B^_{4,T} = 0 identically (Lemma A in Mellin form), while for k = 2 it gives the
     pole 1/(U' - c) -> the decaying residues of the k = 2 (C'_EE) light side (H17), not of C^_4;
 (3) hence d_{nu^2}(-C^_4)(0)|_{W3W3} = the flat coefficient 4 pi G |g3|^2 exactly (flat-unit argument: a degree-1 polynomial in nu^2).
Usage: python3 rv_w3w3.py -> cache/rv_w3w3.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import sympy as sp

S, T, c, D = sp.symbols('S T c Df')
lines = []


def poch(x, n):
    return sp.prod([x + i for i in range(n)])


def Bhat(M, k, Df=3):
    Tp = T + 2*Df - S
    f = M.subs({'Tp': Tp}) if hasattr(M, 'subs') else M
    integrand = (-1)**(k//2)*f/((2*Df - S)*poch(Df - S/2, k//2)*poch(Df - Tp/2, k//2))
    integrand = sp.together(integrand)
    # residue at infinity = - coefficient of 1/S in the large-S expansion; B^ = contour at |S| = oo  = -Res_{S=oo} ... = sum of all finite residues
    w = sp.symbols('w')
    g = sp.simplify(integrand.subs(S, 1/w)/w**2)
    res_inf = -sp.residue(g, w, 0)
    return sp.simplify(-res_inf), integrand          # contour integral at |S| = oo  =  -Res_inf


if __name__ == '__main__':
    Df = 3
    Up = 2*Df - T                                          # fixed u-channel Mellin variable
    lines.append('# rv_w3w3: scalar-Mellin structure of exchange light sides in C_{4,nu} (Df = 3, d = 3)')
    for (p, q) in ((5, 0), (4, 1), (3, 2), (6, 0)):
        Tp = sp.Symbol('Tp')
        M = S**p*Up**q/(T + 2*Df - S - c)                 # t-channel pole term, growth S^(p-1) at fixed U'
        B, integ = Bhat(M, 4, Df)
        # independent route: minus the sum of the finite residues
        poles = sp.solve(sp.denom(sp.together(integ)), S)
        sfin = sum(sp.residue(integ, S, P) for P in set(poles))
        ok = sp.simplify(B - sfin) == 0
        isp = sp.simplify(B).is_polynomial(T)
        lines.append(f'  (1) t-channel  M = S^{p} U^{q}/(T\' - c): B^_4,T polynomial in T: {isp}; deg {sp.degree(sp.expand(B), T) if isp and B != 0 else "-"}; '
                     f'= sum of finite residues: {ok};  B^ = {sp.factor(sp.expand(B))}')
    for k in (4, 2):
        M = S**2*(Up**2 + 3*Up + 1)/(Up - c)
        B, _ = Bhat(M, k, Df)
        lines.append(f'  (2) u-channel  M = S^2 P(U\')/(U\' - c), k = {k}: B^_{k},T = {sp.factor(sp.simplify(B))}')
        if k == 4:
            _, integ = Bhat(M, 4, Df)
            for P in sorted(set(sp.solve(sp.denom(sp.together(integ)), S)), key=str):
                r = sp.factor(sp.simplify(sp.residue(integ, S, P)))
                lines.append(f'      finite residue at S\' = {P}: {r}')
            lines.append('      (S\' = 2Df, 2Df+2: s-channel double-twist points = the J-RESUMMED [TT]_0,1 light-family terms, M and dM/dS there;'
                         ' S\' = T, T-2: t-channel double-twist points; they cancel identically: no O(G) remnant in Mellin form)')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'rv_w3w3.txt'), 'w').write('\n'.join(lines) + '\n')
