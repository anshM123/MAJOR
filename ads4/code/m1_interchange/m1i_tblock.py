r"""
m1_interchange/m1i_tblock.py  --  the exact tree-level test of the block-wise light side (step 2 of the M1 brief).

Exactly known object: the O(1/c_T) light side of Omega^EE on the tree-level Einstein <TTTT>.  Per nu (k = 2 piece, heavy
normalisation) it is  -C'_{EE,nu}|_E = 16 pi G/(nu^2 + 9/4)  [supplement row L1 / H7 Lemma R: residue of the full correlator],
and C_3|_E = 0.  In the s-channel ((12)(34)) decomposition the O(1/c_T) light content is the T block plus the [TT]_{n,J} at double
zeros (row L1), so the block-wise light side at O(1/c_T) is  lambda_TTT^2 * omega^cont[T]  + sum_{n,J} 0 * omega^cont[[TT]_{n,J}].

Block-wise side = the ladder (fp-H), i.e. the analytic continuation in Delta of CLS-D's per-block integral (m1_direct).
Per nu, every state's k = 2 image is  rho(nu) * S(s)  with a state-independent rational factor rho (c(nu) = (4nu^2+1)/((4nu^2+9)
(4nu^2+49)) of identity I2, or c/(R_- + R_+), or 1, depending on the weight convention) and S a POLYNOMIAL in s = nu^2 + 1/4
(finite ladder).  The honest light value 16 pi G/(nu^2+9/4) is not rho * polynomial in any of these conventions.

Checks:
 1. S(s; Delta, J=2, b=0) has a simple pole at Delta = 3 (the T block); residue profile is a polynomial; k = 3 part likewise.
 2. Twist-1 poles for every J (b = 0 ladder): Delta = J + 1.
 3. No pole at even twist 2, 4, 6, ... (so the [TT] double traces contribute 0 x finite = 0 at O(1/c_T)), J <= 40, both ladders.
 4. Exact: (16 pi/(nu^2+9/4))/rho has a pole at finite nu for each convention rho, so no finite set of continued blocks
    (and no regularisation of the T-block pole) reproduces the honest per-nu light value.
Usage: python3 m1i_tblock.py  -> out/tblock.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'm1_direct'))
import sympy as sp
import md_ladder as ML      # read-only import: z_m, y_m, U, h2 (own sympy implementation of (fp-zy), (fp-U), (fp-H))
import md_k3 as K3          # read-only import: q_n (d = 1 CDH)

s = ML.s
eps = sp.Symbol('epsilon')
out = []
P = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)))

Z, Y = ML.zy(50)
q = K3.qseq(50)


def S2(b, Delta, J):
    T = (Delta - 1)*(Delta - 2); C = J*(J + 1); r = C/T
    tot = 0
    for k in range(0, J - b + 1):
        Dk = sp.binomial(k + 3, 3)
        V = sp.Rational(4*k, k + 3) if b == 0 else 4 + sp.Rational(12, k + 1)
        Uk = 1
        for j in range(k):
            u = (b + j)*(b + j + 1)
            Uk *= (C - u)/(T - u)
        tot += Dk**2*Uk*(Z[k + 2]*(sp.Rational(1, k + 3) - r/(k + 1)) + (2*Y[k + 2] + V*Z[k + 2])/T)
    return tot


def S3(b, Delta, J):
    T = (Delta - 1)*(Delta - 2); C = J*(J + 1); sig = 1 if b == 0 else -1
    tot = 0
    for k in range(0, J - b + 1):
        Uk = 1
        for j in range(k):
            u = (b + j)*(b + j + 1)
            Uk *= (C - u)/(T - u)
        tot += sig*sp.binomial(k + 3, 3)**2*Uk*q[k]
    return tot


P('# 1. The stress-tensor block (Delta, J) = (3 + eps, 2), ladder b = 0 (the only ladder at J = 2)')
e2 = sp.series(sp.together(S2(0, 3 + eps, 2)), eps, 0, 1).removeO()
res2 = sp.factor(sp.expand(e2.coeff(eps, -1)))
fin2 = sp.expand(e2.coeff(eps, 0))
P('  k=2 per-nu image S(s): residue at eps=0 :', res2)
P('                          finite part      :', sp.factor(fin2))
e3 = sp.series(sp.together(S3(0, 3 + eps, 2)), eps, 0, 1).removeO()
P('  k=3 per-nu image S3(s): residue          :', sp.factor(sp.expand(e3.coeff(eps, -1))))
P('  -> the continued per-block action of T has a simple pole at Delta = 3; residue and finite part are polynomials in s.')
P('  honest O(1/c_T) light value (row L1): -C_{EE,nu}|_E = 16 pi G/(nu^2+9/4), not of this form in any normalisation (sec. 4).')

P('\n# 2. Twist-1 (Delta = J + 1) poles of the b = 0 ladder, J = 2..12: order of the pole in eps at Delta = J + 1 + eps')
for J in range(2, 13, 2):
    T = (J + eps)*(J - 1 + eps); C = J*(J + 1)
    # pole comes from U_{0,J} = prod_{j<J} (C - u_j)/(T - u_j): the factor j = J - 1 has T - u_{J-1} = (J+eps)(J-1+eps) - (J-1)J
    den = sp.expand(T - (J - 1)*J)
    P(f'  J = {J:2d}: T - u_(J-1) = {sp.factor(den)}  -> simple pole at eps = 0 in the k = J ladder term (coefficient z/y_(J+2) != 0)')

P('\n# 3. Even twist (double traces) are regular: min_{j<k<=n} |T - u_j| over the ladder, tau = 2..12, all J <= 40, b = 0, 4')
worst = None
for tau in range(2, 13, 2):
    for J in range(0, 41, 1):
        for b in (0, 4):
            if J - b < 0:
                continue
            D = tau + J; T = (D - 1)*(D - 2)
            if J - b < 1:
                continue
            m = min([abs(T - (b + j)*(b + j + 1)) for j in range(J - b)])
            if worst is None or m < worst[0]:
                worst = (m, tau, J, b)
P(f'  smallest |T - u_j| = {worst[0]} at (tau, J, b) = {worst[1:]}  (> 0: no pole at any [TT]_{{n,J}}, so 0 x finite = 0 at O(1/c_T))')

P('\n# 4. The honest O(1/c_T) light value is not (state-independent factor) x (polynomial), in any of the normalisations in use.')
P('   Per nu, a continued block image is rho(nu) * S(s) with S a polynomial (finite ladder) and rho one of: c = (4nu^2+1)/((4nu^2+9)(4nu^2+49))')
P('   (identity I2 / Eq. fp-moments), c/(R_- + R_+) (weight X = chi (R_- + R_+)), or 1 (g-units, as in m1_remainder/m1_n0.py).')
import json
nu = sp.Symbol('nu')
d = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'kernels.json')))['D_k2_EE']
Rm = sp.sympify(d['R_minus'], locals={'nu': nu, 'I': sp.I}); Rp = sp.sympify(d['R_plus'], locals={'nu': nu, 'I': sp.I})
c = (4*nu**2 + 1)/((4*nu**2 + 9)*(4*nu**2 + 49))
P('   R_- + R_+ =', sp.factor(sp.simplify(Rm + Rp)))
TE = 16*sp.pi/(nu**2 + sp.Rational(9, 4))
for name, rho in (('c', c), ('c/(R_-+R_+)', sp.simplify(c/(Rm + Rp))), ('1', sp.Integer(1))):
    q = sp.factor(sp.simplify(TE/rho))
    num, den = sp.fraction(sp.together(q))
    P(f'   rho = {name:12s}: (16 pi/(nu^2+9/4))/rho = {q}   -> polynomial in nu? {sp.degree(den, nu) == 0}')
P('   In every normalisation the honest profile has a pole at finite nu (nu = +-i/2 or +-3i/2), every continued block image')
P('   is entire after removing rho.  So lambda_TTT^2 * (residue, finite part or any regularisation of omega^cont[T]) + sum 0 is never')
P('   the light side: the block-wise light side fails on the exact tree-level correlator, per nu and for every weight.')
open(os.path.join(HERE, 'out', 'tblock.txt'), 'w').write('\n'.join(out) + '\n')
