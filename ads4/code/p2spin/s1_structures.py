"""
p2spin/s1_structures.py (P2_SPINNING.md §1): conserved <T T O_{Delta,J}> structures in d = 3 over Q(Delta, J).

Unknowns: coefficients of the 14 even (21 odd) CPPR monomials of p2_core.  Constraints, each imposed as a function identity at random
rational configurations (exact Fractions):
  (C1) conservation of T at x1: d_{x1}^mu D_mu [N e^L] = 0, D the d = 3 Todorov operator (h = 3/2), L = J log V5 - alpha.log X
       (first x1-derivatives exact by dual numbers; the z1-dependence is carried as a symmetric matrix, so D is applied exactly);
       each row is A + Delta B + J C, with A, B, C rational;
  (S)  Bose symmetry f(1<->2) = f; for generic J the factor V5^J -> (-V5)^J is a sign sigma = (-1)^J (sigma = +1: even J, -1: odd J);
  (C2) conservation at x2 is NOT imposed; it is checked on the solution.
The d = 3 kinematic relations (null vectors of the monomial-evaluation matrix) solve everything trivially and are quotiented out.
Null space over the fraction field Q(Delta, J) by sympy DomainMatrix (exact).

Usage: python3 s1_structures.py  -> cache/s1_structures.txt, cache/s1_structures.pkl (exact solution bases, sympy srepr)
"""
import os, sys, pickle, random, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import p2_core as C

Dl, Jl = sp.symbols('Delta J')
h = Fr(3, 2)
OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def basis_vec(k):
    return tuple(Fr(1) if i == k else Fr(0) for i in range(3))


def cons_rows(sector, cfg, monos):
    """conservation-at-x1 row (A, B, C) for each monomial (Delta-coefficient B, J-coefficient C)."""
    x1, z1, x2, z2, x5, z5 = cfg
    x1d = tuple(C.Dual(x1[k], tuple(Fr(1) if i == k else Fr(0) for i in range(3))) for k in range(3))
    Bz = [C.blocks(x1d, basis_vec(k), x2, z2, x5, z5) for k in range(3)]     # linear blocks evaluated on basis vectors
    B0 = Bz[0]                                                               # z1-independent blocks are the same in all three
    dlV5, dl12, dl15 = B0['V5'].dlog(), B0['X12'].dlog(), B0['X15'].dlog()
    gA = [-5*dl12[r] for r in range(3)]
    if sector == 'odd':
        gA = [gA[r] - Fr(1, 2)*(dl12[r] + dl15[r]) for r in range(3)]
    gB = [Fr(1, 2)*(dl12[r] - dl15[r]) for r in range(3)]
    gC = [dlV5[r] + Fr(1, 2)*(dl12[r] - dl15[r]) for r in range(3)]
    zl = [sum(C.GL[m][n]*z1[n] for n in range(3)) for m in range(3)]
    rows = []
    for t in monos:
        lin, sc = C.mono_even_split(t) if sector == 'even' else C.mono_odd_split(t)
        S = sc(B0)
        la = [Bz[k][lin[0]] for k in range(3)]; lb = [Bz[k][lin[1]] for k in range(3)]
        Q = [[S*(la[k]*lb[l] + la[l]*lb[k])*Fr(1, 2) for l in range(3)] for k in range(3)]
        tr = sum((C.GI[a][b]*Q[a][b] for a in range(3) for b in range(3)), C.Dual(0))
        Dv = [2*h*sum((Q[m][n]*z1[n] for n in range(3)), C.Dual(0)) - zl[m]*tr for m in range(3)]
        A = sum(C.GI[m][r]*Dv[m].g[r] for m in range(3) for r in range(3)) + \
            sum(C.GI[m][r]*gA[r]*Dv[m].v for m in range(3) for r in range(3))
        Bc = sum(C.GI[m][r]*gB[r]*Dv[m].v for m in range(3) for r in range(3))
        Cc = sum(C.GI[m][r]*gC[r]*Dv[m].v for m in range(3) for r in range(3))
        rows.append((A, Bc, Cc))
    return rows


def mono_vals(sector, cfg, monos):
    B = C.blocks(*cfg)
    f = C.mono_even if sector == 'even' else C.mono_odd
    return [f(B, t) for t in monos]


def swap(cfg):
    x1, z1, x2, z2, x5, z5 = cfg
    return (x2, z2, x1, z1, x5, z5)


def rcfg(rng):
    while True:
        cfg = (C.rpoint(rng), C.rnull(rng), C.rpoint(rng), C.rnull(rng), C.rpoint(rng), C.rnull(rng))
        try:
            B = C.blocks(*cfg)
        except ZeroDivisionError:
            continue
        if all(B[k] != 0 for k in ('V1', 'V2', 'V5', 'X12', 'X15', 'X25', 'H12', 'H15', 'H25')):
            return cfg


def normalise(v):
    """scale a null vector so its last nonzero entry is 1; factor entries."""
    piv = [x for x in v if x != 0][-1]
    return [sp.factor(sp.cancel(x/piv)) for x in v]


def solve(sector, monos, sigma, Jval=None, nrows=None, seed=1):
    """exact null space of (C1) + (S) over Q(Delta, J) (or Q(Delta) if Jval given); also the relation space."""
    rng = random.Random(seed)
    n = len(monos)
    nrows = nrows or (2*n + 6)
    rows = []
    for _ in range(nrows):
        cfg = rcfg(rng)
        r = cons_rows(sector, cfg, monos)
        J = Jl if Jval is None else Jval
        rows.append([sp.Rational(a.numerator, a.denominator) + Dl*sp.Rational(b.numerator, b.denominator)
                     + J*sp.Rational(c.numerator, c.denominator) for (a, b, c) in r])
    for _ in range(n + 4):
        cfg = rcfg(rng)
        v0 = mono_vals(sector, cfg, monos); v1 = mono_vals(sector, swap(cfg), monos)
        rows.append([sp.Rational((sigma*b - a).numerator, (sigma*b - a).denominator) for a, b in zip(v0, v1)])
    syms = [Dl] if Jval is not None else [Dl, Jl]
    K = sp.QQ.frac_field(*syms)
    M = DomainMatrix.from_Matrix(sp.Matrix(rows)).convert_to(K)
    ns = M.nullspace()
    ns = [list(ns.to_Matrix().row(i)) for i in range(ns.shape[0])]
    # relations
    rel_rows = [[sp.Rational(x.numerator, x.denominator) for x in mono_vals(sector, rcfg(rng), monos)] for _ in range(n + 6)]
    R = DomainMatrix.from_Matrix(sp.Matrix(rel_rows)).convert_to(sp.QQ).nullspace()
    rel = [list(R.to_Matrix().row(i)) for i in range(R.shape[0])]
    return ns, rel, len(rows)


def quotient(ns, rel):
    """basis of span(ns) / span(rel): drop rel directions (rel subset of ns checked)."""
    K = sp.Matrix(ns + rel) if rel else sp.Matrix(ns)
    rk_all = sp.Matrix(ns).rank(simplify=True) if ns else 0
    rk_with = K.rank(simplify=True) if ns else len(rel)
    return rk_all, rk_with


def main():
    t0 = time.time()
    pr('# s1_structures: conserved <T T O_{Delta,J}> in d = 3, exact over Q(Delta, J)')
    res = {}
    for sector, monos in (('even', C.EVEN), ('odd', C.ODD)):
        for sigma, lab in ((1, 'J even'), (-1, 'J odd')):
            ns, rel, nr = solve(sector, monos, sigma)
            rk_ns, rk_union = quotient(ns, rel)
            phys = rk_ns - len(rel)
            pr(f'{sector:4s} sector, generic {lab}: {len(monos)} monomials, {nr} rows; null space dim {len(ns)}; '
               f'kinematic relations {len(rel)} (contained: {rk_union == rk_ns}); PHYSICAL conserved structures = {phys}')
            res[(sector, sigma)] = dict(ns=ns, rel=rel, phys=phys)
            pr(f'   [{time.time() - t0:.0f}s]')
    pickle.dump({k: dict(ns=[[sp.srepr(x) for x in v] for v in d['ns']], rel=[[sp.srepr(x) for x in v] for v in d['rel']],
                         phys=d['phys']) for k, d in res.items()},
                open(os.path.join(HERE, 'cache', 's1_structures.pkl'), 'wb'))
    with open(os.path.join(HERE, 'cache', 's1_structures.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()


# ================================================================ canonical bases and validations

def load_sol(sector, sigma):
    d = pickle.load(open(os.path.join(HERE, 'cache', 's1_structures.pkl'), 'rb'))[(sector, sigma)]
    ev = lambda s: sp.sympify(eval(s, sp.__dict__))
    return [[ev(x) for x in v] for v in d['ns']], [[ev(x) for x in v] for v in d['rel']]


def canonical(ns, rel, order):
    """basis of span(ns)/span(rel): reduce by the relations' RREF (relation pivots set to 0), then RREF in column order 'order'."""
    n = len(ns[0])
    K = sp.QQ.frac_field(Dl, Jl)
    vecs = [list(v) for v in ns]
    if rel:
        Rm = sp.Matrix(rel)
        # pivot the relations on the LAST columns of 'order' (monomials least wanted)
        perm = list(reversed(order))
        Rr, piv = Rm[:, perm].rref()
        pivcols = [perm[p] for p in piv]
        for i, pc in enumerate(pivcols):
            row = [0]*n
            for j, col in enumerate(perm):
                row[col] = Rr[i, j]
            for v in vecs:
                f = v[pc]
                if f != 0:
                    for k in range(n):
                        v[k] = sp.cancel(v[k] - f*row[k])
    M = DomainMatrix.from_Matrix(sp.Matrix(vecs)[:, order]).convert_to(K)
    Rr = M.rref()[0].to_Matrix()
    out = []
    for i in range(Rr.shape[0]):
        if all(x == 0 for x in Rr.row(i)):
            continue
        v = [0]*n
        for j, col in enumerate(order):
            v[col] = sp.factor(Rr[i, j])
        out.append(v)
    return out


def struct_value(sector, coeffs, cfg, Dv, Jv, sigma):
    """sum_mu c_mu N_mu at a configuration (N = monomial / V5^J, no X-prefactor), at numeric Delta, J."""
    monos = C.EVEN if sector == 'even' else C.ODD
    vals = mono_vals(sector, cfg, monos)
    return sum(sp.Rational(v.numerator, v.denominator)*c.subs({Dl: Dv, Jl: Jv}) for v, c in zip(vals, coeffs))
