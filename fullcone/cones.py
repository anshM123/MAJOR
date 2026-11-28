"""Facet lists of the RT holographic entropy cones C_3, C_4, C_5 (orbit representatives -> all instances).

C_3: SA, MMI.  C_4: SA, MMI (Bao et al. 1505.07839).  C_5: SA, MMI, cyclic-5 and the four further
inequalities of Bao et al. 1505.07839 Sec. 4.3; completeness for N=5: Hernandez-Cuenca arXiv:1903.09148
(8 orbits, 372 facets; orbit sizes cross-checked against Avis-Hernandez-Cuenca 2102.07535).  A vector S is indexed by bitmasks 1..2^N-1 over parties 0..N-1; the purifier is party N.
An inequality is a dict {mask over N+1 parties: coeff} meaning sum coeff*S(mask) >= 0, reduced with
S(I) = S(complement of I in N+1 parties)."""
import itertools

def parse(s, letters="ABCDEO"):
    """'2ABC+ABD-AB' style -> {frozenset: coeff}"""
    out = {}
    for tok in s.replace("-", "+-").split("+"):
        tok = tok.strip()
        if not tok: continue
        sg = -1 if tok[0] == "-" else 1
        tok = tok.lstrip("-")
        k = 0
        while tok[k].isdigit(): k += 1
        c = int(tok[:k]) if k else 1
        A = frozenset(letters.index(ch) for ch in tok[k:])
        out[A] = out.get(A, 0) + sg*c
    return out

def lhs_rhs(l, r): 
    d = parse(l); 
    for A, c in parse(r).items(): d[A] = d.get(A, 0) - c
    return d

REPS = {
 3: {"SA": lhs_rhs("A+B", "AB"), "MMI": lhs_rhs("AB+BC+AC", "A+B+C+ABC")},
 4: {"SA": lhs_rhs("A+B", "AB"), "MMI": lhs_rhs("AB+BC+AC", "A+B+C+ABC")},
 5: {"SA": lhs_rhs("A+B", "AB"), "MMI": lhs_rhs("AB+BC+AC", "A+B+C+ABC"),
     "MMI2": lhs_rhs("AB+BCD+ACD", "A+B+CD+ABCD"),
     "cyc5": lhs_rhs("ABC+BCD+CDE+DEA+EAB", "AB+BC+CD+DE+EA+ABCDE"),
     "I2": lhs_rhs("2ABC+ABD+ABE+ACD+ADE+BCE+BDE", "AB+ABCD+ABCE+ABDE+AC+AD+BC+BE+DE"),
     "I3": lhs_rhs("ABE+ABC+ABD+ACD+ACE+ADE+BCE+BDE+CDE", "AB+ABCE+ABDE+AC+ACDE+AD+BCD+BE+CE+DE"),
     "I4": lhs_rhs("ABC+ABD+ABE+ACD+ACE+BC+DE", "AB+ABCD+ABCE+AC+ADE+B+C+D+E"),
     # I5 as transcribed from Bao et al. by pdftotext failed the contraction test; use the representative of
     # Avis-Hernandez-Cuenca 2102.07535 Table 2, no. 5 (digits 1..5 -> A..E):
     "I5": lhs_rhs("3ABC+3ABD+ABE+ACD+3ACE+ADE+BCD+BCE+BDE+CDE",
                   "2AB+2AC+AD+AE+BC+2BD+2CE+DE+2ABCD+2ABCE+ABDE+ACDE")},
}

def reduce(ineq, N):
    """map {frozenset over 0..N (N=purifier)} -> {mask over 0..N-1}"""
    full = frozenset(range(N+1)); out = {}
    for A, c in ineq.items():
        if N in A: A = full - A
        if not A or A == full: continue
        m = sum(1 << i for i in A)
        out[m] = out.get(m, 0) + c
    return tuple(sorted((m, c) for m, c in out.items() if c != 0))

def instances(N):
    """all distinct facet instances, as (name, tuple((mask,coeff)))"""
    res = {}
    for name, rep in REPS[N].items():
        used = set().union(*rep.keys())
        for perm in itertools.permutations(range(N+1)):
            img = {frozenset(perm[i] for i in A): c for A, c in rep.items()}
            key = reduce(img, N)
            if key not in res: res[key] = name
    return [(n, k) for k, n in res.items()]

def evaluate(ineq, S):
    return sum(c*S[m] for m, c in ineq)

def balance(ineq, N, i, j):
    """Delta for a Bell pair between parties i,j (j may be the purifier N): value of the inequality on b_ij"""
    return sum(c*(((m >> i) & 1) != (((m >> j) & 1) if j < N else 0)) for m, c in ineq)

if __name__ == "__main__":
    for N in (3, 4, 5):
        inst = instances(N)
        from collections import Counter
        print(N, len(inst), dict(Counter(n for n, _ in inst)))
