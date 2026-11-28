"""Light single traces of spin 0 and 2 (Letter assumption (i); O(1) amplitudes; b=0 ladder only, k<=J) in Theorem N.
J=0: F=1, constraint W0(Delta)+W_B>0 trivially (W0>=3/2-1/12 for Delta>1).
J=2, any Delta>3: F = 1 + c1*6/T + c2*24/(T(T-2)), c1 = M(q(1)-1), c2 = C(M,2)(q(2)-2q(1)+1).
Sufficient: P(T)=T(T-2)+6c1(T-2)+24c2 > 0 for T>2 (then F>0 and W0+W_B F>0). Exact check for every M in each certificate's range."""
import json, glob, math
from fractions import Fraction as Fr
for f in sorted(glob.glob('cache/certI_*.json')):
    d = json.load(open(f)); one = Fr(1)
    E = d['E']; shift = int(d['shift'])
    # recover q(1), q(2) exactly from S and stored data: recompute Q(1),Q(2)
    from cert_lemmaM import gram_rows
    Mb, n = d['Mb'], d['n']; A = [int(a) for a in d['A']]; es = d['es']
    Q = [0, 0, 0]
    for j, t in gram_rows(Mb, n):
        fj = A[j] << (E - es[j])
        for z in range(3): Q[z] += fj*t[z]
    Q = [v + shift for v in Q]; den = 1 << (60+E)
    q1, q2 = Fr(Q[1], den), Fr(Q[2], den)
    worst = None; okall = True
    for M in range(n, Mb+1):
        c1 = M*(q1-1); c2 = Fr(M*(M-1), 2)*(q2-2*q1+1)
        # P(T) = T^2 + (6c1-2)T + (24c2-12c1); positive for all real T iff disc<0; else check roots <= 2
        B = 6*c1-2; Cc = 24*c2-12*c1; disc = B*B-4*Cc
        if disc < 0: ok = True
        else:
            r = (-B + math.sqrt(float(disc)))/2   # larger root
            ok = r < 2 - 1e-9 and (4 + 2*B + Cc) > 0
        okall &= ok
        ratio = float(disc/(B*B))
        if worst is None or ratio > worst[0]: worst = (ratio, M)
    print(f, 'q1=%.6f q2=%.6f' % (float(q1), float(q2)), 'all M in [n,Mb] pass:', okall, ' max disc/B^2 =', worst)
