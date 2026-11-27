r"""
m1_remainder/m1a_hhat_allJ.py  [NE; ball arithmetic]

Per-state stripped light action hhat_b(Delta, J) of Omega^EE (Delta_g = 1000) on light double-twist families, for EVERY
integer J (both parities of J, although a parity-even theory populates only even J on both ladders, see closed_forms.py),
from the exact ladder formula (sec_fullpos.tex, Eqs. (fp-target), (fp-H)), with the moment balls of ../fullpos (read-only):

   hhat_b = T sum_{k<=J-b} D_k^2 U_{b,k} [ Z_k (1/(k+3) - r/(k+1)) + (2 Y_k + V_{b,k} Z_k + sigma_b Lam^2 Q_k)/T ].

Same formula and code path as ../no_eft/nl1_light_action.py (imported read-only); here every J is scanned and the
single-sequence quantities P_k = D_k^2 L_k (L_k = (Z_k - k Z_{k-1}/(k+3))/(k+3)) and their partial sums S_J are stored,
so that the parity structure in J and the generating-function identities of m1b can be checked.

Output: out/m1a_hhat_<twist>_b<b>.json : {J: [sign, log10|hhat|, log10|hhat/T|]},  out/m1a_Pk.json : P_k, S_k (k <= 1000).
Usage: python3 m1a_hhat_allJ.py   (about 3 min; Q_k exact setup dominates);  M1_DG=300 M1_K=300 python3 m1a_hhat_allJ.py
"""
import sys, os, json, time, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'fullpos'))
sys.path.insert(0, os.path.join(HERE, '..', 'no_eft'))
import fp_common as FC
import fp_poly as FP
from flint import arb, fmpq, ctx

OUT = os.path.join(HERE, 'out')
Dg = int(os.environ.get('M1_DG', '1000'))
K = int(os.environ.get('M1_K', '1000'))
SUF = '' if Dg == 1000 else f'_dg{Dg}'


def lg(x):
    if 0 in x:
        return float('nan')
    return float(abs(x).log()/arb(10).log())


def sgn(x):
    return 0 if 0 in x else (1 if x > 0 else -1)


def main():
    t0 = time.time()
    Z, Y = FP.load_moments(Dg, K)
    ctx.prec = 4000                      # plenty: light states have no cancellation (U ~ 1); radii are tracked
    lam = (Dg - 2)*(Dg - 1)
    Q = [arb(fmpq(q.p, q.q)) if hasattr(q, 'p') else arb(q) for q in FP.Q_exact(K, lam)]
    Lam2 = arb(lam)
    print(f'setup {time.time()-t0:.0f}s', flush=True)
    # single sequence
    from math import comb
    D = [comb(k + 3, 3) for k in range(K + 1)]
    L = [Z[0]/3] + [(Z[k] - arb(k)/(k + 3)*Z[k - 1])/(k + 3) for k in range(1, K + 1)]
    P = [D[k]**2*L[k] for k in range(K + 1)]
    S = []; acc = arb(0)
    for k in range(K + 1):
        acc += P[k]; S.append(acc)
    json.dump({'P': [[sgn(p), lg(p)] for p in P], 'S': [[sgn(s), lg(s)] for s in S],
               'Pstr': [p.mid().str(30, radius=False) for p in P]}, open(os.path.join(OUT, f'm1a_Pk{SUF}.json'), 'w'))
    for tw in (2, 3, 4, 6):
        for b in (0, 4):
            res = {}
            for J in range(max(b, 2), K + b + 1):
                n = J - b
                if n > K:
                    break
                Delta = J + tw
                T = arb((Delta - 1)*(Delta - 2)); C = J*(J + 1); r = arb(C)/T
                sig = 1 if b == 0 else -1
                U = arb(1); tot = arb(0)
                for k in range(n + 1):
                    V = arb(4*k)/(k + 3) if b == 0 else 4 + arb(12)/(k + 1)
                    tot += D[k]**2*U*(Z[k]*(1/arb(k + 3) - r/(k + 1)) + (2*Y[k] + V*Z[k] + sig*Lam2*Q[k])/T)
                    if k < n:
                        u = (b + k)*(b + k + 1)
                        U = U*(C - u)/(T - u)
                h = T*tot
                res[J] = [sgn(h), lg(h), lg(h/T)]
            json.dump(res, open(os.path.join(OUT, f'm1a_hhat_{tw}_b{b}{SUF}.json'), 'w'))
            print(f'twist {tw} b {b} done {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
