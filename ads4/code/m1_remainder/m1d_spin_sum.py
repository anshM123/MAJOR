r"""
m1_remainder/m1d_spin_sum.py  [P for the identity given (fp-H) at light states; NE for the numbers]

The FULL spin sum (all J, no truncation at J = 1000) of the twist-2 light remainder with physical, Bose-symmetric (even-J)
weights, from the generating functions H_b(z) = sum_J hhat_b(J) z^J of m1b (out/m1b_H_*.json):

   W(J) = C (J+1)^{-p} = (C/Gamma(p)) int_0^oo t^{p-1} e^{-t} e^{-J t} dt          (completely monotone)
   sum_{J even} W(J) hhat_b(J) = (C/(2 Gamma(p))) int_0^oo t^{p-1} e^{-t} [ H_b(e^{-t}) + H_b(-e^{-t}) ] dt
                                  - (model terms at J = 0, 2 for b = 0: negligible, |hhat| < 10^5 there).

The "+" term samples H_b on (0,1) (heavy-like, modest), the "-" term samples it on (-1,0) (back-to-back), where it is
~ e^{pi Lam/2}.  C and p are fitted to the a2prime weights c_T|a|^2 of each family (m1c.ctA2) on J in [400, 1600].
Output: out/m1d_spin_sum.txt
"""
import os, sys, json, math, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
mp.mp.dps = 40
OUT = os.path.join(HERE, 'out')


def load_H():
    pts = {}
    for f in glob.glob(os.path.join(OUT, 'm1b_H_*.json')):
        for z, v in json.load(open(f)).items():
            pts[mp.mpf(z)] = {k: mp.mpf(x) for k, x in v.items()}
    return dict(sorted(pts.items()))


def interp_logabs(ts, vals, t):
    """piecewise-linear interpolation of log|H| and of the sign in t (sign must be constant between nodes)."""
    for i in range(len(ts) - 1):
        if ts[i] <= t <= ts[i + 1]:
            a, b = vals[i], vals[i + 1]
            if mp.sign(a) != mp.sign(b):
                return None
            la, lb = mp.log(abs(a)), mp.log(abs(b))
            u = (t - ts[i])/(ts[i + 1] - ts[i])
            return mp.sign(a)*mp.exp(la + u*(lb - la))
    return None


def main():
    import m1c_physical_sum as C
    pts = load_H()
    lines = []
    def pr(s):
        print(s, flush=True); lines.append(s)
    pr('H_b(z) = sum_J hhat_b(J) z^J at the computed points (twist 2, Delta_g = 1000):')
    pr('   z                t=-ln|z|     H_0                      H_4')
    for z, v in pts.items():
        t = -mp.log(abs(z)) if z != 0 else mp.inf
        pr(f'   {mp.nstr(z, 12):>16s}  {mp.nstr(t, 4):>8s}  {mp.nstr(v["H0"], 8):>22s}  {mp.nstr(v["H4"], 8):>22s}')
    neg = [(-mp.log(-z), v) for z, v in pts.items() if z < 0]
    neg.sort(key=lambda x: x[0])
    pos = [(-mp.log(z), v) for z, v in pts.items() if z > 0]
    pos.sort(key=lambda x: x[0])
    # weights
    fam = {'pp': ('pp', 1, 0), 'pm': ('pm', 1, 4), 'N': ('N', 1, 4)}
    pr('\nfits W(J) = c_T|a|^2 = C (J+1)^{-p} on J in [400, 1600]:')
    fits = {}
    for key, (name, h0, b) in fam.items():
        J1, J2 = 400, 1600
        if name == 'N':
            w1, _ = C.ctA2_N(J1, mp.mpf(1)/J1); w2, _ = C.ctA2_N(J2, mp.mpf(1)/J2)
        else:
            w1, _ = C.ctA2(name, h0, J1); w2, _ = C.ctA2(name, h0, J2)
        p = mp.log(w1/w2)/mp.log(mp.mpf(J2 + 1)/(J1 + 1))
        Cc = w1*(J1 + 1)**p
        fits[key] = (Cc, p, b)
        pr(f'   {key:3s}: C = {mp.nstr(Cc, 6)}, p = {mp.nstr(p, 6)}')
    pref = 64*mp.mpf(998*999)**2/9
    total = mp.mpf(0)
    pr('\nfull even-J spin sums  sum_J W(J) hhat_b(J)  =  (+ part)  +  (- part, back-to-back):')
    for key, (Cc, p, b) in fits.items():
        Hk = f'H{b}'
        tn = [x[0] for x in neg]; vn = [x[1][Hk] for x in neg]
        f = lambda t: (t**(p - 1)*mp.exp(-t)*(interp_logabs(tn, vn, t) or 0))
        # integrate on the node intervals (log-linear integrand), then the tail beyond the last node is negligible
        Im = mp.mpf(0)
        for i in range(len(tn) - 1):
            Im += mp.quad(f, [tn[i], tn[i + 1]])
        tp = [x[0] for x in pos]; vp = [x[1][Hk] for x in pos]
        fp_ = lambda t: (t**(p - 1)*mp.exp(-t)*(interp_logabs(tp, vp, t) or 0))
        Ip = mp.mpf(0)
        for i in range(len(tp) - 1):
            Ip += mp.quad(fp_, [tp[i], tp[i + 1]])
        Sm = Cc/(2*mp.gamma(p))*Im; Sp = Cc/(2*mp.gamma(p))*Ip
        total += Sm + Sp
        pr(f'   {key:3s} (ladder b={b}): back-to-back part {mp.nstr(Sm, 6)}  (log10 {float(mp.log10(abs(Sm))):.2f}),'
           f'  heavy-like part (t in [{mp.nstr(tp[0],3)},{mp.nstr(tp[-1],3)}]) {mp.nstr(Sp, 4)}')
    R2 = pref*total
    pr(f'\nR_2 (twist-2 families pp, pm, N; all J; A = kappa = 1) = {mp.nstr(R2, 6)}   log10|R_2| = {float(mp.log10(abs(R2))):.2f}')
    L_A = mp.mpf('14339.59')/16; n0 = mp.mpf('4.53987e-7')
    cstar = mp.pi*abs(R2)/(512*n0*L_A)
    pr(f'threshold c_T* (X = 1, n_0 = {mp.nstr(n0, 4)}): {mp.nstr(cstar, 4)}  log10 = {float(mp.log10(cstar)):.2f}')
    pr(f'sign of R_2: {"negative (weakens the bound)" if R2 < 0 else "positive: the sum rule then needs c_T >= c_T* (Lemma 2 with Sigma_heavy >= 0)"}')
    open(os.path.join(OUT, 'm1d_spin_sum.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
