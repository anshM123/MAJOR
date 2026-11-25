"""
fullpos/fp_norm.py (FULLPOS_CERT.md §2): normalisation of the angular action against T0 in the fixed-impact-parameter limit.
Per-frequency identities (mpmath, 40 digits, series to k = 1500):
 (I2)  (1-r)^3 sum_k D_k^2 r^k z_{k+2}(nu) (1/(k+3) - r/(k+1))
         = -(4x+9)(4x+49)/(4608(4x+1)) [R_-(nu) P_{-1/2+i(nu+i)}(X) + R_+(nu) P_{-1/2+i(nu-i)}(X)],  X = (1+r)/(1-r), x = nu^2
 (I3)  sum_k D_k^2 r^k q_k(nu^2) = P_{-1/2+i nu}(X)/(1-r)^4
With g = -9216(4x+1)/((4x+9)(4x+49)) chi_re theta/Lam, (I2) gives
   lim_{T->oo, r fixed} H_ang = int g (LHS of I2) dnu = (2/Lam) int chi_re theta [R_- P(nu+i) + R_+ P(nu-i)] dnu
                              = (2/Lam) * (realised LR image of the k=2 weight) = 2 int_0^1 psi2(p) P_{-1/2+i Lam p}(cosh b) dp + O(eps_W),
which is T0's lim (B2 + cv2) (u -> 0, D_h -> J0).  (I3) gives the k = 3 analogue  sigma u int psi3 P.
"""
import json, os, sys
sys.dont_write_bytecode = True
import mpmath as mp
import fp_common as FC

mp.mp.dps = 40
kj = json.load(open(os.path.join(FC.DATA, 'kernels.json')))['D_k2_EE']
nu_s = mp.mpf(0)


def Rfun(expr, nu):
    return mp.mpmathify(eval(expr.replace('^', '**'), {'nu': nu, 'I': mp.mpc(0, 1), 'sqrt': mp.sqrt, 'pi': mp.pi}))


def series(nu, r, K=1500):
    s = nu*nu + mp.mpf(1)/4
    f = [mp.mpf(0), mp.mpf(1)]
    for m in range(1, K + 6):
        f.append((-(s - 2*m*m)*f[m] - (m - 1)**2*f[m - 1])/(m + 1)**2)
    q = [mp.mpf(1)]; qm = mp.mpf(0)
    for n in range(K + 1):
        prev = q[n - 1] if n >= 1 else qm
        q.append((-(s + 4 - 2*(n + 2)**2)*q[n] - n*n*prev)/(n + 4)**2)
    S2 = mp.mpf(0); S3 = mp.mpf(0)
    for k in range(K + 1):
        m = k + 2
        z = 12*f[m] - (m - 1)*f[m - 1] + (m + 1)*f[m + 1]
        D2 = mp.binomial(k + 3, 3)**2
        S2 += D2*r**k*z*(mp.mpf(1)/(k + 3) - r/(k + 1))
        S3 += D2*r**k*q[k]
    return (1 - r)**3*S2, S3


def main():
    print('keys', list(kj.keys()))
    Rm, Rp = kj['R_minus'], kj['R_plus']
    out = []
    for nu in ['0.3', '0.7', '2.0', '5.5']:
        for r in ['0.05', '0.3', '0.6']:
            nu_ = mp.mpf(nu); r_ = mp.mpf(r); X = (1 + r_)/(1 - r_); x = nu_**2
            lhs2, lhs3 = series(nu_, r_)
            Pm = mp.legenp(-mp.mpf(1)/2 + 1j*(nu_ + 1j), 0, X, type=3)
            Pp = mp.legenp(-mp.mpf(1)/2 + 1j*(nu_ - 1j), 0, X, type=3)
            rhs2 = -(4*x + 9)*(4*x + 49)/(4608*(4*x + 1))*(Rfun(Rm, nu_)*Pm + Rfun(Rp, nu_)*Pp)
            rhs3 = mp.legenp(-mp.mpf(1)/2 + 1j*nu_, 0, X, type=3)/(1 - r_)**4
            e2 = abs(lhs2 - rhs2)/abs(rhs2); e3 = abs(lhs3 - rhs3)/abs(rhs3)
            line = f'nu={nu:4s} r={r:5s}  I2: lhs={mp.nstr(lhs2, 15):>22s} rhs={mp.nstr(mp.re(rhs2), 15):>22s} (Im {mp.nstr(mp.im(rhs2), 3)}) rel {mp.nstr(e2, 3)}' \
                   f'   I3 rel {mp.nstr(e3, 3)}'
            print(line, flush=True); out.append(line)
    open(os.path.join(FC.CACHE, 'fp_norm.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
