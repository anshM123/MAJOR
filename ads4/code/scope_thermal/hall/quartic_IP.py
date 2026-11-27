"""
Flat 4-graviton contact amplitudes of I^2, P^2 and I P (I = C.C, P = C.C~), same conventions as
../quartic_flat.py (CLPS-normalised Einstein vertex).  Checks that CHLPS's
    (1/4)(alpha4 I^2 + alpha4' P^2 + 2 alpha4~ I P)  ->  g4 = 8piG(alpha4 + alpha4'),  ghat4 = 8piG(alpha4 - alpha4' + i alpha4~)
i.e. that alpha4~ enters the all-plus amplitude with unit weight relative to alpha4 - alpha4'.  [computed]
"""
import itertools, importlib.util, os, sys, io, contextlib
import numpy as np

spec = importlib.util.spec_from_file_location("qf", os.path.join(os.path.dirname(__file__), "..", "quartic_flat.py"))
qf = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(qf)

etaL = qf.etaL
epsmix = np.einsum('abgh,ge,hf->abef', qf.LCu, etaL, etaL)


def sums(Rs):
    Rl = [qf.lower_all(R) for R in Rs]
    out = np.zeros(3, dtype=complex)
    for perm in itertools.permutations(range(4)):
        A, B, C, D = [Rs[k] for k in perm]
        Bl, Dl = Rl[perm[1]], Rl[perm[3]]
        At = 0.5 * np.einsum('abef,efcd->abcd', epsmix, A)
        Ct = 0.5 * np.einsum('abef,efcd->abcd', epsmix, C)
        I_AB = np.einsum('abcd,abcd->', A, Bl); I_CD = np.einsum('abcd,abcd->', C, Dl)
        P_AB = np.einsum('abcd,abcd->', At, Bl); P_CD = np.einsum('abcd,abcd->', Ct, Dl)
        out += [I_AB * I_CD, P_AB * P_CD, I_AB * P_CD]
    return out


for trial in range(2):
    ls, lts = qf.four_kin()
    ps = [qf.vec(ls[i], lts[i]) for i in range(4)]
    rr, rt = qf.rs(), qf.rs()
    es = [qf.pol_plus(ls[i], lts[i], rr) for i in range(4)]
    sq = qf.sq
    ref = (sq(lts[0], lts[1]) ** 2 * sq(lts[2], lts[3]) ** 2 + sq(lts[0], lts[2]) ** 2 * sq(lts[1], lts[3]) ** 2
           + sq(lts[0], lts[3]) ** 2 * sq(lts[1], lts[2]) ** 2) ** 2
    s = sums([qf.riem_upper(ps[i], es[i]) for i in range(4)]) / ref
    print(f"trial {trial} ++++ : [I^2, P^2, I P] sym sums / (sum [ij]^2[kl]^2)^2 =", np.round(s, 10))
    print(f"   all-plus amplitude of (1/4)(a4 I^2 + a4' P^2 + 2 a4~ IP)  ~  {s[0]/4:.3f} a4 + {s[1]/4:.3f} a4' + {2*s[2]/4:.3f} a4~")
    es = [qf.pol_minus(ls[0], lts[0], rt), qf.pol_minus(ls[1], lts[1], rt),
          qf.pol_plus(ls[2], lts[2], rr), qf.pol_plus(ls[3], lts[3], rr)]
    ref = qf.ang(ls[0], ls[1]) ** 4 * qf.sq(lts[2], lts[3]) ** 4
    s = sums([qf.riem_upper(ps[i], es[i]) for i in range(4)]) / ref
    print(f"trial {trial} --++ : [I^2, P^2, I P] sym sums / <12>^4[34]^4 =", np.round(s, 10))

# normalisation of the all-plus structure against CHLPS's h-function: M(++++) = ([12][34]/(<12><34>))^2 h,
# h_low > (1/2) ghat4 (s^2+t^2+u^2)^2.  Ratio below = c in  (sum [ij]^2[kl]^2)^2 = c ([12][34]/<12><34>)^2 (s^2+t^2+u^2)^2
for trial in range(3):
    ls, lts = qf.four_kin()
    ang, sq = qf.ang, qf.sq
    s_ = ang(ls[0], ls[1]) * sq(lts[1], lts[0]); t_ = ang(ls[0], ls[3]) * sq(lts[3], lts[0]); u_ = ang(ls[0], ls[2]) * sq(lts[2], lts[0])
    # check s+t+u = 0 with these sign conventions (Mandelstams up to overall sign convention, squared below)
    S = (sq(lts[0], lts[1]) ** 2 * sq(lts[2], lts[3]) ** 2 + sq(lts[0], lts[2]) ** 2 * sq(lts[1], lts[3]) ** 2
         + sq(lts[0], lts[3]) ** 2 * sq(lts[1], lts[2]) ** 2) ** 2
    pref = (sq(lts[0], lts[1]) * sq(lts[2], lts[3]) / (ang(ls[0], ls[1]) * ang(ls[2], ls[3]))) ** 2
    print(f"c-check trial {trial}: s+t+u = {abs(s_+t_+u_):.1e};  c = {S / (pref * (s_**2 + t_**2 + u_**2) ** 2):.6f}")
