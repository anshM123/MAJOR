r"""
c1_flipcheck (C1_COMPUTE.md §3.2): does the EE functional annihilate every helicity-flip plane-wave configuration at leading Regge?

The O(t4) E x W3 light term, in the Lemma-R / D3 prescription of H7, is a superposition of plane-wave pairs (p, pbar) whose
polarisation content at O(t4) is "flip on one line, preserved on the other" (the d = 3 one-point matrix H = 1 + (t4/4) flip; Q8).
This script evaluates the exact leading-Regge plane-wave image L(a1, a3; a2, a4)(nu, beta) of the de-shifted k = 2 functional with the
EE dual structure (s1_offdiag_kernel.img_vec with A_monos_struct('EE'), imported read-only) for ALL 16 polarisation assignments
(a_i = +-1; line (1,3) preserved iff a3 = -a1, line (2,4) preserved iff a4 = -a2), and the same for C^_3 (o_H q, j' = 1) for reference.
It prints |L(config)| / |L(PRES)| per configuration class.  q^2 is printed as a control (its (2,4)-flip is non-zero, OFFDIAG §1.2).

Usage: OMP_NUM_THREADS=1 python3 c1_flipcheck.py  -> cache/c1_flipcheck.txt     (~1-3 min, 1 core)
"""
import os, sys, io, time, itertools
os.environ['OMP_NUM_THREADS'] = '1'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF); sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
import s1_offdiag_kernel as S

OUT = os.path.join(HERE, 'cache', 'c1_flipcheck.txt')


def cls(p):
    a1, a3, a2, a4 = p
    l13 = 'pres' if a3 == -a1 else 'flip'
    l24 = 'pres' if a4 == -a2 else 'flip'
    return f"(1,3) {l13} x (2,4) {l24}"


def main():
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    pols = list(itertools.product((1, -1), repeat=4))
    nus = [0.3, 1.5, 4.0, 12.0, 40.0]
    betas = np.array([0.05, 0.2, 0.5, 1.0, 1.8])
    pr("# c1_flipcheck: max over beta in %s of |L(config)|/|L(PRES)|, per nu" % list(betas))
    for name, kw in (('EE', dict(kind='C2', Amonos=S.A_monos_struct('EE'))),
                     ('q2 (control)', dict(kind='C2', Amonos=S.A_monos_struct('q2'))),
                     ('C3 o_H q', dict(kind='C3'))):
        pr(f"\n## structure {name}")
        for nu in nus:
            nphi = 512 if nu*betas.max() > 20 else 256
            ref = np.abs(S.img_vec(nu, betas, S.PRES, nphi=nphi, **kw))
            res = {}
            for p in pols:
                v = np.abs(S.img_vec(nu, betas, p, nphi=nphi, **kw))/ref
                c = cls(p)
                res[c] = max(res.get(c, 0.0), float(v.max()))
            pr(f"  nu = {nu:6.2f}: " + "; ".join(f"{c}: {res[c]:.2e}" for c in sorted(res)))
    pr(f"\n# done in {time.time() - t0:.0f} s")
    with open(OUT, 'w') as f:
        f.write(buf.getvalue())


if __name__ == '__main__':
    main()
