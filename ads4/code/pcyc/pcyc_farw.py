"""
pcyc_farw (PCYC_CLOSE.md section 1.4): share of int_{D'} |Phi_scal| (nu-independent on D', A > 0) carried by the far part
max(s_x, s_y) >= r (raw light-cone sums; the 1/Delta series of B converges only for r < 1/3, RD_CLOSE section 2.3), with the D' product
rule of pcyc_dquad at nu = 0 (node mask: a rough estimate, the mask is not smooth; a few % accuracy).
Usage: python3 pcyc_farw.py -> cache/pcyc_farw.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pcyc_dlib as P
import pcyc_dquad as DQ


def run(Delta, ns=(16, 12, 6, 16, 12, 12, 12), rs=(0.2, 0.3, 0.5)):
    R = DQ.rules('S', list(ns), Delta)
    W = [np.real(r.weights([0.0])[:, 0]) for r in R]
    tx, rx, ts, ty, ry, ph, al = (r.x for r in R)
    G = [g.ravel() for g in np.meshgrid(tx, rx, ty, ry, indexing='ij')]
    w4 = np.einsum('a,b,c,d->abcd', W[0], W[1], W[3], W[4]).ravel()
    lw4 = R[0].logW(G[0]) + R[1].logW(G[1]) + R[3].logW(G[2]) + R[4].logW(G[3])
    tot = 0.0; far = np.zeros(len(rs))
    sx, _ = P.sig_map(G[0], G[1]*(1 - G[1]), Delta)          # sigma at s = 0 only for the mask (s_x = sigma)
    for it, tt in enumerate(ts):
        for ip, pp in enumerate(ph):
            for ia, aa in enumerate(al):
                lw = lw4 + R[2].logW(tt) + R[5].logW(pp) + R[6].logW(aa)
                for sg in (1.0, -1.0):
                    n = G[0].size
                    o = P.point(G[0], G[1], np.full(n, sg*np.sqrt(tt)), G[2], G[3], np.full(n, pp), np.full(n, aa), Delta)
                    h = np.exp(o['L'] - lw)*w4*W[2][it]*W[5][ip]*W[6][ia]
                    kx = G[1]*(1 - G[1])*(1 - tt); ky = G[3]*(1 - G[3])*np.cos(pp)**2
                    sx, _ = P.sig_map(G[0], kx, Delta); sy, _ = P.sig_map(G[2], ky, Delta)
                    m = np.maximum(sx, sy)
                    tot += np.sum(h)
                    for k, r in enumerate(rs):
                        far[k] += np.sum(h[m >= r])
    return far/tot


if __name__ == '__main__':
    L = ['# pcyc_farw: share of int_{D\'} |Phi_scal| with max(s_x, s_y) >= r (product rule at nu = 0, node mask)',
         '#  Delta   r = 0.2      0.3      0.5']
    for D in (20, 30, 50, 100, 300):
        f = run(float(D))
        L.append(f'  {D:5d}  ' + '  '.join(f'{v:.2e}' for v in f))
        print(L[-1], flush=True)
    open(os.path.join(HERE, 'cache', 'pcyc_farw.txt'), 'w').write('\n'.join(L) + '\n')
