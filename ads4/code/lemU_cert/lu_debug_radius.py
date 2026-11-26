"""lemU_cert/lu_debug_radius.py [debug]: locate the ball-radius blow-up in the contraction (prints per-m the largest |R| rad(Y) and rad(R)|Y| terms).
Usage: lu_debug_radius.py L prec tgt   (tgt = GL error target; 1e5 for a faster, coarser run)"""
import sys, math, os
sys.path.insert(0, '/home/user/MAJOR/ads4/code/lemU_cert')
from flint import arb, ctx
import lu_core as C, lu_cert1 as C1
L = float(sys.argv[1]); prec = int(sys.argv[2]); tgt = float(sys.argv[3])
ctx.prec = prec
c = arb((1 - arb(-L).exp()).mid()); d = 1 - c
YA, YD, meta = C1.build(c, prec, print, tgt=tgt)
for pc in (prec, meta['pj']):
    ctx.prec = pc
    lr = [float(v.rad().log()) if v.rad() > 0 else -1e9 for v in YA]; lm = [float(abs(v.mid()).log()) if v.mid() != 0 else -1e9 for v in YA]
    print('prec', pc, 'YA: max ln rad %.1f at %d; max ln|mid| %.1f at %d' % (max(lr), lr.index(max(lr)), max(lm), lm.index(max(lm))))
    lrd = [float(v.rad().log()) if v.rad() > 0 else -1e9 for v in YD]; print('YD max ln rad %.1f at %d' % (max(lrd), lrd.index(max(lrd))))
    s = 1; xs = C.xstar(c, s); d2 = d*d; ic = 1/c
    YB = C.dress(YA, YD, c, 0)
    prev = [arb(1)]; worst = None
    for m in range(1, C1.M + 1):
        top = len(prev) - 1
        Bv = [prev[j + 1] - prev[j]*ic for j in range(top)] + [-prev[top]*ic]
        cur = [arb(0)]*(3*m + 1)
        for j in range(3*m):
            acc = d2*(j*(j + s))*Bv[j] if j <= top else arb(0)
            if 1 <= j <= top + 1: acc += d*(j*(2*j + s - 1))*Bv[j - 1]
            if 2 <= j <= top + 2: acc += (j*(j - 1))*Bv[j - 2]
            cur[j + 1] = xs*acc
        t1 = max(range(1, 3*m + 1), key=lambda j: float((abs(cur[j].mid())*YA[j].rad()).log()) if cur[j].mid() != 0 and YA[j].rad() > 0 else -1e9)
        t2 = max(range(1, 3*m + 1), key=lambda j: float((cur[j].rad()*abs(YA[j].mid())).log()) if cur[j].rad() > 0 and YA[j].mid() != 0 else -1e9)
        a1 = float((abs(cur[t1].mid())*YA[t1].rad()).log()); a2 = float((cur[t2].rad()*abs(YA[t2].mid())).log())
        if m in (1, 2, 5, 10, 30, 100, 300, 600, 1000) or max(a1, a2) > 0:
            print(f'  m={m}: max ln |R| rad(Y) = {a1:.1f} (j={t1}), max ln rad(R)|Y| = {a2:.1f} (j={t2}); ln|R_j| {float(abs(cur[t2].mid()).log()) if cur[t2].mid()!=0 else -1e9:.1f} ln rad R_j {float(cur[t2].rad().log()):.1f}', flush=True)
            if max(a1, a2) > 0:
                break
        prev = cur
