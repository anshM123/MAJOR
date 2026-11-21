"""
rv_lower (RESULT_Q_REVISION.md §3.3): the lower-bound functional  L = c_L B4|0 - (dB4|0 - B5|0)  >= 0  on Stage 6's state set (HYB/EXA,
stage6_stress.states, Delta_g = 300, 1000) WITH the measured forward-derivative correction (BELOW_FOLD2/3):  dB4 -> dB4 + 2 u^3 s_d.
  c_L = max_states (dB4 - B5)/B4;  slope 1/(4 c_L) in  g4 Lam^6/(8 pi G) >= slope_W x |g3|^2 Lam^8/(4 c_L)  (slope_W = 1 if the W3W3 light
  side of d_{nu^2}C^_4(0) is the flat 4 pi G |g3|^2, i.e. (P) for the s-/t-channel W3W3 exchanges; RESULT_Q_REVISION §3).
Variants: K = 0 (Stage 3/6); measured s_d on the measured channels only (++ J = 0, 2, 4; +- J = 4); universal sigma_h s_0 on every state
with odd-J +- given the WORSE sign.  Also: implied |g3| closure  |g3|^2 Lam^8 <= 2 L_AdS / (kappa + slope)  for the delta = 0.2 functional.
Usage: python3 rv_lower.py -> cache/rv_lower.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

if __name__ == '__main__':
    import stage6_model as M, stage6_stress as ST
    from stage6_scan import XS, JCAP, BMAX
    name = M.FUNCS[0.2]; F, c2, c3 = M.load(name); a = F['a']
    flat = M.flat_table(name, XS, JCAP)
    lines = ['# rv_lower: lower-bound functional c_L B4 - (dB4 - B5) with measured s_d']
    SDM = {300: {('++', 0): .509, ('++', 2): .509, ('++', 4): .517, ('+-', 4): -.516},
           1000: {('++', 0): .508, ('++', 2): .508, ('++', 4): .510, ('+-', 4): -.510}}
    for Dg in (300, 1000):
        Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
        rows = ST.states(name, Dg, flat, Lam)
        R = []
        for (h, u, J, b2, b3, fw) in rows:
            _, _, _, P = M.hyb(F, h, u, float(J), b2, b3, tab, Lam, fw=fw)
            B4 = float(P['B4']); D = -float(P['FW'])/a          # = dB4 - B5
            R.append((h, u, J, B4, D))
        hs = np.array([r[0] for r in R]); U = np.array([r[1] for r in R]); Js = np.array([r[2] for r in R])
        B4 = np.array([r[3] for r in R]); D = np.array([r[4] for r in R])
        def rep(tag, sd):
            rat = (D + 2*U**3*sd)/B4; i = int(np.argmax(rat)); cL = rat[i]
            lines.append(f"  Dg={Dg} {tag}: c_L = {cL:.4f} at {hs[i]} J={Js[i]} m/Lam={1/np.sqrt(U[i]):.4f}; slope 1/(4 c_L) = {1/(4*cL):.4f}")
            return cL
        cl0 = rep('K=0', 0*U)
        sdm = np.array([SDM[Dg].get((h, int(J)), 0.0) for h, J in zip(hs, Js)])
        clm = rep('measured channels', sdm)
        for s0, g in ((0.509, 'H0'), (0.526, 'SB')):
            sig = np.where(hs == '++', 1.0, -1.0); odd = (hs == '+-') & (Js % 2 == 1)
            r1 = (D + 2*U**3*s0*sig)/B4; r2 = (D + 2*U**3*s0)/B4
            rat = np.where(odd, np.maximum(r1, r2), r1); i = int(np.argmax(rat))
            lines.append(f"  Dg={Dg} universal sigma_h s_0 = {s0} ({g}), odd-J +- worse sign: c_L = {rat[i]:.4f} at {hs[i]} J={Js[i]} "
                         f"m/Lam={1/np.sqrt(U[i]):.4f}; slope {1/(4*rat[i]):.4f}")
            clu = rat[i]
        # implied closure (delta = 0.2): 2 L_AdS = A log(2 Lam/3) + B_eff
        L2 = {300: 72.041, 1000: 92.814}[Dg]; kap = F['kappa']
        for cl, tag in ((cl0, 'K=0'), (clm, 'measured'), (clu, 'universal SB')):
            for sw in (1.0, 0.5, 0.0):
                sl = sw/(4*cl)
                lines.append(f"     closure [{tag}, slope_W = {sw}]: |g3|^2 Lam^8 <= 2L/(kappa + slope) = {L2/(kap + sl):.1f}")
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'rv_lower.txt'), 'w').write('\n'.join(lines) + '\n')
