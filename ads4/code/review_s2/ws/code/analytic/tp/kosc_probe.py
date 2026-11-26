# probe TAIL_CERT Theorem P (tc_power.K_osc) at lower L0 (table rebuilt into analytic/cache, tailcert untouched)
import os, sys, math, json, shutil
L0 = float(sys.argv[1])
os.environ['TC_L0TAB'] = str(int(L0)) if L0 == int(L0) else str(L0)
sys.path.insert(0, '/home/user/MAJOR/ads4/code/tailcert')
import tc_power as P
# redirect the table cache to analytic/cache
P_tab = os.path.join('/home/user/MAJOR/ads4/code/analytic/cache', f'ehat_d_table_L{P.ETAB_L0:g}.json')
def Ehat_table():
    if P._ETAB is None:
        if os.path.exists(P_tab): P._ETAB = json.load(open(P_tab))
        else:
            P._ETAB = [P.Ehat_d(P.ETAB_U0 + P.ETAB_DU*i, 0.92, P.ETAB_DELTA) for i in range(P.ETAB_N)]
            json.dump(P._ETAB, open(P_tab, 'w'))
    return P._ETAB
P.Ehat_table = Ehat_table
for sg in (0.25, 1.25):
    for (t1, t2) in ((0.25, 1), (10, 11), (45, 46), (100, 101), (300, 301), (600, 601), (800, 801), (900, 901), (998, 999), (1050, 1051), (1100, 1101)):
        l = P.logK_osc(sg, t1, t2, L0)
        print(f'L0={L0} sigma={sg} tau in [{t1},{t2}]: log|S^| <= {l:.2f}', flush=True)
