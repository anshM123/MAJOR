"""n3fix/d4_logid.py: which source line calls each recorded log (n1_pw._rlog) in build_spin(kind), J even.  python3 d4_logid.py J Delta x kind"""
import nf_guard  # noqa: F401
import os, sys, traceback
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc')); sys.path.insert(0, SF)
os.environ['HC_CACHE'] = os.path.join(HERE, 'cache'); os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import hc_points as HC, z2_points as Z2
J = int(sys.argv[1]); Dl = float(sys.argv[2]); x = float(sys.argv[3]); kind = sys.argv[4]
HC.setup(odd=(J % 2 == 1))
import n1_pw as PW, n1_direct as D, n3_flow as F
nu = x*HC.mphys(Dl, J)
c, wc, L, _ = F.saddle(nu, Dl, J); Z2.set_qsc(Dl, J, c, wc, L)
ev = F.Evaluator(nu, Dl, J, c, wc, L, kind, nuf=nu)
orig = PW._rlog
calls = []
def rl(a):
    st = traceback.extract_stack()[-4:-1]
    calls.append(' <- '.join(f"{os.path.basename(f.filename)}:{f.lineno}" for f in st[::-1]) + '  | ' + st[-1].line[:90])
    return orig(a)
PW._rlog = rl
ev.run(ev.grav, np.zeros((1, 7), complex), None)
for i, s in enumerate(calls):
    print(i, s)
