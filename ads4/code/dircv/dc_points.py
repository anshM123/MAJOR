"""
dc_points (DIRCV.md §4): NEW independent-seed per-wave + P_J + reference nodes (pw_points, unchanged, run in-process) written ONLY to
dircv/cache (pw_points.C redirected; dc_guard active).  Same production settings as HWIN queue B (PW_KINDS=EE,C3H, PW_ADAPT=1, M=9)
except the seed (PW_SEED) and the replica count (PW_R).
Usage: PW_TAG=_x1 PW_SEED=777 PW_R=16 python3 dc_points.py J Delta nu1,nu2
"""
import dc_guard, os, sys
os.environ.setdefault('PW_KINDS', 'EE,C3H'); os.environ.setdefault('PW_ADAPT', '1'); os.environ.setdefault('PW_M', '9')
os.environ['OMP_NUM_THREADS'] = '1'
import pw_points as PP
PP.C = dc_guard.CACHE
J = int(sys.argv[1]); Ds = float(sys.argv[2]); nus = [float(v) for v in sys.argv[3].split(',')]
print(f"# dc_points J={J} Delta={Ds:g} nus={nus} tag={PP.TAG} seed={PP.SEED} R={PP.R} M={PP.M} -> {PP.fn(J, Ds)}", flush=True)
PP.run(J, Ds, nus)
