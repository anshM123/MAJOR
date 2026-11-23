"""
ef_oddrun (EXT_FOLD.md §3): run highspin/hs_odd.py (unchanged, read-only) on an odd-J own row (Delta = 1000 + J) whose direct BULK row
came from the external job A, with every direct point of the state merged (ef_merge: HIGHSPIN GK15 / HC edge nodes, SMALLNU small-nu
n1 nodes, ext bulk row) and every output redirected to extfold/cache (ef_core sets hs_core.HSC; ef_guard forbids any other write).
This is exactly what smallnu/sn_oddrun.py did for J = 21.
Usage: python3 ef_oddrun.py J [mfl,...]   -> extfold/cache/hs_odd_J<J>.txt / .json
"""
import os, sys, runpy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_core as C

if __name__ == '__main__':
    C.SET[0] = os.environ.get('EF_SET', 'all')
    sys.argv = ['hs_odd.py'] + sys.argv[1:]
    runpy.run_path(os.path.join(HERE, '..', 'highspin', 'hs_odd.py'), run_name='__main__')
