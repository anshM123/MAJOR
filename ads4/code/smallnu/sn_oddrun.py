"""
sn_oddrun (SMALLNU_WINDOW.md): run highspin/hs_odd.py (unchanged, read-only) on a row whose direct BULK row was computed here
(smallnu/cache), with every output redirected to smallnu/cache (hs_core.HSC patched before hs_odd's main runs; the sn_guard forbids
any other write).  Direct points: HIGHSPIN's (highspin/cache: GK15 edge nodes) merged with this note's (bulk + small nu), sn_core.
Usage: python3 sn_oddrun.py J [mfl,...]   -> smallnu/cache/hs_odd_J<J>.txt / .json
"""
import os, sys, runpy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sn_core as S

if __name__ == '__main__':
    S.H.HSC = S.SNC                      # hs_odd writes os.path.join(H.HSC, 'hs_odd_J<J>.*'); reading of points is via sn_core._pts
    sys.argv = ['hs_odd.py'] + sys.argv[1:]
    runpy.run_path(os.path.join(HERE, '..', 'highspin', 'hs_odd.py'), run_name='__main__')
