"""
ef_window (EXT_FOLD.md §3-§4): smallnu/sn_odd.py's window analysis (imported read-only, its main() run unchanged) on the merged
point sets of ef_core.
  * J = 23, 27, 31: extfold/cache/hs_odd_J<J>.json (ef_oddrun: direct bulk row) exists -> sn_odd uses it and the FULL normalisation
    (the grid with all direct points), as SMALLNU did for J = 21.
  * J = 29: HIGHSPIN's hs_odd_J29_final.json and the baseline normalisation (highspin/cache points), as SMALLNU did; with set 'all'
    the n3 nodes at nu = 10 - 17.5 are the POOLED records (31 distinct RQMC streams, ef_merge).
  * set 'old' with an output subdir: reproduction of SMALLNU's own numbers (no ext point used).
Usage: python3 ef_window.py J [set] [outsub]   -> extfold/cache[/outsub]/sn_odd_J<J>.{txt,json}
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_core as C

if __name__ == '__main__':
    J = int(sys.argv[1])
    C.SET[0] = sys.argv[2] if len(sys.argv) > 2 else 'all'
    if len(sys.argv) > 3:
        d = os.path.join(C.EFC, sys.argv[3]); os.makedirs(d, exist_ok=True)
        C.S.SNC = d
    C.H.HSC = C.HS_CACHE                       # sn_odd's fallback: HIGHSPIN's hs_odd_J<J>[_final].json (read-only)
    import sn_odd as SO                        # read-only import
    SO.main(J)
