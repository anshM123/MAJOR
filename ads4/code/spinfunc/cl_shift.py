"""cl_shift: max |lambda_a - lambda_b| over near-gap states (m/Lam <= 1.3) of spin J between two cl_chain outputs (same mode/Dg).
usage: python3 cl_shift.py mode Dg J tag_a tag_b   (prints the number; used by the cl_t*.sh scripts to decide on edge jackknife runs)"""
import os, sys, pickle
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
mode, Dg, J, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
L = lambda t: {(round(r['mfl'], 6), r['DF']): r['lam'] for r in pickle.load(open(os.path.join(C, 'cl_chain_%s_%s.pkl' % (mode, t)), 'rb'))[Dg]
               if r['J'] == J and r['mfl'] <= 1.3 + 1e-9}
A, B = L(a), L(b)
print('%.6f' % max(abs(A[k] - B[k]) for k in A if k in B))
