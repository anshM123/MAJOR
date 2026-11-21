"""Result Q, Stage 4: write the de-shift weights of the forward k = 4 pieces (EE and scalar skeleton; min-norm gauge, dnu = 0.02,
N = 2048) and the exact R_j to functionals_q4_deshift.json (kept, not cache).  Usage: python3 stage4_save.py"""
import os, json, pickle
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
import stage4_hj as HJ
out = {'doc': 'Stage-4 de-shifted forward C4: functional = int_0^oo chi(nu) Psi_{4,nu} dnu, Psi normalised so that its leading-Regge heavy '
              'image is sum_j R_j(nu) P_{-1/2+i(nu - i j)}(cosh beta) (R_j below; h_j = N_com R_j, N_com absorbed in chi). '
              'chi_value: image P_{-1/2}(cosh beta) (= C^_4(0)); chi_deriv: image d_{nu^2}P|_0 (= d_{nu^2} C^_4(0)). Midpoint grid nu_k = (k+1/2) dnu. '
              'min-norm gauge of the chi-normalised Abel operator (3 null directions). See RESULT_Q_STAGE4.md.'}
D = HJ.load()
for key in ('EE3', 'sc3'):
    P = pickle.load(open(os.path.join(HERE, 'cache', f'stage4_deshift_{key}.pkl'), 'rb'))
    out[key] = dict(R_j={str(j): str(v) for j, v in D[key].items()}, dnu=float(P['dnu']), N=len(P['nu']),
                    chi_value=[float(x) for x in P['chi_value']], chi_deriv=[float(x) for x in P['chi_deriv']])
json.dump(out, open(os.path.join(HERE, 'functionals_q4_deshift.json'), 'w'))
print('written', os.path.getsize(os.path.join(HERE, 'functionals_q4_deshift.json')))
