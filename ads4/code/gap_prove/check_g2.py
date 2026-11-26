#!/usr/bin/env python3
"""Rigorous assembly check for Lemma G2 (seed/truncation allowance for 0.56 < c <= 0.9), designer addendum form.
Runs the author's full-band package scripts (copied, unmodified, into PKGDIR) and re-asserts each component with exact rationals:
  (a) old rational-pole seed: cutoff discrepancy 7.1e-24 plus numerical seed widths 1e-40, through the exact imaginary-pole
      action bounds of pole_bound_full_results.json (c<=0.9, h<=998^-2)                         -> < 1.7e-15 per channel
  (b) correlated log seed (raw_seed_budget_full.py, r=0.9/0.99 i.e. c<=0.9)                      -> < 1e-21
  (c) log theta cutoff (raw_seed_budget_full.py)                                                    -> < 1e-400 (asserted there)
  (d) new Gaussian M0 seed (verify_new_seeds.py)                                                    -> < 1e-45 (asserted there)
  (e) spectral-degree truncation M=6000 (spectral_cut.py)                                           -> < 2^-1910 (asserted there)
  total (a)+...+(e) < 1e-14.  (f) ordinary tail: tail(c_hi,K) with envelope 2^2337 (tight_spectral_tail.py) is charged per cell.
Usage: check_g2.py PKGDIR
"""
import sys,os,runpy,json
from fractions import Fraction as F
pk=os.path.abspath(sys.argv[1]);os.chdir(pk);sys.path.insert(0,pk)
g=runpy.run_path('tight_spectral_tail.py');assert g['E']==2337
runpy.run_path('pole_bound_full.py')
rr=json.load(open('pole_bound_full_results.json'))
coeff={1:F(-96),9:F(5184,5),49:F(-4704,5)}
pole={}
for b in (0,4):
 pole[b]=sum(abs(coeff[r['a']])*F(r['action_bound'])*(F('7.1e-24')/r['a']+F('1e-40'))/998 for r in rr if r['b']==b)
 assert pole[b]<F('1.7e-15')
g2=runpy.run_path('raw_seed_budget_full.py')
logs={}
for b in (0,4):
 logs[b]=g2['E']+sum(abs(g2['coef'][x['a']])*F(x['action_bound'])*F('1e-40')*(1+F(3*32*(2+240),x['a']))/998 for x in g2['rr'] if x['b']==b and x['a'] in g2['coef'])
 assert logs[b]<F('1e-21')
assert g2['cut']<F('1e-400')
runpy.run_path('verify_new_seeds.py')   # asserts < 1e-45 for both b
g3=runpy.run_path('spectral_cut.py');assert g3['E']<F(1,2**1910) and g3['c']==F(9,10)
for b in (0,4):
 tot=pole[b]+logs[b]+F('1e-400')+F('1e-45')+F(1,2**1910)
 assert tot<F('1e-14');print('G2 b',b,'pole %.4e log %.3e total(non-tail) %.4e < 1e-14'%(float(pole[b]),float(logs[b]),float(tot)))
tail=g['tail']
for c,K in [('.56',3200),('.6',3316),('.7',4700),('.777',6544),('.9',17199)]:
 print('example tail(',c,',',K,') = %.3e'%float(tail(F(c),K)))
print('G2 PASS: allowance(c_hi<=0.9) = tail(c_hi,K) + 1e-14 (tail charged per cell, exact)')
