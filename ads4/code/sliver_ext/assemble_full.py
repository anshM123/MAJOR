#!/usr/bin/env python3
"""Exact end-to-end sign comparisons, conditional on retained inputs."""
import json
from fractions import Fraction as F
from pathlib import Path
from check_lowfreq import exp_u
from check_norm import parameters
from base_bounds import hi

def check_d(rows,Llo,Lhi):
    d,dd,*_=parameters(rows[0]['config'])
    assert d-dd <= 1/exp_u(F(Lhi)) and d+dd >= exp_u(-F(Llo))
    assert all(parameters(r['config'])[:2]==[d,dd] for r in rows)
    return d-dd,d+dd

def main():
    res=[]
    for L,H,name,stem,floor in [('2.9','3.2','resid_combo_2.9_3.2','2.9_3.2','0.0000025'),('3.2','3.35','resid_3.2_3.35','3.2_3.35','0.00000187'),('3.35','3.5','resid_3.35_3.5','3.35_3.5','0.00000187')]:
        dl,du=check_d(json.load(open(name+'_results.json')),L,H)
        for r in json.load(open(name+'_results_budget.json')):
            b=r['b'];assert F(r['dlo'])==dl and F(r['dhi'])==du
            contour=F(0)
            for p in (3,4):
                c=json.load(open(f'ct_{stem}_{p}_{b}.json'))[0]
                assert c['p']==p and c['b']==b and F(c['dlo'])==dl and F(c['dhi'])==du
                contour+=F(c['total'])
            E=F(r['total_exact'])+contour;floorq=F(floor)
            row={'L':[L,H],'b':b,'total_exact':str(E),'floor_exact':str(floorq),'contour_approx':float(contour),'error_ratio_upper':float(hi(E/floorq)),'lower_approx':float(floorq-E),'status':'PASS_CONDITIONAL' if E<floorq else 'FAIL'}
            print(json.dumps(row));res.append(row)
    check_d(json.load(open('resid_3p5_3p65_results.json')),'3.5','3.65')
    for r in json.load(open('assembly35_results.json')):
        E=F(r['total_exact']);floor=F('1.98e-6');assert E<floor
        row={'L':['3.5','3.65'],'b':r['b'],'total_exact':str(E),'floor_exact':str(floor),'error_ratio_upper':float(hi(E/floor)),'lower_approx':float(floor-E),'status':'PASS_CONDITIONAL'}
        res.append(row);print(json.dumps(row))
    Path('full_ledger.json').write_text(json.dumps(res,indent=2)+'\n')
if __name__=='__main__':main()
