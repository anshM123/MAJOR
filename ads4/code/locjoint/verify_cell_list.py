#!/usr/bin/env python3
"""Verify the emitted cell partition and recorded maxima; stdlib only."""
import csv,json,sys
from pathlib import Path
p=Path(sys.argv[1] if len(sys.argv)>1 else Path(__file__).parent/'output')
s=json.loads((p/'summary.json').read_text()); S=int(s['scale']); total=0
for run in s['runs']:
    vs=run['V']; slabs=run['slabs']; count=0; worst=[None]*len(vs); prev=None
    with (p/f"cells_{run['run']}.csv").open(newline='') as f:
        rows=iter(csv.DictReader(f))
        for sl in slabs:
            tl,th=int(sl['t_lo']),int(sl['t_hi']); assert prev is None or prev==tl; prev=th
            n=sl['beta_cells']; m=sl['angular_cells']; vcuts=[]; xcuts=[]; local=[None]*len(vs)
            for i in range(n):
                for j in range(m):
                    r=next(rows); assert (int(r['slab']),int(r['beta_index']),int(r['angle_index']))==(sl['slab'],i,j)
                    assert (int(r['t_lo']),int(r['t_hi']))==(tl,th)
                    v=(int(r['v_lo']),int(r['v_hi'])); x=(int(r['x_lo']),int(r['x_hi']))
                    if j==0: vcuts.append(v)
                    assert v==vcuts[i]
                    if i==0: xcuts.append(x)
                    assert x==xcuts[j]
                    assert int(r['log_beta_upper'])<=0 and int(r['log_angle_upper'])<=0
                    for h,V in enumerate(vs):
                        ft=int(r[f'log_FT_{V}']); z=r[f'log_tail_{V}']; best=ft if z=='NA' else min(ft,int(z))
                        w=int(r[f'log_weighted_{V}']); assert w==int(r['log_beta_upper'])+int(r['log_angle_upper'])+best
                        local[h]=w if local[h] is None else max(local[h],w)
                    count+=1
            for cuts in [vcuts,xcuts]:
                assert cuts[0][0]==0 and cuts[-1][1]==S
                assert all(a<b for a,b in cuts) and all(cuts[i][1]==cuts[i+1][0] for i in range(len(cuts)-1))
            assert local==list(map(int,sl['max_cell_scaled']))
            worst=[w if old is None else max(w,old) for w,old in zip(local,worst)]
        assert next(rows,None) is None
    assert count==run['joint_cells']; assert worst==list(map(int,run['max_cell_log_scaled']))
    total+=count; print(f"PASS {run['run']}: {count:,} cells, no missing cells, maxima agree")
assert total==s['total_joint_cells']==190031
print(f'PASS: complete archived partition ({total:,} joint cells)')
