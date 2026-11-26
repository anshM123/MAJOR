# Sampling prediction (box_cheb_h.py) for every (cell, h-piece, b) of a plan json: centre, L1/centre.
import sys,json,subprocess
P=json.load(open(sys.argv[1]))
for cell in P['cells']:
    cl,cu=cell['c_lo'],cell['c_hi'];m=cell['m'];K=cell['K']
    for b in (0,4):
        for i in range(m):
            fl,fu=i/m,(i+1)/m
            out=subprocess.run(['python3','box_cheb_h.py',str(float(eval(cl))),str(float(eval(cu))),str(b),str(K),'16','12',str(fl),str(fu)],capture_output=True,text=True).stdout.strip()
            print(out,flush=True)
