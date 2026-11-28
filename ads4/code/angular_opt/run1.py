from lp_primal import *
for g,n in [(10,10),(20,20)]:
    c2,res = solve(g,n,Jmax=6*g+10,verbose=True)
    print(g,n,c2, 4*np.sqrt(c2)*g**3)
    # dual
    print('dual eq marg', res.eqlin.marginals[:12])
