import sys; sys.dont_write_bytecode=True
sys.argv=['x','1000']
import importlib.util, numpy as np
src=open('ls_curvbound.py').read(); src=src[:src.index("out = []")]
g={}; exec(compile(src,'lsc','exec'),g)
for b in (0.05,0.09,0.15,0.3):
    print(b, 2*g['Bk'](2,b), 2*g['Bk_env'](2,b))
x=np.array([50.,100.,300.,1000.])
print(g['Ek'](2,0,x), 190/x**2, g['Ek'](2,1,x))
print(g['K'])
