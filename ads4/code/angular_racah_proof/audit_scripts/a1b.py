import sympy as S
exec(open('a1.py').read().split('# 2: recurrence')[0].split("print('P_k")[0])
nu,Y=S.symbols('nu Y')
def yfun(xx,ss,nuv=nu,YY=Y):
    tot=0; c=S.Integer(1)
    for kk in range(0,xx+1):
        if kk>0:
            i=kk-1
            c = c*((i**2-nuv)/S.Integer(kk**2))*(xx*(xx+ss)-i*(i+ss))/(YY-i*(i+ss))
        tot+=c
    return tot
ss=1
for xx in range(0,10):
    Xx=xx*(xx+ss); q=(2*xx+ss)*nu/(Y-Xx)
    lhs = p(xx,ss)*(yfun(xx+1,ss)-yfun(xx,ss)) + q*yfun(xx,ss)
    if xx>0: lhs -= p(xx-1,ss)*(yfun(xx,ss)-yfun(xx-1,ss))
    num,den=S.fraction(S.cancel(S.together(lhs)))
    print(xx, S.factor(num) if num!=0 else 0)
