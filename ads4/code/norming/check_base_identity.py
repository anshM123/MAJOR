#!/usr/bin/env python3
"""Exact algebraic regression for the s=1 spin recurrence and s=9 shift identity."""
from fractions import Fraction as F
from math import comb
from check_identity import C,cast

def ladder(n,N,s):
    D=N*(N+s);v=[F(1)]
    for j in range(n):v.append(v[-1]*(n*(n+s)-j*(j+s))/(D-j*(j+s)))
    return v

def power(n,N,s,beta):
    v=ladder(n,N,s);c=cast(1);out=cast(1)
    for j in range(1,n+1):c=c*(beta+j-1)/j;out+=v[j]*c
    return out

def main():
    count=0;hel=0
    for N in [F(29,2),F(29),F(1001,1),F(2021,2)]:
        for beta in [C(F(1,4),F(17,3)),C(F(5,2),F(900)),C(F(-15,4),F(15,7))]:
            old=cast(1);cur=cast(1)
            for n in range(1,12):
                k=n-1;eta=F(2*k+1)/((N-k)*(N+k+1));nxt=old+2*eta*beta*cur
                assert nxt==power(n,N,1,beta);old,cur=cur,nxt;count+=1
            for n in (1,2,5,9):
                N0=N+4;n0=n+4;U=ladder(n0,N0,1)[4]
                assert power(n,N,9,beta)==sum((-1)**i*comb(4,i)*power(n0,N0,1,beta-i) for i in range(5))/U
                hel+=1
    print(f'PASS: {count} exact spin steps and {hel} exact helicity identities')
if __name__=='__main__':main()
