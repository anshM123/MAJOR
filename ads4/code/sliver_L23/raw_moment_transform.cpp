// Exact scaled raw CDH transform. Input: K digits, followed by K+1
// (midpoint,radius) pairs for integral g(nu)*(nu/1000)^(2j).
// Output m, midpoint, radius of integral g f_m, m=1,...,K+1.
#include <gmpxx.h>
#include <iostream>
#include <vector>
using I=mpz_class;
I ceildiv(const I&n,unsigned long d){I q;mpz_cdiv_q_ui(q.get_mpz_t(),n.get_mpz_t(),d);return q;}
int main(){int K,D;if(!(std::cin>>K>>D)||K<0||D<1)return 1;std::vector<I>prev(K+1),cur(K+1),next(K+1);I ep=0,e=0,r;for(int j=0;j<=K;j++){if(!(std::cin>>cur[j]>>r)||r<0)return 1;if(r>e)e=r;}for(int m=1;m<=K+1;m++){std::cout<<m<<' '<<cur[0]<<' '<<e<<'\n';if(m==K+1)break;unsigned long A=8ul*m*m-1,B=4ul*(m-1)*(m-1),C=4000000,Dn=4ul*(m+1)*(m+1);for(int j=0;j<=K-m;j++)next[j]=(A*cur[j]-C*cur[j+1]-B*prev[j])/Dn;I en=ceildiv((A+C)*e+B*ep,Dn)+1;ep=e;e=en;prev.swap(cur);cur.swap(next);}return 0;}
