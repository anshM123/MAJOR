// Global signed integrated conical coefficients via triangular moment transport.
// Input: p K digits, followed by K+1 pairs of scaled integer midpoint/radius.
// Output: k midpoint radius of int g_{p,k}(Lambda*t) psi(t) dt.
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <stdexcept>
using I=mpz_class;
I ceildiv(const I& n, unsigned long d){I q;mpz_cdiv_q_ui(q.get_mpz_t(),n.get_mpz_t(),d);return q;}
int main(){int p,K,digits;if(!(std::cin>>p>>K>>digits)|| (p!=3&&p!=4)||K<0||digits<1) return 1;
std::vector<I> prev(K+1),cur(K+1),next(K+1);I ep=0,e=0,r;
for(int j=0;j<=K;j++){if(!(std::cin>>cur[j]>>r)||r<0)return 1;if(r>e)e=r;}
for(int k=0;k<=K;k++){std::cout<<k<<' '<<cur[0]<<' '<<e<<'\n';if(k==K)break;
unsigned long A=8ul*k*k+8ul*p*k+4*p-1,B=4ul*(k+p-1)*(k+p-1),C=4ul*997002,D=4ul*(k+1)*(k+1);
for(int j=0;j<K-k;j++){next[j]=(A*cur[j]-C*cur[j+1]-B*prev[j])/D;}
I en=ceildiv((A+C)*e+B*ep,D)+1;ep=e;e=en;prev.swap(cur);cur.swap(next);
}
return 0;}
