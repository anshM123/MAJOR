#include <gmpxx.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <stdexcept>
using Z=mpz_class;
int main(int argc,char**argv){std::ifstream f(argv[1]);int N,D;f>>N>>D;std::vector<Z>A(N+1),G(N+1),T[2]; T[0].resize(N);T[1].resize(N); Z r,g;for(int k=0;k<=N;k++)f>>A[k]>>r>>G[k]>>g;for(int i=0;i<2;i++){int b=4*i,q=b*(b+1),sg=i?-1:1;for(int k=0;k<N;k++){long k1=(long)(k+1)*(k+1),k2=-(long)k*k+(8-4*b)*(long)k+9-2*b*(b-1)+2*q;T[i][k]=k1*A[k+1]+k2*A[k]+2*sg*G[k];}}
unsigned long C,H,Den;int b,K;Z v,w,tmp;while(std::cin>>b>>K>>C>>H>>Den){if(Den>(1UL<<62)||K>=N||K<0||C>=Den||H>Den)throw std::runtime_error("bad inputs");int i=b/4,s=2*b+1;v=0;for(long k=K;k>=0;k--){unsigned long p=k*(k+s);if(k<K){if(H>0&&p>(Den-1)/H)throw std::runtime_error("pole");unsigned long z=H*p,rd=Den-z;bool negative=z>C;unsigned long rn=negative?z-C:C-z;mpz_mul_ui(v.get_mpz_t(),v.get_mpz_t(),rn);mpz_tdiv_q_ui(v.get_mpz_t(),v.get_mpz_t(),rd);if(negative)mpz_neg(v.get_mpz_t(),v.get_mpz_t());}mpz_mul_ui(w.get_mpz_t(),T[i][k].get_mpz_t(),H);mpz_tdiv_q_ui(w.get_mpz_t(),w.get_mpz_t(),2*Den);v+=w;v+=A[k];}std::cout<<v<<'\n';}}
